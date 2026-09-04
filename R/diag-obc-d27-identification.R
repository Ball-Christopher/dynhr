## R/diag-obc-d27-identification.R
## --------------------------------------------------------------------------
## Phase H: D27 — OBC/Piecewise-Linear Identification Diagnostic
##
## Assesses how occasionally binding constraints (OBCs) affect parameter
## identification. In a piecewise-linear (OccBin) setting with multiple
## regimes, the identification Jacobian and strength can differ across
## regimes.  Parameters that are well-identified in the slack regime may
## become weakly identified (or vice versa) under binding constraints.
##
## Algorithm:
##   1. Build the slack-regime and binding-regime system matrices from
##      the OBC model (using obc_build_binding_sys, obc_solve_binding).
##   2. For each regime, compute the identification Jacobian of the
##      model-implied moments w.r.t. parameters (as in D1/D19).
##   3. Compare identification strength across regimes — flag parameters
##      whose identification is regime-dependent.
##   4. Report the regime-dependent identification matrix and highlight
##      parameters in danger of losing identification under OBC.
##
## References:
##   Guerrieri & Iacoviello (2015). OccBin: A toolkit for solving dynamic
##     models with occasionally binding constraints easily.
##   Harrison & Waldron (2021). Optimal monetary policy with occasionally
##     binding constraints. J. Econ. Dyn. Control.
##   Iskrev, N. (2010). Local identification in DSGE models.
##   Qu, Z., & Tkachenko, D. (2012). Identification and frequency domain
##     analysis of DSGE models.
## --------------------------------------------------------------------------

#' D27. OBC/Piecewise-Linear Identification Diagnostic
#'
#' Assesses how occasionally binding constraints affect parameter
#' identification. Each OBC regime (slack vs. binding for each constraint)
#' induces a different linear state-space representation.  Parameters that
#' are well-identified in the slack regime may be weakly identified (or
#' entirely unidentified) when a constraint binds, and vice versa.
#'
#' The diagnostic:
#' \enumerate{
#'   \item Builds the system matrices for each OBC regime (slack + all
#'     single-constraint-binding regimes by default, or a user-specified
#'     subset).
#'   \item For each regime, computes the first-order perturbation solution
#'     and the resulting identification Jacobian (moments w.r.t. parameters).
#'   \item Compares identification strength (D20-style Fisher information)
#'     across regimes.
#'   \item Flags parameters whose identification rank or strength changes
#'     materially across regimes.
#' }
#'
#' @param model           A dynhr_mod object (must have OBC tags parsed).
#' @param dr              Decision rules object (slack-regime DecisionRules).
#' @param params          Named parameter vector at the calibration point.
#' @param compiled        Compiled model object (from \code{compile_model()}).
#' @param obc_specs       List of OBC specs (from \code{obc_parse_tags()} or
#'   \code{obc_collect_specs()}).  If NULL, attempts to parse MCP tags from
#'   \code{model}.
#' @param regime_subset   Integer vector of regime indices to evaluate.
#'   Default evaluates regime 0 (all slack) and all single-binding regimes
#'   (1, 2, 4, 8, ...). Set to \code{NULL} to evaluate all \code{2^k}
#'   regimes (may be expensive for large k).
#' @param param_names     Character vector of parameter names to analyse.
#'   If NULL, defaults to \code{names(params)}.
#' @param eps             Step size for finite-difference Jacobian (default 1e-5).
#' @param strength_threshold Threshold for D20-style \code{|t|-ratio} below
#'   which a parameter is flagged as weakly identified (default 1.0).
#' @param strength_ratio_threshold Ratio of max/min strength across regimes
#'   above which a parameter is flagged as regime-dependent (default 3.0).
#' @param verbose         Print progress messages.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{List containing:
#'     \itemize{
#'       \item \code{regime_jacobians} — named list of Jacobian matrices,
#'         one per regime.
#'       \item \code{regime_strength} — matrix (n_param x n_regimes) of
#'         identification strength (|t|-ratios).
#'       \item \code{regime_ranks} — integer vector of rank per regime.
#'       \item \code{regime_dependent_params} — character vector of
#'         parameters flagged as regime-dependent.
#'       \item \code{weak_params_by_regime} — named list of weakly
#'         identified parameters per regime.
#'       \item \code{regime_labels} — character vector describing each regime.
#'     }}
#'   \item{pass}{Logical — TRUE if no parameters are both regime-dependent
#'     AND weakly identified in any regime.}
#'   \item{plots}{List of ggplot2 objects (regime comparison heatmap).}
#'   \item{summary}{Human-readable summary.}
#'
#' @noRd
d27_obc_identification <- function(model,
                                    dr,
                                    params,
                                    compiled = NULL,
                                    obc_specs = NULL,
                                    regime_subset = NULL,
                                    param_names = NULL,
                                    eps = 1e-5,
                                    strength_threshold = 1.0,
                                    strength_ratio_threshold = 3.0,
                                    verbose = FALSE,
                                    meta = NULL) {
  # ---- 1. Validate and defaults ----
  if (!inherits(model, "dynhr_mod")) {
    return(.make_result(
      pass    = NA,
      summary = "D27 OBC Identification: model must be a dynhr_mod object."
    ))
  }
  if (is.null(params)) {
    return(.make_result(
      pass    = NA,
      summary = "D27 OBC Identification: params is required."
    ))
  }
    if (is.null(param_names)) param_names <- names(params)
    n_par <- length(param_names)

    obc_specs <- obc_specs %||% obc_parse_tags(model)
    if (is.null(obc_specs) || length(obc_specs) == 0L) {
      return(.make_result(
        result  = NULL,
        pass    = NA,
        plots   = list(),
        summary = paste(
          "D27 OBC Identification: No OBC constraints detected or provided.",
          "Use obc_parse_tags() or the @#obc MCP annotation mechanism",
          "to define constraints, then re-run."
        ),
        llm_summary = paste(
          "[INFO] D27 OBC Identification status=skipped reason=no_obc_specs",
          "action: define OBC constraints via obc_parse_tags"
        )
      ))
    }

    k <- length(obc_specs)
    n_regimes_total <- 2L ^ k
    if (verbose) cat(sprintf("[d27] %d OBC spec(s) detected, %d total regimes.\n", k, n_regimes_total))

    # ---- 2. Resolve system matrices from DR or compiled ----
    sys <- NULL
    if (!is.null(dr) && !is.null(dr$sys_mat)) {
      sys <- dr$sys_mat
    } else if (!is.null(compiled) && !is.null(compiled$sys_mat)) {
      sys <- compiled$sys_mat
    } else if (!is.null(compiled) && !is.null(dr)) {
      # Build cache from compiled model structure
      sc <- cache_system_structure(compiled)
      ss_vals <- dr$ys %||% rep(0, nrow(dr$ghx))
      sys <- extract_system_matrices_fast(sc, ss_vals, params)
    } else {
      # No system matrices available — fall through to numerical Jacobian
      sys <- NULL
    }

    if (is.null(sys)) {
      # Build a minimal solve_fn for numerical Jacobian instead
      if (verbose) cat("[d27] No system matrices available; using numerical Jacobian approach.\n")
      return(.build_d27_via_solve_fn(model, params, obc_specs, dr, compiled,
                                      regime_subset, param_names, eps,
                                      strength_threshold, strength_ratio_threshold,
                                      verbose))
    }

    # ---- 3. Determine regime subset ----
    if (is.null(regime_subset)) {
      # Default: regime 0 (all slack) + all single-binding regimes
      regime_subset <- c(0L, 2L ^ (seq_len(k) - 1L))
    }
    regime_subset <- unique(as.integer(regime_subset))
    regime_subset <- regime_subset[regime_subset >= 0 & regime_subset < n_regimes_total]
    n_regimes <- length(regime_subset)

    # Build regime labels
    regime_labels <- vapply(regime_subset, function(idx) {
      if (idx == 0L) return("all slack")
      flags <- obc_regime_flags(idx, k)
      binding_names <- sprintf("%s(bind)", obc_specs$name %||% paste0("spec", seq_len(k)))[flags]
      paste(binding_names, collapse = " + ")
    }, character(1))
    names(regime_subset) <- regime_labels

    if (verbose) cat(sprintf("[d27] Evaluating %d regime(s): %s\n",
                             n_regimes, paste(regime_labels, collapse = ", ")))

    # ---- 4. Per-regime identification ----
    regime_jacobians <- vector("list", n_regimes)
    regime_strength  <- matrix(NA_real_, nrow = n_par, ncol = n_regimes,
                               dimnames = list(param_names, regime_labels))
    regime_ranks     <- integer(n_regimes)
    names(regime_ranks) <- regime_labels

    # Build a solve function that returns moments for a given theta
    # under a specific regime
    for (r in seq_len(n_regimes)) {
      regime_idx <- regime_subset[r]
      label <- regime_labels[r]

      if (verbose) cat(sprintf("[d27]   Regime %d (%s)...", regime_idx, label))

      solve_fn <- .make_obc_solve_fn(model, params, obc_specs, dr, compiled,
                                      regime_idx)

      J <- .numerical_jacobian(function(th) {
        solve_fn(th)
      }, params[param_names], eps = eps)

      if (is.null(J)) {
        if (verbose) cat(" FAILED\n")
        regime_ranks[r] <- NA_integer_
        next
      }

      colnames(J) <- param_names
      regime_jacobians[[r]] <- J
      sv_J <- svd(J)
      regime_ranks[r] <- .svd_rank(sv_J$d, dim(J))

      # Strength (D20-style: |t|-ratio = |theta| / se)
      I_raw <- crossprod(J)
      I_reg <- I_raw + diag(1e-10, n_par)
      # Singular information matrix (collinear shock-std parameters) falls back
      # to a pseudo-inverse inside .safe_sym_inv().
      I_inv <- .safe_sym_inv(I_reg)
      se <- sqrt(pmax(diag(I_inv), 0))
      regime_strength[, r] <- abs(params[param_names]) / pmax(se, 1e-16)

      if (verbose) cat(sprintf(" rank=%d\n", regime_ranks[r]))
    }

    # ---- 5. Cross-regime analysis ----
    # Find regime-dependent parameters
    regime_dependent <- character(0)
    weak_by_regime <- list()

    for (i in seq_len(n_par)) {
      pname <- param_names[i]
      s_i <- regime_strength[i, ]
      s_valid <- s_i[is.finite(s_i) & s_i > 0]

      # Check which regimes have weak identification
      weak_regimes <- which(s_i < strength_threshold)
      if (length(weak_regimes) > 0) {
        weak_by_regime[[pname]] <- regime_labels[weak_regimes]
      }

      # Check regime dependence
      if (length(s_valid) >= 2) {
        s_ratio <- max(s_valid) / min(s_valid)
        if (s_ratio > strength_ratio_threshold) {
          regime_dependent <- c(regime_dependent, pname)
        }
      }
    }

    # Pass: no parameter is both regime-dependent AND weakly identified anywhere
    pass <- TRUE
    for (pname in regime_dependent) {
      if (pname %in% names(weak_by_regime)) {
        pass <- FALSE
        break
      }
    }

    # ---- 6. Build result and plots ----
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
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        ) +
        ggplot2::labs(
          title = "D27: OBC regime-dependent identification strength",
          subtitle = sprintf("%d OBC spec(s), %d regime(s)", k, n_regimes),
          x = NULL, y = NULL
        )
      plots$regime_strength <- .apply_meta(p_rs, meta)
    }

    # Build result list
    result <- list(
      regime_jacobians          = regime_jacobians,
      regime_strength           = regime_strength,
      regime_ranks              = regime_ranks,
      regime_dependent_params   = regime_dependent,
      weak_params_by_regime     = weak_by_regime,
      regime_labels             = regime_labels,
      regime_indices            = regime_subset,
      n_obc_specs               = k
    )

    # Summary text
    min_rank <- min(regime_ranks, na.rm = TRUE)
    max_rank <- max(regime_ranks, na.rm = TRUE)
    rank_str <- if (min_rank == max_rank) sprintf("rank=%d", min_rank)
                else sprintf("rank range [%d, %d]", min_rank, max_rank)

    dep_str <- if (length(regime_dependent) > 0)
      sprintf("Regime-dependent: %s", paste(regime_dependent, collapse = ", "))
    else "No regime-dependent parameters detected."

    summary_str <- sprintf(
      "D27 OBC Identification: %d spec(s), %d regime(s). %s. %s",
      k, n_regimes, rank_str, dep_str
    )

    .make_result(
      result  = result,
      pass    = pass,
      plots   = plots,
      summary = summary_str,
      llm_summary = sprintf(
        "[%s] D27 OBC Identification n_obc_specs=%d n_regimes=%d min_rank=%d max_rank=%d n_regime_dependent=%d",
        if (isTRUE(pass)) "PASS" else if (is.na(pass)) "INFO" else "FAIL",
        k, n_regimes, min_rank, max_rank, length(regime_dependent)
      )
    )
}


# ==========================================================================
# Internal helpers for D27
# ==========================================================================

#' Build a regime-specific solve function for D27 numerical Jacobian
#'
#' Creates a closure that, given a parameter vector theta, solves the OBC
#' model under the specified regime and returns a moment vector (variances
#' and selected covariances of observables).  When sys matrices are available
#' (via dr$sys_mat or compiled$sys_mat), uses obc_ensure_policy to obtain
#' the regime-specific decision rules (slack vs binding).  Otherwise falls
#' back to re-solving the perturbation for the slack regime only.
#'
#' @param model      dynhr_mod
#' @param params     Full named parameter vector
#' @param obc_specs  OBC spec list
#' @param dr         Slack-regime DecisionRules
#' @param compiled   Compiled model object
#' @param regime_idx Integer regime index (0 = slack, >0 = binding subset)
#' @return Function: theta -> named numeric moment vector
#' @noRd
.make_obc_solve_fn <- function(model, params, obc_specs, dr, compiled,
                                regime_idx) {
  # Resolve sys matrices once at closure-creation time
  sys <- NULL
  obs_idx <- NULL
  if (!is.null(dr) && !is.null(dr$sys_mat)) {
    sys <- dr$sys_mat
  } else if (!is.null(compiled) && !is.null(compiled$sys_mat)) {
    sys <- compiled$sys_mat
  }
  # Try to build sys from compiled if not pre-stored
  if (is.null(sys) && !is.null(compiled) && !is.null(dr)) {
    sc <- cache_system_structure(compiled)
    ss_vals <- dr$ys %||% rep(0, nrow(dr$ghx))
    sys <- extract_system_matrices_fast(sc, ss_vals, params)
  }
  # obs_idx: all variables for full endo moments, or try dr$obs_idx
  if (!is.null(dr$obs_idx)) {
    obs_idx <- dr$obs_idx
  } else if (!is.null(model$obs_mat)) {
    obs_idx <- seq_len(nrow(model$obs_mat))
  }

  # Cache lives in the enclosing function environment so it persists across
  # Jacobian column evaluations (different theta perturbations of the same regime).
  regime_cache <- new.env(parent = emptyenv(), hash = TRUE)

  function(theta) {
    params_new <- params
    params_new[names(theta)] <- theta

    if (is.null(dr)) return(NULL)
    # Obtain regime-specific decision rules
    if (!is.null(sys) && regime_idx > 0L) {
        # Use occbin binding-regime solver via obc_ensure_policy
        obc_ensure_policy(regime_idx, regime_cache, sys, dr, obc_specs, obs_idx)
        entry <- regime_cache[[as.character(regime_idx)]]
        if (is.null(entry) || is.null(entry$dr)) return(NULL)
        dr_regime <- entry$dr
      } else {
        # Slack regime: re-solve perturbation at new params
        dr_regime <- .stoch_simul_internal_diag(model, params_new, dr_order = 1L)
        if (is.null(dr_regime) || is.null(dr_regime$ghx)) return(NULL)
      }

      # Compute moments from regime-specific DR
      moments <- .moments_from_dr(dr_regime, model = model, params = params_new)
      if (is.null(moments)) return(NULL)
      if (is.null(moments$sigma_y) || is.null(moments$acf_y)) return(NULL)

      # Guard: acf_y may be 2-D or 3-D depending on compute_moments version
      acf_y_diag <- if (length(dim(moments$acf_y)) == 3L) {
        diag(moments$acf_y[, , 1])
      } else {
        diag(as.matrix(moments$acf_y))
      }

      m <- c(
        log_var  = log(pmax(diag(moments$sigma_y), 1e-16)),
        log_acf1 = log(pmax(acf_y_diag, 1e-16))
      )
      if (length(m) == 0) return(NULL)
      names(m) <- paste0("m_", seq_along(m))
      m
  }
}


#' Alternative D27 implementation using numerical solve_fn
#'
#' Used when system matrices are not directly available. Builds a
#' model_solve_fn compatible with the D1/D19 numerical Jacobian approach.
#'
#' @noRd
.build_d27_via_solve_fn <- function(model, params, obc_specs, dr, compiled,
                                     regime_subset, param_names, eps,
                                     strength_threshold, strength_ratio_threshold,
                                     verbose) {
  k <- length(obc_specs)

  if (is.null(regime_subset)) {
    regime_subset <- c(0L, 2L ^ (seq_len(k) - 1L))
  }
  regime_subset <- unique(as.integer(regime_subset))
  n_regimes <- length(regime_subset)

  regime_labels <- vapply(regime_subset, function(idx) {
    if (idx == 0L) return("all slack")
    flags <- obc_regime_flags(idx, k)
    binding_names <- sprintf("spec%d", which(flags))
    paste(binding_names, collapse = " + ")
  }, character(1))

  # Compute identification per regime
  regime_strength <- matrix(NA_real_, nrow = length(param_names), ncol = n_regimes,
                            dimnames = list(param_names, regime_labels))
  regime_ranks <- integer(n_regimes)

  for (r in seq_len(n_regimes)) {
    regime_idx <- regime_subset[r]
    solve_fn <- .make_obc_solve_fn(model, params, obc_specs, dr, compiled,
                                    regime_idx)

    # Test: solve_fn must return a non-empty vector
    test_out <- solve_fn(params[param_names])
    if (is.null(test_out) || length(test_out) == 0L) {
      if (verbose) cat(sprintf("  Regime %d solve_fn returned empty, skipping.\n", regime_idx))
      regime_ranks[r] <- NA_integer_
      next
    }

    J <- .numerical_jacobian(function(th) solve_fn(th), params[param_names], eps = eps)

    if (is.null(J) || nrow(J) == 0L) {
      regime_ranks[r] <- NA_integer_
      next
    }

    colnames(J) <- param_names
    sv_J <- svd(J)
    regime_ranks[r] <- .svd_rank(sv_J$d, dim(J))

    I_raw <- crossprod(J)
    np_ <- length(param_names)
    I_reg <- I_raw + diag(1e-10, np_)
    I_inv <- .safe_sym_inv(I_reg)
    se <- sqrt(pmax(diag(I_inv), 0))
    regime_strength[, r] <- abs(params[param_names]) / pmax(se, 1e-16)
  }

  # Cross-regime analysis (same as main function)
  regime_dependent <- character(0)
  weak_by_regime <- list()

  for (i in seq_along(param_names)) {
    pname <- param_names[i]
    s_i <- regime_strength[i, ]
    s_valid <- s_i[is.finite(s_i) & s_i > 0]

    weak_regimes <- which(s_i < strength_threshold)
    if (length(weak_regimes) > 0) {
      weak_by_regime[[pname]] <- regime_labels[weak_regimes]
    }
    if (length(s_valid) >= 2) {
      if (max(s_valid) / min(s_valid) > strength_ratio_threshold) {
        regime_dependent <- c(regime_dependent, pname)
      }
    }
  }

  pass <- TRUE
  for (pname in regime_dependent) {
    if (pname %in% names(weak_by_regime)) {
      pass <- FALSE
      break
    }
  }

  .make_result(
    result = list(
      regime_strength           = regime_strength,
      regime_ranks              = regime_ranks,
      regime_dependent_params   = regime_dependent,
      weak_params_by_regime     = weak_by_regime,
      regime_labels             = regime_labels,
      regime_indices            = regime_subset
    ),
    pass = pass,
    plots = list(),
    summary = sprintf(
      "D27 OBC Identification (%d regimes): ranks = [%s]%s",
      n_regimes,
      paste(regime_ranks, collapse = ", "),
      if (length(regime_dependent) > 0)
        sprintf("; regime-dependent: %s", paste(regime_dependent, collapse = ", "))
      else ""
    )
  )
}
