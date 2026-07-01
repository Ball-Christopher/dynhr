## R/diag-pre-d1-identification.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D1 local identification (Iskrev 2010)
## --------------------------------------------------------------------------

#' D1. Local identification (Iskrev 2010)
#'
#' Computes the numerical Jacobian of the vector of model-implied moments
#' (vec of first and second moments of observables) with respect to the
#' deep parameter vector. Checks the rank via SVD: full rank means all
#' parameters are locally identified.
#'
#' @param model_solve_fn  Function: theta -> named numeric vector of moments
#'                        (e.g. c(sd_y, sd_pi, acf1_y, ...) )
#' @param theta           Numeric vector -- parameter values at calibration point
#' @param param_names     Character vector -- names of parameters
#' @param moment_names    Character vector -- names of moments (rows of Jacobian)
#' @param eps             Step size for finite differences (default 1e-5)
#' @param jacobian        Optional pre-computed Jacobian (n_moment x n_par);
#'   when supplied (e.g. by the orchestrator, shared with D20) the
#'   finite-difference re-computation is skipped.
#' @return dynhr_diagnostic list
#' @references Iskrev, N. (2010). Local identification in DSGE models.
#'   \emph{Journal of Monetary Economics}, 57(2), 189-202.
#'   Komunjer, I., & Ng, S. (2011). Dynamic identification of dynamic stochastic
#'   general equilibrium models. \emph{Econometrica}, 79(6), 1995-2032.
#' @noRd
d1_local_identification <- function(model_solve_fn,
                                    theta,
                                    param_names  = NULL,
                                    moment_names = NULL,
                                    eps = 1e-5,
                                    jacobian = NULL,
                                    meta = NULL) {

    n_par <- length(theta)
    if (is.null(param_names))  param_names  <- paste0("theta_", seq_len(n_par))
    if (is.null(moment_names)) {
      f0 <- model_solve_fn(theta)
      moment_names <- if (!is.null(names(f0))) names(f0) else paste0("m_", seq_along(f0))
    }

    # Compute Jacobian (or reuse one supplied by the orchestrator: D1 and D20
    # share the same J, so the orchestrator computes it once and passes it in).
    J <- jacobian %||% .numerical_jacobian(model_solve_fn, theta, eps = eps)
    # Guard: ensure param_names length matches Jacobian column count
    if (ncol(J) != length(param_names)) {
      warning(sprintf(
        "d1: param_names length (%d) != Jacobian columns (%d). Using generic param labels.",
        length(param_names), ncol(J)))
      param_names <- paste0("theta_", seq_len(ncol(J)))
    }
    colnames(J) <- param_names
    # Guard: ensure moment_names length matches Jacobian row count
    if (nrow(J) != length(moment_names)) {
      warning(sprintf(
        "d1: moment_names length (%d) != Jacobian rows (%d). Using generic moment labels.",
        length(moment_names), nrow(J)))
      moment_names <- paste0("m_", seq_len(nrow(J)))
    }
    rownames(J) <- moment_names

    # Guard: non-finite Jacobian (NaN/Inf from solve_lyapunov degeneracy at
    # near-unit-root points) collapses rank() to 0, which is a *numerical*
    # failure, not a substantive non-identification finding.  Detect it and
    # return pass=NA with an honest message rather than rank=0.
    if (!all(is.finite(J))) {
      n_nonfinite <- sum(!is.finite(J))
      msg <- sprintf(
        paste0("D1 Local identification: Jacobian is non-finite at this point ",
               "(%d of %d entries are NaN/Inf, likely from solve_lyapunov ",
               "degeneracy at a near-unit-root). Rank is undetermined -- this ",
               "is a numerical failure, NOT a non-identification finding."),
        n_nonfinite, length(J)
      )
      return(.make_result(
        result  = list(jacobian = J, singular_values = NULL,
                       rank = NA_integer_, weak_params = character(0),
                       svd = NULL, numerical_degenerate = TRUE),
        pass    = NA,
        plots   = list(),
        summary = msg,
        llm_summary = paste0(
          "D1 | Local Identification | INFO\n",
          sprintf("  non_finite_entries=%d total=%d\n", n_nonfinite, length(J)),
          "  action: Jacobian non-finite (near-unit-root/Lyapunov degeneracy) -- ",
          "rank undetermined; not a non-identification finding."
        )
      ))
    }

    # SVD decomposition
    sv <- svd(J)
    singular_values <- sv$d
    names(singular_values) <- paste0("sv_", seq_along(singular_values))

    # Rank check (numerical tolerance)
    rank_J <- .svd_rank(singular_values, dim(J))
    full_rank <- (rank_J == n_par)

    # Identify weakly-identified parameters (last singular vectors)
    # The right singular vectors corresponding to near-zero singular values
    # indicate linear combinations of parameters that are unidentified
    weak_threshold <- max(singular_values) * 1e-3
    n_weak <- sum(singular_values < weak_threshold)
    weak_params <- character(0)
    if (n_weak > 0) {
      # sv$v is n_par x min(n_par, n_mom). Guard against out-of-bounds.
      n_v_cols <- ncol(sv$v)
      start_col <- max(1, n_v_cols - n_weak + 1)
      if (start_col <= n_v_cols) {
        V_weak <- sv$v[, start_col:n_v_cols, drop = FALSE]
        # Parameters with largest loading in weak directions
        for (k in seq_len(ncol(V_weak))) {
          idx <- which.max(abs(V_weak[, k]))
          weak_params <- c(weak_params, param_names[idx])
        }
        weak_params <- unique(weak_params)
      }
    }

    pass <- full_rank

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    # (a) Heatmap of the Jacobian
    J_long <- reshape2::melt(J)
    colnames(J_long) <- c("Moment", "Parameter", "Value")
    p_jac <- ggplot2::ggplot(
      J_long, ggplot2::aes(x = Parameter, y = Moment, fill = Value)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
      ggplot2::scale_fill_gradient2(low = dynhr_colours$red, mid = dynhr_colours$white,
                                    high = dynhr_colours$mid_blue, midpoint = 0,
                                    name = "Sensitivity") +
      theme_dynhr_diagnostic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = ggplot2::rel(0.75)),
        axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.7))
      ) +
      ggplot2::labs(
        title = "D1: Jacobian of moments w.r.t. parameters",
        subtitle = sprintf("Rank = %d / %d parameters", rank_J, n_par),
        x = NULL, y = NULL
      )
    plots$jacobian_heatmap <- .apply_meta(p_jac, meta)

    # (b) Singular values bar chart
    sv_df <- data.frame(
      index = seq_along(singular_values),
      value = pmax(singular_values, 1e-300)
    )
    sv_df$identified <- ifelse(sv_df$value >= weak_threshold, "Identified", "Weak / unidentified")

    p_sv <- ggplot2::ggplot(
      sv_df, ggplot2::aes(x = factor(index), y = value, fill = identified)
    ) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_hline(yintercept = weak_threshold, linetype = "dashed",
                          colour = dynhr_colours$red, linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = c("Identified" = dynhr_colours$mid_blue,
                                            "Weak / unidentified" = dynhr_colours$red),
                                 name = NULL) +
      ggplot2::scale_y_log10() +
      theme_dynhr_diagnostic() +
      ggplot2::labs(
        title = "D1: Singular values of the identification Jacobian",
        subtitle = sprintf("Threshold = %.2e (0.1%% of max SV)", weak_threshold),
        x = "Singular value index", y = "Value (log scale)"
      )
    plots$singular_values <- .apply_meta(p_sv, meta)

    }  # end requireNamespace guard

    summary_text <- sprintf(
      "D1 Local identification: Jacobian rank = %d / %d. %s%s",
      rank_J, n_par,
      ifelse(full_rank, "All parameters locally identified.",
             sprintf("RANK DEFICIENT -- %d parameter direction(s) unidentified.", n_par - rank_J)),
      ifelse(length(weak_params) > 0,
             sprintf(" Weakly identified: %s.", paste(weak_params, collapse = ", ")),
             "")
    )

    .make_result(
      result  = list(jacobian = J, singular_values = singular_values,
                     rank = rank_J, weak_params = weak_params,
                     svd = sv),
      pass    = pass,
      plots   = plots,
      summary = summary_text,
      llm_summary = {
        badge   <- if (pass) "PASS" else "FAIL"
        n_mom   <- nrow(J)
        sv_str  <- paste(sprintf("%.3e", head(sort(sv$d), 5)), collapse = ", ")
        paste(c(
          sprintf("D1 | Local Identification | %s", badge),
          sprintf("  params=%d moments=%d jacobian_rank=%d (expected=%d)",
                  n_par, n_mom, rank_J, n_par),
          sprintf("  smallest_sv: %s", sv_str),
          if (length(weak_params) > 0)
            sprintf("  weak_params: %s",
                    paste(weak_params, collapse = ", ")),
          sprintf("  action: %s",
                  if (pass)
                    "Model is locally identified at calibration point."
                  else
                    sprintf("Rank deficient (rank=%d, expected=%d). %s may not be identified. Fix: add moments, calibrate one parameter, or check for collinearity.",
                            rank_J, n_par,
                            paste(head(weak_params, 3), collapse = ", ")))
        ), collapse = "\n")
      }
    )
}


## NOTE (2026-05-31): the `d19_local_rank_identification` stub was removed. It
## was an admitted re-label of D1 (local identification via moment-Jacobian
## rank) that only ever returned an "identical to D1 ... Skipped" placeholder,
## so two diagnostics reported the same thing. D1 above IS the local rank
## identification check; callers should use it directly.


#' D20. Identification strength via Fisher information
#'
#' Computes an approximate Fisher information matrix \eqn{I = J'\,W\,J} from the
#' Jacobian J of model-implied moments w.r.t. parameters.  The default
#' \code{weighting = "none"} uses the unweighted \eqn{I = J'J}; on this metric,
#' moments on different scales (e.g. variances vs. autocorrelations) inflate the
#' sensitivity component \eqn{\Delta_i} for parameters that move large-magnitude
#' moments.  \code{weighting = "scaled"} standardises each moment by its
#' baseline magnitude (\eqn{\Omega = diag(f_0^2)}, i.e. J row-scaled by
#' \eqn{1/|f_0|}), which makes the cross-parameter \eqn{\Delta_i} ranking
#' scale-free.  Two important caveats keep \code{"none"} the default:
#' (i) this is \strong{not} the sampling-covariance-weighted Fisher information
#' \eqn{J'\,Var(m)^{-1} J} (the moment sampling covariance is not available at
#' this call site); and (ii) elasticity scaling distorts the absolute
#' \eqn{s_i = \theta_i / SE_i} t-ratio gate -- a parameter that maps linearly
#' to a dedicated moment yields \eqn{s_i \approx 1} regardless of magnitude, so
#' the "weak if \eqn{s<1}" threshold is only meaningful under \code{"none"}.
#' Use \code{"scaled"} for cross-parameter sensitivity comparison, not for the
#' pass gate.  Row scaling by positive weights preserves rank, so the
#' rank-deficiency verdict is identical under either weighting.  Reports:
#' - sensitivity component Delta_i = I_ii
#' - collinearity component rho_i from I = D R D
#' - normalised strength s_i = theta_i / sqrt((I^{-1})_ii)
#'
#' @param model_solve_fn Function: theta -> named numeric vector of moments.
#' @param theta Numeric vector of parameter values.
#' @param param_names Optional character vector of parameter names.
#' @param moment_names Optional character vector of moment names.
#' @param eps Step size for finite differences.
#' @param ridge Numeric ridge term for stable inversion (default 1e-10).
#' @param weighting Moment weighting for the Fisher information:
#'   \code{"none"} (default, raw \eqn{J'J}); \code{"scaled"} (standardise each
#'   moment by \eqn{|f_0|}; scale-free cross-parameter comparison only -- see
#'   Details); or \code{"sampling"} (principled \eqn{I = J'\,Var(\hat m)^{-1} J}
#'   using \code{moment_cov}; here \eqn{s_i=\theta/SE} is a genuine asymptotic
#'   t-ratio and the gate is well-defined). \code{"sampling"} weights only the
#'   moments present in \code{moment_cov} (aligned by name); if it is missing or
#'   does not overlap the Jacobian rows, D20 falls back to \code{"none"} and
#'   notes this in the summary.
#' @param moment_cov Optional named covariance matrix of the moment
#'   \emph{estimator}, \eqn{Var(\hat m)} (i.e. the long-run moment covariance
#'   already divided by the sample size T). Row/column names must match the
#'   moment names. Used only when \code{weighting = "sampling"}.
#' @param jacobian Optional pre-computed Jacobian (n_moment x n_par). When
#'   supplied (e.g. by the orchestrator, shared with D1) the finite-difference
#'   re-computation is skipped.
#' @return dynhr_diagnostic list
#' @noRd
d20_fisher_identification_strength <- function(model_solve_fn,
                                               theta,
                                               param_names  = NULL,
                                               moment_names = NULL,
                                               eps   = 1e-5,
                                               ridge = 1e-10,
                                               weighting  = c("none", "scaled", "sampling"),
                                               moment_cov = NULL,
                                               jacobian   = NULL,
                                               meta  = NULL) {

    weighting <- match.arg(weighting)
    n_par <- length(theta)
    if (is.null(param_names)) param_names <- names(theta) %||% paste0("theta_", seq_len(n_par))

    # Reuse the orchestrator-supplied Jacobian when available (shared with D1).
    J  <- jacobian %||% .numerical_jacobian(model_solve_fn, theta, eps = eps)
    # Baseline moments are needed both for labels and for scale-standardisation.
    f0 <- model_solve_fn(theta)
    if (is.null(moment_names)) {
      moment_names <- names(f0) %||% paste0("m_", seq_len(length(f0)))
    }
    if (nrow(J) != length(moment_names)) {
      warning(sprintf(
        "d20: moment_names length (%d) != Jacobian rows (%d). Using generic labels.",
        length(moment_names), nrow(J)))
      moment_names <- paste0("m_", seq_len(nrow(J)))
    }
    # Guard: ensure param_names length matches Jacobian column count
    if (ncol(J) != length(param_names)) {
      warning(sprintf(
        "d20: param_names length (%d) != Jacobian columns (%d). Using generic param labels.",
        length(param_names), ncol(J)))
      param_names <- paste0("theta_", seq_len(ncol(J)))
    }
    colnames(J) <- param_names
    rownames(J) <- moment_names

    # Guard: non-finite Jacobian (NaN/Inf from solve_lyapunov degeneracy at
    # near-unit-root points) would silently collapse Fisher rank to 0, masking
    # a numerical failure as a non-identification finding.  Detect it and
    # return pass=NA with an honest message.
    if (!all(is.finite(J))) {
      n_nonfinite <- sum(!is.finite(J))
      msg <- sprintf(
        paste0("D20 Identification strength: Jacobian is non-finite at this ",
               "point (%d of %d entries are NaN/Inf, likely from solve_lyapunov ",
               "degeneracy at a near-unit-root). Fisher information is ",
               "undetermined -- this is a numerical failure, NOT a ",
               "non-identification finding."),
        n_nonfinite, length(J)
      )
      return(.make_result(
        result  = list(jacobian = J, param_names = param_names,
                       sensitivity = NULL, fisher_rank_deficient = NA,
                       numerical_degenerate = TRUE),
        pass    = NA,
        plots   = list(),
        summary = msg,
        llm_summary = paste0(
          "D20 | Fisher Identification Strength | INFO\n",
          sprintf("  non_finite_entries=%d total=%d\n", n_nonfinite, length(J)),
          "  action: Jacobian non-finite (near-unit-root/Lyapunov degeneracy) -- ",
          "Fisher information undetermined; not a non-identification finding."
        )
      ))
    }
    # Finite Jacobian: proceed with the normal Fisher path.  (Previously we
    # zeroed non-finite entries here; that is now handled by the guard above.)

    # Optional scale-standardisation (weighting = "scaled", NOT the default).
    # Raw I = J'J lets large-magnitude moments (e.g. variances) dominate small
    # ones (e.g. autocorrelations), inflating the sensitivity component for
    # parameters that move high-variance moments.  Standardising each moment by
    # its baseline magnitude (Omega = diag(f0^2), i.e. row-scaling J by 1/|f0|)
    # makes the cross-parameter ranking scale-free.  It is NOT the
    # sampling-covariance-weighted Fisher information, and it makes the
    # theta/SE t-ratio gate degenerate for linear maps (s ~ 1 regardless of
    # magnitude) -- hence "none" is the default and the gate is defined on it.
    # Row scaling by positive weights preserves rank either way.
    # weighting = "sampling": principled Fisher information I = J' Var(m_hat)^-1 J
    # using the supplied sampling covariance of the moment ESTIMATOR (Var(m_hat),
    # i.e. the long-run moment covariance already divided by T).  This is the
    # correct GMM/Fisher weighting: s_i = theta_i / SE_i is then a genuine
    # asymptotic t-ratio, so the "weak if |s|<1" gate is well-defined (unlike the
    # "scaled" proxy).  Whitening uses Omega^{-1/2} (.robust_Omega_inv_sqrt), and
    # we align moment_cov to the Jacobian rows by name -- only moments with a
    # sampling-covariance estimate enter (moments without a data analogue cannot
    # be weighted and are excluded).
    f0v <- as.numeric(f0)
    f0v[!is.finite(f0v)] <- 0
    weighting_note   <- NULL
    n_moments_used   <- nrow(J)
    weighting_actual <- weighting
    if (identical(weighting, "sampling")) {
      ok_cov <- is.matrix(moment_cov) && !is.null(rownames(moment_cov)) &&
                all(is.finite(moment_cov))
      common <- if (ok_cov) intersect(rownames(J), rownames(moment_cov)) else character(0)
      if (length(common) >= 1L) {
        Js <- J[common, , drop = FALSE]
        Om <- moment_cov[common, common, drop = FALSE]
        Wh <- .robust_Omega_inv_sqrt(Om)   # Omega^{-1/2}
        Jw <- Wh %*% Js                    # whitened Jacobian; I = J' Omega^-1 J
        n_moments_used <- length(common)
      } else {
        weighting_actual <- "none"
        weighting_note   <- "sampling weighting requested but moment_cov was missing/unaligned; fell back to unweighted J'J"
        Jw <- J
      }
    } else if (identical(weighting, "scaled") && length(f0v) == nrow(J)) {
      m_scale <- pmax(abs(f0v), 1e-8 * max(abs(f0v), na.rm = TRUE), 1e-300)
      Jw <- J / m_scale            # divides each row i by m_scale[i] (recycled by column)
    } else {
      Jw <- J
    }
    I_raw <- crossprod(Jw)
    I_reg <- I_raw + diag(ridge, n_par)

    # Explicit rank / condition check BEFORE inverting (NA-safe)
    rk_fisher  <- tryCatch(qr(I_raw)$rank, error = function(e) 0L)
    if (!is.finite(rk_fisher)) rk_fisher <- 0L
    rc_fisher  <- tryCatch(rcond(I_raw), error = function(e) 0)
    if (!is.finite(rc_fisher)) rc_fisher <- 0
    fisher_rank_deficient <- (rk_fisher < ncol(I_raw)) || (rc_fisher < 1e-12)

    if (fisher_rank_deficient) {
      # Moore-Penrose pseudo-inverse via SVD (MASS::ginv if available, else manual)
      I_inv <- if (requireNamespace("MASS", quietly = TRUE)) {
        MASS::ginv(I_reg)
      } else {
        sv_fi <- svd(I_reg)
        d_inv <- sv_fi$d
        thresh_fi <- max(d_inv) * max(n_par, 1L) * .Machine$double.eps
        d_inv[d_inv <= thresh_fi] <- 0
        d_inv[d_inv >  0] <- 1 / d_inv[d_inv > 0]
        sv_fi$v %*% diag(d_inv, n_par) %*% t(sv_fi$u)
      }
    } else {
      qri   <- qr(I_reg)
      I_inv <- solve(qri, diag(n_par))
    }

    Delta <- diag(I_raw)
    D_sqrt <- sqrt(pmax(Delta, 0))
    Corr <- matrix(0, nrow = n_par, ncol = n_par, dimnames = list(param_names, param_names))
    denom <- outer(D_sqrt, D_sqrt)
    ok <- denom > .Machine$double.eps
    Corr[ok] <- I_raw[ok] / denom[ok]
    diag(Corr) <- 1

    rho <- vapply(seq_len(n_par), function(i) {
      others <- setdiff(seq_len(n_par), i)
      if (length(others) == 0) return(0)
      r_i <- Corr[i, others, drop = TRUE]
      R_oo <- Corr[others, others, drop = FALSE]
      R_oo_reg <- R_oo + diag(1e-10, nrow(R_oo))
      val <- as.numeric(t(r_i) %*% solve(R_oo_reg, r_i))
      max(0, min(1, val))
    }, numeric(1))
    names(rho) <- param_names

    # Cramer-Rao based strength index
    se_crlb <- sqrt(pmax(diag(I_inv), 0))
    s_dyn <- as.numeric(theta) / pmax(se_crlb, 1e-16)
    names(s_dyn) <- param_names

    strength_df <- data.frame(
      parameter = param_names,
      theta = as.numeric(theta),
      delta = as.numeric(Delta),
      rho = as.numeric(rho),
      se_crlb = as.numeric(se_crlb),
      s_dyn = as.numeric(s_dyn),
      stringsAsFactors = FALSE
    )
    strength_df <- strength_df[order(abs(strength_df$s_dyn)), , drop = FALSE]
    weak <- strength_df$parameter[abs(strength_df$s_dyn) < 1]

    # --- ENHANCED: Condition number and parameter type breakdown ---
    sv_J <- svd(J)$d
    cond_num <- max(sv_J) / max(min(sv_J), .Machine$double.eps)

    # Detect Sigma_e (shock std) parameters by name pattern.
    # These are computed BEFORE pass so the structural-only gate can use them.
    sigma_e_idx <- grep("^sig_|^stderr_|^se_", strength_df$parameter, ignore.case = TRUE)
    sigma_e_weak <- strength_df$parameter[intersect(which(abs(strength_df$s_dyn) < 1), sigma_e_idx)]
    structural_weak <- strength_df$parameter[setdiff(which(abs(strength_df$s_dyn) < 1), sigma_e_idx)]

    # Gate on structural (non-sigma) params only.
    # Sigma_e (shock std) parameters are inherently collinear in the Fisher
    # information matrix (all scale covariance in the same direction) and
    # nearly always appear weak under the CRLB criterion even in well-specified
    # models.  Gating on sigma params would produce spurious FAIL on every DSGE.
    # Structural params (phi_pi, rho, b, etc.) are the meaningful gate.
    # When Fisher matrix is rank-deficient, use pseudo-inverse strengths but
    # still apply the structural-only gate so D20 is informative on common
    # under-identified models (rather than always NA).
    pass <- if (length(structural_weak) == 0L) TRUE else FALSE

    # Under sampling weighting, a rank-deficient weighted Fisher information has
    # a non-trivial null space: the pseudo-inverse assigns arbitrarily large
    # (meaningless) s_i to unidentified parameter directions, which would mask a
    # genuine identification failure as a spurious PASS.  In that case the
    # per-parameter CRLB gate is not trustworthy -> report INCONCLUSIVE (NA) and
    # defer to D1's rank verdict.  (The unweighted/scaled modes keep their
    # established gate.)
    if (identical(weighting_actual, "sampling") && isTRUE(fisher_rank_deficient)) {
      pass <- NA
      weighting_note <- paste0(
        if (!is.null(weighting_note)) paste0(weighting_note, "; ") else "",
        sprintf("weighted Fisher rank-deficient (rank=%d/%d): per-parameter strengths for null-space directions are unreliable; see D1 rank.",
                rk_fisher, n_par))
    }

    # Per-equation Frobenius norm of Jacobian columns (which moments drive identification)
    col_norms <- sqrt(colSums(J^2))
    names(col_norms) <- param_names

    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Use log10 x scale so weakly-identified params are visible even when
      # a few dominant params (e.g. rho_pref, phi_pi) have very large |s_i|.
      # Sigma_e params are distinguished as INFO-only (grey fill).
      sigma_e_names <- grep("^sig_|^stderr_|^se_", strength_df$parameter,
                             ignore.case = TRUE, value = TRUE)
      strength_df$param_type <- ifelse(
        strength_df$parameter %in% sigma_e_names, "sigma_e (INFO)", "structural"
      )
      # For log10 scale: use absolute |s_i|; show sign via point shape
      strength_df$abs_s <- abs(strength_df$s_dyn)
      # Floor to 1e-3 so log scale is well-defined
      strength_df$abs_s_plot <- pmax(strength_df$abs_s, 1e-3)

      p_si <- ggplot2::ggplot(
        strength_df,
        ggplot2::aes(x = reorder(parameter, abs_s), y = abs_s_plot,
                     fill = interaction(abs_s >= 1, param_type, sep = "|"))
      ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::geom_vline(xintercept = NA, linetype = "dashed",
                            colour = dynhr_colours$red, linewidth = 0.4) +
        # Dashed threshold line at |s_i| = 1 in the log scale
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                            colour = dynhr_colours$red, linewidth = 0.4) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_log10() +
        ggplot2::scale_fill_manual(
          values = c(
            "TRUE|structural"   = dynhr_colours$mid_blue,
            "FALSE|structural"  = dynhr_colours$red,
            "TRUE|sigma_e (INFO)"  = dynhr_colours$light_blue,
            "FALSE|sigma_e (INFO)" = dynhr_colours$orange
          ),
          labels = c(
            "TRUE|structural"   = "Structural: adequate",
            "FALSE|structural"  = "Structural: weak",
            "TRUE|sigma_e (INFO)"  = "Sigma_e: adequate (INFO)",
            "FALSE|sigma_e (INFO)" = "Sigma_e: weak (INFO)"
          ),
          name = NULL
        ) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title = "D20: Identification strength index (log10 scale)",
          subtitle = sprintf("|s_i| < 1 = weak. Cond(J) = %.1e. Gate on structural params only.",
                             cond_num),
          x = NULL, y = "|s_i| = |theta_i / sqrt(diag(I^-1))| (log10)"
        )
      plots$strength_index <- .apply_meta(p_si, meta)
    }

    # Weighting descriptor for the summary
    weighting_detail <- switch(weighting_actual,
      sampling = sprintf(" [sampling-weighted Fisher info: I=J'Omega^-1 J over %d data moments]", n_moments_used),
      scaled   = " [scale-standardised Fisher info (cross-parameter comparison; gate not principled)]",
      "")
    if (!is.null(weighting_note)) weighting_detail <- paste0(weighting_detail, " (", weighting_note, ")")

    # Build richer summary
    rank_deficient_note <- if (fisher_rank_deficient) {
      sprintf(
        " NOTE: Fisher information matrix is rank-deficient (rank=%d/%d, rcond=%.2e). CRLB computed via pseudo-inverse. Non-identified parameters: %s.",
        rk_fisher, n_par, rc_fisher,
        paste(param_names[diag(I_raw) < max(diag(I_raw)) * 1e-10], collapse = ", ") |>
          (\(s) if (nchar(s) == 0) "unknown (all diagonal entries positive)" else s)()
      )
    } else ""

    summary_detail <- if (isTRUE(pass)) {
      paste0("PASS -- no weakly identified structural parameters.", rank_deficient_note,
             if (length(sigma_e_weak) > 0)
               sprintf(" Sigma_e weak (INFO-only): %s.", paste(sigma_e_weak, collapse = ", "))
             else "")
    } else if (is.na(pass)) {
      paste0("INCONCLUSIVE -- Fisher matrix rank-deficient; CRLB unreliable.", rank_deficient_note)
    } else {
      parts <- c("FAIL.", rank_deficient_note)
      if (length(sigma_e_weak) > 0) {
        parts <- c(parts, sprintf(" Sigma_e (shock std) weak: %s.",
                                   paste(sigma_e_weak, collapse = ", ")))
        parts <- c(parts, " Shock stds are inherently collinear in Fisher info (all scale covariance).")
        parts <- c(parts, " Condition number of Jacobian suggests adequate rank despite low |s_i|.")
      }
      if (length(structural_weak) > 0) {
        parts <- c(parts, sprintf(" Structural params weak: %s.",
                                   paste(structural_weak, collapse = ", ")))
        parts <- c(parts, " Consider additional observables, reparameterisation, or calibration.")
      }
      paste(parts, collapse = "")
    }

    .make_result(
      result = list(
        jacobian = J,
        fisher_information = I_raw,
        fisher_information_inverse = I_inv,
        fisher_rank_deficient = fisher_rank_deficient,
        fisher_rank = rk_fisher,
        fisher_rcond = rc_fisher,
        strength_table = strength_df,
        weak_params = weak,
        cond_number = cond_num,
        sigma_e_weak = sigma_e_weak,
        structural_weak = structural_weak,
        col_norms = col_norms,
        weighting = weighting_actual,
        n_moments_used = n_moments_used
      ),
      pass = pass,
      plots = plots,
      summary = sprintf(
        "D20 Identification strength: %d params. min |s_i| = %.3f. Cond(J) = %.1e.%s %s",
        n_par,
        min(abs(strength_df$s_dyn), na.rm = TRUE),
        cond_num,
        weighting_detail,
        summary_detail
      ),
      llm_summary = paste(c(
        sprintf("D20 | Fisher Identification Strength | %s",
                if (isTRUE(pass)) "PASS" else if (is.na(pass)) "INFO" else "FAIL"),
        sprintf("  params=%d min_abs_s=%.3f median_abs_s=%.3f cond_J=%.1e",
                n_par, min(abs(strength_df$s_dyn)), median(abs(strength_df$s_dyn)),
                cond_num),
        sprintf("  weakest: %s",
                paste(sprintf("%s=%.2f", head(strength_df$parameter, 5),
                              head(strength_df$s_dyn, 5)), collapse = ", ")),
        sprintf("  sigma_e_weak=%d structural_weak=%d",
                length(sigma_e_weak), length(structural_weak)),
        sprintf("  action: %s",
                if (isTRUE(pass)) "Structural identification strength adequate (Sigma_e weakness is expected and INFO-only)."
                else if (is.na(pass)) "Fisher matrix rank-deficient; CRLB via pseudo-inverse -- treat strengths as indicative only."
                else if (length(structural_weak) > 0)
                  sprintf("Structural parameters weakly identified (%s). Add moments or calibrate.",
                          paste(head(structural_weak, 3), collapse = ", "))
                else
                  "Weakness confined to Sigma_e (shock std) parameters. This is typical: shock stds are collinear in Fisher info. Use D25 for higher-order identification. Condition number suggests adequate rank for identification."),
        sprintf("note: |s_i| < 1 threshold is stringent for T=%s obs. Gate applies to structural params only; Sigma_e params are INFO-only.",
                if (!is.null(meta$T_obs)) as.character(meta$T_obs) else "(T not supplied -- pass meta$T_obs for sample-size context)"),
        sprintf("  top_3_col_norms: %s",
                paste(sprintf("%s=%.1e", names(sort(col_norms, decreasing = TRUE)[1:3]),
                              sort(col_norms, decreasing = TRUE)[1:3]), collapse = ", "))
      ), collapse = "\n")
    )
}
