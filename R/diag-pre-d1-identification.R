## R/diag-pre-d1-identification.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D1 local identification (Iskrev 2010)
## --------------------------------------------------------------------------

#' D1. Local identification (Iskrev 2010)
#'
#' Computes the numerical Jacobian of the vector of model-implied moments
#' with respect to the deep parameter vector and checks its column rank via
#' SVD: full column rank means all parameters are locally identified from
#' those moments (Iskrev 2010).
#'
#' The rank is taken on an \emph{equilibrated} Jacobian: each row is divided
#' by its largest absolute entry (Dynare's \code{normalize_jacobians}) and each
#' column is then scaled to unit Euclidean norm. Both scalings are
#' nonsingular, so the exact rank is unchanged, but the singular values become
#' free of moment and parameter units. All-zero rows are dropped; an
#' all-zero column (a parameter that moves no moment) is kept and shows up as
#' an exactly zero singular value.
#'
#' \strong{Rank tolerance.} The Jacobian is a central finite difference, so its
#' entries carry truncation and round-off error far above machine precision
#' (about 5e-7 relative on a small RBC model). A machine-epsilon rank
#' tolerance therefore reports an exactly non-identified model (e.g. two
#' parameters that enter only through their product) as full rank. With
#' \code{tol_rank = NULL} (default) D1 estimates that error by re-computing the
#' Jacobian at step \code{2 * eps} and uses
#' \eqn{tol = \max(\max(\dim J)\,\sigma_1\,\epsilon_{mach},\ 10\,\|J_e(h) - J_e(2h)\|_2)};
#' a singular value below it is indistinguishable from zero at the
#' precision of the derivatives. A numeric \code{tol_rank} is used instead as a
#' tolerance relative to the largest singular value.
#'
#' A parameter is reported in \code{unidentified_params} when it loads on a
#' null-space direction (singular value <= tol) with at least half the
#' largest loading of that direction, and in \code{weak_params} when it does
#' so for an unidentified \emph{or} a weak direction (singular value below
#' \code{weak_rel} times the largest).
#'
#' @param model_solve_fn  Function: theta -> named numeric vector of moments
#'                        (e.g. c(sd_y, sd_pi, acf1_y, ...) ). May be NULL
#'                        when \code{jacobian} is supplied (the rank tolerance
#'                        then falls back to machine precision).
#' @param theta           Numeric vector -- parameter values at calibration point
#' @param param_names     Character vector -- names of parameters (default
#'                        \code{names(theta)})
#' @param moment_names    Character vector -- names of moments (rows of Jacobian)
#' @param eps             Step size for finite differences (default 1e-5); must
#'                        match the step of a supplied \code{jacobian}.
#' @param jacobian        Optional pre-computed Jacobian (n_moment x n_par);
#'   when supplied (e.g. by the orchestrator, shared with D20) the
#'   finite-difference re-computation is skipped.
#' @param tol_rank        NULL (adaptive, see Details) or a relative rank
#'   tolerance on the equilibrated Jacobian.
#' @param weak_rel        Relative singular-value threshold for "weak"
#'   directions (default 1e-3).
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
                                    meta = NULL,
                                    tol_rank = NULL,
                                    weak_rel = 1e-3) {

    if (is.null(param_names)) param_names <- names(theta)
    if (is.null(param_names)) param_names <- paste0("theta_", seq_along(theta))
    if (is.null(moment_names) && is.null(colnames(jacobian)) && !is.null(model_solve_fn)) {
      f0 <- model_solve_fn(theta)
      moment_names <- if (!is.null(names(f0))) names(f0) else paste0("m_", seq_along(f0))
    }

    # Compute Jacobian (or reuse one supplied by the orchestrator: D1 and D20
    # share the same J, so the orchestrator computes it once and passes it in).
    J <- jacobian %||% .numerical_jacobian(model_solve_fn, theta, eps = eps)
    J <- as.matrix(J)
    if (is.null(moment_names)) moment_names <- rownames(J) %||% paste0("m_", seq_len(nrow(J)))
    # Guard: ensure param_names length matches Jacobian column count
    if (ncol(J) != length(param_names)) {
      .dynhr_warn(sprintf(
        "d1: param_names length (%d) != Jacobian columns (%d). Using generic param labels.",
        length(param_names), ncol(J)))
      param_names <- paste0("theta_", seq_len(ncol(J)))
    }
    colnames(J) <- param_names
    # Guard: ensure moment_names length matches Jacobian row count
    if (nrow(J) != length(moment_names)) {
      .dynhr_warn(sprintf(
        "d1: moment_names length (%d) != Jacobian rows (%d). Using generic moment labels.",
        length(moment_names), nrow(J)))
      moment_names <- paste0("m_", seq_len(nrow(J)))
    }
    rownames(J) <- moment_names
    n_par <- ncol(J)

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
                       unidentified_params = character(0),
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

    J2 <- if (is.null(tol_rank) && !is.null(model_solve_fn))
      .numerical_jacobian(model_solve_fn, theta, eps = 2 * eps)
    rk <- .ident_equilibrated_rank(J, J2, tol_rank = tol_rank, weak_rel = weak_rel)
    Je <- rk$Je; row_scale <- rk$row_scale; sv <- rk$svd; V <- rk$V
    singular_values <- rk$singular_values; sv_max <- rk$sv_max
    tol <- rk$tol; fd_noise <- rk$fd_noise; tol_source <- rk$tol_source
    weak_threshold <- rk$weak_threshold; rank_J <- rk$rank
    sv_class <- rk$sv_class; unid_dirs <- rk$unid_dirs
    unidentified_params <- rk$unidentified_params
    weak_only <- rk$weak_only
    full_rank <- (rank_J == n_par)
    weak_params <- c(unidentified_params, weak_only)

    pass <- full_rank

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    # (a) Heatmap of the row-normalised Jacobian (each row scaled to max |.|=1,
    # so moments on different scales are all visible)
    Jrow <- J / ifelse(row_scale > 0, row_scale, 1)
    J_long <- data.frame(
      Moment    = factor(rep(moment_names, times = n_par), levels = rev(moment_names)),
      Parameter = factor(rep(param_names, each = nrow(J)), levels = param_names),
      Value     = as.vector(Jrow)
    )
    p_jac <- ggplot2::ggplot(
      J_long, ggplot2::aes(x = Parameter, y = Moment, fill = Value)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
      ggplot2::scale_fill_gradient2(low = dynhr_colours$red, mid = dynhr_colours$white,
                                    high = dynhr_colours$mid_blue, midpoint = 0,
                                    limits = c(-1, 1),
                                    name = "d moment / d param\n(row-normalised)") +
      theme_dynhr_diagnostic() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = ggplot2::rel(0.75)),
        axis.text.y = ggplot2::element_text(size = ggplot2::rel(0.7))
      ) +
      ggplot2::labs(
        title = "D1: Jacobian of moments w.r.t. parameters",
        subtitle = sprintf("Each row divided by its largest absolute entry. Rank = %d / %d parameters",
                           rank_J, n_par),
        x = "Parameter", y = "Moment"
      )
    plots$jacobian_heatmap <- .apply_meta(p_jac, meta)

    # (b) Singular values (lollipop on a log axis: bars on a log axis would
    # hang from 1 and misrepresent values < 1). Exact zeros sit at the floor.
    pos <- singular_values[singular_values > 0]
    cand <- c(pos, tol, weak_threshold)
    cand <- cand[cand > 0]
    floor_val <- 10^(floor(log10(if (length(cand)) min(cand) else 1e-16)) - 1)
    lead <- vapply(seq_len(n_par), function(k) {
      v <- abs(V[, k]); if (max(v) <= 0) "" else param_names[which.max(v)]
    }, character(1))
    sv_df <- data.frame(
      index = factor(seq_len(n_par),
                     labels = sprintf("%d\n%s", seq_len(n_par), lead)),
      value = pmax(singular_values, floor_val),
      class = factor(sv_class, levels = c("Identified", "Weak", "Unidentified")),
      zero  = singular_values <= 0
    )
    ref_df <- data.frame(
      y = c(tol, weak_threshold),
      what = c(sprintf("rank tolerance = %.1e x largest (%s)", tol / max(sv_max, 1e-300), tol_source),
               sprintf("weak threshold = %.0e x largest", weak_rel))
    )
    p_sv <- ggplot2::ggplot(sv_df, ggplot2::aes(x = index, y = value, colour = class)) +
      ggplot2::geom_segment(ggplot2::aes(xend = index, y = floor_val, yend = value),
                            linewidth = 0.8) +
      ggplot2::geom_point(ggplot2::aes(shape = zero), size = 2.6) +
      ggplot2::geom_hline(data = ref_df,
                          ggplot2::aes(yintercept = y, linetype = what),
                          colour = dynhr_colours$grey, linewidth = 0.5) +
      ggplot2::scale_colour_manual(values = c("Identified" = dynhr_colours$mid_blue,
                                              "Weak" = dynhr_colours$orange,
                                              "Unidentified" = dynhr_colours$red),
                                   name = NULL) +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                                  breaks = "TRUE",
                                  labels = c(`TRUE` = "exactly 0 (drawn at floor)"),
                                  name = NULL) +
      ggplot2::scale_linetype_manual(values = c("dashed", "dotted"), name = NULL) +
      ggplot2::scale_y_log10() +
      theme_dynhr_diagnostic() +
      ggplot2::theme(legend.box = "vertical") +
      ggplot2::labs(
        title = "D1: Singular values of the equilibrated identification Jacobian",
        subtitle = sprintf("Rank = %d / %d (count above the rank tolerance). Axis label: parameter with the largest loading",
                           rank_J, n_par),
        x = "Singular value (index / dominant parameter)",
        y = "Singular value (log scale)"
      )
    plots$singular_values <- .apply_meta(p_sv, meta)

    }  # end requireNamespace guard

    tol_txt <- sprintf("rank tolerance %.1e relative (%s%s)",
                       tol / max(sv_max, 1e-300), tol_source,
                       if (is.finite(fd_noise)) sprintf(", FD error %.1e", fd_noise) else "")
    summary_text <- paste0(
      sprintf("D1 Local identification: Jacobian rank = %d / %d (%d moments; %s). ",
              rank_J, n_par, nrow(J), tol_txt),
      if (full_rank) "All parameters locally identified."
      else sprintf("RANK DEFICIENT -- %d parameter direction(s) unidentified%s.",
                   n_par - rank_J,
                   if (length(unidentified_params))
                     sprintf(" (involving: %s)", paste(unidentified_params, collapse = ", "))
                   else ""),
      if (length(weak_only))
        sprintf(" Weakly identified: %s.", paste(weak_only, collapse = ", "))
      else ""
    )

    .make_result(
      result  = list(jacobian = J, jacobian_equilibrated = Je,
                     singular_values = singular_values,
                     rank = rank_J, rank_tolerance = tol / max(sv_max, 1e-300),
                     tolerance_source = tol_source, fd_error = fd_noise,
                     weak_threshold = weak_rel,
                     sv_class = stats::setNames(sv_class, names(singular_values)),
                     null_space = V[, unid_dirs, drop = FALSE],
                     unidentified_params = unidentified_params,
                     weak_params = weak_params,
                     svd = sv),
      pass    = pass,
      plots   = plots,
      summary = summary_text,
      llm_summary = {
        badge   <- if (pass) "PASS" else "FAIL"
        sv_str  <- paste(sprintf("%.3e", head(sort(singular_values), 5)), collapse = ", ")
        paste(c(
          sprintf("D1 | Local Identification | %s", badge),
          sprintf("  params=%d moments=%d jacobian_rank=%d (expected=%d) %s",
                  n_par, nrow(J), rank_J, n_par, tol_txt),
          sprintf("  smallest_sv (equilibrated): %s", sv_str),
          if (length(unidentified_params) > 0)
            sprintf("  unidentified_params: %s", paste(unidentified_params, collapse = ", ")),
          if (length(weak_only) > 0)
            sprintf("  weak_params: %s", paste(weak_only, collapse = ", ")),
          sprintf("  action: %s",
                  if (pass)
                    "Model is locally identified at calibration point."
                  else
                    sprintf("Rank deficient (rank=%d, expected=%d). %s may not be identified. Fix: add moments, calibrate one parameter, or check for collinearity.",
                            rank_J, n_par,
                            paste(head(unidentified_params, 3), collapse = ", ")))
        ), collapse = "\n")
      }
    )
}




#' Rank analysis of an equilibrated identification Jacobian (D1 / D20)
#'
#' Row max-abs then column 2-norm equilibration (rank-preserving), full SVD,
#' and a rank tolerance that is aware of finite-difference error: with
#' \code{J2} (the same Jacobian at step 2h) the tolerance is
#' \code{max(machine, 10 * ||Je(h) - Je(2h)||_2)}, and rows no larger than
#' 10x their own FD error are dropped as numerically zero; with a numeric
#' \code{tol_rank} it is relative to the largest singular value.
#' @return list: Je, row_scale, col_scale, svd, V, singular_values, sv_max,
#'   tol, fd_noise, tol_source, weak_threshold, rank, sv_class, unid_dirs,
#'   null_space, unidentified_params, weak_only.
#' @noRd
.ident_equilibrated_rank <- function(J, J2 = NULL, tol_rank = NULL, weak_rel = 1e-3) {
  param_names <- colnames(J) %||% paste0("theta_", seq_len(ncol(J)))
  n_par <- ncol(J)
  row_scale <- apply(abs(J), 1, max)
  use_J2 <- is.null(tol_rank) && !is.null(J2) && identical(dim(J2), dim(J)) &&
            all(is.finite(J2))
  # A row whose entries are all within 10x of their own finite-difference
  # error is numerically zero (e.g. a whitened moment combination the
  # parameters do not move); max-abs scaling would blow that noise up to 1.
  row_noise <- if (use_J2) apply(abs(J - as.matrix(J2)), 1, max) else 0
  keep_rows <- row_scale > 10 * row_noise & row_scale > 0
  Je <- J[keep_rows, , drop = FALSE] / row_scale[keep_rows]
  col_scale <- sqrt(colSums(Je^2))
  zero_cols <- param_names[col_scale <= .Machine$double.eps]
  col_scale[col_scale <= .Machine$double.eps] <- 1
  Je <- sweep(Je, 2, col_scale, "/")

  # Full SVD (nv = n_par so the null space is complete even when there are
  # fewer (non-zero) moments than parameters).
  if (nrow(Je) > 0L) {
    sv <- svd(Je, nu = 0, nv = n_par)
  } else {
    sv <- list(d = numeric(0), v = diag(n_par))
  }
  singular_values <- c(sv$d, rep(0, n_par - length(sv$d)))[seq_len(n_par)]
  names(singular_values) <- paste0("sv_", seq_len(n_par))
  V <- sv$v
  dimnames(V) <- list(param_names, paste0("sv_", seq_len(n_par)))
  sv_max <- max(c(singular_values, 0))

  tol_machine <- max(dim(Je), 1L) * sv_max * .Machine$double.eps
  fd_noise <- NA_real_
  tol_source <- "machine"
  if (!is.null(tol_rank)) {
    tol <- max(tol_machine, tol_rank * sv_max)
    tol_source <- "user"
  } else {
    if (use_J2 && sv_max > 0) {
      J2 <- as.matrix(J2)
      J2e <- sweep(J2[keep_rows, , drop = FALSE] / row_scale[keep_rows], 2, col_scale, "/")
      fd_noise <- max(svd(Je - J2e, nu = 0, nv = 0)$d)
      tol_source <- "finite-difference"
    }
    tol <- max(tol_machine, 10 * fd_noise, na.rm = TRUE)
  }
  weak_threshold <- max(weak_rel * sv_max, tol)
  sv_class <- ifelse(singular_values <= tol, "Unidentified",
              ifelse(singular_values < weak_threshold, "Weak", "Identified"))

  # Dominant parameters of each flagged right singular vector.
  loaders <- function(k) {
    v <- abs(V[, k])
    if (max(v) <= 0) return(character(0))
    param_names[v >= 0.5 * max(v)]
  }
  unid_dirs <- which(sv_class == "Unidentified")
  unidentified_params <- unique(c(zero_cols,
                                  as.character(unlist(lapply(unid_dirs, loaders)))))
  weak_only <- setdiff(unique(as.character(unlist(lapply(which(sv_class == "Weak"), loaders)))),
                       unidentified_params)
  list(Je = Je, row_scale = row_scale, col_scale = col_scale, svd = sv, V = V,
       singular_values = singular_values, sv_max = sv_max, tol = tol,
       fd_noise = fd_noise, tol_source = tol_source,
       weak_threshold = weak_threshold, rank = sum(singular_values > tol),
       sv_class = sv_class, unid_dirs = unid_dirs,
       null_space = V[, unid_dirs, drop = FALSE],
       unidentified_params = unidentified_params, weak_only = weak_only)
}


## NOTE (2026-05-31): the `d19_local_rank_identification` stub was removed. It
## was an admitted re-label of D1 (local identification via moment-Jacobian
## rank) that only ever returned an "identical to D1 ... Skipped" placeholder,
## so two diagnostics reported the same thing. D1 above IS the local rank
## identification check; callers should use it directly.


#' D20. Identification strength via Fisher information
#'
#' Computes the Fisher information of a moment vector, \eqn{I = J'\,W\,J},
#' from the Jacobian J of model-implied moments w.r.t. parameters, and reports
#' per parameter
#' - the sensitivity component \eqn{\Delta_i = I_{ii}},
#' - the collinearity component \eqn{\rho_i} (squared multiple correlation of
#'   parameter i on the others, from \eqn{I = D R D}),
#' - the strength \eqn{s_i = \theta_i / SE_i} with
#'   \eqn{SE_i = \sqrt{(I^{-1})_{ii}}}.
#'
#' \strong{Weighting.} \code{"none"} (the default) uses \eqn{W = I}, the
#' unweighted \eqn{J'J} of Iskrev (2010, J. Monetary Econ. 57(2):189-202),
#' whose own a-priori identification-strength measure is the diagonal of
#' \eqn{(\partial m/\partial\theta)'(\partial m/\partial\theta)} with columns
#' normalised for cross-parameter comparability -- no efficient/sampling
#' weighting. Ratto & Iskrev (2011, "Identification analysis of DSGE models
#' with Dynare") and Dynare's \code{identification} command follow the same
#' convention (identity or model-implied likelihood metric; Dynare's
#' \code{normalize_jacobians} rescaling is for numerical conditioning only).
#' Under \eqn{W = I} the strength \eqn{s_i} depends on the units of the
#' moments and is \emph{not} a t-ratio.
#'
#' \code{"scaled"} uses \eqn{W = diag(f_0^{-2})} (J row-scaled by
#' \eqn{1/|f_0|}), which makes the cross-parameter \eqn{\Delta_i} ranking
#' scale-free but makes an \eqn{s_i < 1} gate degenerate (a parameter mapped
#' linearly to a dedicated moment gets \eqn{s_i \approx 1}).
#'
#' \code{"sampling"} is the GMM / efficient-minimum-distance form
#' \eqn{W = Var(\hat m)^{-1}}: \eqn{SE_i} is then the asymptotic standard
#' error of an efficient minimum-distance estimator and \eqn{s_i} is a genuine
#' t-ratio. It is \strong{not} the default and no reference implementation
#' (Dynare, MacroModelling.jl, RISE) offers it: a sample-estimated
#' \eqn{\Omega} is itself noisy at macro sample sizes, so making an
#' identification verdict depend on it would make the diagnostic
#' sampling-noise-dependent (the same finite-sample concern that keeps D29's
#' Stock-Wright S test informational).
#'
#' \strong{Rank.} Rank deficiency is decided exactly as in D1 (shared helper):
#' on the row/column-equilibrated weighted Jacobian, with a tolerance that
#' accounts for the finite-difference error of J (estimated from a second
#' Jacobian at step \code{2 * eps}). A parameter that loads on a null-space
#' direction is not identified: its SE is \code{Inf} and \eqn{s_i = 0}. For
#' the remaining parameters SE comes from the inverse of I on its identified
#' subspace. A rank-deficient Fisher information always FAILs.
#'
#' \strong{What gates the badge.} The \emph{rank} test always gates. The
#' strength table (\eqn{s_i}, the weak list) gates the badge \emph{only} under
#' \code{weighting = "sampling"}, where \eqn{s_i} is an actual t-ratio and
#' \eqn{|s_i| \ge 1} is a meaningful bar. Under \code{"none"} and
#' \code{"scaled"} the strengths are \strong{informational}: they carry the
#' units of the moments, so a "weak" flag there is a ranking hint, not a
#' verdict, and a full-rank model PASSes with an INFO note listing the weak
#' parameters. Shock standard deviations (names starting \code{sig_},
#' \code{sigma_}, \code{stderr_}, \code{se_}, \code{std_}) are reported as
#' INFO only under every weighting. The \eqn{|s_i| \ge 1} bar itself is a
#' package choice with no literature source; Iskrev and Dynare report the
#' strength bars without a numeric cutoff.
#'
#' @param model_solve_fn Function: theta -> named numeric vector of moments.
#' @param theta Numeric vector of parameter values.
#' @param param_names Optional character vector of parameter names.
#' @param moment_names Optional character vector of moment names.
#' @param eps Step size for finite differences (must match a supplied
#'   \code{jacobian}).
#' @param weighting \code{"none"} (default), \code{"scaled"} or
#'   \code{"sampling"}; see Details. \code{"sampling"} weights only the moments
#'   present in \code{moment_cov} (aligned by name); if it is missing or does
#'   not overlap the Jacobian rows, D20 falls back to \code{"none"} and notes
#'   this in the summary.
#' @param moment_cov Optional named covariance matrix of the moment
#'   \emph{estimator}, \eqn{Var(\hat m)} (the long-run moment covariance
#'   already divided by the sample size T; see
#'   \code{.d20_moment_sampling_cov}). Used only when
#'   \code{weighting = "sampling"}.
#' @param jacobian Optional pre-computed Jacobian (n_moment x n_par). When
#'   supplied (e.g. by the orchestrator, shared with D1) the step-h
#'   finite-difference re-computation is skipped.
#' @param tol_rank NULL (finite-difference-aware tolerance, as D1) or a rank
#'   tolerance relative to the largest equilibrated singular value.
#' @return dynhr_diagnostic list
#' @noRd
d20_fisher_identification_strength <- function(model_solve_fn,
                                               theta,
                                               param_names  = NULL,
                                               moment_names = NULL,
                                               eps   = 1e-5,
                                               weighting  = c("none", "scaled", "sampling"),
                                               moment_cov = NULL,
                                               jacobian   = NULL,
                                               tol_rank   = NULL,
                                               meta  = NULL) {

    weighting <- match.arg(weighting)
    n_par <- length(theta)
    if (is.null(param_names)) param_names <- names(theta) %||% paste0("theta_", seq_len(n_par))

    # Reuse the orchestrator-supplied Jacobian when available (shared with D1).
    J  <- as.matrix(jacobian %||% .numerical_jacobian(model_solve_fn, theta, eps = eps))
    f0 <- model_solve_fn(theta)
    if (is.null(moment_names)) {
      moment_names <- names(f0) %||% paste0("m_", seq_len(length(f0)))
    }
    if (nrow(J) != length(moment_names)) {
      .dynhr_warn(sprintf(
        "d20: moment_names length (%d) != Jacobian rows (%d). Using generic labels.",
        length(moment_names), nrow(J)))
      moment_names <- paste0("m_", seq_len(nrow(J)))
    }
    if (ncol(J) != length(param_names)) {
      .dynhr_warn(sprintf(
        "d20: param_names length (%d) != Jacobian columns (%d). Using generic param labels.",
        length(param_names), ncol(J)))
      param_names <- paste0("theta_", seq_len(ncol(J)))
    }
    n_par <- ncol(J)
    dimnames(J) <- list(moment_names, param_names)

    # Guard: non-finite Jacobian (NaN/Inf from solve_lyapunov degeneracy at
    # near-unit-root points) is a numerical failure, not a finding.
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

    # --- Weighting: Jw such that I = crossprod(Jw) ---
    # sampling: Jw = Omega^{-1/2} J over the moments that have a sampling
    #   covariance (aligned by name), so I = J' Var(m_hat)^{-1} J.
    # scaled:   Jw = diag(1/|f0|) J.
    f0v <- as.numeric(f0)
    f0v[!is.finite(f0v)] <- 0
    weighting_note   <- NULL
    weighting_actual <- weighting
    rows_used <- moment_names
    Wh <- NULL
    m_scale <- NULL
    if (identical(weighting, "sampling")) {
      ok_cov <- is.matrix(moment_cov) && !is.null(rownames(moment_cov)) &&
                all(is.finite(moment_cov))
      common <- if (ok_cov) intersect(moment_names, rownames(moment_cov)) else character(0)
      if (length(common) >= 1L) {
        rows_used <- common
        Wh <- .robust_Omega_inv_sqrt(moment_cov[common, common, drop = FALSE])
      } else {
        weighting_actual <- "none"
        weighting_note   <- "sampling weighting requested but moment_cov was missing/unaligned; fell back to unweighted J'J"
      }
    } else if (identical(weighting, "scaled") && length(f0v) == nrow(J)) {
      m_scale <- pmax(abs(f0v), 1e-8 * max(abs(f0v)), 1e-300)
    }
    weigh <- function(M) {
      M <- as.matrix(M)
      if (!is.null(m_scale)) M <- M / m_scale      # row i divided by m_scale[i]
      dimnames(M) <- list(moment_names, param_names)
      M <- M[rows_used, , drop = FALSE]
      if (!is.null(Wh)) M <- Wh %*% M
      M
    }
    Jw <- weigh(J)
    n_moments_used <- nrow(Jw)

    # --- Rank: same equilibration + FD-noise-aware tolerance as D1 ---
    J2 <- if (is.null(tol_rank)) .numerical_jacobian(model_solve_fn, theta, eps = 2 * eps)
    J2w <- if (!is.null(J2) && identical(dim(J2), dim(J)) && all(is.finite(J2))) weigh(J2)
    rk <- .ident_equilibrated_rank(Jw, J2w, tol_rank = tol_rank)
    fisher_rank <- rk$rank
    fisher_rank_deficient <- fisher_rank < n_par
    unidentified <- union(rk$unidentified_params,
                          param_names[rowSums(rk$null_space^2) >= 0.01])

    # --- Fisher information and its inverse on the identified subspace ---
    # Column scaling keeps the SVD well conditioned; row weights are part of I.
    I_raw <- crossprod(Jw)
    cs <- sqrt(colSums(Jw^2))
    cs[cs <= 0] <- 1
    sj <- svd(sweep(Jw, 2, cs, "/"), nu = 0, nv = n_par)
    k  <- fisher_rank
    Vk <- sj$v[, seq_len(k), drop = FALSE]
    I_inv <- (Vk %*% diag(1 / sj$d[seq_len(k)]^2, k) %*% t(Vk)) / outer(cs, cs)
    dimnames(I_inv) <- list(param_names, param_names)

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

    # Cramer-Rao standard errors; unidentified parameters have SE = Inf.
    se_crlb <- sqrt(pmax(diag(I_inv), 0))
    se_crlb[param_names %in% unidentified] <- Inf
    s_dyn <- as.numeric(theta) / se_crlb
    s_dyn[!is.finite(s_dyn)] <- 0
    names(se_crlb) <- names(s_dyn) <- param_names

    strength_df <- data.frame(
      parameter = param_names,
      theta = as.numeric(theta),
      delta = as.numeric(Delta),
      rho = as.numeric(rho),
      se_crlb = as.numeric(se_crlb),
      s_dyn = as.numeric(s_dyn),
      identified = !(param_names %in% unidentified),
      stringsAsFactors = FALSE
    )
    strength_df <- strength_df[order(abs(strength_df$s_dyn)), , drop = FALSE]
    rownames(strength_df) <- NULL
    weak <- strength_df$parameter[abs(strength_df$s_dyn) < 1]

    # Condition number of the equilibrated (unit-free) weighted Jacobian.
    sv_e <- rk$singular_values
    cond_num <- if (min(sv_e) > 0) rk$sv_max / min(sv_e) else Inf

    # Shock standard deviations are INFO-only in the gate.
    sigma_pat <- "^(sig|sigma|stderr|se|std)(_|$)"
    is_sigma <- grepl(sigma_pat, strength_df$parameter, ignore.case = TRUE)
    is_weak  <- abs(strength_df$s_dyn) < 1
    sigma_e_weak    <- strength_df$parameter[is_weak & is_sigma]
    structural_weak <- strength_df$parameter[is_weak & !is_sigma]

    # The rank test always gates. The strength table gates only under
    # "sampling", where s_i is a genuine t-ratio; under "none"/"scaled" the
    # s_i carry the moments' units and are informational only (see Details).
    strength_gates <- identical(weighting_actual, "sampling")
    pass <- !fisher_rank_deficient &&
            (!strength_gates || length(structural_weak) == 0L)

    col_norms <- sqrt(colSums(Jw^2))
    names(col_norms) <- param_names

    weight_lab <- switch(weighting_actual,
      sampling = "I = J' Var(m)^-1 J",
      scaled   = "I = J' diag(1/m^2) J",
      "I = J'J (unit moment variances)")

    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      pdf_ <- strength_df
      pdf_$abs_s <- abs(pdf_$s_dyn)
      pos <- pdf_$abs_s[pdf_$abs_s > 0]
      floor_val <- 10^(floor(log10(min(c(pos, 1)))) - 1)
      pdf_$abs_s_plot <- pmax(pdf_$abs_s, floor_val)
      pdf_$status <- ifelse(!pdf_$identified, "Not identified (s = 0, at floor)",
                     paste0(ifelse(grepl(sigma_pat, pdf_$parameter, ignore.case = TRUE),
                                   "Shock s.d. (INFO): ", "Structural: "),
                            ifelse(pdf_$abs_s >= 1, "adequate", "weak")))
      pdf_$parameter <- factor(pdf_$parameter, levels = pdf_$parameter)
      p_si <- ggplot2::ggplot(
        pdf_, ggplot2::aes(x = parameter, y = abs_s_plot, colour = status)
      ) +
        ggplot2::geom_segment(ggplot2::aes(xend = parameter, y = floor_val, yend = abs_s_plot),
                              linewidth = 0.8) +
        ggplot2::geom_point(size = 2.6) +
        ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                            colour = dynhr_colours$grey, linewidth = 0.5) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_log10() +
        ggplot2::scale_colour_manual(
          values = c("Structural: adequate" = dynhr_colours$mid_blue,
                     "Structural: weak" = dynhr_colours$orange,
                     "Shock s.d. (INFO): adequate" = dynhr_colours$light_blue,
                     "Shock s.d. (INFO): weak" = dynhr_colours$grey,
                     "Not identified (s = 0, at floor)" = dynhr_colours$red),
          name = NULL) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title = "D20: Identification strength, abs(theta_i) / SE_i",
          subtitle = sprintf("%s; Fisher rank %d / %d. Dashed line at 1: weak below (%s).",
                             weight_lab, fisher_rank, n_par,
                             if (identical(weighting_actual, "sampling")) "gates the badge"
                             else "informational; the rank test gates the badge"),
          x = NULL,
          y = if (identical(weighting_actual, "sampling"))
                "Strength = absolute asymptotic t-ratio (log scale)"
              else "Strength under unit moment weights, not a t-ratio (log scale)"
        )
      plots$strength_index <- .apply_meta(p_si, meta)
    }

    weighting_detail <- switch(weighting_actual,
      sampling = sprintf(" [sampling-weighted Fisher info I=J'Var(m)^-1 J over %d data moments; s_i is an asymptotic t-ratio]", n_moments_used),
      scaled   = " [scale-standardised Fisher info (cross-parameter comparison; gate not principled)]",
      " [unweighted J'J: s_i assumes unit moment variances and is not a t-ratio]")
    if (!is.null(weighting_note)) weighting_detail <- paste0(weighting_detail, " (", weighting_note, ")")

    tol_txt <- sprintf("rank tolerance %.1e relative (%s%s)",
                       rk$tol / max(rk$sv_max, 1e-300), rk$tol_source,
                       if (is.finite(rk$fd_noise)) sprintf(", FD error %.1e", rk$fd_noise) else "")
    rank_deficient_note <- if (fisher_rank_deficient) {
      sprintf(" Fisher information is RANK DEFICIENT (rank %d/%d; %s). Not identified (SE = Inf, s = 0): %s.",
              fisher_rank, n_par, tol_txt,
              if (length(unidentified)) paste(unidentified, collapse = ", ") else "no single dominant parameter")
    } else ""

    # Under "none"/"scaled" a weak s_i is an INFO note, not part of the verdict.
    structural_weak_txt <- if (length(structural_weak) == 0L) "" else if (strength_gates) {
      sprintf(" Structural params weak (|s_i| < 1): %s. Consider additional observables, reparameterisation, or calibration.",
              paste(structural_weak, collapse = ", "))
    } else {
      sprintf(paste0(" INFO (does not affect the badge): structural params with |s_i| < 1 under ",
                     "weighting=\"%s\": %s. These s_i carry the units of the moments and are not ",
                     "t-ratios; re-run with weighting=\"sampling\" and a moment_cov for a gating ",
                     "strength measure."),
              weighting_actual, paste(structural_weak, collapse = ", "))
    }
    summary_detail <- if (pass) {
      paste0(if (length(structural_weak) == 0L)
               "PASS -- full Fisher rank and no weakly identified structural parameters."
             else "PASS -- full Fisher rank (the strength table is informational here).",
             structural_weak_txt,
             if (length(sigma_e_weak) > 0)
               sprintf(" Shock s.d. weak (INFO-only): %s.", paste(sigma_e_weak, collapse = ", "))
             else "")
    } else {
      paste0("FAIL.", rank_deficient_note, structural_weak_txt,
             if (length(sigma_e_weak) > 0)
               sprintf(" Shock s.d. weak (INFO-only): %s.", paste(sigma_e_weak, collapse = ", "))
             else "")
    }

    top_cn <- sort(col_norms, decreasing = TRUE)[seq_len(min(3L, n_par))]
    .make_result(
      result = list(
        jacobian = J,
        jacobian_weighted = Jw,
        fisher_information = I_raw,
        fisher_information_inverse = I_inv,
        fisher_rank_deficient = fisher_rank_deficient,
        fisher_rank = fisher_rank,
        rank_tolerance = rk$tol / max(rk$sv_max, 1e-300),
        tolerance_source = rk$tol_source,
        fd_error = rk$fd_noise,
        singular_values = sv_e,
        unidentified_params = unidentified,
        strength_table = strength_df,
        weak_params = weak,
        cond_number = cond_num,
        sigma_e_weak = sigma_e_weak,
        structural_weak = structural_weak,
        col_norms = col_norms,
        weighting = weighting_actual,
        strength_gates_badge = strength_gates,
        n_moments_used = n_moments_used
      ),
      pass = pass,
      plots = plots,
      summary = sprintf(
        "D20 Identification strength: %d params, Fisher rank %d. min |s_i| = %.3f. Cond(equilibrated J) = %.1e.%s %s",
        n_par, fisher_rank, min(abs(strength_df$s_dyn)), cond_num,
        weighting_detail, summary_detail
      ),
      llm_summary = paste(c(
        sprintf("D20 | Fisher Identification Strength | %s", if (pass) "PASS" else "FAIL"),
        sprintf("  params=%d fisher_rank=%d min_abs_s=%.3f median_abs_s=%.3f cond_Je=%.1e weighting=%s",
                n_par, fisher_rank, min(abs(strength_df$s_dyn)),
                stats::median(abs(strength_df$s_dyn)), cond_num, weighting_actual),
        sprintf("  %s", tol_txt),
        sprintf("  weakest: %s",
                paste(sprintf("%s=%.2f", head(strength_df$parameter, 5),
                              head(strength_df$s_dyn, 5)), collapse = ", ")),
        if (length(unidentified))
          sprintf("  unidentified_params: %s", paste(unidentified, collapse = ", ")),
        sprintf("  sigma_e_weak=%d structural_weak=%d",
                length(sigma_e_weak), length(structural_weak)),
        sprintf("  strength_gates_badge=%s", strength_gates),
        sprintf("  action: %s",
                if (pass && length(structural_weak) == 0L)
                  "Structural identification strength adequate (shock s.d. weakness is INFO-only)."
                else if (pass)
                  sprintf(paste0("Fisher rank is full -- PASS. INFO only: %s have |s_i| < 1 under ",
                                 "weighting=\"%s\", where s_i carries the moments' units and is not a ",
                                 "t-ratio. Re-run with weighting=\"sampling\" plus a moment_cov to get ",
                                 "a gating strength measure."),
                          paste(head(structural_weak, 3), collapse = ", "), weighting_actual)
                else if (fisher_rank_deficient)
                  "Fisher information rank deficient: some parameter directions are not identified (see D1). Add moments or calibrate."
                else sprintf("Structural parameters weakly identified (%s, sampling-weighted t-ratios). Add moments or calibrate.",
                             paste(head(structural_weak, 3), collapse = ", "))),
        sprintf("  note: s_i gate |s_i| < 1 (package choice, no literature source)%s; T=%s obs.",
                if (strength_gates) "" else " -- INFORMATIONAL under this weighting",
                if (!is.null(meta$T_obs)) as.character(meta$T_obs) else "(not supplied)"),
        sprintf("  top_col_norms: %s",
                paste(sprintf("%s=%.1e", names(top_cn), top_cn), collapse = ", "))
      ), collapse = "\n")
    )
}


#' Model-implied var/autocovariance moments as a function of theta (D20)
#'
#' Re-solves the model (steady state + first-order perturbation) at each
#' theta, so structural parameters AND shock standard deviations move the
#' moments. Names match \code{.compute_data_moments} /
#' \code{.d20_moment_sampling_cov}. A theta at which the model does not solve
#' returns NA moments (D20 then reports a non-finite Jacobian).
#' @noRd
.d20_data_moment_fn <- function(model, compiled, params, obs_names, max_lag = 4L) {
  sys_cache <- cache_system_structure(compiled)
  state <- new.env(parent = emptyenv())
  nm <- c(paste0("var_", obs_names),
          unlist(lapply(seq_len(max_lag), function(l) paste0("acv", l, "_", obs_names))))
  function(theta) {
    pp  <- .apply_theta_to_params(model, theta, params)
    sol <- .solve_dr_for_theta(model, compiled, sys_cache, pp, state)
    if (is.null(sol)) return(stats::setNames(rep(NA_real_, length(nm)), nm))
    mm <- compute_moments(sol$dr, model, n_ar = max_lag, params = sol$params)
    v  <- diag(mm$var_cov)[obs_names]
    acv <- lapply(seq_len(max_lag), function(l)
      vapply(obs_names, function(o) mm$autocorr[o, o, l], numeric(1)) * v)
    stats::setNames(c(v, unlist(acv)), nm)
  }
}


#' Sampling covariance of the var/autocovariance moment estimator (D20)
#'
#' \eqn{Var(\hat m) = \hat\Omega / T}, where \eqn{\hat\Omega} is the Bartlett
#' (Newey-West) long-run covariance of the per-period moment contributions
#' \eqn{y_{t}^2} and \eqn{y_t y_{t-k}} (k = 1..max_lag) of the demeaned data.
#' The truncation lag uses Andrews' (1991) AR(1) plug-in bandwidth
#' \eqn{S_T = 1.1447 (\hat\alpha(1) T)^{1/3}} (lags \eqn{j < S_T} get weight
#' \eqn{1 - j/S_T}); a fixed small lag badly under-states the long-run
#' variance of squared persistent series.
#' @param Y T x n_obs data matrix with column names.
#' @param max_lag Largest autocovariance lag in the moment vector.
#' @return named n_mom x n_mom matrix (names as \code{.compute_data_moments}),
#'   with attribute \code{"bandwidth"}.
#' @noRd
.d20_moment_sampling_cov <- function(Y, max_lag = 4L) {
  Y <- scale(as.matrix(Y), scale = FALSE)
  T_obs <- nrow(Y)
  nm <- colnames(Y) %||% seq_len(ncol(Y))
  cols <- list()
  for (lag in 0:max_lag) {
    lead <- Y[(lag + 1):T_obs, , drop = FALSE]
    lagd <- Y[1:(T_obs - lag), , drop = FALSE]
    prod_ <- rbind(matrix(0, lag, ncol(Y)), lead * lagd)
    colnames(prod_) <- paste0(if (lag == 0) "var" else paste0("acv", lag), "_", nm)
    cols[[lag + 1L]] <- prod_
  }
  u <- scale(do.call(cbind, cols), scale = FALSE)

  # Andrews (1991) AR(1) plug-in, unit weights, Bartlett kernel.
  num <- 0; den <- 0
  for (a in seq_len(ncol(u))) {
    x1 <- u[-1, a]; x0 <- u[-T_obs, a]
    if (sum(x0^2) <= 0) next
    r  <- min(max(sum(x1 * x0) / sum(x0^2), -0.99), 0.99)
    s2 <- mean((x1 - r * x0)^2)
    num <- num + 4 * r^2 * s2^2 / ((1 - r)^6 * (1 + r)^2)
    den <- den + s2^2 / (1 - r)^4
  }
  S_T <- if (den > 0) 1.1447 * (num / den * T_obs)^(1 / 3) else 1
  L <- min(max(ceiling(S_T) - 1L, 0L), T_obs - 1L)
  # .newey_west weights lag j by 1 - j/(L+1); L+1 = ceiling(S_T).
  Om <- .newey_west(u, max_lag = L) / T_obs
  dimnames(Om) <- list(colnames(u), colnames(u))
  attr(Om, "bandwidth") <- L + 1L
  Om
}
