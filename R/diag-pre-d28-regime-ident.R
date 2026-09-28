## R/diag-pre-d28-regime-ident.R
## --------------------------------------------------------------------------
## Phase H: D28 -- Regime-Switching Identification Diagnostic
##
## Assesses parameter identification in Markov-switching DSGE models.
##
## Conventions (checked, 0.9.4 refresh):
##   * transition_matrix is ROW-stochastic:
##       P[i, j] = Pr(s_t = j | s_{t-1} = i),  rowSums(P) = 1.
##   * the ergodic distribution is the LEFT unit eigenvector, pi' P = pi',
##     sum(pi) = 1, obtained by a direct linear solve (not power iteration,
##     which silently stops short for persistent chains).
##   * the ergodic moment map is the ergodic mixture of the regime-conditional
##     moments, m_bar = sum_r pi_r m_r(theta_r).  This is EXACT for moments that
##     are linear in the distribution (raw moments, E[y y'], ...) when each
##     regime-conditional distribution is that regime's own stationary law
##     (e.g. iid regime switching of a static model); for persistent regimes
##     with within-regime dynamics it is the standard "regime-conditional
##     moments" approximation.  Moments on a log scale must not be mixed, so
##     the generic solver returns raw variances and lag-1 autocovariances.
##   * the Markov-switching parameter vector is stacked: a parameter whose
##     value is common to all regimes is ONE column (sum_r pi_r J_r[, p]); a
##     parameter whose value differs across regimes gets one column per regime
##     (pi_r J_r[, p], named "p[regime]"); with estimate_transition = TRUE the
##     free ergodic weights pi_1..pi_{R-1} are columns too (m_r - m_R).  Mixture
##     moments depend on P only through pi, so they can identify at most R - 1
##     of the R (R - 1) free transition probabilities.
##   * label switching: m_bar is invariant to a joint relabelling of regimes, so
##     regime-specific parameters are identified only up to a permutation
##     (global, not local).  Regimes with identical moments are locally
##     unidentified (the weight columns vanish).
##
## References:
##   Davig & Leeper (2007). Generalizing the Taylor Principle. AER.
##   Foerster et al. (2016). Markov-switching DSGE models: a solution
##     algorithm. J. Econ. Dyn. Control.
##   Iskrev, N. (2010). Local identification in DSGE models.
##   Hamilton, J. D. (1994). Time Series Analysis, ch. 22 (ergodic probs).
##   Stephens, M. (2000). Dealing with label switching in mixture models.
## --------------------------------------------------------------------------

#' D28. Regime-Switching Identification Diagnostic
#'
#' Assesses how parameter identification varies across regimes in a
#' Markov-switching DSGE model, and whether the stacked Markov-switching
#' parameter vector is locally identified from the ergodic mixture of the
#' regime-conditional moments.
#'
#' The diagnostic:
#' \enumerate{
#'   \item For each regime, computes the moment Jacobian at the
#'     regime-specific parameter values (central differences at steps
#'     \code{eps} and \code{2 eps}) and its rank with the shared
#'     finite-difference-aware equilibrated rank test of D1.
#'   \item Computes the ergodic distribution \eqn{\pi' P = \pi'} of the
#'     row-stochastic transition matrix and the Jacobian of the ergodic
#'     mixture moments \eqn{\bar m = \sum_r \pi_r m_r(\theta_r)} with respect
#'     to the stacked parameter vector: common parameters (one column),
#'     regime-switching parameters (one column per regime, named
#'     \code{"p[regime]"}) and, if \code{estimate_transition}, the free
#'     ergodic weights (\code{"pi[regime]"}).
#'   \item Computes D20-style identification strength (\eqn{|t|}-ratios with
#'     identity moment weighting) per regime and for the stacked system.
#'   \item Reports label switching: regimes with indistinguishable moments,
#'     and switching parameters whose values are distinct across all regimes
#'     (an ordering constraint on one of them fixes the labels).
#' }
#'
#' \strong{What gates the badge.} The \emph{rank} tests only: the badge is
#' FAIL when any evaluated regime's moment Jacobian is rank deficient, or the
#' stacked ergodic Jacobian is. The \eqn{|t|}-type strengths
#' (\code{strength_threshold}, \code{params_weak_in_all_regimes},
#' \code{weak_params_by_regime}) are \strong{informational}: they are computed
#' under identity moment weighting, so they carry the units of the moments and
#' are not t-ratios in any sampling sense. More importantly, no
#' Markov-switching DSGE identification reference establishes a numeric
#' \eqn{|t|} cutoff: Foerster, Rubio-Ramirez, Waggoner & Zha (2016) is a
#' solution-method paper with no identification test, and RISE (Maih) checks
#' regime identification by reusing the Iskrev / Ratto-Iskrev Jacobian-rank
#' machinery per regime rather than by a strength threshold. The default
#' \code{strength_threshold = 1.0} is therefore a package convention with no
#' literature source, kept as a ranking hint.
#'
#' @param regime_defs       A named list of regime definitions.  Each element is
#'   a list with:
#'   \describe{
#'     \item{\code{params}}{Named numeric vector: parameter values in this
#'       regime (overriding \code{params}).}
#'     \item{\code{model_solve_fn}}{Function: theta -> named numeric moments for
#'       this regime.  Should return moments that are linear in the
#'       distribution (raw moments) so that the ergodic mixture is meaningful.
#'       If omitted, variances and lag-1 autocovariances from \code{model} are
#'       used.  A failed solve must return NULL or non-finite values.}
#'   }
#' @param transition_matrix Square ROW-stochastic Markov transition matrix,
#'   \code{P[i, j] = Pr(s_t = j | s_{t-1} = i)}.  If it has dimnames they must
#'   be the regime names (it is reordered to \code{regime_defs} order).  If
#'   NULL only the per-regime analysis is run.
#' @param model             A dynhr_mod object (generic solver).
#' @param dr                Baseline DecisionRules object (optional).
#' @param params            Named baseline parameter vector.
#' @param param_names       Parameters to analyse (default \code{names(params)}).
#' @param eps               Finite-difference step (default 1e-5).
#' @param strength_threshold \eqn{|t|}-ratio below which a parameter is
#'   reported as weak (default 1.0). INFORMATIONAL only -- it never changes
#'   the badge; see Details.
#' @param tol_rank          Optional relative rank tolerance (see D1); NULL uses
#'   the finite-difference-noise tolerance.
#' @param estimate_transition Logical: are the transition probabilities
#'   estimated (add the ergodic weights to the stacked system)? Default TRUE.
#' @param verbose           Print progress messages.
#' @param meta              Plot metadata.
#'
#' @return A \code{dynhr_diagnostic}; \code{pass} is FALSE when the stacked
#'   ergodic Jacobian is rank deficient or any evaluated regime's Jacobian is;
#'   NA when nothing could be assessed. Weak \eqn{|t|} strengths do not make
#'   it FALSE (see Details).
#' @noRd
d28_regime_switching_identification <- function(regime_defs,
                                                 transition_matrix,
                                                 model = NULL,
                                                 dr = NULL,
                                                 params = NULL,
                                                 param_names = NULL,
                                                 eps = 1e-5,
                                                 strength_threshold = 1.0,
                                                 tol_rank = NULL,
                                                 estimate_transition = TRUE,
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
  if (!is.list(regime_defs)) .dynhr_abort("D28: regime_defs must be a list.")

  regime_names <- names(regime_defs)
  if (is.null(regime_names)) {
    regime_names <- vapply(seq_len(n_regimes), function(r) {
      nm <- regime_defs[[r]]$name
      if (is.character(nm) && length(nm) == 1L) nm else paste0("regime", r)
    }, character(1))
  }
  if (anyDuplicated(regime_names) || any(!nzchar(regime_names))) {
    .dynhr_abort("D28: regime names must be unique and non-empty.")
  }
  names(regime_defs) <- regime_names

  if (!is.null(transition_matrix)) {
    transition_matrix <- .d28_check_transition(transition_matrix, regime_names)
  }

  # ---- 2. Regime parameter vectors ----
  theta_by_regime <- lapply(regime_defs, function(rd) {
    th <- params
    if (!is.null(rd$params)) {
      if (is.null(names(rd$params)) || any(!nzchar(names(rd$params)))) {
        .dynhr_abort("D28: regime_defs[[.]]$params must be a named numeric vector.")
      }
      th[names(rd$params)] <- rd$params
    }
    th
  })

  if (is.null(param_names)) {
    param_names <- names(params)
    if (is.null(param_names)) param_names <- names(theta_by_regime[[1]])
  }
  if (is.null(param_names) || length(param_names) == 0L) {
    return(.make_result(
      pass    = NA,
      summary = "D28 Regime-Switching Identification: No parameter names available."
    ))
  }
  n_par <- length(param_names)
  for (r in seq_len(n_regimes)) {
    miss <- setdiff(param_names, names(theta_by_regime[[r]]))
    if (length(miss) > 0L) {
      .dynhr_abort(sprintf("D28: regime '%s' has no value for parameter(s): %s.",
                           regime_names[r], paste(miss, collapse = ", ")))
    }
  }
  theta_mat <- vapply(theta_by_regime, function(th) as.numeric(th[param_names]),
                      numeric(n_par))
  theta_mat <- matrix(theta_mat, nrow = n_par,
                      dimnames = list(param_names, regime_names))

  # ---- 3. Ergodic distribution ----
  ergodic_dist <- NULL
  if (!is.null(transition_matrix)) {
    ergodic_dist <- .compute_ergodic_dist_ms(transition_matrix)
    if (verbose) {
      .dynhr_cat(sprintf("[d28] Ergodic distribution: %s\n",
                         paste(sprintf("%s=%.4f", regime_names, ergodic_dist),
                               collapse = ", ")))
    }
  }

  # ---- 4. Per-regime identification ----
  regime_jacobians <- stats::setNames(vector("list", n_regimes), regime_names)
  regime_jacobians2 <- regime_jacobians
  regime_moments <- regime_jacobians
  regime_unidentified <- regime_jacobians
  regime_strength <- matrix(NA_real_, nrow = n_par, ncol = n_regimes,
                            dimnames = list(param_names, regime_names))
  regime_ranks <- stats::setNames(rep(NA_integer_, n_regimes), regime_names)
  regime_fail <- stats::setNames(rep(NA_character_, n_regimes), regime_names)

  for (r in seq_len(n_regimes)) {
    rd <- regime_defs[[r]]
    r_name <- regime_names[r]
    lbl <- sprintf("[d28]   Regime %s:", r_name)

    solve_fn <- rd$model_solve_fn
    if (is.null(solve_fn)) {
      solve_fn <- .make_regime_solve_fn(model, dr, params, rd$params, r_name)
    }
    if (!is.function(solve_fn)) {
      regime_fail[r] <- "no model_solve_fn and no model"
      if (verbose) .dynhr_cat(paste0(lbl, " no solve_fn available, skipping.\n"))
      next
    }

    th_r <- stats::setNames(theta_mat[, r], param_names)
    ## the moment function sees the FULL regime vector; only param_names move
    fn_r <- local({
      full <- theta_by_regime[[r]]
      user_fn <- solve_fn
      function(th) { full[names(th)] <- th; user_fn(full) }
    })
    f0 <- fn_r(th_r)
    if (!is.numeric(f0) || length(f0) == 0L || any(!is.finite(f0))) {
      regime_fail[r] <- "moment function failed at the regime parameters"
      if (verbose) .dynhr_cat(paste0(lbl, " moment evaluation failed.\n"))
      next
    }
    J  <- .numerical_jacobian(fn_r, th_r, eps = eps)
    J2 <- .numerical_jacobian(fn_r, th_r, eps = 2 * eps)
    if (any(!is.finite(J)) || any(!is.finite(J2))) {
      regime_fail[r] <- "moment function failed inside the finite-difference stencil"
      if (verbose) .dynhr_cat(paste0(lbl, " Jacobian failed.\n"))
      next
    }
    m_names <- names(f0) %||% paste0("m_", seq_along(f0))
    dimnames(J) <- dimnames(J2) <- list(m_names, param_names)
    f0 <- stats::setNames(as.numeric(f0), m_names)

    regime_moments[[r]] <- f0
    regime_jacobians[[r]] <- J
    regime_jacobians2[[r]] <- J2

    rk <- .ident_equilibrated_rank(J, J2, tol_rank = tol_rank)
    regime_ranks[r] <- as.integer(rk$rank)
    regime_unidentified[[r]] <- rk$unidentified_params
    regime_strength[, r] <- .d28_strength(J, th_r)

    if (verbose)
      .dynhr_cat(sprintf("%s rank=%d of %d\n", lbl, regime_ranks[r], n_par))
  }
  valid <- !is.na(regime_ranks)

  # ---- 5. Stacked ergodic (Markov-switching) identification ----
  switching_params <- param_names[apply(theta_mat, 1, function(v) diff(range(v)) > 0)]
  ergodic_jacobian <- NULL
  ergodic_strength <- NULL
  ergodic_moments <- NULL
  ergodic_rank <- NA_integer_
  ergodic_unidentified <- character(0)
  ergodic_col_info <- NULL
  ergodic_note <- NULL

  if (is.null(ergodic_dist)) {
    ergodic_note <- "no transition_matrix: ergodic analysis not run"
  } else if (!all(valid)) {
    ergodic_note <- sprintf("regime(s) %s failed: ergodic analysis not run",
                            paste(regime_names[!valid], collapse = ", "))
  } else {
    ref_names <- rownames(regime_jacobians[[1]])
    for (r in seq_len(n_regimes)) {
      if (!identical(rownames(regime_jacobians[[r]]), ref_names)) {
        .dynhr_abort(sprintf(paste0(
          "D28: regime '%s' returns moments [%s] but regime '%s' returns [%s]; ",
          "the ergodic mixture needs the same moments in the same order."),
          regime_names[r], paste(rownames(regime_jacobians[[r]]), collapse = ", "),
          regime_names[1], paste(ref_names, collapse = ", ")))
      }
    }
    built  <- .d28_stack_jacobian(regime_jacobians, regime_moments, ergodic_dist,
                                  theta_mat, switching_params, estimate_transition)
    built2 <- .d28_stack_jacobian(regime_jacobians2, regime_moments, ergodic_dist,
                                  theta_mat, switching_params, estimate_transition)
    ergodic_jacobian <- built$J
    ergodic_col_info <- built$info
    ergodic_moments <- Reduce(`+`, Map(function(m, w) w * m, regime_moments,
                                       as.list(ergodic_dist)))
    rk_e <- .ident_equilibrated_rank(ergodic_jacobian, built2$J, tol_rank = tol_rank)
    ergodic_rank <- as.integer(rk_e$rank)
    ergodic_unidentified <- rk_e$unidentified_params
    ergodic_strength <- .d28_strength(ergodic_jacobian, built$info$value)
    names(ergodic_strength) <- built$info$column
  }

  # ---- 6. Label switching ----
  label_switching <- .d28_label_switching(regime_moments[valid], theta_mat[, valid, drop = FALSE],
                                          switching_params)
  if (!all(valid)) {
    label_switching$note <- sprintf("not fully assessed (regime(s) %s failed)",
                                    paste(regime_names[!valid], collapse = ", "))
  }

  # ---- 7. Cross-regime dependence ----
  regime_dependent <- character(0)
  weak_by_regime <- list()
  params_weak_in_all <- character(0)
  for (i in seq_len(n_par)) {
    pname <- param_names[i]
    s_i <- regime_strength[i, ]
    fin <- is.finite(s_i)
    weak_r <- which(fin & s_i < strength_threshold)
    if (length(weak_r) > 0L) weak_by_regime[[pname]] <- regime_names[weak_r]
    if (any(fin) && length(weak_r) == sum(fin)) {
      params_weak_in_all <- c(params_weak_in_all, pname)
    }
    s_valid <- s_i[fin & s_i > 0]
    if (length(s_valid) >= 2L) {
      cv <- stats::sd(s_valid) / max(mean(s_valid), 1e-16)
      if (cv > 1.0) regime_dependent <- c(regime_dependent, pname)
    }
  }
  dep_and_weak <- intersect(regime_dependent, names(weak_by_regime))
  ergodic_ok <- if (is.na(ergodic_rank)) NA else ergodic_rank == ncol(ergodic_jacobian)

  # The badge is driven by the RANK tests only (per-regime and ergodic). The
  # |t|-type strengths stay informational: no MS-DSGE identification paper
  # certifies a numeric |t| cutoff, so strength_threshold is a package
  # convention, not a verdict. See the Details section.
  rank_deficient_regimes <- regime_names[valid][regime_ranks[valid] < n_par]

  pass <- if (!any(valid)) {
    NA
  } else if (length(rank_deficient_regimes) > 0L || isFALSE(ergodic_ok)) {
    FALSE
  } else if (isTRUE(ergodic_ok)) {
    TRUE
  } else {
    NA  # per-regime checks clean, but the Markov-switching system was not assessed
  }

  # ---- 8. Plots ----
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE) && any(valid)) {
    plots$regime_strength <- .apply_meta(
      .d28_plot_regime_strength(regime_strength, strength_threshold,
                                regime_ranks, n_par, ergodic_dist),
      meta)
    if (!is.null(ergodic_strength)) {
      plots$ergodic_comparison <- .apply_meta(
        .d28_plot_ergodic(ergodic_strength, ergodic_col_info, ergodic_unidentified,
                          strength_threshold, ergodic_rank),
        meta)
    }
  }

  # ---- 9. Result ----
  result <- list(
    regime_jacobians           = regime_jacobians,
    regime_moments             = regime_moments,
    regime_params              = theta_mat,
    ergodic_jacobian           = ergodic_jacobian,
    ergodic_columns            = ergodic_col_info,
    ergodic_moments            = ergodic_moments,
    per_regime_strength        = regime_strength,
    ergodic_strength           = ergodic_strength,
    regime_ranks               = regime_ranks,
    rank_deficient_regimes     = rank_deficient_regimes,
    regime_failures            = regime_fail[!valid],
    regime_unidentified_params = regime_unidentified,
    ergodic_rank               = ergodic_rank,
    ergodic_n_params           = if (is.null(ergodic_jacobian)) NA_integer_ else ncol(ergodic_jacobian),
    ergodic_unidentified_params = ergodic_unidentified,
    ergodic_note               = ergodic_note,
    switching_params           = switching_params,
    regime_dependent_params    = regime_dependent,
    weak_params_by_regime      = weak_by_regime,
    params_weak_in_all_regimes = params_weak_in_all,
    label_switching            = label_switching,
    transition_matrix          = transition_matrix,
    ergodic_dist               = ergodic_dist,
    n_regimes                  = n_regimes,
    regime_names               = regime_names
  )

  rank_str <- if (is.na(ergodic_rank)) {
    sprintf("per-regime ranks [%s] of %d; %s",
            paste(sprintf("%s=%s", regime_names, regime_ranks), collapse = ", "),
            n_par, ergodic_note)
  } else {
    sprintf("Markov-switching system rank = %d of %d%s", ergodic_rank,
            ncol(ergodic_jacobian),
            if (length(ergodic_unidentified) > 0L)
              sprintf(" (unidentified: %s)", paste(ergodic_unidentified, collapse = ", "))
            else "")
  }
  pi_str <- if (is.null(ergodic_dist)) "" else
    sprintf(" Ergodic probs: %s.", paste(sprintf("%s=%.3f", regime_names, ergodic_dist),
                                         collapse = ", "))
  dep_str <- if (length(regime_dependent) > 0L)
    sprintf(" Regime-dependent strength: %s.", paste(regime_dependent, collapse = ", "))
  else ""
  weak_str <- if (length(params_weak_in_all) > 0L)
    sprintf(" INFO (does not affect the badge): |t| < %.3g in ALL regimes: %s.",
            strength_threshold, paste(params_weak_in_all, collapse = ", "))
  else ""
  depweak_str <- if (length(dep_and_weak) > 0L)
    sprintf(" INFO: |t| < %.3g in some regime: %s.",
            strength_threshold, paste(dep_and_weak, collapse = ", "))
  else ""
  rankdef_str <- if (length(rank_deficient_regimes) > 0L)
    sprintf(" RANK DEFICIENT in regime(s): %s.",
            paste(rank_deficient_regimes, collapse = ", "))
  else ""
  label_str <- sprintf(" Labels: %s.", label_switching$note)

  summary_text <- sprintf("D28 Regime-Switching Identification: %d regimes, %s.%s%s%s%s%s%s",
                          n_regimes, rank_str, pi_str, rankdef_str, dep_str, weak_str,
                          depweak_str, label_str)

  badge <- if (isTRUE(pass)) "PASS" else if (is.na(pass)) "INFO" else "FAIL"
  action <- if (!any(valid)) {
    "No regime could be evaluated; supply model_solve_fn per regime or a model."
  } else if (isFALSE(ergodic_ok)) {
    sprintf("Markov-switching parameters %s are not locally identified from the ergodic moments -- add moments, calibrate, or fix the transition probabilities.",
            paste(head(ergodic_unidentified, 3), collapse = ", "))
  } else if (length(rank_deficient_regimes) > 0L) {
    sprintf("The moment Jacobian is rank deficient in regime(s) %s -- add moments or calibrate the affected parameters.",
            paste(head(rank_deficient_regimes, 3), collapse = ", "))
  } else if (isTRUE(pass) && length(params_weak_in_all) > 0L) {
    sprintf(paste0("Rank tests pass -- PASS. INFO only: %s have |t| < %.3g in every regime under ",
                   "identity moment weighting. No MS-DSGE identification paper certifies a |t| ",
                   "cutoff, so this is a ranking hint, not a verdict."),
            paste(head(params_weak_in_all, 3), collapse = ", "), strength_threshold)
  } else if (isTRUE(pass) && length(dep_and_weak) > 0L) {
    sprintf(paste0("Rank tests pass -- PASS. INFO only: regime-dependent parameters %s have ",
                   "|t| < %.3g in some regime; consider richer measurement if that matters."),
            paste(head(dep_and_weak, 3), collapse = ", "), strength_threshold)
  } else if (isTRUE(pass)) {
    "Regime-switching identification adequate (up to regime relabelling)."
  } else {
    "Per-regime identification adequate; supply transition_matrix to assess the Markov-switching system."
  }
  llm <- paste(c(
    sprintf("D28 | Regime-Switching Identification | %s", badge),
    sprintf("  n_regimes=%d ergodic_rank=%s/%s n_switching=%d n_regime_dependent=%d",
            n_regimes,
            if (is.na(ergodic_rank)) "NA" else as.character(ergodic_rank),
            if (is.null(ergodic_jacobian)) "NA" else as.character(ncol(ergodic_jacobian)),
            length(switching_params), length(regime_dependent)),
    if (!is.null(ergodic_dist))
      sprintf("  ergodic_dist: %s", paste(sprintf("%s=%.4f", regime_names, ergodic_dist),
                                         collapse = " ")),
    sprintf("  rank_deficient_regimes: %s",
            if (length(rank_deficient_regimes)) paste(rank_deficient_regimes, collapse = ", ") else "none"),
    if (length(params_weak_in_all) > 0L)
      sprintf("  globally_weak (all regimes, INFO-only, |t| < %.3g): %s",
              strength_threshold, paste(params_weak_in_all, collapse = ", ")),
    sprintf("  label_switching: %s", label_switching$note),
    sprintf("  action: %s", action)
  ), collapse = "\n")

  .make_result(
    result      = result,
    pass        = pass,
    plots       = plots,
    summary     = summary_text,
    llm_summary = llm
  )
}


# ==========================================================================
# Internal helpers for D28
# ==========================================================================

#' Validate / reorder a row-stochastic transition matrix
#' @noRd
.d28_check_transition <- function(P, regime_names) {
  n <- length(regime_names)
  if (!is.matrix(P) || !is.numeric(P) || nrow(P) != n || ncol(P) != n) {
    .dynhr_abort(sprintf("D28: transition_matrix must be a numeric %d x %d matrix.", n, n))
  }
  if (any(!is.finite(P)) || any(P < 0)) {
    .dynhr_abort("D28: transition_matrix entries must be finite and non-negative.")
  }
  if (any(abs(rowSums(P) - 1) > 1e-8)) {
    .dynhr_abort(paste0(
      "D28: transition_matrix rows must sum to 1 (row-stochastic convention ",
      "P[i, j] = Pr(s_t = j | s_{t-1} = i)); row sums are ",
      paste(format(rowSums(P), digits = 6), collapse = ", "),
      if (all(abs(colSums(P) - 1) <= 1e-8)) " -- it looks column-stochastic; pass t(P)." else "."))
  }
  rn <- rownames(P); cn <- colnames(P)
  if (!is.null(rn) || !is.null(cn)) {
    if (is.null(rn) || is.null(cn) || !identical(rn, cn) ||
        !setequal(rn, regime_names) || anyDuplicated(rn)) {
      .dynhr_abort(sprintf(
        "D28: transition_matrix dimnames must both be the regime names (%s).",
        paste(regime_names, collapse = ", ")))
    }
    P <- P[regime_names, regime_names, drop = FALSE]
  } else {
    dimnames(P) <- list(regime_names, regime_names)
  }
  P
}


#' |t|-ratio strength with identity moment weighting (D20 style)
#' @noRd
.d28_strength <- function(J, theta) {
  k <- ncol(J)
  I_inv <- .safe_sym_inv(crossprod(J) + diag(1e-10, k))
  se <- sqrt(pmax(diag(I_inv), 0))
  abs(as.numeric(theta)) / pmax(se, 1e-16)
}


#' Stacked Markov-switching Jacobian of the ergodic mixture moments
#' @return list(J, info) with info a data.frame(column, param, regime, value)
#' @noRd
.d28_stack_jacobian <- function(Js, moments, pi, theta_mat, switching,
                                estimate_transition) {
  R <- length(Js)
  rn <- names(Js)
  cols <- list(); info <- list()
  for (p in rownames(theta_mat)) {
    if (p %in% switching) {
      for (r in seq_len(R)) {
        nm <- sprintf("%s[%s]", p, rn[r])
        cols[[nm]] <- pi[r] * Js[[r]][, p]
        info[[nm]] <- data.frame(column = nm, param = p, regime = rn[r],
                                 value = theta_mat[p, r], stringsAsFactors = FALSE)
      }
    } else {
      cols[[p]] <- Reduce(`+`, lapply(seq_len(R), function(r) pi[r] * Js[[r]][, p]))
      info[[p]] <- data.frame(column = p, param = p, regime = "all regimes",
                              value = theta_mat[p, 1], stringsAsFactors = FALSE)
    }
  }
  if (isTRUE(estimate_transition) && R >= 2L) {
    for (r in seq_len(R - 1L)) {
      nm <- sprintf("pi[%s]", rn[r])
      ## pi_R = 1 - sum_{r<R} pi_r
      cols[[nm]] <- moments[[r]] - moments[[R]]
      info[[nm]] <- data.frame(column = nm, param = "ergodic prob.", regime = rn[r],
                               value = pi[r], stringsAsFactors = FALSE)
    }
  }
  J <- do.call(cbind, cols)
  dimnames(J) <- list(rownames(Js[[1]]), names(cols))
  info <- do.call(rbind, info)
  rownames(info) <- NULL
  list(J = J, info = info)
}


#' Label-switching report
#' @noRd
.d28_label_switching <- function(moments, theta_mat, switching) {
  rn <- names(moments)
  R <- length(rn)
  pairs <- character(0)
  if (R >= 2L) {
    for (a in seq_len(R - 1L)) for (b in seq(a + 1L, R)) {
      ma <- moments[[a]]; mb <- moments[[b]]
      if (length(ma) == length(mb) &&
          max(abs(ma - mb)) <= 1e-10 * max(1, abs(ma), abs(mb))) {
        pairs <- c(pairs, paste(rn[a], rn[b], sep = "~"))
      }
    }
  }
  ordering <- switching[vapply(switching, function(p)
    !anyDuplicated(theta_mat[p, ]), logical(1))]
  note <- if (R < 2L) {
    "single regime, no relabelling"
  } else if (length(pairs) > 0L) {
    sprintf("regimes %s have identical moments (weights/labels not identified)",
            paste(pairs, collapse = ", "))
  } else if (length(ordering) > 0L) {
    sprintf("identified only up to relabelling; an ordering constraint on %s fixes it",
            ordering[1])
  } else if (length(switching) > 0L) {
    "identified only up to relabelling; no switching parameter is distinct across all regimes"
  } else {
    "regimes differ only through their moment functions"
  }
  list(indistinguishable_pairs = pairs, ordering_params = ordering, note = note)
}


#' log10-scale |t| fill with |t| labels
#' @noRd
.d28_fill_scale <- function() {
  scale_fill_dynhr_cividis(
    name = "|t| (log scale)",
    breaks = function(l) unique(floor(l[1]):ceiling(l[2])),
    labels = function(b) formatC(10^b, digits = 2, format = "g"))
}


#' @noRd
.d28_plot_regime_strength <- function(S, threshold, ranks, n_par, pi) {
  df <- as.data.frame.table(S, responseName = "Strength", stringsAsFactors = FALSE)
  colnames(df) <- c("Parameter", "Regime", "Strength")
  df$Parameter <- factor(df$Parameter, levels = rev(rownames(S)))
  xlab <- sprintf("%s\nrank %s/%d%s", colnames(S),
                  ifelse(is.na(ranks), "NA", ranks), n_par,
                  if (is.null(pi)) "" else sprintf(", \u03c0=%.2f", pi))
  df$Regime <- factor(df$Regime, levels = colnames(S), labels = xlab)
  df$log10S <- log10(pmax(df$Strength, 1e-12))
  df$label <- ifelse(is.finite(df$Strength),
                     paste0(formatC(df$Strength, digits = 2, format = "g"),
                            ifelse(df$Strength < threshold, " (weak)", "")),
                     "failed")
  ggplot2::ggplot(df, ggplot2::aes(x = Regime, y = Parameter, fill = log10S)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = label,
                                    colour = log10S > stats::median(log10S, na.rm = TRUE)),
                       size = 3.5, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "white"),
                                 na.value = "black") +
    .d28_fill_scale() +
    theme_dynhr_diagnostic() +
    ggplot2::labs(
      title = "D28: Identification strength by regime",
      subtitle = sprintf("|t| = |theta_r| / se_r at regime-specific values; weak if |t| < %g",
                         threshold),
      x = NULL, y = NULL
    )
}


#' @noRd
.d28_plot_ergodic <- function(strength, info, unid, threshold, rank) {
  df <- info
  df$Strength <- as.numeric(strength)
  df$log10S <- log10(pmax(df$Strength, 1e-12))
  df$unid <- df$column %in% unid
  df$label <- paste0(formatC(df$Strength, digits = 2, format = "g"),
                     ifelse(df$unid, " (unid.)",
                            ifelse(df$Strength < threshold, " (weak)", "")))
  df$param <- factor(df$param, levels = rev(unique(df$param)))
  df$regime <- factor(df$regime, levels = unique(c("all regimes", df$regime)))
  ggplot2::ggplot(df, ggplot2::aes(x = regime, y = param, fill = log10S)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_tile(data = df[df$unid, , drop = FALSE], fill = NA,
                       colour = unname(tol_vibrant["red"]), linewidth = 1.2) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.5,
                       colour = ifelse(df$log10S > stats::median(df$log10S), "black", "white")) +
    .d28_fill_scale() +
    theme_dynhr_diagnostic() +
    ggplot2::labs(
      title = "D28: Markov-switching system identification (ergodic moments)",
      subtitle = sprintf(paste0("Rank %d of %d stacked parameters; 'all regimes' = common ",
                                "parameter; red outline = loads on an unidentified direction"),
                         rank, nrow(df)),
      x = NULL, y = NULL
    )
}


#' Build a regime-specific solve function (generic model path)
#'
#' Returns raw variances and lag-1 autocovariances of all model variables
#' (linear in the distribution, so their ergodic mixture is meaningful).
#' @return Function: theta -> named moments (NULL on failure), or NULL.
#' @noRd
.make_regime_solve_fn <- function(model, dr, params, regime_params, regime_name) {
  if (is.null(model)) return(NULL)

  theta_base <- params
  if (!is.null(regime_params)) {
    theta_base[names(regime_params)] <- regime_params
  }

  function(theta) {
    theta_full <- theta_base
    theta_full[names(theta)] <- theta
    dr_new <- .stoch_simul_internal_diag(model, theta_full, dr_order = 1L)
    if (is.null(dr_new) || is.null(dr_new$ghx)) return(NULL)
    moments <- .moments_from_dr(dr_new, model = model, params = theta_full)
    if (is.null(moments) || is.null(moments$acf_y)) return(NULL)
    sig <- as.matrix(moments$sigma_y)
    ac <- moments$acf_y
    k <- nrow(sig)
    ii <- seq_len(k)
    ## diag() of a 1x1 drop()ed array is a 0x0 matrix -- index explicitly.
    v <- sig[cbind(ii, ii)]
    rho1 <- if (length(dim(ac)) == 3L) ac[cbind(ii, ii, 1L)] else as.matrix(ac)[cbind(ii, ii)]
    vn <- rownames(sig) %||% paste0("y", seq_len(k))
    stats::setNames(c(v, rho1 * v), c(paste0("var_", vn), paste0("acov1_", vn)))
  }
}


#' Ergodic distribution of a row-stochastic transition matrix (D28 wrapper)
#'
#' Thin alias for the shared \code{.ergodic_dist()} (\code{R/ms-spec.R}) so
#' D28's error messages keep their \code{D28:} prefix.  There is exactly ONE
#' implementation in the package (0.9.4); \code{ramsey-regime.R} and
#' \code{ms-spec.R} used to carry power-iteration copies that were wrong for
#' persistent chains.
#'
#' @param P Transition matrix (rows sum to 1)
#' @return Named numeric vector of ergodic probabilities
#' @noRd
.compute_ergodic_dist_ms <- function(P) {
  .ergodic_dist(P, "D28: transition_matrix")
}
