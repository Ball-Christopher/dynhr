## R/diag-pre-d37-komunjer-ng.R
## --------------------------------------------------------------------------
## D37: Komunjer & Ng (2011, Econometrica) dynamic identification rank check.
##
## Implements the rank condition for local identification of theta from the
## spectral density / autocovariance generating function of the ABCD
## state-space form:
##   s_t = A(theta) s_{t-1} + B(theta) eps_t
##   y_t = C(theta) s_{t-1} + D(theta) eps_t
##
## Identification requires:
##   1. (A, B, C, D) is a MINIMAL realisation (controllable + observable).
##   2. The stacked Jacobian Delta(theta) = [Delta_Lambda  Delta_T  Delta_U]
##      has full column rank, where:
##        - Delta_Lambda = d vec(A,B,C,D) / d theta'
##        - Delta_T      = tangent directions of similarity transforms
##                          T_sim s -> (E A - A E, E B, -C E, 0)
##        - Delta_U      = tangent directions of orthogonal shock rotations
##                          (square case n_y = n_e only):
##                          (0, B S, 0, D S) for skew-symmetric S
##      Required rank:
##        - square case (n_e == n_y):  n_theta + n_x^2 + n_e(n_e-1)/2
##        - singular case (n_e < n_y): n_theta + n_x^2  (Delta_U dropped)
##        - n_e > n_y: out of scope (KN-S), reported as INFO
## --------------------------------------------------------------------------

#' D37. Dynamic identification rank check (Komunjer & Ng 2011)
#'
#' Computes the Komunjer & Ng (2011) rank condition for local identification
#' of the deep parameter vector \code{theta} from the ABCD state-space
#' representation of the first-order solution. Reuses
#' \code{build_dsge_state_space()} (via \code{model_solve_fn}) to obtain
#' \eqn{A = T}, \eqn{B = R}, \eqn{C = Z}, \eqn{D} for a perturbed
#' \code{theta}, folds \code{Sigma_e} into \code{B} and \code{D} via a
#' Cholesky factor (so the innovations are standardised), and stacks:
#' \itemize{
#'   \item \eqn{\Delta_\Lambda}: numerical Jacobian of
#'     \eqn{vec(A,B,C,D)} w.r.t. \code{theta} (central finite differences).
#'   \item \eqn{\Delta_T}: similarity-transform tangent directions
#'     (\eqn{(EA - AE, EB, -CE, 0)} for each elementary basis matrix
#'     \eqn{E_{ij}}).
#'   \item \eqn{\Delta_U}: orthogonal shock-rotation tangent directions
#'     (\eqn{(0, BS, 0, DS)} for skew-symmetric basis \eqn{S}); included
#'     only in the square case \eqn{n_e = n_y}.
#' }
#' Local identification holds iff
#' \eqn{\Delta = [\Delta_\Lambda\ \Delta_T\ \Delta_U]} has full column rank
#' (rank \eqn{n_\theta + n_x^2 + n_e(n_e-1)/2} in the square case, or
#' \eqn{n_\theta + n_x^2} in the singular case \eqn{n_e < n_y}).
#'
#' A pre-check verifies that \code{(A,B,C,D)} is a minimal realisation
#' (controllable + observable). DSGE state vectors routinely carry redundant
#' (non-minimal) states, so a non-minimal representation does NOT fail the
#' diagnostic outright -- it is reported as \code{pass = NA} (INFO) because
#' the KN rank condition, as stated, applies to minimal realisations only.
#'
#' @param dr             Optional decision-rule object (informational; not
#'                        used by the rank check itself).
#' @param model          Optional parsed model object (informational; not
#'                        used by the rank check itself).
#' @param theta          Numeric vector -- parameter values at the
#'                        evaluation point.
#' @param param_names    Optional character vector of parameter names
#'                        (length = \code{length(theta)}).
#' @param model_solve_fn Function: \code{theta -> list(A, B, C, D)}
#'                        (optionally with \code{$Sigma_e}), re-solving the
#'                        first-order state space at each perturbed
#'                        \code{theta} (e.g. wrapping
#'                        \code{solve_perturbation} +
#'                        \code{build_dsge_state_space}). REQUIRED for the
#'                        rank check: \code{NULL} degrades to
#'                        \code{pass = NA} with an explanatory summary.
#' @param Sigma_e        Optional shock covariance matrix (n_e x n_e). If
#'                        \code{NULL}, looked up on the result of
#'                        \code{model_solve_fn(theta)} (\code{$Sigma_e}), and
#'                        otherwise defaults to the identity.
#' @param obs_vars       Optional character vector of observable names
#'                        (informational only; included in the summary).
#' @param eps            Relative step size for central finite differences
#'                        (default \code{1e-6}; actual step is
#'                        \code{max(eps * |theta_j|, 1e-7)}).
#' @param meta           Plot-provenance metadata (see \code{diag_meta()}).
#' @param ...            Additional arguments (reserved / ignored).
#' @return A \code{dynhr_diagnostic} list.
#' @references
#'   Komunjer, I., & Ng, S. (2011). Dynamic identification of dynamic
#'     stochastic general equilibrium models. \emph{Econometrica}, 79(6),
#'     1995-2032.
#' @noRd
d37_komunjer_ng <- function(dr             = NULL,
                            model          = NULL,
                            theta          = NULL,
                            param_names    = NULL,
                            model_solve_fn = NULL,
                            Sigma_e        = NULL,
                            obs_vars       = NULL,
                            eps            = 1e-6,
                            meta           = NULL,
                            ...) {

  # ---------------------------------------------------------------
  # 0. Input validation / graceful degradation
  # ---------------------------------------------------------------
  if (is.null(theta)) {
    return(.make_result(
      result  = NULL,
      pass    = NA,
      plots   = list(),
      summary = "D37 Komunjer-Ng: theta not provided. Skipped."
    ))
  }

  n_theta <- length(theta)
  if (is.null(param_names)) {
    param_names <- names(theta) %||% paste0("theta_", seq_len(n_theta))
  }
  if (length(param_names) != n_theta) {
    warning(sprintf(
      "d37: param_names length (%d) != length(theta) (%d). Using generic labels.",
      length(param_names), n_theta))
    param_names <- paste0("theta_", seq_len(n_theta))
  }

  # ---------------------------------------------------------------
  # 1. Resolve model_solve_fn: theta -> list(A, B, C, D[, Sigma_e])
  #    Delta_Lambda needs the state space re-solved at perturbed theta, so a
  #    user-supplied closure is REQUIRED (same injection pattern as D23's
  #    dr_solve_fn). Without one the check cannot run; report INFO, not FAIL.
  # ---------------------------------------------------------------
  if (is.null(model_solve_fn)) {
    return(.make_result(
      result  = NULL,
      pass    = NA,
      plots   = list(),
      summary = paste(
        "D37 Komunjer-Ng: the rank check differentiates the ABCD state space",
        "with respect to theta, which requires a re-solve closure.",
        "Provide `model_solve_fn = function(theta) list(A,B,C,D[,Sigma_e])`",
        "(e.g. wrapping solve_perturbation + build_dsge_state_space).",
        "Skipped.")
    ))
  }

  # ---------------------------------------------------------------
  # 2. Baseline ABCD + Sigma_e, fold Sigma_e into B and D
  # ---------------------------------------------------------------
  ss0 <- model_solve_fn(theta)
  if (is.null(ss0) || !all(c("A", "B", "C", "D") %in% names(ss0))) {
    return(.make_result(
      result  = NULL,
      pass    = NA,
      plots   = list(),
      summary = "D37 Komunjer-Ng: model_solve_fn(theta) must return a list with A, B, C, D matrices. Skipped."
    ))
  }

  A0 <- as.matrix(ss0$A)
  B0 <- as.matrix(ss0$B)
  C0 <- as.matrix(ss0$C)
  D0 <- as.matrix(ss0$D)

  n_x <- nrow(A0)
  n_e <- ncol(B0)
  n_y <- nrow(C0)

  Sigma_e_use <- Sigma_e %||% ss0$Sigma_e
  if (is.null(Sigma_e_use)) Sigma_e_use <- diag(n_e)
  Sigma_e_use <- as.matrix(Sigma_e_use)
  if (nrow(Sigma_e_use) != n_e || ncol(Sigma_e_use) != n_e) {
    warning(sprintf(
      "d37: Sigma_e dims (%dx%d) != n_e=%d. Falling back to identity.",
      nrow(Sigma_e_use), ncol(Sigma_e_use), n_e))
    Sigma_e_use <- diag(n_e)
  }

  # Standardising rotation: post-multiply B and D by chol(Sigma_e)' so the
  # innovations entering B/D are unit-variance, uncorrelated. chol() returns
  # the upper-triangular factor U with Sigma_e = U'U, so eps_std = U^{-T} eps
  # and B_std = B U', D_std = D U'.
  L_sig <- tryCatch({
    U <- chol(Sigma_e_use)
    t(U)
  }, error = function(e) {
    # Sigma_e not PD (e.g. near-zero shock variance): fall back to a
    # symmetric square root via eigen-decomposition (handles PSD).
    eg <- eigen(Sigma_e_use, symmetric = TRUE)
    vals <- pmax(eg$values, 0)
    eg$vectors %*% diag(sqrt(vals), n_e) %*% t(eg$vectors)
  })

  .standardise <- function(ss) {
    A <- as.matrix(ss$A)
    B <- as.matrix(ss$B) %*% L_sig
    C <- as.matrix(ss$C)
    D <- as.matrix(ss$D) %*% L_sig
    list(A = A, B = B, C = C, D = D)
  }

  ss0_std <- .standardise(ss0)
  A <- ss0_std$A; B <- ss0_std$B; C <- ss0_std$C; D <- ss0_std$D

  # ---------------------------------------------------------------
  # 3. Minimality pre-check: controllability + observability rank
  # ---------------------------------------------------------------
  .ctrb_obsv_rank <- function(A, B, C) {
    n <- nrow(A)
    # Controllability matrix [B, AB, A^2 B, ..., A^{n-1} B]
    ctrb <- B
    Ak <- diag(n)
    if (n > 1) {
      for (k in seq_len(n - 1)) {
        Ak <- Ak %*% A
        ctrb <- cbind(ctrb, Ak %*% B)
      }
    }
    # Observability matrix [C; CA; CA^2; ...; CA^{n-1}]
    obsv <- C
    Ak <- diag(n)
    if (n > 1) {
      for (k in seq_len(n - 1)) {
        Ak <- Ak %*% A
        obsv <- rbind(obsv, C %*% Ak)
      }
    }
    list(
      ctrb_rank = .svd_rank(svd(ctrb)$d, dim(ctrb)),
      obsv_rank = .svd_rank(svd(obsv)$d, dim(obsv)),
      ctrb = ctrb,
      obsv = obsv
    )
  }

  mo <- .ctrb_obsv_rank(A, B, C)
  is_minimal <- (mo$ctrb_rank == n_x) && (mo$obsv_rank == n_x)

  # ---------------------------------------------------------------
  # 4. n_e > n_y: KN-S (singular, more shocks than observables) is
  #    out of scope -- report INFO.
  # ---------------------------------------------------------------
  if (n_e > n_y) {
    return(.make_result(
      result = list(
        state_space = list(A = A0, B = B0, C = C0, D = D0, Sigma_e = Sigma_e_use),
        n_x = n_x, n_y = n_y, n_e = n_e, n_theta = n_theta,
        minimal = is_minimal,
        ctrb_rank = mo$ctrb_rank, obsv_rank = mo$obsv_rank
      ),
      pass = NA,
      plots = list(),
      summary = sprintf(
        paste("D37 Komunjer-Ng: n_e (%d) > n_y (%d). The KN (2011) conditions",
              "implemented here cover the square (n_e = n_y) and singular",
              "(n_e < n_y) cases; the n_e > n_y case (KN-S with more shocks",
              "than observables) is out of scope. INFO only."),
        n_e, n_y)
    ))
  }

  # ---------------------------------------------------------------
  # 5. Minimality gate: report pass = NA if non-minimal, but continue
  #    to compute Delta for informational purposes.
  # ---------------------------------------------------------------
  minimality_note <- if (!is_minimal) {
    sprintf(
      paste("Representation is NON-MINIMAL: controllability rank = %d/%d,",
            "observability rank = %d/%d. The Komunjer-Ng (2011) rank",
            "condition applies to minimal (controllable + observable)",
            "realisations; reducing to a minimal realisation is out of",
            "scope for this diagnostic. Reporting pass = NA (INFO);",
            "the rank numbers below are still computed but not directly",
            "interpretable as a pass/fail of the KN condition."),
      mo$ctrb_rank, n_x, mo$obsv_rank, n_x)
  } else ""

  # ---------------------------------------------------------------
  # 6. Delta_Lambda: numerical Jacobian of vec(A,B,C,D) wrt theta
  #    (central finite differences, relative step with floor)
  # ---------------------------------------------------------------
  n_vec <- n_x * n_x + n_x * n_e + n_y * n_x + n_y * n_e

  .vec_abcd <- function(th) {
    ss_th <- model_solve_fn(th)
    ss_std <- .standardise(ss_th)
    c(as.numeric(ss_std$A), as.numeric(ss_std$B),
      as.numeric(ss_std$C), as.numeric(ss_std$D))
  }

  Delta_Lambda <- matrix(NA_real_, nrow = n_vec, ncol = n_theta)
  for (j in seq_len(n_theta)) {
    h <- max(eps * abs(theta[j]), 1e-7)
    th_plus  <- theta; th_plus[j]  <- theta[j] + h
    th_minus <- theta; th_minus[j] <- theta[j] - h
    f_plus  <- tryCatch(.vec_abcd(th_plus),  error = function(e) rep(NA_real_, n_vec))
    f_minus <- tryCatch(.vec_abcd(th_minus), error = function(e) rep(NA_real_, n_vec))
    if (length(f_plus) != n_vec || length(f_minus) != n_vec) {
      Delta_Lambda[, j] <- NA_real_
    } else {
      Delta_Lambda[, j] <- (f_plus - f_minus) / (2 * h)
    }
  }
  if (any(!is.finite(Delta_Lambda))) {
    Delta_Lambda[!is.finite(Delta_Lambda)] <- 0
  }
  colnames(Delta_Lambda) <- param_names

  # ---------------------------------------------------------------
  # 7. Delta_T: similarity-transform tangent directions
  #    For each elementary E_ij (n_x x n_x): vec(EA - AE, EB, -CE, 0)
  # ---------------------------------------------------------------
  n_simil <- n_x * n_x
  Delta_T <- matrix(0, nrow = n_vec, ncol = n_simil)
  col_k <- 0L
  if (n_x > 0) {
    for (jj in seq_len(n_x)) {
      for (ii in seq_len(n_x)) {
        col_k <- col_k + 1L
        E <- matrix(0, n_x, n_x)
        E[ii, jj] <- 1
        dA <- E %*% A - A %*% E
        dB <- E %*% B
        dC <- -C %*% E
        dD <- matrix(0, n_y, n_e)
        Delta_T[, col_k] <- c(as.numeric(dA), as.numeric(dB),
                              as.numeric(dC), as.numeric(dD))
      }
    }
  }

  # ---------------------------------------------------------------
  # 8. Delta_U: orthogonal shock-rotation tangent directions
  #    (square case only, n_e == n_y): for each skew-symmetric basis
  #    S_kl = e_k e_l' - e_l e_k' (k < l): vec(0, B S, 0, D S)
  # ---------------------------------------------------------------
  square_case <- (n_e == n_y) && (n_e >= 2)
  n_rot <- if (n_e >= 2) n_e * (n_e - 1L) / 2L else 0L
  Delta_U <- matrix(0, nrow = n_vec, ncol = 0)

  if (square_case) {
    Delta_U <- matrix(0, nrow = n_vec, ncol = n_rot)
    col_k <- 0L
    for (l in seq_len(n_e)) {
      for (k in seq_len(n_e)) {
        if (k >= l) next
        col_k <- col_k + 1L
        S <- matrix(0, n_e, n_e)
        S[k, l] <- 1
        S[l, k] <- -1
        dB <- B %*% S
        dD <- D %*% S
        Delta_U[, col_k] <- c(rep(0, n_x * n_x), as.numeric(dB),
                              rep(0, n_y * n_x), as.numeric(dD))
      }
    }
  }

  # ---------------------------------------------------------------
  # 9. Stack Delta and required rank
  # ---------------------------------------------------------------
  if (n_e == n_y) {
    Delta <- cbind(Delta_Lambda, Delta_T, Delta_U)
    required_rank <- n_theta + n_x * n_x + n_rot
    case_label <- "square (n_e = n_y)"
  } else {
    # n_e < n_y: singular case, drop Delta_U
    Delta <- cbind(Delta_Lambda, Delta_T)
    required_rank <- n_theta + n_x * n_x
    case_label <- "singular (n_e < n_y)"
  }

  sv <- svd(Delta)
  singular_values <- sv$d
  names(singular_values) <- paste0("sv_", seq_along(singular_values))

  tol <- max(dim(Delta)) * max(singular_values, 0) * .Machine$double.eps * 100
  rank_Delta <- sum(singular_values > tol)
  full_rank  <- (rank_Delta >= required_rank) && (rank_Delta <= min(dim(Delta)))

  pass <- if (!is_minimal) NA else full_rank

  # ---------------------------------------------------------------
  # 10. Deficient-direction loadings: when rank-deficient, report the
  #     parameter loadings (theta-block, first n_theta rows of v) of the
  #     right singular vectors corresponding to the smallest singular
  #     values, sorted by |loading|, full vector for each direction.
  # ---------------------------------------------------------------
  deficient_loadings <- NULL
  n_deficient <- max(0L, required_rank - rank_Delta)
  if (n_deficient > 0 && n_theta > 0) {
    n_v_cols <- ncol(sv$v)
    start_col <- max(1L, n_v_cols - n_deficient + 1L)
    if (start_col <= n_v_cols) {
      V_def <- sv$v[seq_len(n_theta), start_col:n_v_cols, drop = FALSE]
      deficient_loadings <- list()
      for (k in seq_len(ncol(V_def))) {
        loadings <- V_def[, k]
        ord <- order(abs(loadings), decreasing = TRUE)
        df_k <- data.frame(
          parameter = param_names[ord],
          loading   = loadings[ord],
          abs_loading = abs(loadings[ord]),
          stringsAsFactors = FALSE
        )
        sv_idx <- n_v_cols - ncol(V_def) + k
        df_k$direction <- sv_idx
        df_k$singular_value <- singular_values[sv_idx]
        deficient_loadings[[k]] <- df_k
      }
      deficient_loadings <- do.call(rbind, deficient_loadings)
      rownames(deficient_loadings) <- NULL
    }
  }

  # ---------------------------------------------------------------
  # 11. Plots: singular-value bar chart in house style
  # ---------------------------------------------------------------
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    sv_df <- data.frame(
      index = seq_along(singular_values),
      value = pmax(singular_values, 1e-300)
    )
    sv_df$status <- ifelse(seq_along(singular_values) <= rank_Delta,
                           "Full rank", "Deficient")

    p_sv <- ggplot2::ggplot(
      sv_df, ggplot2::aes(x = factor(index), y = value, fill = status)
    ) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_hline(yintercept = tol, linetype = "dashed",
                          colour = dynhr_colours$red, linewidth = 0.5) +
      ggplot2::scale_fill_manual(values = c("Full rank" = dynhr_colours$mid_blue,
                                            "Deficient" = dynhr_colours$red),
                                 name = NULL) +
      ggplot2::scale_y_log10() +
      theme_dynhr_diagnostic() +
      ggplot2::labs(
        title = "D37: Singular values of the Komunjer-Ng identification matrix",
        subtitle = sprintf("rank(Delta) = %d, required = %d (%s)",
                           rank_Delta, required_rank, case_label),
        x = "Singular value index", y = "Value (log scale)"
      )
    plots$singular_values <- .apply_meta(p_sv, meta)
  }

  # ---------------------------------------------------------------
  # 12. Summary text
  # ---------------------------------------------------------------
  obs_note <- if (!is.null(obs_vars)) sprintf(" obs = [%s].", paste(obs_vars, collapse = ", ")) else ""

  rank_status <- if (is.na(pass)) {
    "INFO -- minimality fails, KN rank condition not directly applicable"
  } else if (pass) {
    "PASS -- locally identified"
  } else {
    sprintf("FAIL -- rank deficient by %d direction(s)", required_rank - rank_Delta)
  }

  deficient_note <- if (!is.null(deficient_loadings) && nrow(deficient_loadings) > 0) {
    top_params <- unique(deficient_loadings$parameter)
    sprintf(" Deficient direction(s) load most heavily on: %s.",
            paste(head(top_params, 5), collapse = ", "))
  } else ""

  summary_text <- sprintf(
    paste0("D37 Komunjer-Ng (2011) dynamic identification: %s. ",
           "Case = %s. n_x=%d, n_y=%d, n_e=%d, n_theta=%d. ",
           "rank(Delta) = %d, required = %d (n_theta + n_x^2%s).%s%s%s"),
    rank_status, case_label, n_x, n_y, n_e, n_theta,
    rank_Delta, required_rank,
    if (n_e == n_y) sprintf(" + n_e(n_e-1)/2 = %d", n_rot) else "",
    obs_note,
    if (nchar(minimality_note) > 0) paste0(" ", minimality_note) else "",
    deficient_note
  )

  .make_result(
    result = list(
      A = A0, B = B0, C = C0, D = D0, Sigma_e = Sigma_e_use,
      A_std = A, B_std = B, C_std = C, D_std = D,
      Delta_Lambda = Delta_Lambda,
      Delta_T = Delta_T,
      Delta_U = Delta_U,
      Delta = Delta,
      singular_values = singular_values,
      svd = sv,
      rank = rank_Delta,
      required_rank = required_rank,
      n_x = n_x, n_y = n_y, n_e = n_e, n_theta = n_theta,
      case = case_label,
      minimal = is_minimal,
      ctrb_rank = mo$ctrb_rank,
      obsv_rank = mo$obsv_rank,
      deficient_loadings = deficient_loadings
    ),
    pass = pass,
    plots = plots,
    summary = summary_text,
    llm_summary = {
      badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
      sv_str <- paste(sprintf("%.3e", head(sort(singular_values), 5)), collapse = ", ")
      paste(c(
        sprintf("D37 | Komunjer-Ng Dynamic Identification | %s", badge),
        sprintf("  case=%s n_x=%d n_y=%d n_e=%d n_theta=%d", case_label, n_x, n_y, n_e, n_theta),
        sprintf("  rank(Delta)=%d required=%d minimal=%s (ctrb_rank=%d, obsv_rank=%d, n_x=%d)",
                rank_Delta, required_rank, as.character(is_minimal), mo$ctrb_rank, mo$obsv_rank, n_x),
        sprintf("  smallest_sv: %s", sv_str),
        if (!is.null(deficient_loadings) && nrow(deficient_loadings) > 0)
          sprintf("  deficient_loadings_top: %s",
                  paste(sprintf("%s=%.3f", head(deficient_loadings$parameter, 5),
                                head(deficient_loadings$loading, 5)), collapse = ", ")),
        sprintf("  action: %s",
                if (is.na(pass))
                  "Representation non-minimal -- KN rank condition not directly applicable. If states are intentionally redundant (e.g. ME blocks), this is expected; otherwise consider removing redundant states."
                else if (pass)
                  "Model is dynamically (locally) identified per Komunjer & Ng (2011)."
                else
                  sprintf("Rank deficient (rank=%d, expected=%d). Parameters %s are not locally identified from the spectral density / autocovariances. Consider reparameterisation, calibration, or additional observables.",
                          rank_Delta, required_rank,
                          paste(head(unique(deficient_loadings$parameter), 3), collapse = ", ")))
      ), collapse = "\n")
    }
  )
}
