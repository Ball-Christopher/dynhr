## R/diag-pre-d37-komunjer-ng.R
## --------------------------------------------------------------------------
## D37: Komunjer & Ng (2011, Econometrica) dynamic identification rank check.
##
## State space in KN's timing (the lagged convention shared with D24):
##   s_t = A(theta) s_{t-1} + B(theta) eps_t,   eps_t ~ WN(0, Sigma(theta))
##   y_t = C(theta) s_{t-1} + D(theta) eps_t
##
## Two cases (KN Propositions 2-NS and 3-S):
##
##  * n_e >= n_y ("nonsingular"): the spectral density of y determines the
##    minimal innovations representation (A, K, C, Sigma_a) up to a
##    similarity transform only. Lambda = (vec A, vec K, vec C, vech Sigma_a),
##    Delta_T columns (E A - A E, E K, -C E, 0), and identification holds iff
##    rank [Delta_Lambda, Delta_T] = n_theta + n_x^2.  When n_e = n_y, D is
##    invertible and A - B D^{-1} C is stable (KN's Assumption 4-NS) the
##    innovations form is exactly K = B D^{-1}, Sigma_a = D Sigma D' (KN's
##    own construction); otherwise (n_e > n_y, D singular, or a
##    non-invertible system) it is the steady-state Kalman filter, which is
##    the same object and always minimum-phase.
##
##  * n_e < n_y ("singular", stochastic singularity): under minimality and
##    left-invertibility (KN Assumption 5-S), (A, B, C, D, Sigma) is
##    determined up to similarity T AND the full shock transformation U in
##    GL(n_e): (T A T^-1, T B U, C T^-1, D U, U^-1 Sigma U^-T).
##    Lambda = (vec A, vec B, vec C, vec D, vech Sigma), Delta_U has n_e^2
##    columns (0, B F, 0, D F, vech(-F Sigma - Sigma F')), and
##    identification holds iff rank Delta = n_theta + n_x^2 + n_e^2.
##
## Delta_Lambda is a central finite difference; its rank is judged with the
## shared FD-noise-aware helper .ident_equilibrated_rank() (Jacobian at h and
## 2h), with the analytic Delta_T / Delta_U blocks identical in both.
## --------------------------------------------------------------------------

#' D37. Dynamic identification rank check (Komunjer & Ng 2011)
#'
#' Computes the Komunjer & Ng (2011) rank condition for local identification
#' of \code{theta} from the spectral density (all autocovariances) of the
#' observables, using the state space
#' \eqn{s_t = A s_{t-1} + B\epsilon_t}, \eqn{y_t = C s_{t-1} + D\epsilon_t},
#' \eqn{\epsilon_t \sim (0, \Sigma)} returned by \code{model_solve_fn}.
#'
#' \strong{Nonsingular case} (\eqn{n_e \ge n_y}): the check uses the
#' innovations representation \eqn{(A, K, C, \Sigma_a)} (KN Proposition
#' 2-NS). \eqn{\Delta = [\Delta_\Lambda\ \Delta_T]} with
#' \eqn{\Delta_\Lambda = \partial(vec A, vec K, vec C, vech \Sigma_a)/\partial\theta'}
#' and similarity columns \eqn{(EA - AE, EK, -CE, 0)}; the required rank is
#' \eqn{n_\theta + n_x^2}. With \eqn{n_e = n_y}, \eqn{D} invertible and
#' \eqn{A - BD^{-1}C} stable, \eqn{K = BD^{-1}} and
#' \eqn{\Sigma_a = D\Sigma D'} exactly; otherwise \eqn{(K, \Sigma_a)} come
#' from the steady-state Kalman filter.
#'
#' \strong{Singular case} (\eqn{n_e < n_y}): \eqn{\Lambda = (vec A, vec B,
#' vec C, vec D, vech \Sigma)}, \eqn{\Delta = [\Delta_\Lambda\ \Delta_T\
#' \Delta_U]} with shock-transformation columns
#' \eqn{(0, BF, 0, DF, vech(-F\Sigma - \Sigma F'))} for every
#' \eqn{F = E_{kl}} (\eqn{n_e^2} of them); the required rank is
#' \eqn{n_\theta + n_x^2 + n_e^2} (KN Proposition 3-S). Left-invertibility
#' (no finite zeros of \eqn{[zI - A, B; -C, D]}, KN Assumption 5-S) is
#' checked when \eqn{D} has full column rank.
#'
#' The rank condition applies to minimal, stationary representations.
#' Minimality is checked with the PBH test on \eqn{(A, K, C)} (nonsingular)
#' or \eqn{(A, B, C)} (singular); a non-minimal or non-stationary
#' representation, a failed 5-S check, or a numerical failure of the
#' re-solve gives \code{pass = NA} (INFO) with the reason in the summary.
#'
#' \strong{Shock covariance.} \code{$Sigma_e} on the closure's result is
#' used at every \code{theta} (so shock-std parameters are identified
#' through it); otherwise \code{Sigma_e}; otherwise the identity, with a
#' warning.
#'
#' @param dr             Unused (kept for the orchestrator's call).
#' @param model          Unused (kept for the orchestrator's call).
#' @param theta          Numeric vector -- parameter values at the
#'                        evaluation point.
#' @param param_names    Optional character vector of parameter names
#'                        (length = \code{length(theta)}).
#' @param model_solve_fn Function: \code{theta -> list(A, B, C, D[, Sigma_e])}
#'                        in the timing above (the same contract as D24's
#'                        \code{abcd_solve_fn}). A failed solve should be
#'                        signalled by returning non-finite matrices.
#'                        \code{NULL} degrades to \code{pass = NA}.
#' @param Sigma_e        Optional \eqn{n_e \times n_e} shock covariance used
#'                        when the closure's result carries no
#'                        \code{$Sigma_e}.
#' @param obs_vars       Optional observable names (summary only).
#' @param eps            Relative finite-difference step (default
#'                        \code{1e-5}; the step for \code{theta[j]} is
#'                        \code{eps * max(1, |theta[j]|)}, and the
#'                        Jacobian is also taken at twice that step to
#'                        size the rank tolerance).
#' @param meta           Plot-provenance metadata (see \code{diag_meta()}).
#' @param ...            Reserved / ignored.
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
                            eps            = 1e-5,
                            meta           = NULL,
                            ...) {

  info <- function(msg, result = NULL) {
    .make_result(result = result, pass = NA, plots = list(),
                 summary = paste("D37 Komunjer-Ng:", msg),
                 llm_summary = paste0("D37 | Komunjer-Ng Dynamic Identification | INFO\n  ", msg))
  }

  # ---------------------------------------------------------------
  # 0. Inputs
  # ---------------------------------------------------------------
  if (is.null(theta)) return(info("theta not provided. Skipped."))
  if (!is.numeric(theta) || length(theta) == 0L || !all(is.finite(theta)))
    .dynhr_abort("d37_komunjer_ng: `theta` must be a non-empty finite numeric vector.")
  n_theta <- length(theta)
  if (is.null(param_names)) param_names <- names(theta)
  if (length(param_names) != n_theta) {
    if (!is.null(param_names))
      .dynhr_warn(sprintf(
        "d37: param_names length (%d) != length(theta) (%d); using names(theta).",
        length(param_names), n_theta))
    param_names <- names(theta)
  }
  if (is.null(param_names)) param_names <- paste0("theta_", seq_len(n_theta))
  param_names <- as.character(param_names)

  if (!is.function(model_solve_fn)) {
    return(info(paste(
      "the rank check differentiates the state space with respect to theta,",
      "which requires a re-solve closure.",
      "Provide `model_solve_fn = function(theta) list(A,B,C,D[,Sigma_e])`",
      "(e.g. wrapping solve_perturbation + build_dsge_state_space). Skipped.")))
  }
  if (!is.numeric(eps) || length(eps) != 1L || !is.finite(eps) || eps <= 0)
    .dynhr_abort("d37_komunjer_ng: `eps` must be a positive number.")

  # ---------------------------------------------------------------
  # 1. Baseline state space and Sigma_e
  # ---------------------------------------------------------------
  ss0 <- model_solve_fn(theta)
  if (!is.list(ss0) || !all(c("A", "B", "C", "D") %in% names(ss0)))
    .dynhr_abort("d37_komunjer_ng: model_solve_fn(theta) must return list(A, B, C, D[, Sigma_e]).")
  A0 <- as.matrix(ss0$A); B0 <- as.matrix(ss0$B)
  C0 <- as.matrix(ss0$C); D0 <- as.matrix(ss0$D)
  n_x <- nrow(A0); n_e <- ncol(B0); n_y <- nrow(C0)
  if (ncol(A0) != n_x || nrow(B0) != n_x || ncol(C0) != n_x ||
      nrow(D0) != n_y || ncol(D0) != n_e || n_x < 1L || n_e < 1L || n_y < 1L)
    .dynhr_abort(sprintf(paste(
      "d37_komunjer_ng: non-conformable state space (A %dx%d, B %dx%d,",
      "C %dx%d, D %dx%d); need A n_x x n_x, B n_x x n_e, C n_y x n_x, D n_y x n_e."),
      nrow(A0), ncol(A0), nrow(B0), ncol(B0), nrow(C0), ncol(C0), nrow(D0), ncol(D0)))

  check_sigma <- function(S, what) {
    S <- as.matrix(S)
    if (!identical(dim(S), c(n_e, n_e)))
      .dynhr_abort(sprintf("d37_komunjer_ng: %s is %dx%d but the state space has n_e = %d shocks.",
                           what, nrow(S), ncol(S), n_e))
    S
  }
  sigma_fixed <- if (!is.null(Sigma_e)) check_sigma(Sigma_e, "`Sigma_e`") else NULL
  if (!is.null(ss0$Sigma_e)) {
    sigma_source <- "model_solve_fn(theta)$Sigma_e (re-evaluated at each theta)"
    S0 <- check_sigma(ss0$Sigma_e, "model_solve_fn(theta)$Sigma_e")
    if (!is.null(sigma_fixed) && all(is.finite(S0)) && all(is.finite(sigma_fixed)) &&
        max(abs(S0 - sigma_fixed)) > 1e-8 * max(abs(S0)))
      .dynhr_warn(paste(
        "d37: `Sigma_e` differs from model_solve_fn(theta)$Sigma_e;",
        "the closure's (theta-dependent) Sigma_e is used."))
  } else if (!is.null(sigma_fixed)) {
    sigma_source <- "Sigma_e argument (held fixed across theta)"
    S0 <- sigma_fixed
  } else {
    .dynhr_warn(paste(
      "d37: no shock covariance (neither model_solve_fn(theta)$Sigma_e nor",
      "`Sigma_e`); assuming Sigma_e = I, which is only right for unit-variance shocks."))
    sigma_source <- "identity (DEFAULT -- no Sigma_e supplied)"
    sigma_fixed <- diag(n_e)
    S0 <- sigma_fixed
  }
  sigma_at <- function(ss) {
    S <- if (!is.null(ss$Sigma_e)) check_sigma(ss$Sigma_e, "model_solve_fn(theta)$Sigma_e")
         else sigma_fixed %||% diag(n_e)
    (S + t(S)) / 2
  }
  S0 <- (S0 + t(S0)) / 2

  base_result <- list(A = A0, B = B0, C = C0, D = D0, Sigma_e = S0,
                      sigma_source = sigma_source,
                      n_x = n_x, n_y = n_y, n_e = n_e, n_theta = n_theta)
  if (!all(is.finite(c(A0, B0, C0, D0, S0))))
    return(info("model_solve_fn(theta) returned non-finite matrices. Skipped.",
                base_result))

  S0_eig <- eigen(S0, symmetric = TRUE, only.values = TRUE)$values
  ## relative to Sigma_e's own scale (no max(1, .) floor)
  if (min(S0_eig) < -1e-10 * max(abs(S0_eig)))
    return(info(sprintf("Sigma_e is not positive semi-definite (smallest eigenvalue %.3g). Skipped.",
                        min(S0_eig)), base_result))

  rho_A <- max(Mod(eigen(A0, only.values = TRUE)$values))
  base_result$spectral_radius <- rho_A
  if (rho_A >= 1 - 1e-8)
    return(info(sprintf(paste(
      "A has spectral radius %.6f >= 1. The KN condition concerns the spectral",
      "density of a stationary process. INFO only."), rho_A), base_result))

  singular_case <- n_e < n_y
  case_label <- if (singular_case) "singular (n_e < n_y)" else
    if (n_e == n_y) "nonsingular (n_e = n_y)" else "nonsingular (n_e > n_y)"
  base_result$case <- case_label

  # ---------------------------------------------------------------
  # 2. Case-specific Lambda map
  # ---------------------------------------------------------------
  vech <- function(M) M[lower.tri(M, diag = TRUE)]
  representation <- NULL
  left_invertible <- NA
  notes <- character(0)

  if (singular_case) {
    representation <- "structural (A, B, C, D, Sigma)"
    if (min(S0_eig) <= 1e-12 * max(abs(S0_eig)))
      return(info(paste(
        "Sigma_e is singular. The singular-case KN condition (Proposition 3-S)",
        "needs a nonsingular shock covariance: drop the degenerate shocks. INFO only."),
        base_result))
    lambda_at <- function(ss) {
      S <- sigma_at(ss)
      c(as.numeric(ss$A), as.numeric(ss$B), as.numeric(ss$C),
        as.numeric(ss$D), vech(S))
    }
    left_invertible <- .d37_left_invertible(A0, B0, C0, D0)
    if (is.na(left_invertible))
      notes <- c(notes, paste(
        "Left-invertibility (KN 5-S) not verified: D does not have full column",
        "rank."))
    lam_names <- c(sprintf("A[%d,%d]", row(A0), col(A0)),
                   sprintf("B[%d,%d]", row(B0), col(B0)),
                   sprintf("C[%d,%d]", row(C0), col(C0)),
                   sprintf("D[%d,%d]", row(D0), col(D0)),
                   sprintf("Sigma[%d,%d]", row(S0)[lower.tri(S0, diag = TRUE)],
                           col(S0)[lower.tri(S0, diag = TRUE)]))
    K0 <- NULL; Sa0 <- NULL
  } else {
    # Nonsingular: innovations form. Choose KN's closed form when it is valid
    # at the baseline (decided once, so every FD evaluation uses one map).
    closed_form <- FALSE
    if (n_e == n_y && rcond(D0) > 1e-10) {
      K_try <- B0 %*% solve(D0)
      rho_min <- max(Mod(eigen(A0 - K_try %*% C0, only.values = TRUE)$values))
      closed_form <- rho_min < 1 - 1e-6
    }
    # Stochastic singularity guard: Var(y) must be positive definite.
    V0 <- solve_lyapunov(A0, B0 %*% S0 %*% t(B0))
    Om <- C0 %*% V0 %*% t(C0) + D0 %*% S0 %*% t(D0)
    om_eig <- eigen((Om + t(Om)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (min(om_eig) <= 1e-12 * max(abs(om_eig), 1e-300))
      return(info(paste(
        "Var(y_t) is singular (stochastic singularity), so there is no",
        "nonsingular innovations representation. INFO only."), base_result))

    innov_at <- function(ss) {
      A <- as.matrix(ss$A); B <- as.matrix(ss$B)
      C <- as.matrix(ss$C); D <- as.matrix(ss$D); S <- sigma_at(ss)
      if (closed_form)
        return(list(K = B %*% solve(D), Sa = D %*% S %*% t(D), ok = TRUE))
      Q <- B %*% S %*% t(B)
      H <- D %*% S %*% t(D)
      G <- B %*% S %*% t(D)
      ## Stop RELATIVE to the size of P and Q (kalman_filter()'s DARE rule):
      ## 1e-14 x max(1, |Q|, |H|) was absolute below unit scale, so
      ## with every shock std x 1e-4 the fixed point stopped ~1e-6-relative
      ## early and K / Sigma_a entered the finite differences with that error.
      dare <- .solve_dare(A, C, Q, H, G, tol = .LYAP_TOL,
                          max_iter = 20000L, relative = TRUE)
      list(K = dare$K, Sa = dare$F, ok = isTRUE(dare$converged))
    }
    in0 <- innov_at(ss0)
    if (!in0$ok)
      return(info("the steady-state Kalman filter did not converge at theta. INFO only.",
                  base_result))
    K0 <- in0$K; Sa0 <- (in0$Sa + t(in0$Sa)) / 2
    representation <- if (closed_form)
      "innovations, K = B D^-1, Sigma_a = D Sigma D' (KN 4-NS holds)" else
      "innovations, steady-state Kalman filter"
    if (n_e == n_y && !closed_form)
      notes <- c(notes, paste(
        "D is singular or A - B D^-1 C is not stable (KN Assumption 4-NS);",
        "the invertible (Kalman) innovations representation is used instead."))
    lambda_at <- function(ss) {
      inn <- innov_at(ss)
      if (!inn$ok) return(NULL)
      c(as.numeric(ss$A), as.numeric(inn$K), as.numeric(ss$C),
        vech((inn$Sa + t(inn$Sa)) / 2))
    }
    lam_names <- c(sprintf("A[%d,%d]", row(A0), col(A0)),
                   sprintf("K[%d,%d]", row(K0), col(K0)),
                   sprintf("C[%d,%d]", row(C0), col(C0)),
                   sprintf("Sigma_a[%d,%d]", row(Sa0)[lower.tri(Sa0, diag = TRUE)],
                           col(Sa0)[lower.tri(Sa0, diag = TRUE)]))
  }
  base_result$representation <- representation
  base_result$K <- K0
  base_result$Sigma_a <- Sa0
  base_result$left_invertible <- left_invertible

  # ---------------------------------------------------------------
  # 3. Minimality (PBH test)
  # ---------------------------------------------------------------
  M_in <- if (singular_case) B0 else K0
  ctrb_rank <- .d37_pbh_rank(A0, M_in, controllability = TRUE)
  obsv_rank <- .d37_pbh_rank(A0, C0, controllability = FALSE)
  is_minimal <- ctrb_rank == n_x && obsv_rank == n_x
  base_result$minimal <- is_minimal
  base_result$ctrb_rank <- ctrb_rank
  base_result$obsv_rank <- obsv_rank

  # ---------------------------------------------------------------
  # 4. Delta_Lambda at steps h and 2h
  # ---------------------------------------------------------------
  lambda_theta <- function(th) {
    ss <- model_solve_fn(th)
    if (!is.list(ss) || !all(c("A", "B", "C", "D") %in% names(ss))) return(NULL)
    ss$A <- as.matrix(ss$A); ss$B <- as.matrix(ss$B)
    ss$C <- as.matrix(ss$C); ss$D <- as.matrix(ss$D)
    if (!identical(dim(ss$A), dim(A0)) || !identical(dim(ss$B), dim(B0)) ||
        !identical(dim(ss$C), dim(C0)) || !identical(dim(ss$D), dim(D0)))
      return(NULL)
    if (!all(is.finite(c(ss$A, ss$B, ss$C, ss$D, sigma_at(ss))))) return(NULL)
    if (!singular_case &&
        max(Mod(eigen(ss$A, only.values = TRUE)$values)) >= 1) return(NULL)
    out <- lambda_at(ss)
    if (is.null(out) || !all(is.finite(out))) NULL else out
  }
  n_lam <- length(lam_names)
  fd_jac <- function(mult) {
    J <- matrix(NA_real_, n_lam, n_theta)
    for (j in seq_len(n_theta)) {
      h <- mult * eps * max(1, abs(theta[j]))
      tp <- theta; tp[j] <- theta[j] + h
      tm <- theta; tm[j] <- theta[j] - h
      fp <- lambda_theta(tp)
      fm <- lambda_theta(tm)
      if (is.null(fp) || is.null(fm)) return(NULL)
      J[, j] <- (fp - fm) / (2 * h)
    }
    dimnames(J) <- list(lam_names, param_names)
    J
  }
  Delta_Lambda <- fd_jac(1)
  Delta_Lambda2 <- if (!is.null(Delta_Lambda)) fd_jac(2)
  if (is.null(Delta_Lambda) || is.null(Delta_Lambda2))
    return(info(paste(
      "model_solve_fn failed (non-finite, non-conformable or non-stationary",
      "result) at a finite-difference point, so Delta_Lambda is undetermined.",
      "This is a numerical failure, not an identification finding."),
      base_result))

  # ---------------------------------------------------------------
  # 5. Delta_T (similarity) and Delta_U (shock transformation, singular)
  # ---------------------------------------------------------------
  M_sim <- if (singular_case) B0 else K0
  n_tail <- n_lam - n_x * n_x - length(M_sim) - n_y * n_x   # D + Sigma or Sigma_a
  Delta_T <- matrix(0, n_lam, n_x * n_x,
                    dimnames = list(lam_names, character(n_x * n_x)))
  k <- 0L
  for (jj in seq_len(n_x)) for (ii in seq_len(n_x)) {
    k <- k + 1L
    E <- matrix(0, n_x, n_x); E[ii, jj] <- 1
    Delta_T[, k] <- c(as.numeric(E %*% A0 - A0 %*% E), as.numeric(E %*% M_sim),
                      as.numeric(-C0 %*% E), rep(0, n_tail))
    colnames(Delta_T)[k] <- sprintf("T[%d,%d]", ii, jj)
  }
  Delta_U <- matrix(0, n_lam, 0, dimnames = list(lam_names, NULL))
  if (singular_case) {
    Delta_U <- matrix(0, n_lam, n_e * n_e,
                      dimnames = list(lam_names, character(n_e * n_e)))
    k <- 0L
    for (ll in seq_len(n_e)) for (kk in seq_len(n_e)) {
      k <- k + 1L
      Fm <- matrix(0, n_e, n_e); Fm[kk, ll] <- 1
      Delta_U[, k] <- c(rep(0, n_x * n_x), as.numeric(B0 %*% Fm),
                        rep(0, n_y * n_x), as.numeric(D0 %*% Fm),
                        vech(-(Fm %*% S0 + S0 %*% t(Fm))))
      colnames(Delta_U)[k] <- sprintf("U[%d,%d]", kk, ll)
    }
  }

  Delta  <- cbind(Delta_Lambda,  Delta_T, Delta_U)
  Delta2 <- cbind(Delta_Lambda2, Delta_T, Delta_U)
  required_rank <- ncol(Delta)
  order_ok <- n_lam >= required_rank

  rk <- .ident_equilibrated_rank(Delta, Delta2)
  singular_values <- rk$singular_values
  rank_Delta <- rk$rank
  full_rank <- rank_Delta == required_rank
  n_deficient <- required_rank - rank_Delta

  assumption_fail <- !is_minimal || identical(left_invertible, FALSE)
  pass <- if (assumption_fail) NA else full_rank

  # ---------------------------------------------------------------
  # 6. Deficient directions: theta-block loadings of the null space
  # ---------------------------------------------------------------
  deficient_loadings <- NULL
  if (n_deficient > 0L) {
    dirs <- rk$unid_dirs
    rows <- lapply(dirs, function(d) {
      l <- rk$V[seq_len(n_theta), d]
      ord <- order(abs(l), decreasing = TRUE)
      data.frame(parameter = param_names[ord], loading = unname(l[ord]),
                 abs_loading = unname(abs(l[ord])), direction = d,
                 singular_value = unname(singular_values[d]),
                 stringsAsFactors = FALSE)
    })
    deficient_loadings <- do.call(rbind, rows)
    rownames(deficient_loadings) <- NULL
  }
  # Parameters that dominate some null direction (theta block only).
  unid_params <- if (n_deficient > 0L) {
    unique(unlist(lapply(rk$unid_dirs, function(d) {
      v <- abs(rk$V[seq_len(n_theta), d])
      if (max(v) <= 1e-8) character(0) else param_names[v >= 0.5 * max(v)]
    })))
  } else character(0)

  # ---------------------------------------------------------------
  # 7. Plot: equilibrated singular values vs the rank tolerance
  # ---------------------------------------------------------------
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plots$singular_values <- .apply_meta(
      .d37_sv_plot(rk, colnames(Delta), n_theta, rank_Delta, required_rank,
                   case_label, pass),
      meta)
  }

  # ---------------------------------------------------------------
  # 8. Text
  # ---------------------------------------------------------------
  req_txt <- if (singular_case)
    sprintf("n_theta + n_x^2 + n_e^2 = %d + %d + %d", n_theta, n_x^2, n_e^2) else
    sprintf("n_theta + n_x^2 = %d + %d", n_theta, n_x^2)
  tol_txt <- sprintf("rank tolerance %.1e relative (%s%s)",
                     rk$tol / max(rk$sv_max, 1e-300), rk$tol_source,
                     if (is.finite(rk$fd_noise)) sprintf(", FD error %.1e", rk$fd_noise) else "")
  why_na <- c(
    if (!is_minimal) sprintf(paste(
      "Representation is NON-MINIMAL (PBH controllability rank %d/%d,",
      "observability rank %d/%d): the KN condition applies to minimal",
      "realisations, so the rank below is not a pass/fail verdict."),
      ctrb_rank, n_x, obsv_rank, n_x),
    if (identical(left_invertible, FALSE)) paste(
      "The system has a finite zero (KN Assumption 5-S fails), so the",
      "singular-case rank condition is not applicable."))
  rank_status <- if (is.na(pass)) "INFO -- KN assumptions fail, rank not a verdict"
    else if (pass) "PASS -- locally identified"
    else sprintf("FAIL -- rank deficient by %d direction(s)", n_deficient)
  summary_text <- paste0(
    sprintf("D37 Komunjer-Ng (2011) dynamic identification: %s. ", rank_status),
    sprintf("Case = %s; %s. n_x=%d, n_y=%d, n_e=%d, n_theta=%d. ",
            case_label, representation, n_x, n_y, n_e, n_theta),
    sprintf("rank(Delta) = %d, required = %d (%s); %s.",
            rank_Delta, required_rank, req_txt, tol_txt),
    if (!order_ok) sprintf(" Order condition fails: Lambda has only %d entries.", n_lam) else "",
    sprintf(" Sigma_e: %s.", sigma_source),
    if (!is.null(obs_vars)) sprintf(" obs = [%s].", paste(obs_vars, collapse = ", ")) else "",
    if (length(why_na)) paste0(" ", paste(why_na, collapse = " ")) else "",
    if (length(notes)) paste0(" Note: ", paste(notes, collapse = " ")) else "",
    if (length(unid_params))
      sprintf(" Deficient direction(s) load most heavily on: %s.",
              paste(head(unid_params, 5), collapse = ", ")) else
    if (n_deficient > 0L)
      " The deficient direction(s) involve only the similarity/shock-transformation block (the representation is not minimal / not canonical)." else ""
  )

  badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
  llm <- paste(c(
    sprintf("D37 | Komunjer-Ng Dynamic Identification | %s", badge),
    sprintf("  case=%s representation=%s n_x=%d n_y=%d n_e=%d n_theta=%d",
            case_label, representation, n_x, n_y, n_e, n_theta),
    sprintf("  rank(Delta)=%d required=%d minimal=%s (ctrb_rank=%d, obsv_rank=%d) %s",
            rank_Delta, required_rank, is_minimal, ctrb_rank, obsv_rank, tol_txt),
    sprintf("  smallest_sv (equilibrated): %s",
            paste(sprintf("%.3e", head(sort(singular_values), 5)), collapse = ", ")),
    if (length(unid_params))
      sprintf("  unidentified_params: %s", paste(unid_params, collapse = ", ")),
    sprintf("  action: %s",
            if (is.na(pass))
              "KN assumptions fail (non-minimal or 5-S) -- reduce to a minimal representation (drop redundant states / shocks) and re-run."
            else if (pass)
              "Model is locally identified from the spectral density (Komunjer & Ng 2011)."
            else
              sprintf("Rank deficient (rank=%d, required=%d). %s not locally identified from the autocovariances. Consider reparameterisation, calibration, or additional observables.",
                      rank_Delta, required_rank,
                      if (length(unid_params)) paste(head(unid_params, 3), collapse = ", ") else "Some parameter combination is"))
  ), collapse = "\n")

  .make_result(
    result = c(base_result, list(
      Delta_Lambda = Delta_Lambda,
      Delta_T = Delta_T,
      Delta_U = Delta_U,
      Delta = Delta,
      Delta_equilibrated = rk$Je,
      singular_values = singular_values,
      sv_class = stats::setNames(rk$sv_class, names(singular_values)),
      svd = rk$svd,
      rank = rank_Delta,
      required_rank = required_rank,
      rank_tolerance = rk$tol / max(rk$sv_max, 1e-300),
      tolerance_source = rk$tol_source,
      fd_error = rk$fd_noise,
      null_space = rk$null_space,
      unidentified_params = unid_params,
      deficient_loadings = deficient_loadings
    )),
    pass = pass,
    plots = plots,
    summary = summary_text,
    llm_summary = llm
  )
}


#' PBH rank of a pair: n minus the dimension of the uncontrollable
#' (unobservable) part, from rank [A - lambda I, M] ([A - lambda I; M]) at
#' each distinct eigenvalue. Both blocks are normalised to unit 2-norm first
#' (the test is invariant to that), unlike a Krylov matrix whose high powers
#' of A underflow for n_x > ~15 and falsely report non-minimality.
#' @noRd
.d37_pbh_rank <- function(A, M, controllability = TRUE, tol_rel = 1e-8) {
  n <- nrow(A)
  a_n <- max(svd(A, nu = 0, nv = 0)$d)
  m_n <- max(svd(M, nu = 0, nv = 0)$d)
  if (m_n <= 0) return(0L)
  if (a_n <= 0) a_n <- 1
  As <- A / a_n
  Ms <- M / m_n
  ev <- eigen(As, only.values = TRUE)$values
  distinct <- ev[1]
  for (lam in ev[-1]) if (all(Mod(lam - distinct) > 1e-7)) distinct <- c(distinct, lam)
  deficit <- 0L
  for (lam in distinct) {
    Al <- As - lam * diag(n)
    P <- if (controllability) cbind(Al, Ms) else rbind(Al, Ms)
    d <- svd(P, nu = 0, nv = 0)$d
    deficit <- deficit + (n - sum(d > tol_rel))
  }
  as.integer(max(0L, n - deficit))
}


#' KN Assumption 5-S: rank [zI - A, B; -C, D] = n_x + n_e for all finite z.
#' Candidate zeros are the eigenvalues of A - B (D'D)^-1 D'C (the zeros of
#' the square pencil diag(I, D') P(z), a superset of the true zeros); each is
#' then tested on P(z) itself. NA when D lacks full column rank.
#' @noRd
.d37_left_invertible <- function(A, B, C, D, tol_rel = 1e-8) {
  n_x <- nrow(A); n_e <- ncol(D)
  sv_D <- svd(D, nu = 0, nv = 0)$d
  if (length(sv_D) < n_e || min(sv_D) <= 1e-10 * max(sv_D, 1e-300)) return(NA)
  cand <- eigen(A - B %*% solve(crossprod(D), crossprod(D, C)),
                only.values = TRUE)$values
  sys_n <- max(svd(rbind(cbind(A, B), cbind(C, D)), nu = 0, nv = 0)$d)
  for (z in cand) {
    P <- rbind(cbind(z * diag(n_x) - A, B), cbind(-C, D))
    if (min(svd(P, nu = 0, nv = 0)$d) <= tol_rel * sys_n * max(1, Mod(z)))
      return(FALSE)
  }
  TRUE
}


#' Singular-value lollipop for D37 (equilibrated Delta). Each index is
#' labelled with the column that loads most on its right singular vector
#' (a parameter, or a similarity T[i,j] / shock-transformation U[k,l] block).
#' @noRd
.d37_sv_plot <- function(rk, col_names, n_theta, rank_Delta, required_rank,
                         case_label, pass) {
  sv <- rk$singular_values
  n <- length(sv)
  pos <- sv[sv > 0]
  cand <- c(pos, rk$tol)
  cand <- cand[cand > 0]
  floor_val <- 10^(floor(log10(if (length(cand)) min(cand) else 1e-16)) - 1)
  lead <- vapply(seq_len(n), function(k) {
    v <- abs(rk$V[, k])
    if (max(v) <= 0) "" else col_names[which.max(v)]
  }, character(1))
  cls <- ifelse(sv > rk$tol, "Above tolerance (identified direction)",
                "At/below tolerance (unidentified direction)")
  df <- data.frame(
    index = factor(seq_len(n), labels = sprintf(if (n > 12L) "%d: %s" else "%d\n%s",
                                                seq_len(n), lead)),
    value = pmax(sv, floor_val),
    class = factor(cls, levels = c("Above tolerance (identified direction)",
                                   "At/below tolerance (unidentified direction)")),
    zero = sv <= 0
  )
  ref <- data.frame(y = rk$tol, what = sprintf(
    "rank tolerance = %.1e x largest (%s)", rk$tol / max(rk$sv_max, 1e-300), rk$tol_source))
  badge <- if (is.na(pass)) "INFO (KN assumptions fail)" else if (pass) "PASS" else "FAIL"
  p <- ggplot2::ggplot(df, ggplot2::aes(x = index, y = value, colour = class)) +
    ggplot2::geom_segment(ggplot2::aes(xend = index, y = floor_val, yend = value),
                          linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(shape = zero), size = 2.6) +
    ggplot2::geom_hline(data = ref, ggplot2::aes(yintercept = y, linetype = what),
                        colour = dynhr_colours$grey, linewidth = 0.5) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(c(dynhr_colours$mid_blue, dynhr_colours$red), levels(df$class)),
      name = NULL) +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), breaks = "TRUE",
                                labels = c(`TRUE` = "exactly 0 (drawn at floor)"),
                                name = NULL) +
    ggplot2::scale_linetype_manual(values = "dashed", name = NULL) +
    ggplot2::scale_y_log10() +
    theme_dynhr_diagnostic() +
    ggplot2::theme(legend.box = "vertical") +
    ggplot2::labs(
      title = sprintf("D37: Komunjer-Ng rank condition -- %s", badge),
      subtitle = sprintf(paste0(
        "rank(Delta) = %d of %d required; %s.\n",
        "Columns: %d parameters + similarity T%s. Label: column loading most on each direction"),
        rank_Delta, required_rank, case_label, n_theta,
        if (grepl("^singular", case_label)) " + shock transformation U" else ""),
      x = "Singular value (index / dominant column)",
      y = "Singular value (log)"
    )
  if (n > 12L)
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                                hjust = 1, size = 8))
  p
}
