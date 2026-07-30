## R/solve-perturbation-order4-5-sigma.R
## --------------------------------------------------------------------------
## Sigma-correction (augmented-state) terms for fourth- and fifth-order
## perturbation solutions.
##
## Extends Mutschler (2022) / Levintal (2017) sigma-cross methodology to
## orders 4 and 5.  For Gaussian shocks (the default), only even powers of
## σ appear.  The non-zero sigma-correction blocks are:
##
##   Order 2:  ghss              (g_0,2   — σ² at SS)
##   Order 3:  ghxss, ghuss     (g_1,2   — x·σ² and u·σ² cross terms)
##   Order 4:  ghxxss            (g_2,2   — x⊗x·σ² cross)
##             ghss2             (g_0,4   — σ⁴ at SS)
##   Order 5:  ghxxxss           (g_3,2   — x⊗x⊗x·σ² cross)
##             ghxss2            (g_1,4   — x·σ⁴ cross)
##
## The key computational pattern (same as order-3 sigma-cross):
##   1. Build compound-derivative matrices involving σ (T_ss, W_σ, etc.)
##   2. Contract with H2 (Hessian), H3 (3rd deriv), and F_k (higher derivs)
##   3. Solve Sylvester equations for mixed x-σ terms
##   4. Solve linear systems for pure σ terms
##
## References:
##   Mutschler (2022), perturbation_solver_nonsymmetric_order3.m
##   Levintal (2017), JEDC 80, 1-16
##   Andreasen, Fernandez-Villaverde & Rubio-Ramirez (2018), ReStud 85(1)
## --------------------------------------------------------------------------


# =====================================================================
# FD-forcing infrastructure for sigma-correction terms (TEST CROSS-CHECK ONLY)
# =====================================================================
#
# The sigma-correction term g_{x^k σ^{2m}} satisfies the SAME linear /
# Sylvester operator as the deterministic k-th state term:
#
#     A_L · g + fp · g · hx^{⊗k} = -Φ
#
# where Φ is the mixed derivative  ∂^k_x ∂^{2m}_σ  of the *expected* dynamic
# residual, evaluated at the deterministic point (x = 0, σ = 0), using the
# policy that already contains all LOWER-order terms (deterministic and
# sigma).  The expectation is over the future shock u' ~ N(0, σ²·Σ_e),
# computed by Gauss–Hermite quadrature.
#
# NOTE: the PRODUCTION forcing is now the analytic moment-based assembler in
# solve-perturbation-sigma-analytic.R (.build_phi_sigma_full /
# .extract_sigma_block2) — exact for Gaussian shocks, ~1000x faster, and
# machine-precision accurate where this FD path was 4%-200%+ off on strongly
# nonlinear models (see test-sigma-exp-golden.R).  `.build_phi_fd_sigma` and
# `.gh_tensor` below are retained ONLY as an independent second method for the
# test suite (test-sigma-order45.R cross-checks ghss/ghxss); they are no longer
# on any user-facing code path.  They are calibrated against the
# Dynare-validated ghss (k=0,m=1) and ghxss (k=1,m=1).

#' Tensor-product Gauss–Hermite nodes/weights for E[f(ε)], ε ~ N(0, Σ).
#' Returns a list with `nodes` (n_pts × n_u matrix of ε-values, BEFORE the σ
#' scaling) and `wts` (length n_pts, summing to 1).
#' @noRd
.gh_tensor <- function(Sigma_e, n_u, n_pt = 5L) {
  if (n_u == 0L) return(list(nodes = matrix(0, 1L, 0L), wts = 1))
  # 5-point physicists' Gauss–Hermite (weight e^{-x^2})
  t5 <- c(-2.020182870456086, -0.958572464613819, 0,
           0.958572464613819,  2.020182870456086)
  w5 <- c(0.019953242059046, 0.393619323152241, 0.945308720482942,
          0.393619323152241, 0.019953242059046)
  sqrtpi <- sqrt(pi)
  # Cholesky: Σ = L Lᵀ ; ε = L · (sqrt(2)·t)
  L <- tryCatch(t(chol(Sigma_e)), error = function(e) {
    # fall back to symmetric sqrt for PSD/near-singular Σ
    es <- eigen(Sigma_e, symmetric = TRUE)
    es$vectors %*% diag(sqrt(pmax(es$values, 0)), n_u) %*% t(es$vectors)
  })
  grid_idx <- expand.grid(rep(list(seq_len(n_pt)), n_u))
  n_pts <- nrow(grid_idx)
  nodes <- matrix(0, n_pts, n_u)
  wts   <- numeric(n_pts)
  for (p in seq_len(n_pts)) {
    ii <- as.integer(grid_idx[p, ])
    z  <- sqrt(2) * t5[ii]                 # standard-normal GH points
    nodes[p, ] <- as.numeric(L %*% z)      # ε ~ N(0, Σ)
    wts[p]     <- prod(w5[ii]) / sqrtpi^n_u
  }
  list(nodes = nodes, wts = wts)
}

#' Build the FD-forcing Φ for a sigma-correction term g_{x^k σ^{2m}}.
#'
#' @param policy function(x, u, sig) -> n-vector deviation from SS, containing
#'        ALL known lower-order terms (NOT the target term).
#' @param k_x   number of state (x) derivative indices (0..3).
#' @param p_sig sigma derivative order: 2 (σ²) or 4 (σ⁴).
#' @param h     base finite-difference step.
#' @param gh    output of .gh_tensor().
#' @return Φ  (n × n_s^k_x), declaration-variable order.
#' @noRd
.build_phi_fd_sigma <- function(dyn, policy, ss_endo, params, ss_full,
                                 state_idx, endo_names, exo_names,
                                 n_s, n_u, n, res_perm, gh,
                                 k_x, p_sig, h) {
  dcm <- dyn$dyn_col_map
  zero_u <- numeric(max(n_u, 1L))
  build_dy <- function(y_lag, y_now, y_lead, u_now) {
    dy <- numeric(dyn$total_cols); keys <- character(dyn$total_cols)
    for (kc in seq_len(nrow(dcm))) {
      c <- dcm$col[kc]; nm <- dcm$name[kc]; ll <- dcm$lead_lag[kc]
      sfx <- if (ll == 0L) "__0" else if (ll > 0L) paste0("__p", ll)
             else paste0("__m", abs(ll))
      keys[c] <- paste0(nm, sfx)
      if (nm %in% exo_names) {
        dy[c] <- u_now[which(exo_names == nm)]
      } else {
        idx <- which(endo_names == nm); if (length(idx) != 1L) next
        dy[c] <- if (ll == -1L) y_lag[idx] else if (ll == 0L) y_now[idx]
                 else if (ll == 1L) y_lead[idx] else NA_real_
      }
    }
    names(dy) <- keys; dy
  }

  cache <- new.env(parent = emptyenv())
  # Expected dynamic residual at state deviation xv, future-shock scale sig.
  Eresid <- function(xv, sig) {
    key <- paste(c(round(xv, 12), round(sig, 12)), collapse = "|")
    hit <- cache[[key]]
    if (!is.null(hit)) return(hit)
    u0     <- numeric(max(n_u, 1L))
    y_now  <- ss_endo + policy(xv, u0, sig)
    x_now  <- y_now[state_idx] - ss_endo[state_idx]
    y_lag  <- ss_endo; y_lag[state_idx] <- ss_endo[state_idx] + xv
    acc <- numeric(dyn$n_eq)
    for (q in seq_along(gh$wts)) {
      eps <- if (n_u > 0L) sig * gh$nodes[q, ] else numeric(0)
      y_lead <- ss_endo + policy(x_now, eps, sig)
      r <- as.numeric(dyn$residuals_fn(
             build_dy(y_lag, y_now, y_lead, eps), params, ss_full))
      acc <- acc + gh$wts[q] * r
    }
    val <- acc[res_perm]
    assign(key, val, envir = cache)
    val
  }

  # sigma central-difference stencil
  if (p_sig == 2L) { s_off <- c(-1L, 0L, 1L);          s_w <- c(1, -2, 1) }
  else             { s_off <- c(-2L, -1L, 0L, 1L, 2L); s_w <- c(1, -4, 6, -4, 1) }

  n_cols <- if (k_x == 0L) 1L else n_s^k_x
  Phi <- matrix(0, n, n_cols)
  x_signs <- if (k_x == 0L) matrix(0, 1L, 0L)
             else as.matrix(expand.grid(rep(list(c(-1L, 1L)), k_x)))

  idx <- rep(1L, max(k_x, 1L))
  for (col in seq_len(n_cols)) {
    val <- numeric(n)
    for (sr in seq_len(nrow(x_signs))) {
      sg <- if (k_x == 0L) integer(0) else as.integer(x_signs[sr, ])
      xv0 <- numeric(n_s)
      if (k_x > 0L) for (ki in seq_len(k_x)) xv0[idx[ki]] <- xv0[idx[ki]] + sg[ki] * h
      psign <- if (k_x == 0L) 1 else prod(sg)
      for (si in seq_along(s_off)) {
        val <- val + psign * s_w[si] * Eresid(xv0, s_off[si] * h)
      }
    }
    Phi[, col] <- val / ((2 * h)^k_x * h^p_sig)
    # advance odometer (last index fastest, matching .build_phi_fd)
    if (k_x > 0L) for (ki in rev(seq_len(k_x))) {
      idx[ki] <- idx[ki] + 1L
      if (idx[ki] <= n_s) break
      idx[ki] <- 1L
    }
  }
  Phi
}

# =====================================================================
# 4th-order sigma-correction entry point
# =====================================================================

#' Add fourth-order sigma-correction terms to a DecisionRules4 object
#'
#' Computes:
#'   \code{ghxxss} (n × n_s²) — mixed 2nd state × sigma^2 correction
#'   \code{ghss2}  (n × 1)     — pure sigma^4 steady-state correction
#'
#' For Gaussian shocks the sigma^3 and sigma^1 cross terms vanish.
#'
#' @param dr4       DecisionRules4 from \code{solve_perturbation_order4()}
#' @param compiled  dynhr_compiled
#' @param ss        Named numeric steady state
#' @param params    Named numeric parameter vector
#' @param Sigma_e   Optional n_u × n_u shock covariance (default: from model)
#' @param verbose   Print progress
#' @return \code{dr4} with fields \code{ghxxss}, \code{ghss2} added
#' @export
solve_sigma_order4 <- function(dr4, compiled, ss, params,
                                Sigma_e = NULL, verbose = FALSE) {
  if (!inherits(dr4, "DecisionRules4")) {
    stop("dr4 must be a DecisionRules4 object.")
  }
  endo_names <- dr4$endo_names
  exo_names  <- dr4$exo_names
  state_idx  <- dr4$state_idx
  n          <- length(endo_names)
  n_s        <- length(state_idx)
  n_u        <- length(exo_names)

  # ---- Sigma_e ---- 
  if (is.null(Sigma_e)) {
    Sigma_e <- dr4$Sigma_e
    if (is.null(Sigma_e)) {
      stderr  <- .get_shock_stderr(compiled$model, exo_names, params)
      Sigma_e <- diag(stderr^2, n_u, n_u)
    }
  }
  if (!is.matrix(Sigma_e) || nrow(Sigma_e) != n_u || ncol(Sigma_e) != n_u) {
    stop(sprintf("Sigma_e must be n_u x n_u (n_u = %d).", n_u))
  }
  SIGMA2 <- as.numeric(Sigma_e)

  dyn <- compiled$dynamic
  sys <- extract_system_matrices(compiled, ss, params)
  has_lead <- sys$is_fwd | sys$is_mixed

  # ---- Operator matrices (the analytic forcing reads policy tensors from dr4) ----
  ghx <- dr4$ghx
  hx  <- ghx[state_idx, , drop = FALSE]

  S <- matrix(0, n, n_s)
  for (s in seq_along(state_idx)) S[state_idx[s], s] <- 1
  A_L <- sys$f_zero + sys$f_plus %*% ghx %*% t(S)
  fp  <- sys$f_plus

  # ---- Analytic moment-based forcing setup ----
  # Reuse sys mapping so forcing rows align with A_L (see eq-to-decl invariant).
  eq_to_decl <- sys$eq_to_decl %||% .build_eq_to_decl(compiled$model)
  res_perm   <- order(eq_to_decl)
  M2 <- .M2_vec(Sigma_e); M4 <- .M4_vec(Sigma_e)

  if (is.null(dr4$ghss))
    stop("dr4$ghss required; run solve_perturbation_order2() first.")

  # Full un-folded forcing tensor over w = (x, sigma, e); built once for K=4,
  # then both order-4 sigma blocks are extracted by folding the e-moments.
  Phi_full <- .build_phi_sigma_full(dyn, dr4, ss, params, state_idx,
                                    endo_names, exo_names, n_s, n_u, n,
                                    4L, res_perm)

  # ---- ghxxss : g_{x² σ²}  (A_L·X + fp·X·hx^{⊗2} = -Φ) ----
  if (verbose) cat("  Computing ghxxss (4th-order x^2-sigma^2, analytic forcing)...\n")
  phi_xxss <- .extract_sigma_block2(Phi_full, 4L, n_s, n_u, 2L, 2L, n, M2, M4)
  ghxxss <- .solve_kron_compact(A_L, fp, hx, 2L, -phi_xxss, verbose = verbose)
  rownames(ghxxss) <- endo_names
  colnames(ghxxss) <- .quad_names_internal(endo_names[state_idx],
                                           endo_names[state_idx])
  for (e in seq_len(n)) {                          # symmetrize two state indices
    Xe <- matrix(ghxxss[e, ], n_s, n_s)
    ghxxss[e, ] <- as.numeric((Xe + t(Xe)) / 2)
  }

  # ---- ghss2 : g_{σ⁴}  ((A_L + fp)·X = -Φ) ----
  if (verbose) cat("  Computing ghss2 (4th-order sigma^4 SS, analytic forcing)...\n")
  phi_ssss <- .extract_sigma_block2(Phi_full, 4L, n_s, n_u, 0L, 4L, n, M2, M4)
  ghss2 <- as.numeric(.solve_equilibrated(A_L + fp, -as.numeric(phi_ssss)))
  names(ghss2) <- endo_names

  # ---- Write back ----
  dr4$ghxxss <- ghxxss
  dr4$ghss2  <- ghss2
  dr4$sigma_order4 <- TRUE
  dr4$Sigma_e_used <- Sigma_e
  dr4
}


# =====================================================================
# 5th-order sigma-correction entry point
# =====================================================================

#' Add fifth-order sigma-correction terms to a DecisionRules5 object
#'
#' Computes:
#'   \code{ghxxxss} (n × n_s³) — mixed 3rd state × sigma^2 correction
#'   \code{ghxss2}  (n × n_s)   — mixed 1st state × sigma^4 correction
#'
#' Requires order-4 sigma terms (\code{ghxxss}, \code{ghss2}) from
#' \code{solve_sigma_order4()}.
#'
#' @param dr5       DecisionRules5 from \code{solve_perturbation_order5()}
#' @param compiled  dynhr_compiled
#' @param ss        Named numeric steady state
#' @param params    Named numeric parameter vector
#' @param Sigma_e   Optional n_u × n_u shock covariance
#' @param verbose   Print progress
#' @return \code{dr5} with \code{ghxxxss}, \code{ghxss2} added
#' @export
solve_sigma_order5 <- function(dr5, compiled, ss, params,
                                Sigma_e = NULL, verbose = FALSE) {
  if (!inherits(dr5, "DecisionRules5")) {
    stop("dr5 must be a DecisionRules5 object.")
  }
  endo_names <- dr5$endo_names
  exo_names  <- dr5$exo_names
  state_idx  <- dr5$state_idx
  n          <- length(endo_names)
  n_s        <- length(state_idx)
  n_u        <- length(exo_names)

  if (is.null(Sigma_e)) {
    Sigma_e <- dr5$Sigma_e
    if (is.null(Sigma_e)) {
      stderr  <- .get_shock_stderr(compiled$model, exo_names, params)
      Sigma_e <- diag(stderr^2, n_u, n_u)
    }
  }
  SIGMA2 <- as.numeric(Sigma_e)

  dyn <- compiled$dynamic
  sys <- extract_system_matrices(compiled, ss, params)
  has_lead <- sys$is_fwd | sys$is_mixed

  # ---- Operator matrices (the analytic forcing reads policy tensors from dr5) ----
  ghx <- dr5$ghx
  hx  <- ghx[state_idx, , drop = FALSE]

  S <- matrix(0, n, n_s)
  for (s in seq_along(state_idx)) S[state_idx[s], s] <- 1
  A_L <- sys$f_zero + sys$f_plus %*% ghx %*% t(S)
  fp  <- sys$f_plus

  # ---- Analytic moment-based forcing setup ----
  # Reuse sys mapping so forcing rows align with A_L (see eq-to-decl invariant).
  eq_to_decl <- sys$eq_to_decl %||% .build_eq_to_decl(compiled$model)
  res_perm   <- order(eq_to_decl)
  M2 <- .M2_vec(Sigma_e); M4 <- .M4_vec(Sigma_e)

  if (is.null(dr5$ghss))   stop("dr5$ghss required.")
  if (is.null(dr5$ghxss))  stop("dr5$ghxss required (3rd-order sigma-cross).")
  if (is.null(dr5$ghxxss)) stop("dr5$ghxxss required; run solve_sigma_order4() first.")
  if (is.null(dr5$ghss2))  stop("dr5$ghss2 required; run solve_sigma_order4() first.")

  # Full un-folded forcing tensor over w = (x, sigma, e); built once for K=5,
  # carrying the full order-5 deterministic policy plus all order<=4 sigma
  # corrections (ghss, ghxss, ghuss, ghxxss, ghss2).  Both order-5 sigma blocks
  # are extracted by folding the e-moments.
  Phi_full <- .build_phi_sigma_full(dyn, dr5, ss, params, state_idx,
                                    endo_names, exo_names, n_s, n_u, n,
                                    5L, res_perm)

  # ---- ghxxxss : g_{x³ σ²}  (A_L·X + fp·X·hx^{⊗3} = -Φ) ----
  if (verbose) cat("  Computing ghxxxss (5th-order x^3-sigma^2, analytic forcing)...\n")
  phi_xxxss <- .extract_sigma_block2(Phi_full, 5L, n_s, n_u, 3L, 2L, n, M2, M4)
  ghxxxss <- .solve_kron_compact(A_L, fp, hx, 3L, -phi_xxxss, verbose = verbose)
  rownames(ghxxxss) <- endo_names
  colnames(ghxxxss) <- .triple_names_internal(endo_names[state_idx])

  # ---- ghxss2 : g_{x σ⁴}  (A_L·X + fp·X·hx = -Φ) ----
  if (verbose) cat("  Computing ghxss2 (5th-order x-sigma^4, analytic forcing)...\n")
  phi_xss2 <- .extract_sigma_block2(Phi_full, 5L, n_s, n_u, 1L, 4L, n, M2, M4)
  ghxss2 <- .solve_kron_compact(A_L, fp, hx, 1L, -phi_xss2, verbose = verbose)
  rownames(ghxss2) <- endo_names
  colnames(ghxss2) <- endo_names[state_idx]

  # ---- Write back ----
  dr5$ghxxxss <- ghxxxss
  dr5$ghxss2  <- ghxss2
  dr5$sigma_order5 <- TRUE
  dr5$Sigma_e_used <- Sigma_e
  dr5
}


# =====================================================================
# Internal helpers
# =====================================================================

#' Generate column names for 2-fold outer product (quadratic names).
#' @noRd
.quad_names_internal <- function(a, b) {
  out <- character(length(a) * length(b))
  idx <- 1L
  for (j in seq_along(b))
    for (i in seq_along(a)) {
      out[idx] <- paste(a[i], b[j], sep = "__x__")
      idx <- idx + 1L
    }
  out
}

#' Generate column names for 3-fold outer product (triple names).
#' @noRd
.triple_names_internal <- function(a) {
  out <- character(length(a)^3)
  idx <- 1L
  for (k in seq_along(a))
    for (j in seq_along(a))
      for (i in seq_along(a)) {
        out[idx] <- paste(a[i], a[j], a[k], sep = "__x__")
        idx <- idx + 1L
      }
  out
}
