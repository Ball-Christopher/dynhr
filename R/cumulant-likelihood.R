## R/cumulant-likelihood.R
## --------------------------------------------------------------------------
## Phase F2 — Polyspectra / Cumulant Estimation (Mutschler 2015)
##
## Implements closed-form 3rd/4th cumulants of the pruned state-space,
## and hooks into the estimation API as likelihood = "cumulant".
##
## References:
##   Mutschler, W. (2015). Identification of DSGE models — The effect of
##     higher-order approximation and pruning. Journal of Economic Dynamics
##     and Control, 56, 1-33. [Henceforth "M2015"]
##   Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez, J. F.
##     (2018). The pruned state-space system for non-linear DSGE models.
##     Review of Economic Studies, 85(1), 1-49.
##
## Notation matches M2015 Sections 3.1-3.2:
##   x_t     — first-order (linear) state layer, Gaussian
##   x_t^(2) — second-order correction layer
##   s_t     — stacked pruned state [x_t; x_t^(2)] (dimension 2*n_s)
##   y_t     — observables
##   κ_k(y)  — k-th cumulant tensor of y (vec form)
## --------------------------------------------------------------------------


# ============================================================================
# Helper: perturbation-order detection
# ============================================================================

## Package-private store for one-time warnings (avoids MCMC warning spam).
.cumulant_warn_env <- new.env(parent = emptyenv())

#' Emit a warning at most once per session, keyed by `key`.
#' @noRd
.cumulant_warn_once <- function(key, msg) {
  if (isTRUE(.cumulant_warn_env[[key]])) return(invisible(NULL))
  .cumulant_warn_env[[key]] <- TRUE
  warning(msg, call. = FALSE)
}

#' Robustly detect the perturbation order of a decision-rule object
#'
#' \code{solve_perturbation()} records the order as a list element
#' (\code{dr$order}) and in the S3 class hierarchy ("DecisionRules2", etc.);
#' it does NOT set an \code{attr(dr, "order")} attribute. Earlier code gated
#' the order-3/4 cumulant terms on that attribute alone, so a standard
#' order-2 decision rule was silently treated as first order and the
#' skewness/kurtosis terms were dropped with no warning.
#'
#' Detection precedence:
#'   1. an explicit \code{attr(dr, "order")} (legacy / manual override),
#'   2. the \code{dr$order} list element set by the solvers,
#'   3. structural fallback: presence of the second-order blocks
#'      (\code{ghxx}/\code{ghss}) implies order >= 2.
#'
#' @param dr A DecisionRules (or higher) object.
#' @return Integer perturbation order (>= 1).
#' @noRd
.dr_perturbation_order <- function(dr) {
  a <- attr(dr, "order")
  if (!is.null(a) && length(a) == 1L && is.finite(a)) return(as.integer(a))
  if (!is.null(dr$order) && length(dr$order) == 1L && is.finite(dr$order))
    return(as.integer(dr$order))
  if (!is.null(dr$ghxx) || !is.null(dr$ghss)) return(2L)
  1L
}


# ============================================================================
# Helper: stationary moments of the Gaussian first-order layer
# ============================================================================

#' Compute stationary covariance of the first-order state layer
#'
#' Solves Σ_x = h_x·Σ_x·h_x' + h_u·Σ_ε·h_u'.
#'
#' @param hx n_s × n_s state transition matrix
#' @param hu n_s × n_u shock impact matrix
#' @param Sigma_e n_u × n_u shock covariance
#' @return n_s × n_s stationary covariance matrix
#' @noRd
.state_covariance <- function(hx, hu, Sigma_e) {
  B <- hu %*% Sigma_e %*% t(hu)
  solve_lyapunov(hx, B)
}


#' Stationary mean of the second-order correction layer
#'
#' μ_x2 = ½·(I - h_x)^{-1}·h_xx·vec(Σ_x)
#'
#' This is the unconditional mean of x_t^(2) under Gaussian shocks.
#'
#' @param hx  n_s × n_s state transition
#' @param hxx n_s × n_s^2 second-order state terms (full, expanded)
#' @param Sigma_x n_s × n_s state covariance
#' @return n_s vector: E[x_t^(2)]
#' @noRd
.second_order_mean <- function(hx, hxx, Sigma_x) {
  n_s <- nrow(hx)
  vec_Sigma <- as.numeric(Sigma_x)  # n_s^2, col-major
  rhs <- hxx %*% vec_Sigma          # n_s
  mu  <- 0.5 * solve(diag(n_s) - hx, rhs)
  as.numeric(mu)
}


# ============================================================================
# Third cumulant of pruned state-space  (M2015 eqs. 12-16)
# ============================================================================

#' Build the full RHS for the third cross-cumulant Lyapunov equation
#'
#' Computes RHS ∈ ℝ^{n_s × n_s^2} for solving C_3^{211} via
#'   (I - h_x^{⊗3}) · vec(C_3^{211}) = ½ · vec(RHS)
#'
#' Three driving-term families (each contributes to the stationary Lyapunov):
#'
#' hxx term (main, from the quadratic-in-x1 part of the x2 equation):
#'   RHS_hxx[i,(j,k)] = 2·(Σ_x·H_i·Σ_x)[j,k]
#'   where H_i = reshape(hxx[i,], n_s, n_s).
#'   Derivation: E[hxx_i*(x1⊗x1 - Σ_x) * x1_j * x1_k] by Isserlis
#'   = Σ_x[r,j]Σ_x[s,k] + Σ_x[r,k]Σ_x[s,j] contracted with H_i[r,s]
#'   = 2*(Σ_x*H_i*Σ_x)[j,k]  (symmetric H_i, symmetric Σ_x).
#'   NOTE: the trace term tr(H_i·Σ_x)·Σ_x[j,k] that appears in some
#'   references corresponds to the UN-centered 4th moment and is zero
#'   after centering (the Σ_x[r,s]·Σ_x[j,k] Isserlis term cancels the
#'   Σ_x[r,s] shift in the centered quadratic).
#'
#' hxu term (from the e⊗x1 cross-term in the x2 equation):
#'   When x2_i gets driven by sum_{u,p} Hxu_i[u,p]*e_u*x1_{t-1,p},
#'   and x1_{t,j} = hx*x1_{t-1} + hu*e_t, the current e_t creates a
#'   non-zero cross-cumulant through the hu*e part of x1_t:
#'   RHS_hxu[i,(j,k)] = (t(M_i)%*%N + t(N)%*%M_i)[j,k]
#'   where M_i = Hxu_i %*% Σ_x %*% t(hx)  [n_exo × n_s]
#'         N   = Σ_e %*% t(hu)              [n_exo × n_s]
#'   Hxu_i = reshape(hxu[i,], n_exo, n_s, byrow=TRUE).
#'
#' huu term (from the e⊗e quadratic in the x2 equation):
#'   RHS_huu[i,(j,k)] = 2*(hu·Σ_e·Uu_i·Σ_e·t(hu))[j,k]
#'   where Uu_i = reshape(huu[i,], n_exo, n_exo, byrow=TRUE).
#'   (The 2 comes from Isserlis on the centered e⊗e product.)
#'
#' @param hxx     n_s × n_s^2 second-order Hessian rows (state eqs only)
#' @param Sigma_x n_s × n_s stationary state covariance
#' @param hxu     Optional n_s × (n_exo*n_s) cross-term rows (state eqs)
#' @param huu     Optional n_s × n_exo^2 quadratic-in-e rows (state eqs)
#' @param hu      Optional n_s × n_exo shock impact matrix
#' @param hx      Optional n_s × n_s state transition (needed for hxu term)
#' @param Sigma_e Optional n_exo × n_exo shock covariance
#' @return n_s × n_s^2 matrix: full RHS for the third-cumulant Lyapunov
#' @noRd
.third_cumulant_rhs <- function(hxx, Sigma_x,
                                hxu = NULL, huu = NULL,
                                hu = NULL, hx = NULL, Sigma_e = NULL) {
  ## Stationary third cross-cumulant of the pruned order-2 state:
  ##   C211[p,(r,s)] = cum(x2_{t,p}, x1_{t,r}, x1_{t,s})
  ## satisfies the discrete Lyapunov equation
  ##   (I - hx^{⊗3}) vec(C211) = vec(R)
  ## where R[p] = sum over the hxx/hxu/huu driving terms of x2_p of
  ##   cum(driving_p, x1_{t,r}, x1_{t,s}).
  ##
  ## Derivation (e ⟂ x1_{t-1}; Gaussian Isserlis).  Two key facts that the
  ## old code got wrong:
  ##  (1) The linear factors x1_{t,r} = hx x1_{t-1} + hu e_t contribute to the
  ##      hxx-driving term ONLY through their hx x1_{t-1} part (the hu e_t part
  ##      is independent of the lagged-state quadratic), so the contraction
  ##      carries hx Σ_x, NOT Σ_x.  Hence the hxx driving is
  ##        R_hxx[p] = (hx Σ_x) Hxx_p (hx Σ_x)^T          (was Σ_x Hxx_p Σ_x).
  ##  (2) For the hxu cross term, the current shock e_t connects to the hu e_t
  ##      part of one linear factor while x1_{t-1} connects to the hx x1_{t-1}
  ##      part of the other:
  ##        R_hxu[p] = W_p + W_p^T,  W_p = (hu Σ_e) Xu_p (hx Σ_x)^T.
  ##  (3) The huu quadratic-in-e term connects to both hu e_t factors:
  ##        R_huu[p] = (hu Σ_e) Uu_p (hu Σ_e)^T.
  ## where Hxx_p = mat(hxx[p,], n_s, n_s), Xu_p = mat(hxu[p,], n_exo, n_s,
  ## byrow), Uu_p = mat(huu[p,], n_exo, n_exo, byrow).
  ##
  ## The solver .solve_third_cross_cumulant() multiplies the returned `rhs` by
  ## 1/2 (legacy convention), so we return rhs = 2 * R to recover R exactly.
  ## The old code applied that implicit 1/2 inconsistently across the three
  ## families (correct for huu, half-strength for hxu, wrong matrix for hxx).
  n_s <- nrow(Sigma_x)
  rhs <- matrix(0, n_s, n_s * n_s)

  # hxx contribution: R_hxx[p] = (hx Σ_x) Hxx_p (hx Σ_x)^T  (needs hx)
  hxS <- if (!is.null(hx)) hx %*% Sigma_x else Sigma_x
  for (i in seq_len(n_s)) {
    Hi  <- matrix(hxx[i, ], n_s, n_s)
    rhs[i, ] <- rhs[i, ] + as.numeric(2 * (hxS %*% Hi %*% t(hxS)))
  }

  # hxu contribution: R_hxu[p] = W_p + W_p^T,  W_p = (hu Σ_e) Xu_p (hx Σ_x)^T
  if (!is.null(hxu) && !is.null(hu) && !is.null(hx) && !is.null(Sigma_e)) {
    n_exo <- nrow(Sigma_e)
    huSe  <- hu %*% Sigma_e   # n_s × n_exo
    for (i in seq_len(n_s)) {
      Xu_i <- matrix(hxu[i, ], n_exo, n_s, byrow = TRUE)
      W_i  <- huSe %*% Xu_i %*% t(hxS)        # n_s × n_s
      rhs[i, ] <- rhs[i, ] + as.numeric(2 * (W_i + t(W_i)))
    }
  }

  # huu contribution: R_huu[p] = (hu Σ_e) Uu_p (hu Σ_e)^T
  if (!is.null(huu) && !is.null(hu) && !is.null(Sigma_e)) {
    n_exo  <- nrow(Sigma_e)
    huSe   <- hu %*% Sigma_e                  # n_s × n_exo
    for (i in seq_len(n_s)) {
      Uu_i <- matrix(huu[i, ], n_exo, n_exo, byrow = TRUE)
      rhs[i, ] <- rhs[i, ] + as.numeric(2 * (huSe %*% Uu_i %*% t(huSe)))
    }
  }

  rhs
}


#' Solve for the third cross-cumulant of the pruned state
#'
#' Solves:  vec(C_3^{211}) = (I - h_x^{⊗3})^{-1} · ½ · vec(RHS_3)
#' where C_3^{211}[i, j, k] = Cum[x_t^(2)[i], x_t[j], x_t[k]].
#'
#' This is M2015 eq. (13): κ(x_t^(2), x_t, x_t) at stationary.
#'
#' @param hx   n_s × n_s state transition
#' @param rhs  n_s × n_s^2 RHS from .third_cumulant_rhs()
#' @return n_s × n_s^2 matrix: the third cross-cumulant
#' @noRd
.solve_third_cross_cumulant <- function(hx, rhs) {
  n_s <- nrow(hx)

  # Solve (I - hx ⊗ hx ⊗ hx) · vec(C3) = ½ · vec(rhs)
  # via eigen-decomposition of hx for efficiency (like the compact Sylvester)
  eig   <- eigen(hx)
  V     <- eig$vectors
  lam   <- eig$values

  # Transform RHS into eigenbasis
  # For tensor T of dims (n_s, n_s, n_s), apply V^{-1} to mode 1,
  # and V^{-1} ⊗ V^{-1} to modes 2,3.
  # Mode 2,3 joint transformation: vec(V^{-1} · X · V^{-T})
  Vi      <- solve(V)
  rhs_tfm <- Vi %*% rhs                     # mode 1: V^{-1} · rhs
  rhs_tfm <- .apply_kron2(Vi, rhs_tfm)      # modes 2,3: V^{-1} · X · V^{-T}

  # Element-wise solve in eigenbasis
  c3_tfm <- matrix(0, n_s, n_s * n_s)
  for (i in seq_len(n_s)) {
    # For row i (mode-1 eigenvalue λ_i), the denominator is 1 - λ_i·λ_j·λ_k
    row_rhs <- matrix(rhs_tfm[i, ], n_s, n_s)
    for (j in seq_len(n_s)) {
      for (k in seq_len(n_s)) {
        denom <- 1 - lam[i] * lam[j] * lam[k]
        if (abs(denom) > 1e-14) {
          row_rhs[j, k] <- 0.5 * row_rhs[j, k] / denom
        } else {
          row_rhs[j, k] <- 0
        }
      }
    }
    c3_tfm[i, ] <- as.numeric(row_rhs)
  }

  # Transform back: V on mode 1, V⊗V on modes 2,3
  # Mode 2,3: V · X · V'  (apply V on left, V' on right)
  c3 <- V %*% c3_tfm
  c3 <- .apply_kron2(V, c3)

  # Symmetrize: C_3 should be symmetric in the last two indices (j,k)
  # C_3[i, (j,k)] = C_3[i, (k,j)] by definition of joint cumulant
  for (i in seq_len(n_s)) {
    Ci <- matrix(c3[i, ], n_s, n_s)
    Ci <- (Ci + t(Ci)) * 0.5
    c3[i, ] <- as.numeric(Ci)
  }

  c3
}


#' Apply (M ⊗ M) to rows of an n_s × n_s^2 matrix
#'
#' Each row of X is reshaped as n_s × n_s, then X_new = M·X·M'.
#'
#' @param M n_s × n_s transformation matrix
#' @param X n_s × n_s^2 matrix
#' @return Transformed n_s × n_s^2 matrix
#' @noRd
.apply_kron2 <- function(M, X) {
  n_s <- nrow(M)
  out <- matrix(0, n_s, n_s * n_s)
  for (i in seq_len(n_s)) {
    Xi  <- matrix(X[i, ], n_s, n_s)
    out[i, ] <- as.numeric(M %*% Xi %*% t(M))
  }
  out
}


#' Compute the unconditional third cumulant of observables
#'
#' Full analytic third cumulant of the pruned order-2 state-space under
#' Gaussian shocks. The pruned observable deviation from its mean is:
#'   dy = Z*x1 + D*e + Z*x2_cent + (1/2)*Ghxx*(x1⊗x1 - vec(Σ_x))
#'        + Ghxu*(e⊗x1) + (1/2)*Ghuu*(e⊗e - vec(Σ_e))
#' where (x1, x2_cent) are lagged states independent of the current shock e.
#'
#' The fully symmetric third cumulant kappa_3(y_i, y_j, y_k) receives four
#' families of contributions (all three permutations included):
#'
#' \describe{
#'   \item{B-term (x2 state)}{
#'     Three permutations of \code{cum(Z*x2_cent, Z*x1, Z*x1)} projected
#'     through the state cross-cumulant
#'     \eqn{C_3^{211}[p,r,s] = \mathrm{cum}(x2_p, x1_r, x1_s)}.
#'   }
#'   \item{C-term (ghxx quadratic)}{
#'     Three permutations of
#'     \eqn{E[Cx_i \cdot (Z_j x1) \cdot (Z_k x1)]}
#'     where \eqn{Cx_i = (1/2)\, ghxx_i (x1 \otimes x1 - \Sigma_x)}.
#'     By Isserlis: \eqn{Z_j \Sigma_x H_i \Sigma_x Z_k^T}.
#'   }
#'   \item{D-term (ghxu cross)}{
#'     Three permutations of
#'     \eqn{E[D_i \cdot (Z_j x1) \cdot (D_k e)]}
#'     where \eqn{D_i = ghxu_i (e \otimes x1)}, using independence of e and x1.
#'     Evaluates to \eqn{2 \, D_j \Sigma_e Xu_i \Sigma_x Z_k^T}
#'     where \eqn{Xu_i = \mathrm{mat}(ghxu_i, n_u, n_s)}.
#'   }
#'   \item{E-term (ghuu quadratic)}{
#'     Three permutations of
#'     \eqn{E[Eu_i \cdot (D_j e) \cdot (D_k e)]}
#'     where \eqn{Eu_i = (1/2)\, ghuu_i (e \otimes e - \Sigma_e)}.
#'     Evaluates to \eqn{(D_j \Sigma_e) Uu_i (D_k \Sigma_e)^T}
#'     where \eqn{Uu_i = \mathrm{mat}(ghuu_i, n_u, n_u)}.
#'   }
#' }
#'
#' @param dr    DecisionRules2 (or higher) with ghx, ghu, ghxx, ghuu, ghxu
#' @param model dynhr_mod (for shock variances)
#' @param params Named numeric parameter vector
#' @return List:
#'   \item{c3_obs}{n_obs x n_obs^2 matrix: vec(kappa_3(y)) with row=var1,
#'     cols=flattened (var2, var3)}
#'   \item{skewness}{Named numeric vector of marginal skewness (gamma_1)}
#'   \item{mean_obs}{Named numeric: unconditional mean of observables}
#'
#' @references
#'   Mutschler, W. (2015). Identification of DSGE models -- The effect of
#'     higher-order approximation and pruning. \emph{Journal of Economic
#'     Dynamics and Control}, 56, 1-33.
#'   Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez, J. F.
#'     (2018). The pruned state-space system for non-linear DSGE models.
#'     \emph{Review of Economic Studies}, 85(1), 1-49.
#' @export
compute_third_cumulant <- function(dr, model, params = NULL) {
  # ---- Extract structure ----
  endo      <- dr$endo_names
  exo       <- dr$exo_names
  state_idx <- dr$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)
  n_s       <- length(state_idx)

  if (is.null(params)) params <- model$param_values

  ghx  <- dr$ghx;   ghu  <- dr$ghu
  ghxx <- dr$ghxx %||% NULL
  ghuu <- dr$ghuu %||% NULL
  ghxu <- dr$ghxu %||% NULL
  ghss <- dr$ghss %||% rep(0, n_endo)

  # Shock covariance
  shock_stderr <- .get_shock_stderr(model, exo, params)
  Sigma_e <- diag(shock_stderr^2, n_exo)
  dimnames(Sigma_e) <- list(exo, exo)

  # ---- State-row submatrices ----
  hx  <- ghx [state_idx, , drop = FALSE]   # n_s × n_s
  hu  <- ghu [state_idx, , drop = FALSE]   # n_s × n_exo
  hxx <- if (!is.null(ghxx)) ghxx[state_idx, , drop = FALSE] else NULL
  hxu <- if (!is.null(ghxu)) ghxu[state_idx, , drop = FALSE] else NULL
  huu <- if (!is.null(ghuu)) ghuu[state_idx, , drop = FALSE] else NULL

  # ---- Stationary moments of the first-order layer ----
  Sigma_x <- .state_covariance(hx, hu, Sigma_e)

  # ---- Unconditional mean ----
  # E[y] = ys + ½·ghss (for order 2+)
  # For order 1, ghss is absent; E[y] = ys
  ys_obs <- dr$ys[endo]
  mean_obs <- ys_obs
  if (!is.null(dr$ghss)) {
    mean_obs <- mean_obs + 0.5 * dr$ghss
  }
  names(mean_obs) <- endo

  # ---- Third cumulant ----
  c3_obs <- matrix(0, n_endo, n_endo * n_endo)
  c3_211 <- NULL

  if (n_s > 0L && !is.null(hxx)) {
    # Observation matrix Z maps states to all endogenous variables.
    # The state columns of ghx correspond to state_idx; there are n_s states.
    Z   <- ghx[, seq_len(n_s), drop = FALSE]   # n_endo × n_s
    Zu  <- ghu                                  # n_endo × n_exo (shock loadings)

    # --- B-term: contribution from the second-order state x2_cent ---
    # C_3^{211}[p,r,s] = cum(x2_cent_p, x1_r, x1_s) solves the Lyapunov
    # (I - hx^{⊗3}) vec(C3) = ½ vec(RHS)
    # where RHS has three contributions:
    #   RHS_hxx: from the 0.5*hxx*(x1⊗x1 - Σ_x) driving term
    #   RHS_hxu: from the hxu*(e⊗x1) driving term — non-zero because e_t and
    #             x1_t = hx*x1_{t-1} + hu*e_t share the current shock e_t
    #   RHS_huu: from the 0.5*huu*(e⊗e - Σ_e) driving term — non-zero
    #             because x1_t contains hu*e_t
    # See .third_cumulant_rhs() documentation for derivations.
    rhs    <- .third_cumulant_rhs(hxx, Sigma_x,
                                  hxu = hxu, huu = huu,
                                  hu = hu, hx = hx, Sigma_e = Sigma_e)
    c3_211 <- .solve_third_cross_cumulant(hx, rhs)  # n_s × n_s^2

    # Project C3^{211} to observables with ALL THREE permutations of (i,j,k):
    # Perm 1: x2 in position i  → c3_obs[i,(j,k)] += Z[i,p]*Z[j,r]*Z[k,s]*C211[p,r,s]
    # Perm 2: x2 in position j  → c3_obs[i,(j,k)] += Z[j,p]*Z[i,r]*Z[k,s]*C211[p,r,s]
    # Perm 3: x2 in position k  → c3_obs[i,(j,k)] += Z[k,p]*Z[i,r]*Z[j,s]*C211[p,r,s]
    #
    # Pre-compute "ZC3" = Z %*% c3_211  [n_endo × n_s^2]: row m contains
    # sum_p Z[m,p]*C211[p,r,s] (flattened over r,s).
    ZC3 <- Z %*% c3_211   # n_endo × n_s^2

    for (i in seq_len(n_endo)) {
      # Perm 1: x2 in position i (current row)
      c3_obs[i, ] <- c3_obs[i, ] +
        as.numeric(Z %*% matrix(ZC3[i, ], n_s, n_s) %*% t(Z))

      # Perms 2 & 3: x2 in positions j and k.
      # For each pair (j,k): c3_obs[i,(j,k)] += ZC3[j, r_s]*Z[i,r]*Z[k,s] + ZC3[k,r_s]*Z[i,r]*Z[j,s]
      # Pre-compute SigmaTerm_i = Sigma_x^{contracted with Z[i,]}: for each state r,
      # Z[i,r] weights the C211 rows. We vectorise over j,k by building
      # the n_endo × n_endo matrices for each permutation.
      #
      # Perm 2: kappa(y_i, y_j, y_k) += Z_j * C211_flat * (Z_i ⊗ Z_k)
      # = sum_p ZC3[j,p_rs]*Z[i,r]*Z[k,s]
      # Build the n_endo x n_endo matrix M2 where M2[j,k] = perm2 contribution:
      # M2 = ZC3 %*% apply_kron2_to_rows(c3_211, Z_i, Z)
      # Simplified: for each j, ZC3[j,] is 1 x n_s^2; reshape as n_s x n_s;
      # then M2[j, k] = Z[i,] %*% reshape(ZC3[j,]) %*% Z[k,]^T.
      zi <- Z[i, , drop = FALSE]   # 1 × n_s
      M2 <- matrix(0, n_endo, n_endo)
      M3 <- matrix(0, n_endo, n_endo)
      for (jj in seq_len(n_endo)) {
        M2_jj <- matrix(ZC3[jj, ], n_s, n_s)  # C211 projected through Z[j,]
        # M2[jj, k] = Z[i,] %*% M2_jj %*% Z[k,]^T
        M2[jj, ] <- as.numeric(zi %*% M2_jj %*% t(Z))
      }
      # Perm 3: kappa(y_i, y_j, y_k) += ZC3[k,] projected through (Z_i ⊗ Z_j)
      # M3[j,k] = perm3 = Z[i,] %*% reshape(ZC3[k,]) %*% Z[j,]^T = t(M2)[j,k]
      M3 <- t(M2)

      # Add both permutations into c3_obs[i, (j-1)*n_endo+k]:
      c3_obs[i, ] <- c3_obs[i, ] + as.numeric(M2) + as.numeric(M3)
    }

    # --- C-term: contribution from the ghxx quadratic in x1 ---
    # For each endo variable m, the centered quadratic is:
    #   Cx_m = (1/2) ghxx[m,] (x1⊗x1 - vec(Σ_x))
    # Its cross-cumulant with (ghx[j,]*x1, ghx[k,]*x1) is (by Isserlis):
    #   E[Cx_m * (Z_j x1) * (Z_k x1)] = Z_j * Σ_x * H_m * Σ_x * Z_k^T
    # where H_m = mat(ghxx[m,], n_s, n_s).
    #
    # The three permutations of where Cx appears in kappa_3(y_i, y_j, y_k):
    # Perm 1: Cx_i contributes Z_j * Σ_x * H_i * Σ_x * Z_k^T to c3_obs[i,(j,k)]
    # Perm 2: Cx_j contributes Z_i * Σ_x * H_j * Σ_x * Z_k^T to c3_obs[i,(j,k)]
    # Perm 3: Cx_k contributes Z_i * Σ_x * H_k * Σ_x * Z_j^T to c3_obs[i,(j,k)]
    #
    # Pre-compute the "Sigma_x H_m Sigma_x projected" matrix for each m:
    # ZSigHSigZT_m[j,k] = Z_j * Σ_x * H_m * Σ_x * Z_k^T  — this is an n_endo×n_endo matrix.
    # Pre-build: A_m = Z %*% Sigma_x %*% H_m %*% Sigma_x %*% t(Z)  [n_endo × n_endo].
    # Then: Perm1 adds A_i to c3_obs[i, (j,k)].
    #       Perm2 adds A_j[i,k] to c3_obs[i, (j,k)].
    #       Perm3 adds A_k[i,j] to c3_obs[i, (j,k)].
    #
    # Perm1: for row i, add the n_endo×n_endo matrix A_i flattened.
    # Perms 2&3: for each (j,k), c3_obs[i,(j,k)] += A_j[i,k] + A_k[i,j].
    #
    # Build A[m] = Z %*% Sigma_x %*% H_m %*% Sigma_x %*% t(Z) for all m.
    SigZ <- Sigma_x %*% t(Z)   # n_s × n_endo  (Σ_x · Z^T)
    # A[m, j, k] = Z[j,] %*% Sigma_x %*% H_m %*% Sigma_x %*% Z[k,]^T
    #            = (Z %*% Sigma_x)[j,] %*% H_m %*% (Z %*% Sigma_x)[k,]^T
    ZSig <- Z %*% Sigma_x    # n_endo × n_s  (Z · Σ_x)

    # Materialise n_endo matrices A_m, one per observable:
    A_m_list <- vector("list", n_endo)
    for (m in seq_len(n_endo)) {
      H_m <- matrix(ghxx[m, ], n_s, n_s)
      # A_m = ZSig %*% H_m %*% t(ZSig)  [n_endo × n_endo]
      A_m_list[[m]] <- ZSig %*% H_m %*% t(ZSig)
    }

    for (i in seq_len(n_endo)) {
      # Perm 1: Cx in position i → add A_i to c3_obs[i, :]
      c3_obs[i, ] <- c3_obs[i, ] + as.numeric(A_m_list[[i]])

      # Perms 2 & 3: Cx in positions j and k.
      # c3_obs[i, (j-1)*n_endo+k] += A_j[i,k] + A_k[i,j]
      # Build n_endo × n_endo matrices P2 and P3:
      P2 <- matrix(0, n_endo, n_endo)
      P3 <- matrix(0, n_endo, n_endo)
      for (jj in seq_len(n_endo)) {
        # Perm 2: Cx in position j=jj; contribution to (i, jj, k):
        #   P2[jj, k] = A_jj[i, k]
        P2[jj, ] <- A_m_list[[jj]][i, ]
        # Perm 3: Cx in position k=jj; contribution to (i, j, jj):
        #   P3[j, jj] = A_jj[i, j] → column jj of P3 = row i of A_jj
        P3[, jj] <- A_m_list[[jj]][i, ]
      }
      c3_obs[i, ] <- c3_obs[i, ] + as.numeric(P2) + as.numeric(P3)
    }

    # --- D-term: contribution from ghxu*(e⊗x1) cross term ---
    # D_m = ghxu[m,]*(e⊗x1) where Xu_m = mat(ghxu[m,], n_exo, n_s).
    # e (current shock) is independent of (x1, x2) (lagged).
    # E[D_m * (Z_j x1) * (Zu_k e)]:
    #   = sum_{u,p,r,v} Xu_m[u,p]*Z[j,r]*Zu[k,v] * E[e_u*x1_p*x1_r*e_v]
    #   = sum_{u,p,r,v} Xu_m[u,p]*Z[j,r]*Zu[k,v] * Σ_e[u,v]*Σ_x[p,r]
    #   = Zu[k,]*Σ_e * Xu_m * Σ_x * Z[j,]^T   (scalar for fixed m,j,k)
    # E[D_m * (Zu_j e) * (Z_k x1)] = Zu[j,]*Σ_e * Xu_m * Σ_x * Z[k,]^T (swap j,k labels)
    # E[D_m * (Z_j x1) * (Z_k x1)] = 0 (e indep x1 so E[e_u*x1_p*x1_r*x1_s]=E[e_u]*...=0)
    # E[D_m * (Zu_j e) * (Zu_k e)] = 0 (E[e_u*x1_p*e_v*e_w]=E[x1_p]*...=0)
    #
    # So for each (m, j, k) where e appears in exactly ONE of the two linear partners:
    # E[D_m * T_j * T_k] where T_j = Z_j*x1 + Zu_j*e:
    #   = E[D_m*(Z_j x1)*(Zu_k e)] + E[D_m*(Zu_j e)*(Z_k x1)]
    #   = Zu[k,]*Σ_e*Xu_m*Σ_x*Z[j,]^T + Zu[j,]*Σ_e*Xu_m*Σ_x*Z[k,]^T
    #   = (SeXuSig_m)^T * Zu[k,]^T at Z[j,] + swap
    # where SeXuSig_m = Zu * Σ_e * Xu_m * Σ_x  [n_endo × n_s] ... actually n_endo × n_endo:
    # Define B_m = Zu %*% Sigma_e %*% Xu_m %*% Sigma_x %*% t(Z)  [n_endo × n_endo]
    # Then E[D_m * T_j * T_k] = B_m[k, j] + B_m[j, k]  (both sides contribute)
    #
    # But note: kappa_3(y_i, y_j, y_k) involves ALL THREE permutations of which variable is "D":
    # Perm 1: D in position i → c3_obs[i,(j,k)] += E[D_i*(T_j)*(T_k)] = B_i[k,j] + B_i[j,k]
    # Perm 2: D in position j → c3_obs[i,(j,k)] += E[D_j*(T_i)*(T_k)] = B_j[k,i] + B_j[i,k]
    # Perm 3: D in position k → c3_obs[i,(j,k)] += E[D_k*(T_i)*(T_j)] = B_k[j,i] + B_k[i,j]
    if (!is.null(ghxu)) {
      SeZ <- Sigma_e %*% t(Zu)   # n_exo × n_endo: (Σ_e · Zu^T)

      # Build B_m = Zu %*% Sigma_e %*% Xu_m %*% Sigma_x %*% t(Z)  for each m.
      # Xu_m = mat(ghxu[m,], n_exo, n_s)
      B_m_list <- vector("list", n_endo)
      for (m in seq_len(n_endo)) {
        # byrow=TRUE: ghxu[m,] is laid out as (exo SLOW, state FAST) per the
        # Kronecker product e⊗x1; reshape to [u, p] = [exo, state] row-major.
        Xu_m <- matrix(ghxu[m, ], n_exo, n_s, byrow = TRUE)
        B_m_list[[m]] <- Zu %*% Sigma_e %*% Xu_m %*% Sigma_x %*% t(Z)  # n_endo × n_endo
      }

      for (i in seq_len(n_endo)) {
        # Perm 1: D in position i
        # c3_obs[i,(j,k)] += B_i[k,j] + B_i[j,k]  = (B_i + t(B_i))[j,k]
        c3_obs[i, ] <- c3_obs[i, ] + as.numeric(B_m_list[[i]] + t(B_m_list[[i]]))

        # Perms 2 & 3: D in positions j and k
        # Perm 2: c3_obs[i,(j,k)] += E[D_j * T_i * T_k] = B_j[k,i] + B_j[i,k]
        # Perm 3: c3_obs[i,(j,k)] += E[D_k * T_i * T_j] = B_k[j,i] + B_k[i,j]
        # Build n_endo × n_endo matrices D2 and D3:
        D2 <- matrix(0, n_endo, n_endo)
        D3 <- matrix(0, n_endo, n_endo)
        for (jj in seq_len(n_endo)) {
          # Perm 2: D in position j=jj; c3_obs[i,(jj,k)] += B_jj[k,i] + B_jj[i,k]
          # D2[jj, k] = B_jj[k,i] + B_jj[i,k]
          D2[jj, ] <- B_m_list[[jj]][, i] + B_m_list[[jj]][i, ]
          # Perm 3: D in position k=jj; c3_obs[i,(j,jj)] += B_jj[j,i] + B_jj[i,j]
          # D3[j, jj] = B_jj[j,i] + B_jj[i,j]
          D3[, jj] <- B_m_list[[jj]][, i] + B_m_list[[jj]][i, ]
        }
        c3_obs[i, ] <- c3_obs[i, ] + as.numeric(D2) + as.numeric(D3)
      }
    }

    # --- E-term: contribution from ghuu*(e⊗e - Σ_e) quadratic in e ---
    # Eu_m = (1/2)*ghuu[m,]*(e⊗e - Σ_e), Uu_m = mat(ghuu[m,], n_exo, n_exo).
    # E[Eu_m * (Zu_j e) * (Zu_k e)]:
    #   = (1/2)*sum_{u,v,w,z} Uu_m[u,v]*Zu[j,w]*Zu[k,z]*E[(e_u e_v-Σ_e[u,v])*e_w*e_z]
    # By Isserlis: E[(e_u e_v - Σ_e[u,v])*e_w*e_z] = Σ_e[u,w]Σ_e[v,z] + Σ_e[u,z]Σ_e[v,w]
    # E[Eu_m * (Zu_j e) * (Zu_k e)] = (1/2)*(ZuSe_j * Uu_m * ZuSe_k^T + ZuSe_j * Uu_m^T * ZuSe_k^T)
    #   where ZuSe_j = Zu[j,] %*% Σ_e  [1 × n_exo]
    # Since Uu_m is symmetric: = ZuSe_j %*% Uu_m %*% ZuSe_k^T
    # E[Eu_m * (Z_j x1) * (Zu_k e)] = 0 (x1 indep of e, E[x1]=0)
    # E[Eu_m * (Z_j x1) * (Z_k x1)] = 0 (Eu_m is centered, indep of x1)
    # E[Eu_m * (Zu_j e) * (Z_k x1)] = 0 (same reason)
    #
    # So: E[Eu_m * T_j * T_k] = ZuSe_j %*% Uu_m %*% ZuSe_k^T   (scalar for fixed m,j,k)
    # where ZuSe_m = Zu[m,] %*% Σ_e  (1 × n_exo).
    #
    # Three permutations of kappa_3(y_i, y_j, y_k):
    # Perm 1: Eu in position i → c3_obs[i,(j,k)] += ZuSe_j %*% Uu_i %*% ZuSe_k^T
    # Perm 2: Eu in position j → c3_obs[i,(j,k)] += ZuSe_i %*% Uu_j %*% ZuSe_k^T
    # Perm 3: Eu in position k → c3_obs[i,(j,k)] += ZuSe_i %*% Uu_k %*% ZuSe_j^T
    if (!is.null(ghuu)) {
      # Pre-compute ZuSe[m,] = Zu[m,] %*% Sigma_e  [n_endo × n_exo]
      ZuSe <- Zu %*% Sigma_e   # n_endo × n_exo

      # For each m: E_m_mat[j,k] = ZuSe[j,] %*% Uu_m %*% ZuSe[k,]^T
      #           = (ZuSe %*% Uu_m %*% t(ZuSe))[j,k]
      E_m_list <- vector("list", n_endo)
      for (m in seq_len(n_endo)) {
        # byrow=TRUE: ghuu[m,] is laid out (exo SLOW, exo FAST) per e⊗e;
        # reshape to [u, v] = [exo, exo] row-major.
        Uu_m <- matrix(ghuu[m, ], n_exo, n_exo, byrow = TRUE)
        E_m_list[[m]] <- ZuSe %*% Uu_m %*% t(ZuSe)   # n_endo × n_endo
      }

      for (i in seq_len(n_endo)) {
        # Perm 1: Eu in position i → add E_i_mat to c3_obs[i, :]
        c3_obs[i, ] <- c3_obs[i, ] + as.numeric(E_m_list[[i]])

        # Perms 2 & 3: Eu in positions j and k
        # Perm 2: c3_obs[i,(j,k)] += E_j_mat[i,k] = ZuSe[i,]*Uu_j*ZuSe[k,]^T
        # Perm 3: c3_obs[i,(j,k)] += E_k_mat[i,j] = ZuSe[i,]*Uu_k*ZuSe[j,]^T
        E2 <- matrix(0, n_endo, n_endo)
        E3 <- matrix(0, n_endo, n_endo)
        for (jj in seq_len(n_endo)) {
          # Perm 2: E in j=jj; E2[jj,k] = E_jj_mat[i,k]
          E2[jj, ] <- E_m_list[[jj]][i, ]
          # Perm 3: E in k=jj; E3[j,jj] = E_jj_mat[i,j]
          E3[, jj] <- E_m_list[[jj]][i, ]
        }
        c3_obs[i, ] <- c3_obs[i, ] + as.numeric(E2) + as.numeric(E3)
      }
    }
  }

  # ---- Marginal skewness ----
  # γ_1 = κ_3(y_i) / σ_i^3  where κ_3(y_i) = c3_obs[i, (i,i)]
  skewness <- rep(0, n_endo)
  names(skewness) <- endo

  # First need variance (first-order is fine for normalisation)
  var_y <- compute_moments(dr, model, params = params)$var_cov
  sd_y <- sqrt(pmax(diag(var_y), 0))

  for (i in seq_len(n_endo)) {
    # c3_obs[i, :] is n_obs^2 vector; extract element at (i,i)
    idx <- (i - 1) * n_endo + i
    if (sd_y[i] > 1e-14) {
      skewness[i] <- c3_obs[i, idx] / sd_y[i]^3
    }
  }

  dimnames(c3_obs) <- list(endo, paste0(rep(endo, each = n_endo), "_", endo))

  list(
    c3_obs   = c3_obs,
    skewness = skewness,
    mean_obs = mean_obs,
    Sigma_x  = Sigma_x,
    c3_211   = c3_211
  )
}


# ============================================================================
# Fourth cumulant of pruned state-space  (M2015 eqs. 17-20)
# ============================================================================

#' Compute fourth-order Gaussian moment E[x⊗x⊗x⊗x]
#'
#' For zero-mean Gaussian x with covariance Σ, the fourth moment is:
#'   E[x_i x_j x_k x_l] = Σ_{ij}·Σ_{kl} + Σ_{ik}·Σ_{jl} + Σ_{il}·Σ_{jk}
#'
#' Returns an (n_s^2 × n_s^2) matrix where column (k,l) indexes the
#' n_s^2 flattened vector for pair (i,j).
#'
#' @param Sigma_x n_s × n_s state covariance
#' @return n_s^2 × n_s^2 matrix: E[(x⊗x) · (x⊗x)']
#' @noRd
.fourth_moment_gaussian <- function(Sigma_x) {
  n_s <- nrow(Sigma_x)
  n_s2 <- n_s * n_s
  M4 <- matrix(0, n_s2, n_s2)

  # Vectorized computation
  vec_Sigma <- as.numeric(Sigma_x)  # n_s^2

  for (ij in seq_len(n_s2)) {
    i <- ((ij - 1L) %/% n_s) + 1L
    j <- ((ij - 1L) %% n_s) + 1L

    for (kl in seq_len(n_s2)) {
      k <- ((kl - 1L) %/% n_s) + 1L
      l <- ((kl - 1L) %% n_s) + 1L

      M4[ij, kl] <- Sigma_x[i, j] * Sigma_x[k, l] +
                    Sigma_x[i, k] * Sigma_x[j, l] +
                    Sigma_x[i, l] * Sigma_x[j, k]
    }
  }
  M4
}


#' Marginal 4th cumulant of pruned observables via a linear+quadratic form
#'
#' The centered pruned observable for a single variable i is exactly a
#' linear-plus-quadratic form in the Gaussian innovation history
#' \eqn{w = (e_t, e_{t-1}, \dots)}:
#'   \deqn{y_{i,t} - \mathrm{E}[y_i] = a' w + (w' M w - \mathrm{tr}(M S))}
#' where \eqn{S} is the (block-diagonal) covariance of \eqn{w}.  This holds
#' because (a) the first-order state is an MA(\eqn{\infty}) of \eqn{e};
#' (b) the second-order state \eqn{x^{(2)}} is the \eqn{h_x}-discounted
#' accumulation of the quadratic driving terms; and (c) the contemporaneous
#' \eqn{g_{hxx}/g_{hxu}/g_{huu}} families are quadratic in \eqn{(x_{t-1}, e_t)}.
#'
#' For such a form the cumulants are closed-form (generalized chi-square):
#'   \deqn{\kappa_2 = a'Sa + 2\,\mathrm{tr}((MS)^2)}
#'   \deqn{\kappa_4 = 48\,\mathrm{tr}((MS)^4) + 48\, a' S M S M S a}
#' (the odd cross terms \eqn{\kappa(L,L,L,Q)} and \eqn{\kappa(L,Q,Q,Q)} vanish
#' by Gaussian parity; \eqn{\kappa_4(L)=0} as \eqn{L} is exactly Gaussian).
#' This was verified against the univariate oracle \eqn{Y=g+\tfrac12 g^2}:
#' \eqn{\kappa_4 = 12 v^3 + 3 v^4}, \eqn{v=\mathrm{Var}(g)}.
#'
#' The innovation history is truncated at \code{n_lag} lags; the contribution
#' of \eqn{x^{(2)}} decays like \eqn{\rho(h_x)^{\text{lag}}}, so \code{n_lag}
#' is chosen so the truncation error sits far below Monte-Carlo error.
#'
#' @param Z,Zu    n_endo × n_s / n_endo × n_exo observation loadings (= ghx
#'   state columns / ghu).
#' @param hx,hu   n_s × n_s / n_s × n_exo first-order state recursion.
#' @param hxx,hxu,huu state-row second-order blocks (n_s × n_s^2 etc.).
#' @param ghxx,ghxu,ghuu full endo-row second-order blocks.
#' @param Sigma_x,Sigma_e stationary state / shock covariance.
#' @param n_lag   innovation-history truncation depth.
#' @return List with marginal \code{kappa2} and \code{kappa4} (length n_endo).
#' @noRd
.fourth_cumulant_qform_marginal <- function(Z, Zu, hx, hu,
                                            hxx, hxu, huu,
                                            ghxx, ghxu, ghuu,
                                            Sigma_x, Sigma_e, n_lag) {
  n_s   <- nrow(hx)
  n_exo <- nrow(Sigma_e)
  n_endo <- nrow(Z)

  # Innovation window w = (e_t, e_{t-1}, ..., e_{t-n_lag}); S = blockdiag(Sigma_e).
  L  <- n_lag
  nw <- (L + 1L) * n_exo
  S  <- matrix(0, nw, nw)
  for (j in 0:L) S[j * n_exo + seq_len(n_exo), j * n_exo + seq_len(n_exo)] <- Sigma_e
  eblk <- function(k) k * n_exo + seq_len(n_exo)

  # Powers of hx: hxpow[[p+1]] = hx^p, p = 0..L+1.
  hxpow <- vector("list", L + 2L)
  hxpow[[1]] <- diag(n_s)
  for (p in seq_len(L + 1L)) hxpow[[p + 1L]] <- hxpow[[p]] %*% hx

  # x1_{t-q} as a linear map from w (n_s × nw):
  #   x1_{t-q} = sum_{m>=0} hx^m hu e_{t-q-m}; pick e_{t-k}, k = q + m.
  x1map <- function(q) {
    Lm <- matrix(0, n_s, nw)
    if (L >= q) for (m in 0:(L - q)) {
      k <- q + m
      Lm[, eblk(k)] <- hxpow[[m + 1L]] %*% hu
    }
    Lm
  }
  X1lags <- lapply(seq_len(L + 1L), x1map)  # X1lags[[q]] = x1_{t-q}
  E0 <- matrix(0, n_exo, nw); E0[, eblk(0)] <- diag(n_exo)  # e_t selector

  # Closed-form cumulants of y = a'w + w'Mw (M symmetric).
  k2f <- function(a, M) {
    MS <- M %*% S
    as.numeric(crossprod(a, S %*% a)) + 2 * sum(MS * t(MS))
  }
  ## kappa4 = 48 tr((MS)^4) [trace] + 48 a'SMSMSa [chain].  Return both so the
  ## caller can swap in an exact closed-form chain and keep only the (small,
  ## fast-converging) trace from this truncated-window quadrature.
  k4f <- function(a, M) {
    MS  <- M %*% S
    MS2 <- MS %*% MS
    SMSMSa <- S %*% (M %*% (S %*% (M %*% (S %*% a))))
    c(trace = 48 * sum(MS2 * t(MS2)),
      chain = 48 * as.numeric(crossprod(a, SMSMSa)))
  }

  # Accumulate a bilinear contribution  Arows' Coef Brows.  The caller
  # symmetrises M exactly once at the end (M <- (M + t(M))/2), so the raw
  # (unsymmetrised) crossprod is accumulated here: sum_b 0.5*(P_b + t(P_b)) =
  # 0.5*(sum_b P_b + t(sum_b P_b)).  Dropping the per-call symmetrisation avoids
  # ~3*n_lag transposes of the dense n_w x n_w matrix per observable (the
  # dominant cost -- t.default was ~35% of this routine's runtime).
  addquad <- function(M, Arows, Brows, Coef) {
    M + crossprod(Arows, Coef %*% Brows)   # t(Arows) %*% Coef %*% Brows
  }

  kappa2       <- numeric(n_endo)
  kappa4       <- numeric(n_endo)
  kappa4_trace <- numeric(n_endo)
  kappa4_chain <- numeric(n_endo)

  X1tm1 <- X1lags[[1]]  # x1_{t-1}
  for (i in seq_len(n_endo)) {
    # Linear part: a' w = Z_i x1_{t-1} + Zu_i e_t.
    a <- as.numeric(Z[i, ] %*% X1tm1 + Zu[i, ] %*% E0)

    M <- matrix(0, nw, nw)
    Hi  <- matrix(ghxx[i, ], n_s, n_s)
    Xui <- matrix(ghxu[i, ], n_exo, n_s, byrow = TRUE)
    Uui <- matrix(ghuu[i, ], n_exo, n_exo, byrow = TRUE)

    # Contemporaneous quadratic families (use x1_{t-1}, e_t):
    M <- addquad(M, X1tm1, X1tm1, 0.5 * Hi)   # 0.5 ghxx (x1_{t-1} ⊗ x1_{t-1})
    M <- addquad(M, E0,    X1tm1, Xui)        # ghxu (e_t ⊗ x1_{t-1})
    M <- addquad(M, E0,    E0,    0.5 * Uui)  # 0.5 ghuu (e_t ⊗ e_t)

    # Second-order state term  Z_i x2c_{t-1} = sum_{j>=1} (Z_i hx^{j-1}) driver_{t-j}
    # driver_{t-j} = 0.5 hxx(x1_{t-1-j}⊗x1_{t-1-j}) + hxu(e_{t-j}⊗x1_{t-1-j})
    #               + 0.5 huu(e_{t-j}⊗e_{t-j}).
    Zi <- Z[i, ]
    for (j in seq_len(L)) {
      coefv <- as.numeric(Zi %*% hxpow[[j]])   # Z_i hx^{j-1}
      Xq <- X1lags[[j + 1L]]                   # x1_{t-1-j}
      Ek <- matrix(0, n_exo, nw); Ek[, eblk(j)] <- diag(n_exo)  # e_{t-j}
      Hc <- matrix(0, n_s, n_s)
      Uc <- matrix(0, n_exo, n_exo)
      Xc <- matrix(0, n_exo, n_s)
      for (p in seq_len(n_s)) {
        Hc <- Hc + coefv[p] * matrix(hxx[p, ], n_s, n_s)
        Uc <- Uc + coefv[p] * matrix(huu[p, ], n_exo, n_exo, byrow = TRUE)
        Xc <- Xc + coefv[p] * matrix(hxu[p, ], n_exo, n_s, byrow = TRUE)
      }
      M <- addquad(M, Xq, Xq, 0.5 * Hc)
      M <- addquad(M, Ek, Xq, Xc)
      M <- addquad(M, Ek, Ek, 0.5 * Uc)
    }
    M <- (M + t(M)) / 2

    kappa2[i]   <- k2f(a, M)
    k4          <- k4f(a, M)
    kappa4_trace[i] <- k4[["trace"]]
    kappa4_chain[i] <- k4[["chain"]]
    kappa4[i]   <- k4[["trace"]] + k4[["chain"]]
  }

  list(kappa2 = kappa2, kappa4 = kappa4,
       kappa4_trace = kappa4_trace, kappa4_chain = kappa4_chain)
}


#' Solve the discrete Lyapunov tensor equation (I - hx^{⊗4}) vec(C) = vec(R).
#'
#' Uses the eigendecomposition of hx so the equation decouples element-wise in
#' the eigenbasis (the same device as \code{.solve_third_cross_cumulant}).
#'
#' @param hx n_s × n_s state transition.
#' @param R  n_s × n_s × n_s × n_s driving tensor.
#' @return n_s × n_s × n_s × n_s solution tensor.
#' @noRd
.solve_lyap4 <- function(hx, R) {
  n_s <- nrow(hx)
  eg  <- eigen(hx); V <- eg$vectors; lam <- eg$values; Vi <- solve(V)
  applymode <- function(A, M, mode) {
    d  <- dim(A); rest <- setdiff(1:4, mode)
    A2 <- matrix(aperm(A, c(mode, rest)), d[mode])
    A2 <- array(M %*% A2, c(nrow(M), d[rest]))
    aperm(A2, order(c(mode, rest)))
  }
  Rt <- applymode(applymode(applymode(applymode(R, Vi, 1), Vi, 2), Vi, 3), Vi, 4)
  Ct <- array(0 + 0i, dim(Rt))
  for (p in seq_len(n_s)) for (q in seq_len(n_s))
    for (r in seq_len(n_s)) for (s in seq_len(n_s))
      Ct[p, q, r, s] <- Rt[p, q, r, s] /
        (1 - lam[p] * lam[q] * lam[r] * lam[s])
  Re(applymode(applymode(applymode(applymode(Ct, V, 1), V, 2), V, 3), V, 4))
}

#' Fourth-order state cross-cumulant C2211[p,q,r,s] = cum(x2_p, x2_q, x1_r, x1_s).
#'
#' The stationary 4th-order cross-cumulant of the pruned order-2 state solves
#' (I - hx^{⊗4}) vec(C2211) = vec(RHS4), where RHS4 collects the driving terms
#' from substituting the x2 recursion x2 = hx x2_{-1} + W (W the quadratic
#' driving) and x1 = hx x1_{-1} + hu e into the cumulant.  Because x2 appears at
#' most twice and its only nonzero joint cumulants with x1 are Cov(x2,x2),
#' cum(x2,x1,x1)=C211 and C2211 itself, RHS4 reduces to (using e ⟂ x1_{-1}):
#'   1-x2 terms (use C211 + Sigma_x), and 0-x2 terms = cum(W_p, W_q, x1_r, x1_s)
#'   which is a pure-Gaussian 2-quadratic-2-linear cumulant (Wick necklace).
#' Validated to Monte-Carlo error against simulated state cross-cumulants.
#'
#' @param hx,hu state transition / shock-impact (n_s × n_s, n_s × n_exo).
#' @param hxx,hxu,huu state second-order rows (n_s × ...).
#' @param Sigma_x,Sigma_e stationary state / shock covariance.
#' @param C211 n_s × n_s^2 third cross-cumulant (col r + (s-1)n_s).
#' @return n_s × n_s × n_s × n_s tensor.
#' @noRd
.fourth_cross_cumulant <- function(hx, hu, hxx, hxu, huu, Sigma_x, Sigma_e, C211) {
  n_s   <- nrow(hx); n_exo <- nrow(Sigma_e)
  Hxx <- lapply(seq_len(n_s), function(p) matrix(hxx[p, ], n_s, n_s))
  Xu  <- lapply(seq_len(n_s), function(p) matrix(hxu[p, ], n_exo, n_s, byrow = TRUE))
  Uu  <- lapply(seq_len(n_s), function(p) matrix(huu[p, ], n_exo, n_exo, byrow = TRUE))
  C3  <- function(i, j, k) C211[i, j + (k - 1L) * n_s]

  ## ---- 1-x2 terms (Pa in slot1, W in slot2) + p<->q transpose ----
  ## (hxx family) cum(hx_p a, 0.5 Hxx_q(b⊗b), hx_r b, hx_s b)
  T1 <- array(0, c(n_s, n_s, n_s, n_s))
  for (p in seq_len(n_s)) for (q in seq_len(n_s))
    for (r in seq_len(n_s)) for (s in seq_len(n_s)) {
      acc <- 0
      for (i in seq_len(n_s)) for (k in seq_len(n_s)) for (l in seq_len(n_s)) {
        hh <- hx[p, i] * hx[r, k] * hx[s, l]; if (hh == 0) next
        for (mm in seq_len(n_s)) for (nn in seq_len(n_s)) {
          Hq <- Hxx[[q]][mm, nn]; if (Hq == 0) next
          K <- C3(i, mm, k) * Sigma_x[nn, l] + C3(i, nn, l) * Sigma_x[mm, k] +
               C3(i, mm, l) * Sigma_x[nn, k] + C3(i, nn, k) * Sigma_x[mm, l]
          acc <- acc + 0.5 * hh * Hq * K
        }
      }
      T1[p, q, r, s] <- acc
    }
  ## (hxu family) cross of e_t with one linear, x1_{-1} with the other
  T2 <- array(0, c(n_s, n_s, n_s, n_s))
  for (p in seq_len(n_s)) for (q in seq_len(n_s))
    for (r in seq_len(n_s)) for (s in seq_len(n_s)) {
      acc <- 0
      for (i in seq_len(n_s)) for (k in seq_len(n_s)) for (l in seq_len(n_s)) {
        cc <- C3(i, k, l); if (cc == 0) next
        for (u in seq_len(n_exo)) for (v in seq_len(n_exo)) {
          xu <- Xu[[q]][u, k]; if (xu == 0) next
          se <- Sigma_e[u, v]; if (se == 0) next
          acc <- acc + hx[p, i] * xu * hu[r, v] * hx[s, l] * se * cc
          acc <- acc + hx[p, i] * xu * hx[r, l] * hu[s, v] * se * cc
        }
      }
      T2[p, q, r, s] <- acc
    }
  oneX2 <- T1 + T2
  oneX2 <- oneX2 + aperm(oneX2, c(2, 1, 3, 4))

  ## ---- 0-x2 terms = cum(W_p, W_q, x1_r, x1_s) pure-Gaussian Wick ----
  d  <- n_s + n_exo
  Sg <- matrix(0, d, d); Sg[seq_len(n_s), seq_len(n_s)] <- Sigma_x
  Sg[(n_s + 1L):d, (n_s + 1L):d] <- Sigma_e
  bi <- seq_len(n_s); ei <- (n_s + 1L):d; Sgi <- solve(Sg)
  Wm <- vector("list", n_s); ar <- matrix(0, d, n_s)
  for (p in seq_len(n_s)) {
    M <- matrix(0, d, d)
    M[bi, bi] <- 0.5 * Hxx[[p]]
    M[bi, ei] <- 0.5 * t(Xu[[p]]); M[ei, bi] <- 0.5 * Xu[[p]]
    M[ei, ei] <- 0.5 * Uu[[p]]
    Wm[[p]] <- M
  }
  for (r in seq_len(n_s)) { ar[bi, r] <- hx[r, ]; ar[ei, r] <- hu[r, ] }
  zeroX2 <- array(0, c(n_s, n_s, n_s, n_s))
  for (p in seq_len(n_s)) for (q in seq_len(n_s)) {
    Cp <- Sg %*% Wm[[p]] %*% Sg; Cq <- Sg %*% Wm[[q]] %*% Sg
    CpSq <- Cp %*% Sgi %*% Cq; CqSp <- Cq %*% Sgi %*% Cp
    for (r in seq_len(n_s)) for (s in seq_len(n_s))
      zeroX2[p, q, r, s] <- 4 * as.numeric(ar[, r] %*% CpSq %*% ar[, s]) +
                            4 * as.numeric(ar[, r] %*% CqSp %*% ar[, s])
  }
  .solve_lyap4(hx, oneX2 + zeroX2)
}

#' Exact closed-form "chain" part of the marginal 4th cumulant.
#'
#' For the observable y_i (centered) = a'g + g'Mc g + Z_i x2 (g = (x1_{-1}, e_t)
#' Gaussian, Mc the contemporaneous quadratic, Z_i x2 the second-order state),
#' the marginal kappa4 = 48 a'S M S M S a (chain) + 48 tr((MS)^4) (trace), with
#' M = Mc + Mx the full quadratic.  The dominant chain expands exactly into
#'   contemp = 48 a'S Mc S Mc S a,
#'   domB    = 6 Z_i^{⊗4} : C2211   (= 48 a'S Mx S Mx S a),
#'   cross   = 48 sum_p Z_i[p] sum_{j,k} a_x1[j] C211(p,j,k) (Mc Sg a)_x1[k]
#'             (= 48 (a'S Mc S Mx S a + a'S Mx S Mc S a)).
#' No innovation-window truncation -> exact at any rho(hx) < 1.
#'
#' @return n_endo numeric vector of the chain contribution to kappa4.
#' @noRd
.fourth_cumulant_chain_closed <- function(ghx, ghu, ghxx, ghxu, ghuu,
                                          hx, Sigma_x, Sigma_e, C211, C2211) {
  n_endo <- nrow(ghx); n_s <- nrow(hx); n_exo <- nrow(Sigma_e)
  C3 <- function(i, j, k) C211[i, j + (k - 1L) * n_s]
  d  <- n_s + n_exo
  Sg <- matrix(0, d, d); Sg[seq_len(n_s), seq_len(n_s)] <- Sigma_x
  Sg[(n_s + 1L):d, (n_s + 1L):d] <- Sigma_e
  bi <- seq_len(n_s); ei <- (n_s + 1L):d
  out <- numeric(n_endo)
  for (i in seq_len(n_endo)) {
    Zi <- ghx[i, seq_len(n_s)]
    a  <- numeric(d); a[bi] <- Zi; a[ei] <- ghu[i, ]
    Mc <- matrix(0, d, d)
    Mc[bi, bi] <- 0.5 * matrix(ghxx[i, ], n_s, n_s)
    Xui <- matrix(ghxu[i, ], n_exo, n_s, byrow = TRUE)
    Mc[bi, ei] <- 0.5 * t(Xui); Mc[ei, bi] <- 0.5 * Xui
    Mc[ei, ei] <- 0.5 * matrix(ghuu[i, ], n_exo, n_exo, byrow = TRUE)
    Mc <- (Mc + t(Mc)) / 2
    contemp <- 48 * as.numeric(a %*% (Sg %*% Mc %*% Sg %*% Mc %*% Sg) %*% a)
    domB <- 0
    for (p in seq_len(n_s)) for (q in seq_len(n_s))
      for (r in seq_len(n_s)) for (s in seq_len(n_s))
        domB <- domB + Zi[p] * Zi[q] * Zi[r] * Zi[s] * C2211[p, q, r, s]
    domB <- 6 * domB
    MSga <- as.numeric(Mc %*% Sg %*% a)
    cross <- 0
    for (p in seq_len(n_s)) {
      acc <- 0
      for (j in seq_len(n_s)) for (k in seq_len(n_s))
        acc <- acc + a[j] * C3(p, j, k) * MSga[k]
      cross <- cross + Zi[p] * acc
    }
    cross <- 48 * cross
    out[i] <- contemp + domB + cross
  }
  out
}

#' Compute the unconditional fourth cumulant (excess kurtosis) of observables
#'
#' Following Mutschler (2015) eqs. (17)-(20). The fourth cumulant captures
#' unconditional excess kurtosis induced by the second-order layer of the
#' pruned state-space.
#'
#' For Gaussian shocks, the first-order layer has zero excess kurtosis.
#' Non-zero kurtosis arises from the quadratic (second-order) layer.
#'
#' @param dr    DecisionRules2 (or higher)
#' @param model dynhr_mod
#' @param params Named numeric parameter vector
#' @return List:
#'   \item{kurtosis_obs}{Named numeric vector of marginal excess kurtosis γ_2}
#'   \item{c4_obs}{n_obs × n_obs^3 matrix: fourth cumulant of observables}
#' @export
compute_fourth_cumulant <- function(dr, model, params = NULL) {
  endo      <- dr$endo_names
  exo       <- dr$exo_names
  state_idx <- dr$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)
  n_s       <- length(state_idx)

  if (is.null(params)) params <- model$param_values

  ghx  <- dr$ghx;   ghu  <- dr$ghu
  ghxx <- dr$ghxx %||% NULL
  ghxu <- dr$ghxu %||% NULL
  ghuu <- dr$ghuu %||% NULL

  shock_stderr <- .get_shock_stderr(model, exo, params)
  Sigma_e <- diag(shock_stderr^2, n_exo)

  hx  <- ghx [state_idx, , drop = FALSE]
  hu  <- ghu [state_idx, , drop = FALSE]
  hxx <- if (!is.null(ghxx)) ghxx[state_idx, , drop = FALSE] else NULL
  hxu <- if (!is.null(ghxu)) ghxu[state_idx, , drop = FALSE] else NULL
  huu <- if (!is.null(ghuu)) ghuu[state_idx, , drop = FALSE] else NULL

  # First-order variance
  Sigma_x <- .state_covariance(hx, hu, Sigma_e)

  # Fourth cumulant marginal estimates
  kurtosis <- rep(0, n_endo)
  names(kurtosis) <- endo

  # c4_obs: n_endo × n_endo^3 (flattened)
  c4_obs <- matrix(0, n_endo, n_endo * n_endo * n_endo)

  if (n_s > 0L && !is.null(hxx) && !is.null(hxu) && !is.null(huu)) {
    Z  <- ghx[, seq_len(n_s), drop = FALSE]  # n_endo × n_s
    Zu <- ghu                                 # n_endo × n_exo

    # Innovation-history truncation depth: the x^(2) accumulation enters the
    # observable with weight ~ rho(hx)^lag, and it enters the (normalized)
    # kurtosis quadratically, so rho^n_lag <= 1e-4 already puts the truncation
    # error far below any Monte-Carlo error (verified: kurtosis at n_lag=120 vs
    # 480 agrees to 4 decimals even at rho=0.95).  Capped at 150 to keep the
    # dense quadratic-form tractable for the GMM / gradient path.
    rho   <- max(abs(eigen(hx, only.values = TRUE)$values))
    rho   <- min(max(rho, 1e-6), 1 - 1e-8)

    ## kappa4 = chain (48 a'SMSMSa, exact closed-form, no truncation) + trace
    ## (48 tr((MS)^4)).  The chain is the dominant ~99.5% and is the piece that
    ## decays only as rho^{2*lag}; computing it in closed form removes the
    ## innovation-window truncation entirely (exact at any rho < 1, fixing
    ## near-unit-root models).  The trace decays as rho^{4*lag} -- twice as fast
    ## -- so it converges at half the n_lag, keeping the dense quadratic-form
    ## build small.
    C211_4 <- .solve_third_cross_cumulant(
      hx, .third_cumulant_rhs(hxx, Sigma_x, hxu = hxu, huu = huu,
                              hu = hu, hx = hx, Sigma_e = Sigma_e))
    C2211  <- .fourth_cross_cumulant(hx, hu, hxx, hxu, huu,
                                     Sigma_x, Sigma_e, C211_4)
    chain  <- .fourth_cumulant_chain_closed(
      ghx, ghu, ghxx, ghxu, ghuu, hx, Sigma_x, Sigma_e, C211_4, C2211)

    n_lag <- as.integer(min(max(ceiling(log(1e-7) / (4 * log(rho))), 20L), 150L))
    qf <- .fourth_cumulant_qform_marginal(
      Z, Zu, hx, hu, hxx, hxu, huu, ghxx, ghxu, ghuu,
      Sigma_x, Sigma_e, n_lag = n_lag)
    kappa4_vec <- chain + qf$kappa4_trace

    ## Exact marginal variance for the kurtosis denominator.  The q-form's
    ## kappa2 = a'Sa + 2 tr((MS)^2) converges only as rho^{2*lag}, so the
    ## reduced trace-n_lag would under-resolve it; compute_moments() gives the
    ## exact pruned order-2 variance with no truncation.
    var_total <- diag(compute_moments(dr, model, params = params)$var_cov)

    for (i in seq_len(n_endo)) {
      k2 <- var_total[i]
      if (k2 > 1e-30) {
        kurtosis[i] <- kappa4_vec[i] / k2^2     # excess kurtosis γ_2 = κ_4 / σ^4

        # Fill the marginal (i,i,i,i) entry of the c4 tensor (κ_4, not γ_2):
        # col-major flatten of (i2,i3,i4) into the row-i column index.
        col_idx <- 1 + (i - 1L) + n_endo * (i - 1L) + n_endo^2 * (i - 1L)
        if (col_idx <= ncol(c4_obs)) c4_obs[i, col_idx] <- kappa4_vec[i]
      }
    }
  }

  list(
    kurtosis_obs = kurtosis,
    c4_obs       = c4_obs
  )
}


# ============================================================================
# Sample cumulant computation
# ============================================================================

#' Compute sample (K-statistics) cumulants from observed data
#'
#' Unbiased k-statistics for cumulants up to order 4, following the
#' classical definition (Fisher 1929).
#'
#' @param Y T × n data matrix
#' @param max_order Maximum cumulant order (2, 3, or 4; default 4)
#' @return List:
#'   \item{mean}{n-vector of sample mean}
#'   \item{var_cov}{n × n sample variance-covariance}
#'   \item{c3}{n × n^2 sample third cumulant (K_3)}
#'   \item{c4}{n × n^3 sample fourth cumulant (K_4)}
#'   \item{n_obs}{number of observations}
#' @export
sample_cumulants <- function(Y, max_order = 4L) {
  if (is.null(dim(Y))) Y <- matrix(Y, ncol = 1L)
  T_obs <- nrow(Y)
  n     <- ncol(Y)

  if (T_obs < 4L) stop("Need at least 4 observations for cumulant estimates.")

  # Helper: de-mean
  Yc <- scale(Y, center = TRUE, scale = FALSE)  # T × n

  # 1. Sample mean
  mu_hat <- colMeans(Y)
  names(mu_hat) <- colnames(Y)

  # 2. Sample variance (unbiased: divide by T-1)
  var_hat <- crossprod(Yc) / (T_obs - 1)  # n × n

  result <- list(
    mean    = mu_hat,
    var_cov = var_hat,
    n_obs   = T_obs
  )

  if (max_order >= 3L) {
    # 3. Third cumulant (k-statistic k_3)
    # k_3[i,j,k] = T/((T-1)(T-2)) * Σ_t Yc[t,i]·Yc[t,j]·Yc[t,k]
    c3 <- array(0, dim = c(n, n, n))
    factor3 <- T_obs / ((T_obs - 1) * (T_obs - 2))
    for (t in seq_len(T_obs)) {
      yt <- Yc[t, ]
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          for (k in seq_len(n)) {
            c3[i, j, k] <- c3[i, j, k] + yt[i] * yt[j] * yt[k]
          }
        }
      }
    }
    c3 <- c3 * factor3

    # Flatten to n × n^2 (row = first index, col = (j,k) col-major)
    c3_flat <- matrix(0, n, n * n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        for (k in seq_len(n)) {
          col_idx <- (j - 1L) * n + k
          c3_flat[i, col_idx] <- c3[i, j, k]
        }
      }
    }
    dimnames(c3_flat) <- list(colnames(Y), NULL)
    result$c3    <- c3_flat
    result$c3_arr <- c3
  }

  if (max_order >= 4L) {
    # 4. Fourth cumulant (k-statistic k_4)
    # k_4[i,j,k,l] = T^2/((T-1)(T-2)(T-3)) * ...
    #   [ Σ(y_i y_j y_k y_l) - (T-1)/(T(T+1))·Σ(y_i y_j)·Σ(y_k y_l)·(all pairings) ]
    # Simplified: compute the fourth central moment and subtract the
    # variance-pairing contribution.
    factor4 <- T_obs^2 / ((T_obs - 1) * (T_obs - 2) * (T_obs - 3))

    # Compute the fourth product moment
    m4 <- array(0, dim = c(n, n, n, n))
    for (t in seq_len(T_obs)) {
      yt <- Yc[t, ]
      for (i in seq_len(n)) {
        for (j in seq_len(n)) {
          for (k in seq_len(n)) {
            for (l in seq_len(n)) {
              m4[i, j, k, l] <- m4[i, j, k, l] +
                yt[i] * yt[j] * yt[k] * yt[l]
            }
          }
        }
      }
    }
    m4 <- m4 / T_obs  # fourth moment (raw)

    # Subtract variance-pairing for the cumulant
    # For Gaussian: κ_4 = m_4 - 3·vec(Σ)'s pairing
    c4 <- array(0, dim = c(n, n, n, n))
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        for (k in seq_len(n)) {
          for (l in seq_len(n)) {
            c4[i, j, k, l] <- m4[i, j, k, l] -
              (var_hat[i, j] * var_hat[k, l] +
               var_hat[i, k] * var_hat[j, l] +
               var_hat[i, l] * var_hat[j, k]) *
              (T_obs - 1) / (T_obs + 1)
          }
        }
      }
    }
    c4 <- c4 * factor4

    # Flatten to n × n^3 (row = first index, col = (j,k,l) col-major)
    c4_flat <- matrix(0, n, n * n * n)
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        for (k in seq_len(n)) {
          for (l in seq_len(n)) {
            col_idx <- (j - 1L) * n * n + (k - 1L) * n + l
            c4_flat[i, col_idx] <- c4[i, j, k, l]
          }
        }
      }
    }
    dimnames(c4_flat) <- list(colnames(Y), NULL)
    result$c4    <- c4_flat
    result$c4_arr <- c4
  }

  result
}


# ============================================================================
# Cumulant-based log-likelihood
# ============================================================================

#' Cumulant-based log-likelihood
#'
#' Computes log p(Y|θ) by matching empirical and model-implied cumulants
#' of the observables (Mutschler 2015, Section 3). The cumulant vector
#' is asymptotically normal, giving a quadratic-form log-likelihood:
#'
#'   log L(θ|Y) = -½ · (m̂ - m(θ))' · W · (m̂ - m(θ))
#'
#' where m̂ is the vector of sample cumulants, m(θ) is the model-implied
#' cumulants, and W is a weight matrix (identity or precision-weighted).
#'
#' @param data T × n_obs data matrix with column names matching obs_vars
#' @param dr Decision rules object (order ≥ 1)
#' @param model dynhr_mod
#' @param params Named numeric parameter vector
#' @param obs_vars Character vector of observable variable names
#' @param orders Integer vector: which cumulant orders to match
#'   (default 1:4, i.e. mean, variance, skewness, kurtosis)
#' @param weight_method How to weight cumulant discrepancies:
#'   "identity" — equal weight (default)
#'   "precision" — inverse of estimated asymptotic variance (stub; use
#'     \code{weight_matrix} instead)
#' @param weight_matrix Optional p × p numeric matrix.  When non-NULL, the
#'   GMM objective \code{-0.5 * T * t(delta) W delta} is used directly,
#'   overriding \code{weight_method}.  Build \code{W} via
#'   \code{estimate_gmm_weight_matrix()}.
#' @param me_variance Measurement error variance added to model variance
#' @return Log-likelihood value (scalar, -Inf on failure)
#' @noRd
.cumulant_loglik <- function(data, dr, model, params, obs_vars,
                              orders = 1:4, weight_method = "identity",
                              weight_matrix = NULL,
                              me_variance = 0) {
  n_obs <- length(obs_vars)
  T_obs <- nrow(data)

  # ---- 1. Compute sample cumulants ----
  max_order <- max(orders)
  sc <- sample_cumulants(data, max_order = max_order)

  # ---- 2. Compute model-implied cumulants ----
  # First-order variance (always needed for moments)
  moments <- compute_moments(dr, model, params = params)
  Sigma_y_model <- moments$var_cov[obs_vars, obs_vars, drop = FALSE]
  mean_model    <- dr$ys[obs_vars]

  # Account for second-order mean correction
  if (!is.null(dr$ghss)) {
    mean_model <- mean_model + 0.5 * dr$ghss[obs_vars]
  }

  # Add measurement error to model variance
  if (me_variance > 0) {
    diag(Sigma_y_model) <- diag(Sigma_y_model) + me_variance
  }

  # Third and fourth model cumulants
  c3_model <- NULL
  c4_model <- NULL

  dr_order <- .dr_perturbation_order(dr)
  if (any(orders >= 3L) && dr_order < 2L) {
    .cumulant_warn_once(
      "cumulant_orders34_needs_order2",
      paste0(
        "Cumulant orders 3-4 were requested but the decision rule is ",
        "first-order (no ghxx/ghss): only orders 1-2 contribute to the ",
        "cumulant log-likelihood. Re-solve with ",
        "solve_perturbation(order = 2) to activate the skewness/kurtosis ",
        "terms."))
  }
  if (any(orders >= 3L) && dr_order >= 2L) {
    c3_result <- compute_third_cumulant(dr, model, params)
    # Select only observables
    obs_idx <- match(obs_vars, dr$endo_names)
    c3_model <- c3_result$c3_obs[obs_idx, , drop = FALSE]
    # Keep only columns corresponding to (obs, obs) pairs
    cols_keep <- rep(obs_idx, each = n_obs) * 0  # need to subset n_obs^2 cols
    c3_obs_only <- matrix(0, n_obs, n_obs * n_obs)
    for (a in seq_len(n_obs)) {
      for (b in seq_len(n_obs)) {
        src_col <- (obs_idx[a] - 1L) * length(dr$endo_names) + obs_idx[b]
        dst_col <- (a - 1L) * n_obs + b
        c3_obs_only[a, dst_col] <- c3_model[a, src_col]
      }
    }
    c3_model <- c3_obs_only
  }

  if (any(orders >= 4L) && dr_order >= 2L) {
    c4_result <- compute_fourth_cumulant(dr, model, params)
    obs_idx <- match(obs_vars, dr$endo_names)
    c4_model <- c4_result$c4_obs[obs_idx, , drop = FALSE]
  }

  # ---- 3. Build the moment vector ----
  # m = [mean; vec(var); vec(c3); vec(c4)]
  m_model <- numeric(0)
  m_emp   <- numeric(0)

  if (1L %in% orders) {
    m_model <- c(m_model, mean_model)
    m_emp   <- c(m_emp, sc$mean[obs_vars])
  }

  if (2L %in% orders) {
    m_model <- c(m_model, as.numeric(Sigma_y_model))
    m_emp   <- c(m_emp, as.numeric(sc$var_cov[obs_vars, obs_vars, drop = FALSE]))
  }

  if (3L %in% orders && !is.null(c3_model)) {
    m_model <- c(m_model, as.numeric(c3_model))
    sc_c3 <- sc$c3
    if (!is.null(sc_c3)) {
      # Subset to observables
      sc_c3_obs <- matrix(0, n_obs, n_obs * n_obs)
      # c3 is already in obs_vars order if colnames match
      sc_c3_obs <- sc_c3
      m_emp <- c(m_emp, as.numeric(sc_c3_obs))
    }
  }

  if (4L %in% orders && !is.null(c4_model)) {
    m_model <- c(m_model, as.numeric(c4_model))
    sc_c4 <- sc$c4
    if (!is.null(sc_c4)) {
      m_emp <- c(m_emp, as.numeric(sc_c4))
    }
  }

  # ---- 4. Check finite values ----
  if (any(!is.finite(m_model))) return(-Inf)

  # ---- 5. Compute log-likelihood ----
  delta <- m_emp - m_model

  if (!is.null(weight_matrix)) {
    # GMM weighted objective: -T/2 * delta' W delta
    # With W = diag(p)/p this equals the identity path EXACTLY (Oracle 1).
    loglik <- -0.5 * T_obs * as.numeric(t(delta) %*% weight_matrix %*% delta)
  } else if (weight_method == "identity") {
    # Equal weight on all moments: -T/(2p) * sum(delta^2)
    loglik <- -0.5 * sum(delta^2) / length(delta) * T_obs
  } else if (weight_method == "precision") {
    # Precision weighting not yet implemented without a pre-estimated W.
    # Build W via estimate_gmm_weight_matrix() and pass as weight_matrix.
    stop("weight_method = \"precision\" requires a pre-estimated weight matrix. ",
         "Call estimate_gmm_weight_matrix() to build W and pass it as ",
         "weight_matrix to .cumulant_loglik() or make_log_posterior_cumulant().")
  } else {
    loglik <- -0.5 * sum(delta^2) / length(delta) * T_obs
  }

  loglik
}


# ============================================================================
# log-posterior factory (cumulant version)
# ============================================================================

#' Create a cached log-posterior evaluator using cumulant-based likelihood
#'
#' Mirrors \code{make_log_posterior()} but uses the cumulant-matching
#' likelihood instead of the Gaussian Kalman filter. Works with any
#' perturbation order >= 1; higher orders provide more cumulant content.
#'
#' At order 1 (Gaussian), the cumulant likelihood matches only the mean
#' and variance (identical to Gaussian MLE). At order >= 2, the model
#' generates non-zero third and fourth cumulants which are matched
#' against the data.
#'
#' @param model       dynhr_mod from \code{\link{parse_mod}}
#' @param data        Observation matrix (T × n_obs), columns = obs_vars
#' @param prior_spec  Prior spec data.frame from \code{\link{prior_spec}}
#' @param obs_vars    Character vector of observed variable names
#' @param compiled    dynhr_compiled from \code{\link{compile_model}}
#' @param me_variance Measurement error variance (default 0)
#' @param cumulant_orders Integer vector: which cumulant orders to match
#'   (default 1:4, i.e. mean, variance, skewness, kurtosis)
#' @param cumulant_weight Weight method: "identity" or "precision"
#' @param weight_matrix Optional p × p numeric matrix pre-computed by
#'   \code{\link{estimate_gmm_weight_matrix}}.  When non-\code{NULL},
#'   the GMM objective \code{-T/2 * t(delta) W delta} is used at every
#'   evaluation, overriding \code{cumulant_weight}.  \code{NULL} (default)
#'   falls back to \code{cumulant_weight}.
#' @param system_priors  Optional \code{system_prior_spec} object (from
#'   \code{\link{system_prior_spec}}) providing system-prior log-density
#'   contributions evaluated at each draw.  \code{NULL} disables system priors.
#' @param ...         Additional arguments forwarded to the perturbation
#'   solve.  \code{order} selects the perturbation order (1 or 2); when
#'   omitted it defaults to 2 if \code{cumulant_orders} requests an
#'   order-3/4 cumulant (so the skewness/kurtosis terms actually activate)
#'   and 1 otherwise.  \code{h} sets the order-2 finite-difference step.
#' @return A function \code{function(theta)} returning a named list
#'   \code{list(logpost, loglik, logprior)}
#' @export
make_log_posterior_cumulant <- function(model, data, prior_spec, obs_vars,
                                         compiled,
                                         me_variance = 0,
                                         cumulant_orders = 1:4,
                                         cumulant_weight = "identity",
                                         weight_matrix = NULL,
                                         system_priors = NULL,
                                         ...) {
  ## data/prior_spec/obs_vars/me_variance/cumulant_* are only referenced inside
  ## the returned closure, so without forcing they remain unevaluated promises
  ## pointing at the caller's frame. A mirai daemon that ships this closure
  ## before it has been called once would fail to resolve them.
  force(data); force(prior_spec); force(obs_vars); force(me_variance)
  force(cumulant_orders); force(cumulant_weight); force(weight_matrix)
  force(system_priors)

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence

  sys_cache <- cache_system_structure(compiled)

  # Capture extra args for solve_perturbation (order, h, etc.)
  solver_args <- list(...)

  ## Resolve the perturbation order the closure will solve at. The cumulant
  ## skewness/kurtosis terms (orders 3-4) are functions of ghxx/ghss, which
  ## only the order-2 (or higher) decision rule provides; the bare
  ## .solve_from_system() core is strictly first-order. Honour an explicit
  ## order = passed through `...`; otherwise default to 2 whenever
  ## cumulant_orders requests an order-3/4 cumulant, so the default
  ## cumulant_orders = 1:4 actually activates the documented terms.
  solve_order <- solver_args$order
  if (is.null(solve_order))
    solve_order <- if (any(cumulant_orders >= 3L)) 2L else 1L
  solve_order <- as.integer(solve_order)
  if (!solve_order %in% c(1L, 2L))
    stop("make_log_posterior_cumulant() supports order 1 or 2 (got ",
         solve_order, "); cumulants above the fourth are not implemented.")
  solver_h <- solver_args$h %||% 1e-4

  function(theta) {
    lp <- log_prior(theta, prior_spec)
    if (!is.finite(lp))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    params <- .apply_theta_to_params(model, theta)

    ss_result <- solve_steady_state(model, compiled, params, verbose = FALSE)
    if (is.null(ss_result) || !isTRUE(ss_result$converged))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## Re-derive SSM-computed params for a consistent linearization point
    ## (no-op for non-SSM-parameter models; Tier 13 #1).
    params <- ss_result$params %||% params
    sys <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)
    dr  <- .solve_from_system(sys, model, compiled, ss_result$ss, params, FALSE)
    if (is.null(dr) || !isTRUE(dr$bk_satisfied))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ghx_state <- dr$ghx[dr$state_idx, , drop = FALSE]
    if (max(Mod(eigen(ghx_state, only.values = TRUE)$values)) >= 1)
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## Higher-order solve: lift the first-order rule to order 2 so that
    ## ghxx/ghss are available for the cumulant orders 3-4. Stationarity is
    ## already checked on the first-order block above; the order-2 solve only
    ## adds the quadratic terms. A failed quadratic solve falls back to the
    ## first-order rule, which trips the order-detection warning downstream.
    if (solve_order >= 2L) {
      Sigma_e <- .get_shock_cov(model, model$varexo_names, params)
      dr2 <- tryCatch(
        solve_perturbation_order2(model, compiled, ss_result$ss, params,
                                  dr1 = dr, Sigma_e = Sigma_e,
                                  h = solver_h, verbose = FALSE),
        error = function(e) NULL)
      if (!is.null(dr2)) dr <- dr2
    }

    # Compute cumulant-based log-likelihood
    loglik <- .cumulant_loglik(data, dr, model, params, obs_vars,
                               orders = cumulant_orders,
                               weight_method = cumulant_weight,
                               weight_matrix = weight_matrix,
                               me_variance = me_variance)
    if (!is.finite(loglik))
      return(list(logpost = -Inf, loglik = -Inf, logprior = lp))

    ## System priors: penalty on model features evaluated from the solved dr.
    if (!is.null(system_priors)) {
      sp_lp <- .eval_system_priors(
        system_priors,
        list(theta   = theta,
             model   = model,
             dr      = dr,
             Sigma_e = .get_shock_cov(model, model$varexo_names, params),
             params  = params))
      if (!is.finite(sp_lp))
        return(list(logpost = -Inf, loglik = loglik, logprior = lp))
      lp <- lp + sp_lp
    }

    list(logpost = .dynhr_opt("power_posterior", default = 1) * loglik + lp,
         loglik = loglik, logprior = lp)
  }
}


# ============================================================================
# Analytic long-run covariance (Omega) for orders 1-2
# ============================================================================

#' Build the commutation matrix K_{n,n}
#'
#' K_{n,n} is the n^2 x n^2 permutation matrix such that
#' K_{n,n} vec(A) = vec(A').
#'
#' @param n Integer dimension
#' @return n^2 × n^2 sparse-in-spirit permutation matrix
#' @noRd
.commutation_matrix <- function(n) {
  n2 <- n * n
  K  <- matrix(0, n2, n2)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      # vec(A) element (i,j) is at row (j-1)*n + i
      # vec(A') element (j,i) is at row (i-1)*n + j
      row_src <- (j - 1L) * n + i   # position of (i,j) in vec(A)
      row_dst <- (i - 1L) * n + j   # position of (j,i) in vec(A') = K vec(A)
      K[row_dst, row_src] <- 1
    }
  }
  K
}


#' Analytic long-run covariance Omega for GMM orders 1-2
#'
#' For a stationary Gaussian linear state-space with observable covariance
#' \eqn{\Sigma_y} and lag-\eqn{h} autocovariance \eqn{\Gamma(h) = Z hx^h \Sigma_s Z'}:
#'
#' \describe{
#'   \item{Omega_11}{T * Sigma_y  — asymptotic variance of the sample mean}
#'   \item{Omega_22}{T * long-run covariance of vec(sample_cov)}
#'   \item{Omega_12}{0  — Gaussian odd-even cumulant cross-independence}
#' }
#'
#' The moment vector for order 2 uses vec(Sigma_y) (full n_obs^2 elements, not
#' vech), matching .build_moment_vector / .cumulant_loglik exactly.  The
#' long-run covariance of vec(yc_t yc_t') is therefore an n_obs^2 x n_obs^2
#' matrix assembled from the Isserlis-theorem lag contributions:
#'
#'   A(h) = Gamma(h) kron Gamma(h) + K_{n,n} (Gamma(h) kron Gamma(-h))
#'        = Gamma(h) kron Gamma(h) + K_{n,n} (Gamma(h) kron Gamma(h)')
#'
#' Omega_22 = sum_{h=-inf}^{inf} A(h)
#'          = A(0) + sum_{h>=1} [A(h) + A(-h)]
#'
#' A(-h) = Gamma(-h) kron Gamma(-h) + K (Gamma(-h) kron Gamma(h))
#'       = Gamma(h)' kron Gamma(h)' + K (Gamma(h)' kron Gamma(h))
#'
#' For h=0 (Sigma_y symmetric):
#'   A(0) = Sigma_y kron Sigma_y + K (Sigma_y kron Sigma_y)
#'        = (I + K)(Sigma_y kron Sigma_y)
#'
#' The lag sum is truncated when the Frobenius norm of A(h) falls below tol.
#'
#' @param dr        Decision rule (order >= 1)
#' @param model     dynhr_mod
#' @param params    Named parameter vector
#' @param obs_vars  Character vector of observable names
#' @param orders    Integer vector (must be subset of 1:2)
#' @param ridge     Ridge fraction for regularization
#' @param max_lags  Maximum lags to sum (default 2000; truncated earlier if tol met)
#' @param tol       Frobenius-norm tolerance for lag truncation (default 1e-12)
#' @return p x p weight matrix W = Omega_reg^{-1}, with attributes
#' @noRd
.analytic_gmm_weight_matrix <- function(dr, model, params, obs_vars,
                                        orders, ridge,
                                        max_lags = 2000L,
                                        tol = 1e-12) {
  n_obs <- length(obs_vars)
  obs_idx <- match(obs_vars, dr$endo_names)
  if (any(is.na(obs_idx)))
    stop(".analytic_gmm_weight_matrix: obs_vars not all found in dr$endo_names")

  # ---- Model-implied observable covariance ----
  moments <- compute_moments(dr, model, params = params)
  Sigma_y <- moments$var_cov[obs_vars, obs_vars, drop = FALSE]   # n_obs x n_obs

  if (!all(is.finite(Sigma_y)))
    stop(".analytic_gmm_weight_matrix: Sigma_y contains non-finite values ",
         "(unit root or non-stationary model?)")

  # ---- Extract state-space matrices for lag autocovariance recursion ----
  ghx       <- dr$ghx
  state_idx <- dr$state_idx
  n_state   <- length(state_idx)

  ghx_state   <- ghx[state_idx, , drop = FALSE]         # n_state x n_state
  Z           <- ghx[obs_idx, seq_len(n_state), drop = FALSE]  # n_obs x n_state
  Sigma_state <- moments$Sigma_state   # n_state x n_state (from compute_moments)

  # ---- Determine total moment dimension p ----
  p_11 <- if (1L %in% orders) n_obs      else 0L
  p_22 <- if (2L %in% orders) n_obs * n_obs else 0L
  p    <- p_11 + p_22

  # ---- Build block-diagonal Omega ----
  Omega <- matrix(0, p, p)

  # -- Block (1,1): Omega_11 = Sigma_y  [not multiplied by T here; T enters
  #    in estimate_gmm_weight_matrix when it inverts to get W] --
  # Convention: the Newey-West Omega is NOT multiplied by T; it equals the
  # large-T limit of T * Var(sqrt(T) sample_moment).  The cumulant loglik
  # uses -T/2 * delta' W delta, so W = Omega^{-1} is correct.
  # The per-observation score for order-1 is y_t, with long-run variance
  # = Omega_11 = sum_{h=-inf}^{inf} Cov(y_t, y_{t-h}) = sum_{h} Gamma(h).
  # For a stationary Gaussian: sum_{h>=0} Gamma(h) = (I - hx)^{-1} Q (I-hx)^{-T}
  # ... but we compute it as sum of lag-h autocovariances.
  #
  # Actually: Omega_11 = sum_{h=-inf}^{inf} Gamma(h) where Gamma(h) = Cov(y_t, y_{t-h}).
  # For h=0: Gamma(0) = Sigma_y.
  # For h>=1: Gamma(h) = Z hx^h Sigma_state Z'.
  # For h<=-1: Gamma(h) = Gamma(-h)'.
  # Sum: Omega_11 = Sigma_y + sum_{h=1}^{inf} [Gamma(h) + Gamma(h)'].

  if (p_11 > 0L) {
    Omega_11 <- Sigma_y   # h=0 contribution

    # Add lag contributions via recursion: Gamma_h = Z * (hx^h * Sigma_state) * Z'
    Gamma_state_h <- Sigma_state   # will be updated as hx^h * Sigma_state * (hx^h)'
    # Actually Gamma(h) = Z hx_state^h Sigma_state Z', so the recursion is:
    # Let S_h = hx_state^h Sigma_state (the propagated covariance);
    # then Gamma(h) = Z S_h Z', and S_{h+1} = hx_state S_h.
    S_h <- Sigma_state
    for (h in seq_len(max_lags)) {
      S_h     <- ghx_state %*% S_h
      Gamma_h <- Z %*% S_h %*% t(Z)
      contrib <- Gamma_h + t(Gamma_h)
      if (max(abs(contrib)) < tol) break
      Omega_11 <- Omega_11 + contrib
    }
    Omega[1L:p_11, 1L:p_11] <- Omega_11
  }

  # -- Block (2,2): Omega_22 = long-run covariance of vec(yc_t yc_t') --
  # Isserlis theorem for Gaussian:
  #   A(h)_{(ij),(kl)} = Gamma(h)[i,k]*Gamma(h)[j,l] + Gamma(h)[i,l]*Gamma(-h)[j,k]
  #                    = Gamma(h)[i,k]*Gamma(h)[j,l] + Gamma(h)[i,l]*Gamma(h)[k,j]
  # (using Gamma(-h) = Gamma(h)')
  # In matrix form: A(h) = Gamma(h) kron Gamma(h) + K (Gamma(h) kron Gamma(h)')
  # where K = K_{n_obs, n_obs} is the commutation matrix.
  #
  # For h=0: A(0) = (I+K)(Sigma_y kron Sigma_y).
  # For h>=1: A(h) + A(-h) =
  #   [Gamma_h kron Gamma_h + K(Gamma_h kron Gamma_h')] +
  #   [Gamma_h' kron Gamma_h' + K(Gamma_h' kron Gamma_h)]
  # Omega_22 = A(0) + sum_{h>=1} [A(h) + A(-h)]

  if (p_22 > 0L) {
    row0 <- p_11 + 1L
    row1 <- p_11 + p_22

    K_nn <- .commutation_matrix(n_obs)   # n_obs^2 x n_obs^2

    # h=0 contribution
    Syky <- kronecker(Sigma_y, Sigma_y)   # n_obs^2 x n_obs^2
    Omega_22 <- Syky + K_nn %*% Syky     # (I + K)(Sigma_y kron Sigma_y)

    # Lag contributions
    S_h <- Sigma_state
    for (h in seq_len(max_lags)) {
      S_h     <- ghx_state %*% S_h
      Gamma_h <- Z %*% S_h %*% t(Z)
      Gamma_ht <- t(Gamma_h)   # Gamma(-h)

      # A(h): Gamma_h kron Gamma_h + K (Gamma_h kron Gamma_h')
      GkG    <- kronecker(Gamma_h, Gamma_h)
      GkGt   <- kronecker(Gamma_h, Gamma_ht)
      Ah     <- GkG + K_nn %*% GkGt

      # A(-h): Gamma_ht kron Gamma_ht + K (Gamma_ht kron Gamma_h)
      GtkGt  <- kronecker(Gamma_ht, Gamma_ht)
      GtkG   <- kronecker(Gamma_ht, Gamma_h)
      Amh    <- GtkGt + K_nn %*% GtkG

      contrib <- Ah + Amh
      if (max(abs(contrib)) < tol) break
      Omega_22 <- Omega_22 + contrib
    }

    # Symmetrize (numerical drift)
    Omega_22 <- (Omega_22 + t(Omega_22)) * 0.5
    Omega[row0:row1, row0:row1] <- Omega_22
  }

  # Cross-block Omega_12 = 0 (Gaussian: odd-even cumulant cross-covariance
  # vanishes — mean and variance are independent for Gaussian).

  # ---- Ridge regularization ----
  diag_max <- max(diag(Omega))
  if (!is.finite(diag_max) || diag_max <= 0) diag_max <- 1
  delta_ridge <- ridge * diag_max
  Omega_reg   <- Omega + delta_ridge * diag(p)

  # ---- Condition number ----
  cond_num <- tryCatch({
    ev <- eigen(Omega_reg, only.values = TRUE, symmetric = TRUE)$values
    max(ev) / min(ev)
  }, error = function(e) Inf)

  if (!is.finite(cond_num) || cond_num > 1e12) {
    warning(
      ".analytic_gmm_weight_matrix: condition number of Omega_reg is ",
      if (is.finite(cond_num)) format(cond_num, scientific = TRUE) else "Inf",
      ".  Increase 'ridge' or check model stationarity.", call. = FALSE)
  }

  # ---- Invert via Cholesky ----
  ch <- tryCatch(chol(Omega_reg), error = function(e) {
    stop(".analytic_gmm_weight_matrix: Cholesky of Omega_reg failed. ",
         "Increase 'ridge'.  Original error: ", e$message)
  })
  W <- chol2inv(ch)
  W <- (W + t(W)) * 0.5

  attr(W, "p")                <- p
  attr(W, "T_obs")            <- NA_integer_
  attr(W, "bandwidth")        <- NA_integer_
  attr(W, "condition_number") <- cond_num
  attr(W, "method")           <- "analytic"
  W
}


# ============================================================================
# GMM optimal weight matrix estimator
# ============================================================================

#' Estimate the GMM optimal weight matrix for cumulant-moment matching
#'
#' Computes the inverse long-run covariance matrix of the empirical moment
#' conditions using a Newey–West HAC estimator.  The returned matrix \code{W}
#' can be passed as \code{weight_matrix} to \code{make_log_posterior_cumulant()}
#' or \code{.cumulant_loglik()} to implement the two-step efficient GMM
#' estimator.
#'
#' @section Two-step GMM workflow:
#' \enumerate{
#'   \item Optimize with the default \code{cumulant_weight = "identity"} to
#'     obtain a first-step estimate \code{theta_hat}.
#'   \item Evaluate \code{W <- estimate_gmm_weight_matrix(data, dr, model,
#'     params, obs_vars)} at \code{theta_hat}.
#'   \item Re-optimize with \code{weight_matrix = W} to get the efficient
#'     second-step estimate.
#' }
#'
#' @section Finite-sample caution:
#' Optimal (efficient) GMM weighting is asymptotically efficient but can
#' \strong{worsen} finite-sample behavior when the number of moments \eqn{p}
#' is large relative to the sample size \eqn{T}.  With \code{orders = 1:4}
#' and \eqn{n_{obs} = 4} observables, \eqn{p \approx 340}, so the sample
#' long-run covariance is severely rank-deficient for typical macro time series
#' (\eqn{T \lesssim 250}).  In this over-identified regime the "optimal"
#' weight is a known finite-sample pathology (the too-many-moments problem;
#' see Donald & Newey 2001).  \strong{Recommendation}: use
#' \code{orders = 1:2} (mean + variance only, \eqn{p \leq n_{obs}^2 +
#' n_{obs}}) with large \eqn{T} when precision weighting is desired, or apply
#' additional moment-selection.  A ridge penalty (\code{ridge}) regularizes
#' the inversion but does not solve the underlying over-identification issue.
#'
#' @param data      T × n_obs numeric matrix of observations.
#' @param dr        Decision rule object (output of \code{solve_perturbation()}).
#' @param model     dynhr_mod from \code{parse_mod()}.
#' @param params    Named numeric vector of structural parameters.
#' @param obs_vars  Character vector of observed variable names
#'   (column names of \code{data}).
#' @param orders    Integer vector of cumulant orders to include (default
#'   \code{1:4}).  Must match the \code{orders} argument used in
#'   \code{.cumulant_loglik()}.
#' @param method    \code{"newey_west"} (default): Newey–West HAC long-run
#'   covariance of per-observation moment contributions.
#'   \code{"analytic"}: block-diagonal analytic long-run covariance for
#'   \code{orders} \eqn{\subseteq \{1, 2\}} (Gaussian state-space formula using
#'   model-implied autocovariances).  Errors for \code{orders} including 3 or 4
#'   (requires order-6/8 cumulants not yet implemented).
#' @param bandwidth Integer or \code{NULL}.  Newey–West lag truncation \eqn{L}.
#'   \code{NULL} (default) uses the rule \eqn{L = \lfloor 4 (T/100)^{2/9} \rfloor}.
#' @param ridge     Ridge regularization fraction (default \code{1e-6}).
#'   The regularized covariance is
#'   \eqn{\hat\Omega_{reg} = \hat\Omega_{NW} + \delta \cdot \max_i(\hat\Omega_{ii}) \cdot I_p},
#'   where \eqn{\delta} = \code{ridge}.
#'
#' @return A \eqn{p \times p} symmetric positive-definite numeric matrix
#'   \eqn{W = \hat\Omega_{reg}^{-1}}.  Attributes: \code{p} (number of
#'   moments), \code{T_obs} (sample size), \code{bandwidth} (lag truncation
#'   used), \code{condition_number} (condition of \eqn{\hat\Omega_{reg}}).
#'
#' @seealso \code{\link{make_log_posterior_cumulant}}, \code{.cumulant_loglik}
#' @export
estimate_gmm_weight_matrix <- function(data, dr, model, params, obs_vars,
                                        orders = 1:4,
                                        method = c("newey_west", "analytic"),
                                        bandwidth = NULL,
                                        ridge = 1e-6) {
  method <- match.arg(method)

  if (method == "analytic" && any(orders >= 3L)) {
    stop("method = \"analytic\" with orders >= 3 is not yet implemented: ",
         "computing the asymptotic covariance of skewness/kurtosis moments ",
         "requires order-6/8 cumulants which dynhr does not yet have.  ",
         "Use method = \"newey_west\", or restrict to orders = 1:2 for the ",
         "analytic path.  For orders 3-4, use method = \"newey_west\" instead.")
  }

  if (method == "analytic") {
    return(.analytic_gmm_weight_matrix(dr, model, params, obs_vars,
                                       orders, ridge))
  }

  # ---- Validate inputs ----
  if (!is.matrix(data) || !is.numeric(data))
    stop("'data' must be a numeric matrix")
  T_obs <- nrow(data)
  n_obs <- length(obs_vars)

  if (T_obs < 2L)
    stop("'data' must have at least 2 rows")
  if (n_obs < 1L)
    stop("'obs_vars' must be non-empty")
  if (!all(obs_vars %in% colnames(data)))
    stop("all 'obs_vars' must be column names of 'data'")

  # ---- Build model moment vector via .build_moment_vector() ----
  # This reuses the EXACT assembly logic from .cumulant_loglik(), guaranteeing
  # the same moment ordering.
  m_model <- .build_moment_vector(dr, model, params, obs_vars, orders,
                                   me_variance = 0)
  if (is.null(m_model))
    stop("Could not compute model moment vector: non-finite values in DR")

  p <- length(m_model)

  # ---- Build T × p matrix G of per-observation raw moment contributions ----
  # For the Newey-West long-run covariance we need the per-obs "influence
  # function" contribution to each sample cumulant.  These match exactly how
  # sample_cumulants() computes the aggregate moments:
  #
  #   order 1: contribution at t is y_t  (raw, not demeaned)
  #   order 2: contribution at t is vec(yc_t yc_t'), yc_t = y_t - sample_mean
  #   order 3: contribution at t is the n_obs × n_obs^2 raw product using yc_t,
  #             flattened in the same col-major order as sc$c3 (i.e. row=i,
  #             col=(j-1)*n_obs+k, value = yc_t[i]*yc_t[j]*yc_t[k])
  #   order 4: analogous n_obs × n_obs^3 raw product using yc_t
  #
  # The moment conditions are g_t = raw_t - m_model (constant in t), where
  # m_model is the same for all t.  Subtracting m_model gives G; the
  # long-run cov of g_t is what Newey-West estimates.

  Y  <- data[, obs_vars, drop = FALSE]   # T × n_obs
  mu <- colMeans(Y)                       # n_obs sample mean
  Yc <- sweep(Y, 2L, mu, FUN = "-")      # T × n_obs demeaned

  G_list <- vector("list", length(orders))

  for (oi in seq_along(orders)) {
    ord <- orders[oi]
    if (ord == 1L) {
      # Raw order-1 contribution at each t: y_t (not demeaned; mean is m_emp)
      G_list[[oi]] <- Y  # T × n_obs
    } else if (ord == 2L) {
      # Raw order-2: vec(yc_t yc_t') — demeaned, matching sc$var_cov
      G2 <- matrix(0, T_obs, n_obs * n_obs)
      for (t in seq_len(T_obs)) {
        yct <- Yc[t, ]
        G2[t, ] <- as.numeric(outer(yct, yct))
      }
      G_list[[oi]] <- G2
    } else if (ord == 3L) {
      # Raw order-3: n_obs × n_obs^2 product using yc_t, vec'd row-major in (i),
      # col-major in (j,k): col = (j-1)*n_obs + k — matches sc$c3 layout
      G3 <- matrix(0, T_obs, n_obs * n_obs * n_obs)
      for (t in seq_len(T_obs)) {
        yct <- Yc[t, ]
        m3t <- matrix(0, n_obs, n_obs * n_obs)
        for (a in seq_len(n_obs)) {
          for (b in seq_len(n_obs)) {
            for (cc in seq_len(n_obs)) {
              m3t[a, (b - 1L) * n_obs + cc] <- yct[a] * yct[b] * yct[cc]
            }
          }
        }
        G3[t, ] <- as.numeric(m3t)
      }
      G_list[[oi]] <- G3
    } else if (ord == 4L) {
      # Raw order-4: n_obs × n_obs^3 product using yc_t, vec'd — matches sc$c4
      G4 <- matrix(0, T_obs, n_obs * n_obs^3)
      for (t in seq_len(T_obs)) {
        yct <- Yc[t, ]
        m4t <- matrix(0, n_obs, n_obs^3)
        for (a in seq_len(n_obs)) {
          col <- 0L
          for (b in seq_len(n_obs)) {
            for (cc in seq_len(n_obs)) {
              for (d in seq_len(n_obs)) {
                col <- col + 1L
                m4t[a, col] <- yct[a] * yct[b] * yct[cc] * yct[d]
              }
            }
          }
        }
        G4[t, ] <- as.numeric(m4t)
      }
      G_list[[oi]] <- G4
    }
  }

  # Bind all order contributions column-wise: T × p
  G_raw <- do.call(cbind, G_list)

  # Moment conditions: g_t = raw_t - m_model (center on model prediction)
  G <- sweep(G_raw, 2L, m_model, FUN = "-")

  # ---- Newey-West bandwidth ----
  if (is.null(bandwidth)) {
    bandwidth <- floor(4 * (T_obs / 100)^(2 / 9))
  }
  bandwidth <- as.integer(bandwidth)

  # ---- Warn when p > T (over-identified, rank-deficient regime) ----
  if (p > T_obs) {
    warning(
      "estimate_gmm_weight_matrix: number of moments (p = ", p, ") exceeds ",
      "sample size (T = ", T_obs, ").  The long-run covariance matrix is ",
      "rank-deficient; the ridge-regularized inverse may be unreliable.  ",
      "Consider using fewer cumulant orders (e.g. orders = 1:2) or a ",
      "larger sample.", call. = FALSE)
  }

  # ---- Newey-West HAC long-run covariance ----
  # Gamma_0 = (1/T) G' G
  # Gamma_j = (1/T) sum_{t=j+1}^{T} g_t g_{t-j}'  for j >= 1
  # Omega_NW = Gamma_0 + sum_{j=1}^{L} (1 - j/(L+1)) (Gamma_j + Gamma_j')

  Gamma0 <- crossprod(G) / T_obs   # p × p

  Omega_NW <- Gamma0
  if (bandwidth > 0L) {
    for (j in seq_len(bandwidth)) {
      w_j <- 1 - j / (bandwidth + 1)
      # Gamma_j: columns of G are demeaned; rows t=j+1,...,T paired with t-j
      Gj_lead  <- G[(j + 1L):T_obs, , drop = FALSE]
      Gj_lag   <- G[1L:(T_obs - j), , drop = FALSE]
      Gamma_j  <- crossprod(Gj_lead, Gj_lag) / T_obs
      Omega_NW <- Omega_NW + w_j * (Gamma_j + t(Gamma_j))
    }
  }

  # Symmetrize numerically
  Omega_NW <- (Omega_NW + t(Omega_NW)) * 0.5

  # ---- Ridge regularization ----
  diag_max <- max(diag(Omega_NW))
  if (!is.finite(diag_max) || diag_max <= 0) {
    # Fall back to a safe diagonal scale if Omega is degenerate
    diag_max <- 1
    warning("estimate_gmm_weight_matrix: Omega_NW has non-positive diagonal; ",
            "using ridge = ", ridge, " * I_p.  Check data and model.",
            call. = FALSE)
  }
  delta_ridge <- ridge * diag_max
  Omega_reg <- Omega_NW + delta_ridge * diag(p)

  # ---- Check condition number and warn if extreme ----
  cond_num <- tryCatch({
    ev <- eigen(Omega_reg, only.values = TRUE, symmetric = TRUE)$values
    max(ev) / min(ev)
  }, error = function(e) Inf)

  if (!is.finite(cond_num) || cond_num > 1e12) {
    warning(
      "estimate_gmm_weight_matrix: condition number of Omega_reg is ",
      if (is.finite(cond_num)) format(cond_num, scientific = TRUE)
      else "Inf",
      ".  The weight matrix W = Omega_reg^{-1} may be numerically unreliable. ",
      "Increase 'ridge' or reduce 'orders'.", call. = FALSE)
  }

  # ---- Invert via Cholesky for stability ----
  ch <- tryCatch(chol(Omega_reg), error = function(e) {
    stop("estimate_gmm_weight_matrix: Cholesky decomposition of Omega_reg ",
         "failed (matrix not positive definite even after ridge). ",
         "Increase 'ridge' or reduce 'orders'.  Original error: ", e$message)
  })
  W <- chol2inv(ch)

  # Symmetrize result
  W <- (W + t(W)) * 0.5

  attr(W, "p")                <- p
  attr(W, "T_obs")            <- T_obs
  attr(W, "bandwidth")        <- bandwidth
  attr(W, "condition_number") <- cond_num
  W
}
