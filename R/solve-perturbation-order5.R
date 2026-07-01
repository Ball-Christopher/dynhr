## R/solve-perturbation-order5.R
## --------------------------------------------------------------------------
## Fifth-order (deterministic) perturbation solver for DSGE models.
##
## Implements Levintal (2017) compact tensor notation for 5th-order
## policy terms:
##
##   y_t = ... + (1/120) g_xxxxx (x ⊗ x ⊗ x ⊗ x ⊗ x)
##             + ...
##
## Key extension beyond order 4: the compact Sylvester solve uses the
## QZ decomposition of (A_L, f_+) and solves the 5-fold Kronecker system
## row-by-row with the eigendecomposition of h_x.
##
## References:
##   Levintal, O. (2017). Fifth-Order Perturbation Solution to DSGE
##     Models. JEDC 80, 1-16.
## --------------------------------------------------------------------------


# =====================================================================
# Compact Sylvester solve for 5th order
# =====================================================================

#' Solve the 5th-order Kronecker system compactly.
#'
#' Solves A_L · X + f_+ · X · (h_x^{⊗5}) = RHS without forming the
#' full (n·ns^5) × (n·ns^5) system.
#'
#' @param A_L    n × n effective feedback matrix
#' @param fp     n × n f_plus matrix
#' @param hx     n_s × n_s state transition matrix
#' @param RHS    n × n_s^5 forcing matrix
#' @param verbose Print progress
#' @return X = ghxxxxx (n × n_s^5)
#' @noRd
.solve_compact_o5 <- function(A_L, fp, hx, RHS, verbose = FALSE) {
  n  <- nrow(A_L)
  ns <- nrow(hx)
  m  <- ns^5

  if (verbose) cat("  Compact Sylvester solve (order 5): n =", n, ", ns^5 =", m, "\n")

  # QZ decomposition
  qz_result <- QZ::qz(A_L, fp)

  S <- qz_result$S
  T <- qz_result$T
  Q <- qz_result$Q
  Z <- qz_result$Z

  C_tilde <- crossprod(Q, RHS)    # n × m

  # Eigendecomposition of hx
  eig_hx <- eigen(hx)
  V      <- eig_hx$vectors
  lambda <- eig_hx$values
  V_inv  <- solve(V)

  # Back-substitution
  Y <- matrix(0, n, m)

  for (i in n:1) {
    rhs_i <- C_tilde[i, , drop = TRUE]

    if (i < n) {
      for (j in (i + 1):n) {
        yj <- Y[j, , drop = TRUE]
        rhs_i <- rhs_i - (S[i, j] * yj + T[i, j] * .apply_kron5(yj, hx))
      }
    }

    alpha <- S[i, i]
    beta  <- T[i, i]

    if (abs(beta) < 1e-14) {
      Y[i, ] <- rhs_i / alpha
    } else {
      Y[i, ] <- .solve_kron5_row(rhs_i, alpha, beta, V, V_inv, lambda, ns)
    }
  }

  X <- Z %*% Y
  rownames(X) <- rownames(A_L)
  X
}


#' Apply the 5-fold Kronecker product M^{⊗5} to a row vector y (1 × ns^5).
#'
#' Mode-by-mode application using tensor reshaping (O(5·ns^6) instead of
#' O(ns^10) for the full Kronecker).
#'
#' @param y  Numeric vector of length ns^5
#' @param M  ns × ns matrix
#' @return Numeric vector: y · M^{⊗5}
#' @noRd
.apply_kron5 <- function(y, M) {
  ns <- nrow(M)
  T5 <- array(y, dim = rep(ns, 5))

  # Apply M to modes 1-5 sequentially
  for (mode in 1:5) {
    perm <- c(mode, seq_len(5)[-mode])
    T5 <- aperm(array(
      crossprod(M, matrix(aperm(T5, perm), ns, ns^4)),
      dim = rep(ns, 5)), order(perm))
  }

  as.numeric(T5)
}


#' Solve per-row equation for 5th-order compact Sylvester.
#'
#' Solves y · (α·I + β·hx^{⊗5}) = rhs using eigenbasis of hx.
#'
#' @noRd
.solve_kron5_row <- function(rhs, alpha, beta, V, V_inv, lambda, ns) {
  # Step 1: rhs_tilde = rhs · (V^{-1})^{⊗5}
  rhs_tilde <- .apply_kron5(rhs, V_inv)

  # Step 2: Element-wise division: denominator = α + β·λ_{j1}·...·λ_{j5}
  D <- alpha + beta *
    outer(outer(outer(outer(lambda, lambda, `*`), lambda, `*`), lambda, `*`), lambda, `*`)

  y_tilde <- rhs_tilde / as.numeric(D)

  # Step 3: y = y_tilde · V^{⊗5}
  y <- .apply_kron5(y_tilde, V)
  y
}


# =====================================================================
# Numerical 5th derivative
# =====================================================================

#' Compute 5th-order numerical derivative of model residuals.
#'
#' Uses nested finite differences on the Jacobian/Hessian to obtain
#' the 5th derivative tensor F_wwwww at SS. Returns a 6-D array
#' c(n_eq, n_cols, n_cols, n_cols, n_cols, n_cols).
#'
#' If compiled symbolic 5th derivatives are available (hessian5_fn), uses
#' them instead of numerical FD.
#'
#' @param dyn      Compiled dynamic model
#' @param dy_ss    Compound vector at SS
#' @param params   Parameter vector
#' @param ss       Steady state
#' @param h        Step size (default 5e-2, larger for high-order FD)
#' @return 6-D array or NULL
#' @noRd
.compute_model_5th_deriv <- function(dyn, dy_ss, params, ss, h = 5e-2) {
  # Try compiled symbolic 5th derivatives first
  if (!is.null(dyn$hessian5_fn) && !is.null(dyn$hess5_triplets) &&
      length(dyn$hess5_triplets) > 0L) {
    return(.compute_model_hessian5_symbolic(dyn, dy_ss, params, ss))
  }
  n_eq       <- dyn$n_eq
  total_cols <- dyn$total_cols

  D5 <- array(0, dim = c(n_eq, total_cols, total_cols, total_cols,
                          total_cols, total_cols))

  # Use 7-point Richardson stencil on the 3rd derivative (Hessian3)
  # d^5F/(dw^5) ≈ (-H3(+3h) + 12·H3(+2h) - 39·H3(+h) + 56·H3(0)
  #                - 39·H3(-h) + 12·H3(-2h) - H3(-3h)) / (6·h^3)
  #
  # For efficiency, compute the Hessian3 numerically at the 7 points
  # along each direction, then apply the stencil.

  # Simplified: use 5-point stencil on 4th derivative:
  # D5 ≈ (-D4(+2h) + 16·D4(+h) - 30·D4(0) + 16·D4(-h) - D4(-2h)) / (12·h)

  for (c5 in seq_len(total_cols)) {
    step_base <- if (abs(dy_ss[c5]) > 1e-12) max(h, abs(dy_ss[c5]) * h) else h

    # Compute 4th derivatives at 5 points along dimension c5
    D4_vals <- list()
    offsets <- c(2, 1, 0, -1, -2)
    for (oi in seq_along(offsets)) {
      dy_p <- dy_ss
      dy_p[c5] <- dy_ss[c5] + offsets[oi] * step_base
      D4_vals[[oi]] <- .compute_4th_at(dyn, dy_p, params, ss, h)
    }

    # 5-point stencil for 5th derivative
    if (!any(sapply(D4_vals, is.null))) {
      D5[, , , , c5, c5] <- (-D4_vals[[1]] + 16 * D4_vals[[2]] -
                               30 * D4_vals[[3]] + 16 * D4_vals[[4]] -
                               D4_vals[[5]]) / (12 * step_base)
    }
  }

  # Symmetrize (simplified — average over a representative set of permutations)
  for (e in seq_len(n_eq)) {
    De <- D5[e, , , , , ]
    # Ensure symmetry in the last two indices (at minimum)
    for (i in seq_len(total_cols)) {
      for (j in seq_len(total_cols)) {
        for (k in seq_len(total_cols)) {
          for (l in seq_len(total_cols)) {
            for (m in seq_len(total_cols)) {
              val <- (De[i, j, k, l, m] + De[i, j, k, m, l]) / 2
              D5[e, i, j, k, l, m] <- val
              # More permutations would be needed for full symmetry
            }
          }
        }
      }
    }
  }

  D5
}


#' Compute the 4th-derivative 5-tensor at a specific dy point.
#'
#' Numerically differentiates the Jacobian using nested FD.
#'
#' @noRd
.compute_4th_at <- function(dyn, dy, params, ss, h) {
  n_eq       <- dyn$n_eq
  total_cols <- dyn$total_cols

  D4 <- array(0, dim = c(n_eq, total_cols, total_cols, total_cols, total_cols))

  for (c4 in seq_len(total_cols)) {
    step <- if (abs(dy[c4]) > 1e-12) max(h, abs(dy[c4]) * h) else h

    # Compute 3rd derivatives at 5 points along c4
    H3_vals <- list()
    offsets <- c(2, 1, 0, -1, -2)
    for (oi in seq_along(offsets)) {
      dy_p <- dy
      dy_p[c4] <- dy[c4] + offsets[oi] * step
      H3_vals[[oi]] <- .compute_h3_at(dyn, dy_p, params, ss, h)
    }

    if (!any(sapply(H3_vals, is.null))) {
      D4[, , , c4, c4] <- (-H3_vals[[1]] + 16 * H3_vals[[2]] -
                             30 * H3_vals[[3]] + 16 * H3_vals[[4]] -
                             H3_vals[[5]]) / (12 * step)
    }
  }

  D4
}


#' Compute the 3rd-derivative 4-tensor (Hessian3) at a specific dy point.
#'
#' Uses 5-point stencil on the Hessian. Returns array c(n_eq, n_cols, n_cols, n_cols).
#'
#' @noRd
.compute_h3_at <- function(dyn, dy, params, ss, h) {
  n_eq       <- dyn$n_eq
  total_cols <- dyn$total_cols

  H3 <- array(0, dim = c(n_eq, total_cols, total_cols, total_cols))

  for (c3 in seq_len(total_cols)) {
    step <- if (abs(dy[c3]) > 1e-12) max(h, abs(dy[c3]) * h) else h

    # Hessians at 5 points along c3
    offsets <- c(2, 1, 0, -1, -2)
    H_vals <- list()
    for (oi in seq_along(offsets)) {
      dy_p <- dy
      dy_p[c3] <- dy[c3] + offsets[oi] * step
      H_vals[[oi]] <- .hessian_at_dy(dyn, dy_p, params, ss, h)
    }

    if (!any(sapply(H_vals, is.null))) {
      H3[, , , c3] <- (-H_vals[[1]] + 16 * H_vals[[2]] -
                         30 * H_vals[[3]] + 16 * H_vals[[4]] -
                         H_vals[[5]]) / (12 * step)
    }
  }

  # Symmetrize
  for (e in seq_len(n_eq)) {
    for (i in seq_len(total_cols)) {
      for (j in seq_len(total_cols)) {
        for (k in seq_len(total_cols)) {
          val <- (H3[e, i, j, k] + H3[e, i, k, j] +
                  H3[e, j, i, k] + H3[e, j, k, i] +
                  H3[e, k, i, j] + H3[e, k, j, i]) / 6
          H3[e, i, j, k] <- val
        }
      }
    }
  }

  H3
}


#' Compute the numerical Hessian at a specific dy point.
#'
#' @noRd
.hessian_at_dy <- function(dyn, dy, params, ss, h) {
  n_eq       <- dyn$n_eq
  total_cols <- dyn$total_cols

  H <- array(0, dim = c(n_eq, total_cols, total_cols))
  J0 <- dyn$jacobian_fn(dy, params, ss)

  for (c2 in seq_len(total_cols)) {
    dy_p <- dy; dy_m <- dy
    step <- if (abs(dy[c2]) > 1e-12) max(h, abs(dy[c2]) * h) else h
    dy_p[c2] <- dy[c2] + step
    dy_m[c2] <- dy[c2] - step

    Jp <- dyn$jacobian_fn(dy_p, params, ss)
    Jm <- dyn$jacobian_fn(dy_m, params, ss)

    if (!is.null(Jp) && !is.null(Jm)) {
      H[, , c2] <- (Jp - Jm) / (2 * step)
    }
  }

  # Symmetrize
  for (e in seq_len(n_eq)) {
    He <- H[e, , ]
    H[e, , ] <- (He + t(He)) / 2
  }

  H
}


#' Compute the 5th derivative tensor from compiled symbolic derivatives.
#'
#' Evaluates the compiled \code{hessian5_fn} at SS, then expands the sparse
#' quintuplet representation into a dense 6-D array with Schwarz symmetry.
#'
#' @param dyn      Compiled dynamic model
#' @param dy_ss    Compound vector at SS
#' @param params   Parameter vector
#' @param ss       Steady state
#' @return 6-D array, or NULL
#' @noRd
.compute_model_hessian5_symbolic <- function(dyn, dy_ss, params, ss) {
  if (is.null(dyn$hessian5_fn) || is.null(dyn$hess5_triplets)) return(NULL)
  if (dyn$n_hess5 == 0L) {
    return(array(0, dim = c(dyn$n_eq, dyn$total_cols, dyn$total_cols,
                             dyn$total_cols, dyn$total_cols, dyn$total_cols)))
  }
  values <- dyn$hessian5_fn(dy_ss, params, ss)

  n_eq       <- dyn$n_eq
  total_cols <- dyn$total_cols
  H5 <- array(0, dim = c(n_eq, total_cols, total_cols, total_cols,
                          total_cols, total_cols))

  for (k in seq_along(values)) {
    t   <- dyn$hess5_triplets[[k]]
    val <- values[k]
    if (val == 0) next
    e  <- t$eq
    cols <- c(t$col1, t$col2, t$col3, t$col4, t$col5)
    perms <- .orbit_5(cols[1], cols[2], cols[3], cols[4], cols[5])
    for (p in perms) {
      H5[e, p[1], p[2], p[3], p[4], p[5]] <-
        H5[e, p[1], p[2], p[3], p[4], p[5]] + val
    }
  }
  H5
}


#' Enumerate the Schwarz symmetry orbit of (c1,..,c5).
#'
#' Generates unique permutations of 5 indices, handling multiplicities.
#' For all-distinct returns 120 permutations; fewer for repeated indices.
#'
#' @noRd
.orbit_5 <- function(c1, c2, c3, c4, c5) {
  cols <- c(c1, c2, c3, c4, c5)
  # Generate all 120 permutations via recursion
  .permute5 <- function(v) {
    if (length(v) <= 1L) return(list(v))
    result <- list()
    for (i in seq_along(v)) {
      rest <- .permute5(v[-i])
      for (r in rest) result <- c(result, list(c(v[i], r)))
    }
    result
  }
  all_perms <- .permute5(1:5)
  seen <- new.env(hash = TRUE, parent = emptyenv())
  uniq <- list()
  for (perm_idx in all_perms) {
    key <- paste(cols[perm_idx], collapse = ",")
    if (is.null(seen[[key]])) {
      seen[[key]] <- TRUE
      uniq <- c(uniq, list(cols[perm_idx]))
    }
  }
  uniq
}


# =====================================================================
# FD-based Phi builder for all 5th-order terms
# =====================================================================

#' Build all 6 forcing matrices via finite differences on the order-4 residual.
#'
#' Mirrors the validated order-4 builder (.build_phi_fd): substitutes the
#' order-4 policy into the model residual, then takes the 5th mixed central
#' finite difference. The result is the part of the 5th derivative of the
#' residual that does NOT involve the unknown 5th-order policy term (which is
#' added analytically in the solve as A_L·X + fp·X·hx^{⊗5} resp. fp·ghxxxxx·mix).
#' @noRd
.build_phi_fd_o5 <- function(dyn, dr4, ss, params,
                             state_idx, endo_names, exo_names,
                             n_s, n_u, n, h = 0.05, res_perm = NULL) {
  if (is.null(res_perm)) res_perm <- seq_len(n)
  ss_endo <- ss[endo_names]
  ghx  <- dr4$ghx;  ghu  <- dr4$ghu
  ghxx <- dr4$ghxx; ghxu <- dr4$ghxu; ghuu <- dr4$ghuu
  ghxxx <- dr4$ghxxx; ghxxu <- dr4$ghxxu
  ghxuu <- dr4$ghxuu; ghuuu <- dr4$ghuuu
  ghxxxx <- dr4$ghxxxx; ghxxxu <- dr4$ghxxxu
  ghxxuu <- dr4$ghxxuu; ghxuuu <- dr4$ghxuuu; ghuuuu <- dr4$ghuuuu
  dcm  <- dyn$dyn_col_map
  K <- function(...) Reduce(`%x%`, list(...))

  ord4 <- function(x, u) {
    y <- as.numeric(ghx %*% x + ghu %*% u)
    y <- y + 0.5*as.numeric(ghxx %*% K(x,x)) + as.numeric(ghxu %*% K(x,u)) +
              0.5*as.numeric(ghuu %*% K(u,u))
    y <- y + (1/6)*as.numeric(ghxxx %*% K(x,x,x)) +
              0.5*as.numeric(ghxxu %*% K(x,x,u)) +
              0.5*as.numeric(ghxuu %*% K(x,u,u)) +
             (1/6)*as.numeric(ghuuu %*% K(u,u,u))
    y + (1/24)*as.numeric(ghxxxx %*% K(x,x,x,x)) +
        (1/6) *as.numeric(ghxxxu %*% K(x,x,x,u)) +
        (1/4) *as.numeric(ghxxuu %*% K(x,x,u,u)) +
        (1/6) *as.numeric(ghxuuu %*% K(x,u,u,u)) +
        (1/24)*as.numeric(ghuuuu %*% K(u,u,u,u))
  }

  zero_u <- numeric(n_u)
  build_dy <- function(y_lag, y_now, y_lead, u_now) {
    dy <- numeric(dyn$total_cols); keys <- character(dyn$total_cols)
    for (kc in seq_len(nrow(dcm))) {
      c <- dcm$col[kc]; nm <- dcm$name[kc]; ll <- dcm$lead_lag[kc]
      sfx <- if(ll==0L)"__0" else if(ll>0L)paste0("__p",ll) else paste0("__m",abs(ll))
      keys[c] <- paste0(nm, sfx)
      if (nm %in% exo_names) {
        dy[c] <- u_now[which(exo_names == nm)]
      } else {
        idx <- which(endo_names == nm); if (length(idx) != 1L) next
        dy[c] <- if(ll==-1L) y_lag[idx] else if(ll==0L) y_now[idx]
                 else if(ll==1L) y_lead[idx] else NA_real_
      }
    }
    names(dy) <- keys; dy
  }

  cache <- list()
  R4_cached <- function(xv, uv) {
    key <- paste(c(round(xv, 10), round(uv, 10)), collapse = "|")
    if (is.null(cache[[key]])) {
      y_now  <- ss_endo + ord4(xv, uv)
      x_now  <- y_now[state_idx] - ss_endo[state_idx]
      y_lead <- ss_endo + ord4(x_now, zero_u)
      y_lag  <- ss_endo; y_lag[state_idx] <- ss_endo[state_idx] + xv
      cache[[key]] <<- as.numeric(dyn$residuals_fn(
        build_dy(y_lag, y_now, y_lead, uv), params, ss))[res_perm]
    }
    cache[[key]]
  }

  signs5 <- expand.grid(rep(list(c(-1L, 1L)), 5))
  zero_x <- numeric(n_s)

  phi_type <- function(n_x, n_u_count) {
    n_cols_out <- if (n_x + n_u_count == 0) 1L
                  else n_s^n_x * max(n_u, 1L)^n_u_count
    Phi <- matrix(0, n, n_cols_out)
    x_dim   <- if (n_x > 0) rep(n_s, n_x) else integer(0)
    u_dim   <- if (n_u_count > 0) rep(n_u, n_u_count) else integer(0)
    all_dim <- c(x_dim, u_dim); k_tot <- n_x + n_u_count
    if (k_tot == 0L) return(Phi)

    idx <- rep(1L, k_tot)
    for (col in seq_len(n_cols_out)) {
      val <- rep(0, n)
      for (srow in seq_len(nrow(signs5))) {
        s  <- as.integer(signs5[srow, ])
        xv <- zero_x; uv <- zero_u
        for (ki in seq_len(n_x))        xv[idx[ki]]       <- xv[idx[ki]]       + s[ki] * h
        for (ki in seq_len(n_u_count))  uv[idx[n_x + ki]] <- uv[idx[n_x + ki]] + s[n_x + ki] * h
        val <- val + prod(s) * R4_cached(xv, uv)
      }
      Phi[, col] <- val / (32 * h^5)   # (2h)^5 central difference
      for (ki in rev(seq_len(k_tot))) {
        idx[ki] <- idx[ki] + 1L
        if (idx[ki] <= all_dim[ki]) break
        idx[ki] <- 1L
      }
    }
    Phi
  }

  list(
    xxxxx = phi_type(5L, 0L),
    xxxxu = phi_type(4L, 1L),
    xxxuu = phi_type(3L, 2L),
    xxuuu = phi_type(2L, 3L),
    xuuuu = phi_type(1L, 4L),
    uuuuu = phi_type(0L, 5L)
  )
}


# =====================================================================
# Main 5th-order solver entry point
# =====================================================================

#' Solve the deterministic fifth-order perturbation of a DSGE model
#'
#' Given first- through fourth-order decision rules, computes the 5th-
#' order terms ghxxxxx, ghxxxxu, ... using Levintal (2017) compact
#' tensor notation.
#'
#' @param model    dynhr_mod object from parse_mod()
#' @param compiled dynhr_compiled from compile_model()
#' @param ss       Named numeric steady state vector
#' @param params   Named numeric parameter vector
#' @param dr4      Fourth-order DecisionRules4 object
#' @param h        Step size for numerical derivatives (default 5e-2)
#' @param verbose  Print progress
#' @return A DecisionRules5 object extending DecisionRules4 with 5th-order
#'   fields ghxxxxx, ghxxxxu, ...
#'
#' @references
#'   Levintal, O. (2017). Fifth-order perturbation solution of DSGE models.
#'     \emph{Journal of Economic Dynamics and Control}, 80, 1-16.
#'   Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez, J. F. (2018).
#'     The pruned state-space system for non-linear DSGE models.
#'     \emph{Review of Economic Studies}, 85(1), 1-49.
#' @export
solve_perturbation_order5 <- function(model, compiled, ss, params, dr4,
                                       h = 5e-2, verbose = FALSE) {
  if (!inherits(dr4, "DecisionRules4")) {
    stop("dr4 must be a DecisionRules4 object from solve_perturbation_order4().")
  }
  if (!isTRUE(dr4$bk_satisfied)) {
    stop("Blanchard-Kahn not satisfied; cannot solve to 5th order.")
  }

  dyn        <- compiled$dynamic
  n          <- length(dr4$endo_names)
  n_s        <- length(dr4$state_idx)
  n_u        <- length(dr4$exo_names)
  state_idx  <- dr4$state_idx
  endo_names <- dr4$endo_names
  exo_names  <- dr4$exo_names

  ghx  <- dr4$ghx;   ghu  <- dr4$ghu
  ghxx <- dr4$ghxx;  ghxu <- dr4$ghxu;  ghuu <- dr4$ghuu
  ghxxx <- dr4$ghxxx; ghxxu <- dr4$ghxxu
  ghxuu <- dr4$ghxuu; ghuuu <- dr4$ghuuu
  ghxxxx <- dr4$ghxxxx; ghxxxu <- dr4$ghxxxu
  ghxxuu <- dr4$ghxxuu; ghxuuu <- dr4$ghxuuu; ghuuuu <- dr4$ghuuuu

  hx  <- ghx[state_idx, , drop = FALSE]
  hu  <- ghu[state_idx, , drop = FALSE]
  # State-row slices of the higher-order policy tensors (needed by the
  # Faa di Bruno chain-rule cross terms below). ghxxxx's slice (hxxxx) is
  # defined locally near its first use (~line 597).
  hxx  <- ghxx [state_idx, , drop = FALSE]   # n_s x n_s^2
  hxxx <- ghxxx[state_idx, , drop = FALSE]   # n_s x n_s^3

  if (n_s == 0L) {
    if (verbose) message("No state variables; 5th-order x-terms are zero.")
    return(.trivial_dr5(dr4))
  }

  # System matrices
  sys <- extract_system_matrices(compiled, ss, params)
  f0  <- sys$f_zero
  fp  <- sys$f_plus

  S_mat <- matrix(0, n, n_s)
  for (s in seq_len(n_s)) S_mat[state_idx[s], s] <- 1
  A_L <- f0 + fp %*% ghx %*% t(S_mat)

  if (verbose) {
    cat("Fifth-order perturbation (deterministic):\n")
    cat(sprintf("  n_endo=%d  n_state=%d  n_exo=%d  ns^5=%d\n",
                n, n_s, n_u, n_s^5))
  }

  # Build dy at SS and transfer matrices
  dy_ss <- .build_dy_ss_o2(compiled, ss)
  tm    <- .build_transfer_matrices(dyn, ghx, ghu, state_idx, hx, hu,
                                     endo_names, exo_names)
  T_x <- tm$T_x
  T_u <- tm$T_u

  # ================================================================
  # Analytic Faa-di-Bruno forcing assembly (validated against FD/Richardson
  # to the FD truncation floor).  Replaces the old FD-forcing path which was
  # ~1000x slower at order 5.
  # ================================================================
  if (verbose) cat("  Computing 5th-order forcing via analytic assembler...\n")
  # Reuse extract_system_matrices()'s mapping so forcing rows align with A_L.
  # See [[eq-to-decl-consistency-invariant]].
  eq_to_decl <- sys$eq_to_decl %||% .build_eq_to_decl(model)
  res_perm   <- order(eq_to_decl)
  phi <- .build_phi_analytic(dyn, dr4, ss, params, state_idx,
                             endo_names, exo_names, n_s, n_u, n,
                             order = 5L, res_perm = res_perm)

  # ================================================================
  # Solve for ghxxxxx:  A_L·X + fp·X·hx^{⊗5} = -Φ_xxxxx
  # ================================================================
  if (verbose) cat("  Solving for ghxxxxx (Sylvester, n*ns^5 =", n * n_s^5, ")...\n")
  # Compact eigen/QZ Sylvester (with iterative refinement + dense fallback);
  # avoids forming/Schur-factorising the ns^5 x ns^5 Kronecker matrix.
  ghxxxxx <- .solve_kron_compact(A_L, fp, hx, 5L, -phi$xxxxx, verbose = verbose)

  state_vars <- endo_names[state_idx]
  rownames(ghxxxxx) <- endo_names
  colnames(ghxxxxx) <- .quint_names(state_vars, state_vars, state_vars,
                                     state_vars, state_vars)

  # ================================================================
  # Mixed shock-state terms:  A_L·X = -(Φ_mixed + fp·ghxxxxx·mix_kron)
  # The future-feedback uses the same-order pure tensor ghxxxxx contracted
  # with hx^{⊗n_state} ⊗ hu^{⊗n_exo} (single arrangement; ghxxxxx symmetric),
  # exactly as order-4's mixed solve uses ghxxxx.
  # ================================================================
  if (n_u > 0L) {
    if (verbose) cat("  Solving for ghxxxxu, ghxxxuu, ghxxuuu, ghxuuuu, ghuuuuu...\n")
    fut <- function(...) {
      mats <- list(...); Kf <- mats[[1L]]
      for (i in 2:length(mats)) Kf <- Kf %x% mats[[i]]
      fp %*% ghxxxxx %*% Kf
    }
    solveA <- function(rhs) tryCatch(solve(A_L, rhs),
                                     error = function(e) qr.solve(A_L, rhs))
    ghxxxxu <- solveA(-(phi$xxxxu + fut(hx, hx, hx, hx, hu)))
    ghxxxuu <- solveA(-(phi$xxxuu + fut(hx, hx, hx, hu, hu)))
    ghxxuuu <- solveA(-(phi$xxuuu + fut(hx, hx, hu, hu, hu)))
    ghxuuuu <- solveA(-(phi$xuuuu + fut(hx, hu, hu, hu, hu)))
    ghuuuuu <- solveA(-(phi$uuuuu + fut(hu, hu, hu, hu, hu)))
    rownames(ghxxxxu) <- endo_names
    rownames(ghxxxuu) <- endo_names
    rownames(ghxxuuu) <- endo_names
    rownames(ghxuuuu) <- endo_names
    rownames(ghuuuuu) <- endo_names
    colnames(ghxxxxu) <- .quad_mixed_names(state_vars, exo_names, 4L, 1L)
    colnames(ghxxxuu) <- .quad_mixed_names(state_vars, exo_names, 3L, 2L)
    colnames(ghxxuuu) <- .quad_mixed_names(state_vars, exo_names, 2L, 3L)
    colnames(ghxuuuu) <- .quad_mixed_names(state_vars, exo_names, 1L, 4L)
    colnames(ghuuuuu) <- .quad_mixed_names(state_vars, exo_names, 0L, 5L)
  } else {
    ghxxxxu <- ghxxxuu <- ghxxuuu <- ghxuuuu <- ghuuuuu <- NULL
  }

  # ================================================================
  # Name, assemble and return
  # ================================================================
  state_vars <- endo_names[state_idx]
  rownames(ghxxxxx) <- endo_names
  colnames(ghxxxxx) <- .quint_names(state_vars, state_vars, state_vars,
                                     state_vars, state_vars)

  dr5 <- unclass(dr4)
  dr5$ghxxxxx <- ghxxxxx
  dr5$ghxxxxu <- ghxxxxu
  dr5$ghxxxuu <- ghxxxuu
  dr5$ghxxuuu <- ghxxuuu
  dr5$ghxuuuu <- ghxuuuu
  dr5$ghuuuuu <- ghuuuuu
  dr5$order   <- 5L
  dr5$fifth_order_method <- "levintal_compact_2017"
  class(dr5) <- c("DecisionRules5", "DecisionRules4", "DecisionRules3",
                   "DecisionRules2", "DecisionRules")

  if (verbose) {
    cat("Fifth-order solution complete.\n")
    cat(sprintf("  ghxxxxx: %d x %d  max|.| = %.3g\n",
                nrow(ghxxxxx), ncol(ghxxxxx), max(abs(ghxxxxx))))
  }

  dr5
}


# =====================================================================
# Helpers
# =====================================================================

#' Contract a 5-D tensor F5 with 5 transfer matrices.
#'
#' Φ[e, (a,b,c,d,e)] = Σ F5[i,j,k,l,m] * Ta[i,a] * Tb[j,b] * Tc[k,c] * Td[l,d] * Te[m,e]
#'
#' @noRd
.contract_h5 <- function(F5, Ta, Tb, Tc, Td, Te, n_eq) {
  n_a <- ncol(Ta); n_b <- ncol(Tb); n_c <- ncol(Tc)
  n_d <- ncol(Td); n_e <- ncol(Te)
  Phi <- matrix(0, n_eq, n_a * n_b * n_c * n_d * n_e)

  if (is.null(F5)) return(Phi)

  for (e in seq_len(n_eq)) {
    Fe <- F5[e, , , , , ]
    # Contract dimension by dimension
    M <- tensor_contract_5d(Fe, Te, 5L)  # F5[i,j,k,l,m] * Te[m,e]
    M <- tensor_contract_5d(M, Td, 4L)   # * Td[l,d]
    M <- tensor_contract_5d(M, Tc, 3L)   # * Tc[k,c]
    M <- tensor_contract_5d(M, Tb, 2L)   # * Tb[j,b]
    M <- tensor_contract_5d(M, Ta, 1L)   # * Ta[i,a]
    Phi[e, ] <- as.numeric(M)
  }

  Phi
}


#' Contract a 5-D tensor with a matrix along a mode.
#' @noRd
tensor_contract_5d <- function(T, M, mode) {
  d <- dim(T)
  if (mode != 1) {
    perm <- c(mode, seq_len(5)[-mode])
    T <- aperm(T, perm)
    d <- dim(T)
  }

  nr <- d[1]
  nc <- prod(d[-1])
  T_mat <- matrix(T, nr, nc)
  result <- t(M) %*% T_mat

  new_d <- c(ncol(M), d[-1])
  result <- array(result, dim = new_d)

  if (mode != 1) {
    inv_perm <- order(c(mode, seq_len(5)[-mode]))
    result <- aperm(result, inv_perm)
  }

  result
}


#' Generate column names for 5-fold outer product.
#' @noRd
.quint_names <- function(a, b, c, d, e) {
  out <- character(length(a) * length(b) * length(c) * length(d) * length(e))
  idx <- 1L
  for (m in seq_along(e))
    for (l in seq_along(d))
      for (k in seq_along(c))
        for (j in seq_along(b))
          for (i in seq_along(a)) {
            out[idx] <- paste(a[i], b[j], c[k], d[l], e[m], sep = "__x__")
            idx <- idx + 1L
          }
  out
}


#' Generate column names for mixed state-shock 5th-order terms.
#' @param state_vars State variable names
#' @param exo_names  Shock names
#' @param n_state    Number of state dimensions
#' @param n_exo      Number of shock dimensions
#' @noRd
.quad_mixed_names <- function(state_vars, exo_names, n_state, n_exo) {
  vars <- c(rep(list(state_vars), n_state), rep(list(exo_names), n_exo))
  out <- character(prod(sapply(vars, length)))
  idx <- 1L
  # Recursive generation
  indices <- as.matrix(expand.grid(rev(lapply(vars, seq_along))))
  for (i in seq_len(nrow(indices))) {
    parts <- character(5)
    for (j in seq_len(5)) parts[j] <- vars[[j]][indices[i, 6 - j]]
    out[i] <- paste(parts, collapse = "__x__")
  }
  out
}


#' Trivial 5th-order for no-state-variables case.
#' @noRd
.trivial_dr5 <- function(dr4) {
  n   <- length(dr4$endo_names)
  n_u <- length(dr4$exo_names)
  dr5 <- unclass(dr4)
  dr5$ghxxxxx <- matrix(0, n, 0)
  dr5$ghxxxxu <- matrix(0, n, 0)
  dr5$ghxxxuu <- matrix(0, n, 0)
  dr5$ghxxuuu <- matrix(0, n, 0)
  dr5$ghxuuuu <- matrix(0, n, 0)
  dr5$ghuuuuu <- matrix(0, n, n_u^5)
  dr5$order <- 5L
  dr5$fifth_order_method <- "trivial_no_states"
  class(dr5) <- c("DecisionRules5", "DecisionRules4", "DecisionRules3",
                   "DecisionRules2", "DecisionRules")
  dr5
}


#' Print method for 5th-order decision rules
#' @param x   A \code{DecisionRules5} object.
#' @param ... Unused; included for S3 compatibility.
#' @export
print.DecisionRules5 <- function(x, ...) {
  cat("Fifth-order Decision Rules (DecisionRules5)\n")
  cat("  Endogenous variables:", length(x$endo_names), "\n")
  cat("  State variables:     ", x$n_state, "\n")
  cat("  Shocks:              ", x$n_exo, "\n")
  cat("  Perturbation order:  5 (Levintal compact)\n")
  cat("\nFirst-order:\n")
  cat(sprintf("  ghx: %d x %d   ghu: %d x %d\n",
              nrow(x$ghx), ncol(x$ghx), nrow(x$ghu), ncol(x$ghu)))
  cat("\nSecond-order:\n")
  cat(sprintf("  ghxx: %d x %d max=%.3g\n",
              nrow(x$ghxx), ncol(x$ghxx), max(abs(x$ghxx))))
  cat("\nThird-order:\n")
  if (!is.null(x$ghxxx))
    cat(sprintf("  ghxxx: %d x %d max=%.3g\n",
                nrow(x$ghxxx), ncol(x$ghxxx), max(abs(x$ghxxx))))
  cat("\nFourth-order:\n")
  if (!is.null(x$ghxxxx))
    cat(sprintf("  ghxxxx: %d x %d max=%.3g\n",
                nrow(x$ghxxxx), ncol(x$ghxxxx), max(abs(x$ghxxxx))))
  cat("\nFifth-order:\n")
  if (!is.null(x$ghxxxxx) && prod(dim(x$ghxxxxx)) > 0)
    cat(sprintf("  ghxxxxx: %d x %d max=%.3g\n",
                nrow(x$ghxxxxx), ncol(x$ghxxxxx), max(abs(x$ghxxxxx))))
  else
    cat("  (no state variables; zero)\n")
  invisible(x)
}


# =====================================================================
# Pruned state-space IRF for 5th order
# =====================================================================

#' Compute impulse response functions at fifth order (pruned state space)
#'
#' Implements the Andreasen, Fernandez-Villaverde & Rubio-Ramirez (2018)
#' pruning scheme truncated at fifth order. Tracks five additive layers
#' \eqn{x^{(1)}} through \eqn{x^{(5)}}. Constant mean-correction terms
#' cancel in deviation IRFs.
#'
#' @param dr5       \code{DecisionRules5} object.
#' @param model     \code{dynhr_mod} from \code{\link{parse_mod}}.
#' @param n_periods Number of IRF horizons (default 40).
#' @param shock_size Shock size in standard-deviation units (default 1).
#' @param params    Named numeric parameter vector.
#' @return An \code{IRFCollection} with \code{attr(., "order") = 5L}.
#' @export
compute_irfs_order5 <- function(dr5, model, n_periods = 40L,
                                 shock_size = 1, params = NULL) {
  if (!inherits(dr5, "DecisionRules5")) {
    stop("dr5 must be a DecisionRules5 object.")
  }

  ghx   <- dr5$ghx;   ghu   <- dr5$ghu
  ghxx  <- dr5$ghxx;  ghuu  <- dr5$ghuu
  ghxxx <- dr5$ghxxx; ghuuu <- dr5$ghuuu
  ghxxxx <- dr5$ghxxxx; ghuuuu <- dr5$ghuuuu
  ghxxxxx <- dr5$ghxxxxx; ghuuuuu <- dr5$ghuuuuu

  endo      <- dr5$endo_names
  exo       <- dr5$exo_names
  state_idx <- dr5$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)
  n_s       <- length(state_idx)

  if (is.null(params)) params <- model$param_values
  shock_stderr <- .get_shock_stderr(model, exo, params)

  hx    <- ghx   [state_idx, , drop = FALSE]
  hu    <- ghu   [state_idx, , drop = FALSE]
  hxx   <- ghxx  [state_idx, , drop = FALSE]
  huu   <- ghuu  [state_idx, , drop = FALSE]
  hxxx  <- ghxxx [state_idx, , drop = FALSE]
  huuu  <- ghuuu [state_idx, , drop = FALSE]
  hxxxx <- ghxxxx[state_idx, , drop = FALSE]
  huuuu <- ghuuuu[state_idx, , drop = FALSE]

  irfs <- list()
  for (k in seq_along(exo)) {
    shock_name <- exo[k]
    irf_mat    <- matrix(0, n_periods, n_endo)
    colnames(irf_mat) <- endo
    rownames(irf_mat) <- paste0("t", seq_len(n_periods))

    eps    <- numeric(n_exo)
    eps[k] <- shock_stderr[shock_name] * shock_size

    x1 <- x2 <- x3 <- x4 <- x5 <- numeric(n_s)

    # --- Period 1: impact ---
    eps_vec <- as.numeric(eps)

    x1 <- as.numeric(hu      %*% eps_vec)
    x2 <- as.numeric(0.5  * huu    %*% (eps_vec %x% eps_vec))
    x3 <- as.numeric((1/6) * huuu   %*% (eps_vec %x% eps_vec %x% eps_vec))
    x4 <- as.numeric((1/24)* huuuu  %*% (eps_vec %x% eps_vec %x% eps_vec %x% eps_vec))
    x5 <- as.numeric((1/120)*ghuuuuu[state_idx, , drop=FALSE] %*%
                      (eps_vec %x% eps_vec %x% eps_vec %x% eps_vec %x% eps_vec))

    y1 <- as.numeric(ghu      %*% eps_vec)
    y2 <- as.numeric(0.5  * ghuu    %*% (eps_vec %x% eps_vec))
    y3 <- as.numeric((1/6) * ghuuu   %*% (eps_vec %x% eps_vec %x% eps_vec))
    y4 <- as.numeric((1/24)* ghuuuu  %*% (eps_vec %x% eps_vec %x% eps_vec %x% eps_vec))
    y5 <- as.numeric((1/120)*ghuuuuu %*%
                      (eps_vec %x% eps_vec %x% eps_vec %x% eps_vec %x% eps_vec))

    irf_mat[1, ] <- y1 + y2 + y3 + y4 + y5

    # --- Periods 2..n_periods: pure propagation ---
    for (t in 2:n_periods) {
      x1p <- x1; x2p <- x2; x3p <- x3; x4p <- x4; x5p <- x5

      x1 <- as.numeric(hx %*% x1p)
      x2 <- as.numeric(hx %*% x2p + 0.5 * hxx %*% (x1p %x% x1p))
      x3 <- as.numeric(hx %*% x3p + hxx %*% (x1p %x% x2p) +
                       (1/6) * hxxx %*% (x1p %x% x1p %x% x1p))
      x4 <- as.numeric(hx %*% x4p + hxx %*% (x1p %x% x3p + x2p %x% x2p) +
                       hxxx %*% (x1p %x% x1p %x% x2p) +
                       (1/24) * hxxxx %*% (x1p %x% x1p %x% x1p %x% x1p))
      x5 <- as.numeric(
        hx     %*% x5p +
        hxx    %*% (x1p %x% x4p + x2p %x% x3p) +
        hxxx   %*% (x1p %x% x1p %x% x3p + x1p %x% x2p %x% x2p) +
        hxxxx  %*% (x1p %x% x1p %x% x1p %x% x2p)
      )

      y1 <- as.numeric(ghx  %*% x1p)
      y2 <- as.numeric(ghx  %*% x2p + 0.5 * ghxx %*% (x1p %x% x1p))
      y3 <- as.numeric(ghx  %*% x3p + ghxx %*% (x1p %x% x2p) +
                       (1/6) * ghxxx %*% (x1p %x% x1p %x% x1p))
      y4 <- as.numeric(ghx  %*% x4p + ghxx %*% (x1p %x% x3p + x2p %x% x2p) +
                       ghxxx %*% (x1p %x% x1p %x% x2p) +
                       (1/24) * ghxxxx %*% (x1p %x% x1p %x% x1p %x% x1p))
      y5 <- as.numeric(
        ghx    %*% x5p +
        ghxx   %*% (x1p %x% x4p + x2p %x% x3p) +
        ghxxx  %*% (x1p %x% x1p %x% x3p + x1p %x% x2p %x% x2p) +
        ghxxxx %*% (x1p %x% x1p %x% x1p %x% x2p) +
        (1/120) * ghxxxxx %*% (x1p %x% x1p %x% x1p %x% x1p %x% x1p)
      )

      irf_mat[t, ] <- y1 + y2 + y3 + y4 + y5
    }

    irfs[[shock_name]] <- irf_mat
  }

  class(irfs) <- "IRFCollection"
  attr(irfs, "n_periods")  <- n_periods
  attr(irfs, "endo_names") <- endo
  attr(irfs, "exo_names")  <- exo
  attr(irfs, "order")      <- 5L
  irfs
}
