## R/solve-perturbation-order3-fd.R
## --------------------------------------------------------------------------
## Finite-difference based third-order perturbation solver.
##
## This is a robust, numerically-derived alternative to the symbolic
## Faà di Bruno solver in solve-perturbation-order3.R. The symbolic
## solver has subtle index-convention bugs in its Phi_xxx forcing term
## (see ORDER3_INVESTIGATION.md) that produce incorrect ghxxx when
## ns ≠ nu. This file extracts the Phi_*** forcings directly by
## finite-differencing the policy-substituted residual function, which
## is provably correct (just slower; O(ns^3 + nu^3) residual evaluations).
##
## Best for small/medium models (n * (ns+nu)^3 < ~10000). For larger
## models, the symbolic solver is needed for performance — but it must
## be fixed.
## --------------------------------------------------------------------------


#' Solve order-3 perturbation using finite-difference Phi forcings.
#'
#' Computes ghxxx, ghxxu, ghxuu, ghuuu by:
#'   1. Setting all order-3 policy terms to zero
#'   2. FD-computing the symmetric 3-tensor d^3 R / d(x,u)^3 at SS
#'      where R is the policy-substituted system residual
#'   3. Solving the appropriate linear/Sylvester system per block
#'
#' Sigma-correction terms (ghxss, ghuss, ghs3) are NOT computed here.
#'
#' @noRd
.solve_perturbation_order3_fd <- function(model, compiled, ss, params,
                                          dr2, verbose = FALSE,
                                          fd_eps = 1e-3) {
  dyn        <- compiled$dynamic
  endo       <- dr2$endo_names
  exo        <- dr2$exo_names
  state_idx  <- dr2$state_idx
  n  <- length(endo)
  ns <- length(state_idx)
  nu <- length(exo)
  ss_vec <- ss[endo]

  if (ns == 0L) {
    if (verbose) message("No state variables; third-order x-terms are zero.")
    return(.trivial_dr3(dr2))
  }

  # Row-permutation: residuals_fn returns in compiled-eq order;
  # we need declaration order to match A_L (which is in decl order).
  # CRITICAL: reuse the SAME mapping extract_system_matrices() used to
  # reorder A_L (sys$eq_to_decl), NOT a freshly recomputed one -- the
  # latter can differ and misalign the residual rows against A_L,
  # injecting spurious order-3 forcing into linear equations.
  # See [[eq-to-decl-consistency-invariant]].
  sys <- extract_system_matrices(compiled, ss, params)
  eq_to_decl <- sys$eq_to_decl %||% .build_eq_to_decl(model)
  res_perm <- order(eq_to_decl)

  # Build dy from y_lag, y_now, y_lead, u_now (state perturbations only;
  # works only for single-period lead/lag models without AUX vars)
  build_dy <- function(y_lag, y_now, y_lead, u_now) {
    dcm <- dyn$dyn_col_map
    dy  <- numeric(dyn$total_cols)
    keys <- character(dyn$total_cols)
    for (kc in seq_len(nrow(dcm))) {
      c  <- dcm$col[kc]; nm <- dcm$name[kc]; ll <- dcm$lead_lag[kc]
      sfx <- if (ll == 0L) "__0" else if (ll > 0L) paste0("__p", ll) else paste0("__m", abs(ll))
      keys[c] <- paste0(nm, sfx)
      if (nm %in% exo) {
        dy[c] <- u_now[which(exo == nm)]
      } else {
        idx <- which(endo == nm); if (length(idx) != 1L) next
        dy[c] <- if (ll == -1L) y_lag[idx]
                 else if (ll == 0L) y_now[idx]
                 else if (ll == 1L) y_lead[idx]
                 else NA_real_
      }
    }
    names(dy) <- keys; dy
  }

  has_aux <- any(grepl("^AUX_", endo))
  has_multi_lead <- any(dyn$dyn_col_map$lead_lag > 1L) ||
                    any(dyn$dyn_col_map$lead_lag < -1L)
  if (has_aux || has_multi_lead) {
    stop("FD order-3 solver requires single-period leads/lags and no AUX ",
         "variables. Use the symbolic solver for this model.")
  }

  # Permute dr2$ghxu from internal convention (state FAST, exo SLOW) to
  # standard Kron (state SLOW, exo FAST) so `ghxu_std %*% (x %x% u)` is
  # correct. For nu==1 this is identity; for nu>1 it matters.
  ghxu_std <- dr2$ghxu
  if (ns > 0L && nu > 0L) {
    .pxu_fd <- integer(ns * nu)
    for (s in seq_len(ns)) for (u in seq_len(nu)) {
      .pxu_fd[(s-1L)*nu + u] <- (u-1L)*ns + s
    }
    ghxu_std <- dr2$ghxu[, .pxu_fd, drop = FALSE]
  }

  # Order-2 policy evaluation (no order-3 terms)
  policy_o2 <- function(x, u) {
    y <- as.numeric(dr2$ghx %*% x + dr2$ghu %*% u)
    y <- y + 0.5*as.numeric(dr2$ghxx %*% (x %x% x)) +
            as.numeric(ghxu_std %*% (x %x% u)) +
            0.5*as.numeric(dr2$ghuu %*% (u %x% u))
    y
  }

  # Residual function R(x_lag, u_now) with order-3 policy = 0 (sigma=0).
  # u_{t+1} = 0 deterministically.
  R_fn <- function(x_lag, u_now) {
    y_now    <- ss_vec + policy_o2(x_lag, u_now)
    x_now    <- y_now[state_idx] - ss_vec[state_idx]
    y_lead   <- ss_vec + policy_o2(x_now, numeric(nu))
    y_lag_lv <- ss_vec
    y_lag_lv[state_idx] <- ss_vec[state_idx] + x_lag
    dy <- build_dy(y_lag_lv, y_now, y_lead, u_now)
    res <- as.numeric(dyn$residuals_fn(dy, params, ss))
    res[res_perm]  # decl-order
  }

  # System matrices in declaration order (sys computed above)
  f0  <- sys$f_zero
  fp  <- sys$f_plus
  S   <- matrix(0, n, ns); for (s in seq_len(ns)) S[state_idx[s], s] <- 1
  A_L <- f0 + fp %*% dr2$ghx %*% t(S)
  hx  <- dr2$ghx[state_idx, , drop = FALSE]
  hu  <- dr2$ghu[state_idx, , drop = FALSE]

  if (verbose) {
    cat("FD-based third-order perturbation:\n")
    cat(sprintf("  n=%d ns=%d nu=%d\n", n, ns, nu))
  }

  # ----------------------------------------------------------------------
  # Helper: FD extract symmetric 3-tensor d^3 R / d(z_i d z_j d z_k) where
  # z is some combination of (x, u). For each canonical (i,j,k) compute the
  # mixed third difference. Universal stencil (works for repeated indices):
  #
  #   d^3R/dz_i dz_j dz_k = (1/(8h^3)) Σ_{s ∈ ±1}^3 s_i s_j s_k
  #                          R(h(s_i e_i + s_j e_j + s_k e_k))
  #
  # The result is symmetric in (i,j,k).
  # ----------------------------------------------------------------------
  fd_third <- function(eval_fn, idx_i, idx_j, idx_k, h) {
    acc <- numeric(n)
    for (si in c(-1, 1)) for (sj in c(-1, 1)) for (sk in c(-1, 1)) {
      acc <- acc + si*sj*sk * eval_fn(idx_i, idx_j, idx_k, si*h, sj*h, sk*h)
    }
    acc / (8 * h^3)
  }

  # ============================================================
  # 1. Compute Phi_xxx by FD (with all order-3 policy = 0)
  # ============================================================
  if (verbose) cat("  FD-computing Phi_xxx...\n")
  eval_x <- function(i, j, k, hi, hj, hk) {
    x <- numeric(ns); x[i] <- x[i]+hi; x[j] <- x[j]+hj; x[k] <- x[k]+hk
    R_fn(x, numeric(nu))
  }
  # Two-step Richardson on h = fd_eps and h/2 to remove O(h^2) bias
  h1 <- fd_eps; h2 <- fd_eps/2
  Phi_xxx <- matrix(0, n, ns^3)
  for (i in 1:ns) for (j in i:ns) for (k in j:ns) {
    d1 <- fd_third(eval_x, i, j, k, h1)
    d2 <- fd_third(eval_x, i, j, k, h2)
    d_rich <- (4*d2 - d1) / 3
    perms <- unique(list(c(i,j,k), c(i,k,j), c(j,i,k), c(j,k,i), c(k,i,j), c(k,j,i)))
    for (p in perms) {
      col <- p[1] + (p[2]-1)*ns + (p[3]-1)*ns^2
      Phi_xxx[, col] <- d_rich
    }
  }

  # ============================================================
  # 2. Solve Sylvester for ghxxx
  #    A_L * X + fp * X * (hx ⊗ hx ⊗ hx) = -Phi_xxx
  # ============================================================
  if (verbose) cat("  Solving Sylvester for ghxxx...\n")
  K_x <- t(hx) %x% t(hx) %x% t(hx)
  Sylv_lhs <- kronecker(diag(ns^3), A_L) + kronecker(K_x, fp)
  ghxxx <- matrix(solve(Sylv_lhs, -as.vector(Phi_xxx)), nrow = n, ncol = ns^3)

  # ============================================================
  # 3. Phi_xxu: FD with ghxxx now known (substituted back into policy)
  #    Form: A_L * ghxxu = -Phi_xxu (no Sylvester recursion since u_{t+1}=0)
  # ============================================================
  if (verbose) cat("  FD-computing Phi_xxu...\n")
  # Re-define R with ghxxx substituted (order-3 x-only part of policy)
  policy_with_xxx <- function(x, u) {
    y <- policy_o2(x, u)
    y <- y + (1/6)*as.numeric(ghxxx %*% (x %x% x %x% x))
    y
  }
  R_with_xxx <- function(x_lag, u_now) {
    y_now    <- ss_vec + policy_with_xxx(x_lag, u_now)
    x_now    <- y_now[state_idx] - ss_vec[state_idx]
    y_lead   <- ss_vec + policy_with_xxx(x_now, numeric(nu))
    y_lag_lv <- ss_vec
    y_lag_lv[state_idx] <- ss_vec[state_idx] + x_lag
    dy <- build_dy(y_lag_lv, y_now, y_lead, u_now)
    res <- as.numeric(dyn$residuals_fn(dy, params, ss))
    res[res_perm]
  }
  # d^3 R / (dx_i dx_j du_k) — symmetric in (i,j) only
  eval_xxu <- function(i, j, k, hi, hj, hk) {
    x <- numeric(ns); x[i] <- x[i]+hi; x[j] <- x[j]+hj
    u <- numeric(nu); u[k] <- u[k]+hk
    R_with_xxx(x, u)
  }
  Phi_xxu <- matrix(0, n, ns^2 * nu)
  for (i in 1:ns) for (j in i:ns) for (k in 1:nu) {
    d1 <- fd_third(eval_xxu, i, j, k, h1)
    d2 <- fd_third(eval_xxu, i, j, k, h2)
    d_rich <- (4*d2 - d1) / 3
    # Standard Kron convention: col = (s1-1)*ns*nu + (s2-1)*nu + k
    for (ord_pair in unique(list(c(i,j), c(j,i)))) {
      s1 <- ord_pair[1]; s2 <- ord_pair[2]
      col <- (s1-1)*ns*nu + (s2-1)*nu + k
      Phi_xxu[, col] <- d_rich
    }
  }
  ghxxu <- solve(A_L, -Phi_xxu)

  # ============================================================
  # 4. Phi_xuu: similar
  # ============================================================
  if (verbose) cat("  FD-computing Phi_xuu...\n")
  eval_xuu <- function(i, j, k, hi, hj, hk) {
    x <- numeric(ns); x[i] <- x[i]+hi
    u <- numeric(nu); u[j] <- u[j]+hj; u[k] <- u[k]+hk
    R_with_xxx(x, u)
  }
  Phi_xuu <- matrix(0, n, ns * nu^2)
  for (i in 1:ns) for (j in 1:nu) for (k in j:nu) {
    d1 <- fd_third(eval_xuu, i, j, k, h1)
    d2 <- fd_third(eval_xuu, i, j, k, h2)
    d_rich <- (4*d2 - d1) / 3
    # Standard Kron convention: col = (i-1)*nu^2 + (u1-1)*nu + u2
    for (ord_pair in unique(list(c(j,k), c(k,j)))) {
      u1 <- ord_pair[1]; u2 <- ord_pair[2]
      col <- (i-1)*nu^2 + (u1-1)*nu + u2
      Phi_xuu[, col] <- d_rich
    }
  }
  ghxuu <- solve(A_L, -Phi_xuu)

  # ============================================================
  # 5. Phi_uuu: fully symmetric in (u1, u2, u3)
  # ============================================================
  if (verbose) cat("  FD-computing Phi_uuu...\n")
  eval_uuu <- function(i, j, k, hi, hj, hk) {
    u <- numeric(nu); u[i] <- u[i]+hi; u[j] <- u[j]+hj; u[k] <- u[k]+hk
    R_with_xxx(numeric(ns), u)
  }
  Phi_uuu <- matrix(0, n, nu^3)
  for (i in 1:nu) for (j in i:nu) for (k in j:nu) {
    d1 <- fd_third(eval_uuu, i, j, k, h1)
    d2 <- fd_third(eval_uuu, i, j, k, h2)
    d_rich <- (4*d2 - d1) / 3
    perms <- unique(list(c(i,j,k), c(i,k,j), c(j,i,k), c(j,k,i), c(k,i,j), c(k,j,i)))
    for (p in perms) {
      col <- p[1] + (p[2]-1)*nu + (p[3]-1)*nu^2
      Phi_uuu[, col] <- d_rich
    }
  }
  ghuuu <- solve(A_L, -Phi_uuu)

  # ============================================================
  # Build DecisionRules3 object
  # ============================================================
  dr3 <- dr2
  state_vars <- endo[state_idx]
  triple_names <- function(a, b, c) {
    out <- character(length(a)*length(b)*length(c)); idx <- 1L
    for (c3 in seq_along(c)) for (c2 in seq_along(b)) for (c1 in seq_along(a)) {
      out[idx] <- paste(a[c1], b[c2], c[c3], sep = "__x__"); idx <- idx + 1L
    }
    out
  }
  rownames(ghxxx) <- endo; rownames(ghxxu) <- endo
  rownames(ghxuu) <- endo; rownames(ghuuu) <- endo
  colnames(ghxxx) <- triple_names(state_vars, state_vars, state_vars)
  colnames(ghxxu) <- triple_names(state_vars, state_vars, exo)
  colnames(ghxuu) <- triple_names(state_vars, exo, exo)
  colnames(ghuuu) <- triple_names(exo, exo, exo)

  # Permute dr2$ghxu to standard Kron convention (state SLOW, exo FAST) on
  # output so `dr3$ghxu %*% (x %x% u)` is correct downstream. For nu==1
  # this is identity; for nu>1 it matters.
  dr3$ghxu  <- ghxu_std
  dr3$ghxxx <- ghxxx
  dr3$ghxxu <- ghxxu
  dr3$ghxuu <- ghxuu
  dr3$ghuuu <- ghuuu
  dr3$order <- 3L
  dr3$sigma_correction <- "not_implemented"
  dr3$solver <- "fd_order3"
  class(dr3) <- c("DecisionRules3", class(dr2))
  dr3
}
