## R/solve-perturbation-order3-sigma.R
## --------------------------------------------------------------------------
## Phase 7c (sigma-correction terms for third-order perturbation).
##
## This file currently implements ONE piece of Phase 7c:
##   ghs3 -- the third-cumulant correction term.
## It is the "smallest beachhead": for Gaussian shocks (SIGMA3 = 0) it
## collapses to zero, and is needed only when the user supplies non-zero
## third-order shock moments E[u_i u_j u_k] (e.g., skewed innovations).
##
## ghxss and ghuss (the time-varying risk-premium / state-sigma cross
## terms) are NOT implemented here; they require the full Levintal-style
## augmented-state recursion with the W_xup, W_uup, W_upup compound-
## derivative matrices and a Sylvester solve for ghxss.  Slated for
## Phase 7c.2 / 7c.3.
##
## Reference (canonical translation):
##   Mutschler (2022), perturbation_solver_nonsymmetric_order3.m,
##   lines 158-164 (the ghs3 block).
##
## Formula:
##   (A + B) ghs3 = -[ f_+ ghuuu(yp) + d2f(yp,yp) (sum_3perm kron(ghu,ghuu)(yp))
##                    + d3f(yp,yp,yp) kron3(ghu(yp)) ] SIGMA3
## where:
##   yp           = jumper variables (those that appear at t+1)
##   A            = first-order LHS matrix = f0 + f_+ ghx S' (= dynhr's A_L)
##   B            = f_+ acting on jumper columns (= dynhr's f_plus)
##   SIGMA3       = vec(E[u (x) u (x) u]) in col-major (k3 fastest) layout
## --------------------------------------------------------------------------


# =====================================================================
# Helper: column-permutation index vectors for the symmetric sum of
# three "single u vs paired uu" arrangements
# =====================================================================

#' Build a column-permutation index for the (n_u^3)-wide symmetric sum.
#'
#' Treat each column of an (anything x n_u^3) matrix as indexed by
#' (k1, k2, k3) in col-major order (k3 fastest, k1 slowest).
#'
#' `perm_dims = c(a, b, c)` means: at output column position (k1, k2, k3),
#' pull the source-matrix column whose (k1', k2', k3') = (k_a, k_b, k_c).
#'
#' Used to translate Mutschler's identities:
#'   P_u1_u2u3 ~ c(1, 2, 3)  -- identity (a=k1, b=k2, c=k3)
#'   P_u2_u1u3 ~ c(2, 1, 3)  -- swap singleton and first paired (a=k2, b=k1)
#'   P_u3_u1u2 ~ c(3, 1, 2)  -- cycle (a=k3, b=k1, c=k2)
#'
#' @noRd
.perm_idx_u3 <- function(perm_dims, n_u) {
  out <- integer(n_u^3)
  i <- 0L
  for (k1 in seq_len(n_u)) {
    for (k2 in seq_len(n_u)) {
      for (k3 in seq_len(n_u)) {
        i <- i + 1L
        ks <- c(k1, k2, k3)
        tk1 <- ks[perm_dims[1]]
        tk2 <- ks[perm_dims[2]]
        tk3 <- ks[perm_dims[3]]
        out[i] <- (tk1 - 1L) * n_u * n_u + (tk2 - 1L) * n_u + tk3
      }
    }
  }
  out
}


#' Reorder dense/sparse derivative rows from compiled equation order to
#' declaration-variable order.
#'
#' `extract_system_matrices()` applies this same row convention to the
#' Jacobian before building `A_L` and `fp`.  Higher-order derivative tensors
#' evaluated directly from the compiled dynamic function still arrive in
#' compiled equation order, so sigma/cumulant RHS terms must be permuted before
#' they are combined with first-order system matrices.
#'
#' @noRd
.reorder_derivatives_to_decl <- function(H2, H3, compiled, n,
                                         verbose = FALSE,
                                         label = "derivatives",
                                         eq_to_decl = NULL) {
  if (is.null(compiled$model$equations)) {
    return(list(H2 = H2, H3 = H3))
  }

  # Reuse the caller's extract_system_matrices() mapping when provided so the
  # forcing rows align with A_L; recomputing it here can silently diverge.
  # See [[eq-to-decl-consistency-invariant]].
  if (is.null(eq_to_decl)) eq_to_decl <- .build_eq_to_decl(compiled$model)
  if (!all(eq_to_decl > 0L) || identical(eq_to_decl, seq_len(n))) {
    return(list(H2 = H2, H3 = H3))
  }

  perm <- order(eq_to_decl)

  if (!is.null(H2)) {
    H2_perm <- array(0, dim = dim(H2))
    for (k in seq_len(n)) H2_perm[k, , ] <- H2[perm[k], , ]
    H2 <- H2_perm
  }

  if (!is.null(H3) && length(H3$triplets) > 0L) {
    inv_perm <- order(perm)
    for (k in seq_along(H3$triplets)) {
      H3$triplets[[k]]$eq <- inv_perm[H3$triplets[[k]]$eq]
    }
  }

  if (verbose) {
    cat(sprintf("  %s rows reordered: compiled -> declaration order.\n",
                label))
  }

  list(H2 = H2, H3 = H3)
}


# =====================================================================
# Helper: identify the jumper compound-columns
# =====================================================================

#' For each jumper (lead) variable in endo order, find the compound
#' column index in dyn_col_map whose (name = var, lead_lag = +1).
#'
#' Returns a list with:
#'   jumper_idx        : integer, indices into endo_names of jumpers (has_lead)
#'   jumper_compound_c : integer, compound-col index of each jumper at ll=+1
#'                       (same length & order as jumper_idx)
#'
#' If a "jumper" by has_lead has no lead column in dyn_col_map (shouldn't
#' happen, but defensive), it is dropped with a warning.
#'
#' @noRd
.identify_jumpers <- function(dyn, endo_names, has_lead) {
  jumper_idx <- which(has_lead)
  dcm        <- dyn$dyn_col_map
  jc         <- integer(length(jumper_idx))
  keep       <- logical(length(jumper_idx))
  for (i in seq_along(jumper_idx)) {
    nm   <- endo_names[jumper_idx[i]]
    hits <- which(dcm$name == nm & dcm$lead_lag == 1L)
    if (length(hits) == 1L) {
      jc[i]   <- dcm$col[hits]
      keep[i] <- TRUE
    } else {
      keep[i] <- FALSE
    }
  }
  if (!all(keep)) {
    warning("Some jumper variables have no lead compound-column; dropping.")
  }
  list(jumper_idx        = jumper_idx[keep],
       jumper_compound_c = jc[keep])
}


# =====================================================================
# Main entry: recompute ghs3 given a third-moment vector
# =====================================================================

#' Recompute the third-cumulant correction term `ghs3`
#'
#' For Gaussian shocks, SIGMA3 = 0 and `ghs3 = 0`. Supply a non-zero
#' `sigma3` (the n_u^3 vector of `E[kron(u, u, u)]`) to obtain the
#' non-Gaussian third-cumulant correction.
#'
#' This is the order-3 analog of Schmitt-Grohe & Uribe's `ghss` for
#' covariance: `ghs3` is the unconditional bias the third moment of the
#' innovations introduces into the policy at order (sigma^3).
#'
#' Layout convention for `sigma3`: column-major, k3 fastest, k1 slowest
#' (matches `as.vector(array(E[u_i u_j u_k], dim = c(n_u, n_u, n_u)))`
#' when the first index varies slowest -- i.e. the same layout as
#' `kron(u, kron(u, u))` in MATLAB / Mutschler 2022).
#'
#' @param dr3      DecisionRules3 from `solve_perturbation_order3()`
#' @param compiled dynhr_compiled from `compile_model()`
#' @param ss       Named numeric steady state vector
#' @param params   Named numeric parameter vector
#' @param sigma3   Numeric vector of length n_u^3 with the third-order
#'                 product moments of the shocks. Default NULL is treated
#'                 as Gaussian (all zero), yielding `ghs3 = 0`.
#' @return The `dr3` object with `$ghs3` (length-n_endo vector) and
#'   `$sigma3` (the moment vector used) added. `$sigma_correction` is
#'   updated to `"ghs3_only"` (when `sigma3` is non-zero) or stays
#'   `"not_implemented"` (when Gaussian -- to signal that the rest of
#'   the sigma terms are still unimplemented).
#' @export
solve_third_cumulant <- function(dr3, compiled, ss, params, sigma3 = NULL) {
  if (!inherits(dr3, "DecisionRules3")) {
    stop("dr3 must be a DecisionRules3 object from solve_perturbation_order3().")
  }

  n         <- length(dr3$endo_names)
  n_u       <- length(dr3$exo_names)

  # Validate sigma3 length BEFORE the all-zero shortcut (so the user
  # gets a clear error for a mis-sized vector, even one of zeros).
  if (!is.null(sigma3) && length(sigma3) != n_u^3) {
    stop(sprintf(
      "sigma3 must have length n_u^3 = %d (got %d).",
      n_u^3, length(sigma3)))
  }

  # -----------------------------------------------------------------
  # Default Gaussian path: ghs3 == 0
  # -----------------------------------------------------------------
  if (is.null(sigma3) || all(sigma3 == 0)) {
    ghs3        <- numeric(n)
    names(ghs3) <- dr3$endo_names
    dr3$ghs3    <- ghs3
    dr3$sigma3  <- numeric(n_u^3)
    return(dr3)
  }

  # No jumpers => no future-shock contribution; ghs3 is zero
  endo_names <- dr3$endo_names
  exo_names  <- dr3$exo_names
  dyn        <- compiled$dynamic
  sys        <- extract_system_matrices(compiled, ss, params)
  has_lead   <- sys$is_fwd | sys$is_mixed
  if (!any(has_lead)) {
    ghs3        <- numeric(n)
    names(ghs3) <- endo_names
    dr3$ghs3    <- ghs3
    dr3$sigma3  <- as.numeric(sigma3)
    dr3$sigma_correction <- "ghs3_only"
    return(dr3)
  }

  # -----------------------------------------------------------------
  # LHS matrix: A + B  (Mutschler 2022 notation)
  #   A = f0 + f+ * ghx * S'  =  A_L  (already used by Phase 7b)
  #   B = f+ restricted to jumper columns; in dynhr f_plus is already
  #       zero on non-jumper columns, so B = f_plus directly.
  # -----------------------------------------------------------------
  state_idx <- dr3$state_idx
  ghx       <- dr3$ghx
  S <- matrix(0, n, length(state_idx))
  for (s in seq_along(state_idx)) S[state_idx[s], s] <- 1
  A_L <- sys$f_zero + sys$f_plus %*% ghx %*% t(S)
  AB  <- A_L + sys$f_plus

  # -----------------------------------------------------------------
  # Identify jumpers and their compound-column indices (for the
  # lead-block of the dense Hessian / Hessian3).
  # -----------------------------------------------------------------
  jm   <- .identify_jumpers(dyn, endo_names, has_lead)
  yp   <- jm$jumper_idx           # endo indices of jumpers
  ypc  <- jm$jumper_compound_c    # compound-col indices, lead block
  n_yp <- length(yp)

  ghu_yp   <- dr3$ghu  [yp, , drop = FALSE]   # n_yp x n_u
  ghuu_yp  <- dr3$ghuu [yp, , drop = FALSE]   # n_yp x n_u^2
  ghuuu_yp <- dr3$ghuuu[yp, , drop = FALSE]   # n_yp x n_u^3

  # -----------------------------------------------------------------
  # Build d2f(yp, yp) and d3f(yp, yp, yp) as compact (n_eq x n_yp^k)
  # matrices by slicing the dense (n_eq x total_cols^k) representations
  # at the lead-block compound-column indices.
  # -----------------------------------------------------------------
  dy_ss <- .build_dy_ss_o2(compiled, ss)
  H2    <- .compute_model_hessian_symbolic(compiled, dy_ss, params, ss)
  if (is.null(H2)) {
    stop("Symbolic Hessian2 unavailable; required for ghs3.")
  }
  H3    <- .compute_model_hessian3_symbolic(compiled, dy_ss, params, ss)
  H_ord <- .reorder_derivatives_to_decl(H2, H3, compiled, n,
                                        label = "ghs3 derivative",
                                        eq_to_decl = sys$eq_to_decl)
  H2    <- H_ord$H2
  H3    <- H_ord$H3

  n_eq  <- dyn$n_eq

  # d2f_yp_yp[e, (j1, j2)] = H2[e, ypc[j1], ypc[j2]],  flattened col-major (j2 fast)
  d2f_yp_yp <- matrix(0, n_eq, n_yp * n_yp)
  for (e in seq_len(n_eq)) {
    d2f_yp_yp[e, ] <- as.numeric(H2[e, ypc, ypc])
  }

  # d3f_yp_yp_yp[e, (j1, j2, j3)] = H3-dense[e, ypc[j1], ypc[j2], ypc[j3]]
  # We assemble the lead-block 3-tensor from H3's sparse triplet form,
  # expanding the Schwarz orbit (orbits of 1, 3, or 6 entries).
  d3f_yp_yp_yp <- matrix(0, n_eq, n_yp * n_yp * n_yp)
  if (!is.null(H3) && length(H3$triplets) > 0L) {
    # Reverse map: compound-col -> jumper index (0 if not a jumper-lead col)
    rev_yp <- integer(dyn$total_cols)
    rev_yp[ypc] <- seq_along(ypc)
    stride_j2 <- n_yp
    stride_j3 <- n_yp * n_yp
    for (k in seq_along(H3$triplets)) {
      t   <- H3$triplets[[k]]
      val <- H3$values[k]
      if (val == 0) next
      # Each triplet (c1<=c2<=c3) is in canonical order; require all three
      # cols to be lead-block jumper cols to contribute to d3f_yp_yp_yp.
      j1 <- rev_yp[t$col1]
      j2 <- rev_yp[t$col2]
      j3 <- rev_yp[t$col3]
      if (j1 == 0L || j2 == 0L || j3 == 0L) next
      orbit <- .orbit_3(j1, j2, j3)
      e <- t$eq
      for (perm in orbit) {
        flat <- (perm[3] - 1L) * stride_j3 + (perm[2] - 1L) * stride_j2 + perm[1]
        d3f_yp_yp_yp[e, flat] <- d3f_yp_yp_yp[e, flat] + val
      }
    }
  }

  # -----------------------------------------------------------------
  # Build the three n_eq x n_u^3 contributions, sum, multiply by SIGMA3
  # -----------------------------------------------------------------
  # f+ * ghuuu(yp)  ->  n_eq x n_u^3
  term_uuu  <- sys$f_plus[, yp, drop = FALSE] %*% ghuuu_yp

  # kron(ghu_yp, ghuu_yp) -> (n_yp^2) x (n_u^3)
  K_pair    <- ghu_yp %x% ghuu_yp

  # The three symmetric-sum permutations
  p1 <- seq_len(n_u^3)
  p2 <- .perm_idx_u3(c(2L, 1L, 3L), n_u)
  p3 <- .perm_idx_u3(c(3L, 1L, 2L), n_u)
  K_pair_sym <- K_pair[, p1, drop = FALSE] +
                K_pair[, p2, drop = FALSE] +
                K_pair[, p3, drop = FALSE]

  term_pair <- d2f_yp_yp %*% K_pair_sym   # n_eq x n_u^3

  # kron3(ghu_yp)  -> (n_yp^3) x (n_u^3)
  K_triple  <- ghu_yp %x% ghu_yp %x% ghu_yp
  term_trip <- d3f_yp_yp_yp %*% K_triple  # n_eq x n_u^3

  RHS_mat <- term_uuu + term_pair + term_trip       # n_eq x n_u^3
  RHS_vec <- as.numeric(RHS_mat %*% as.numeric(sigma3))   # n_eq

  # -----------------------------------------------------------------
  # Solve (A+B) ghs3 = -RHS (equilibrated for badly-scaled models)
  # -----------------------------------------------------------------
  ghs3 <- as.numeric(.solve_equilibrated(AB, -RHS_vec))
  names(ghs3) <- endo_names

  dr3$ghs3   <- ghs3
  dr3$sigma3 <- as.numeric(sigma3)
  dr3$sigma_correction <- "ghs3_only"
  dr3
}


# =====================================================================
# Phase 7c.2 / 7c.3: ghxss + ghuss (time-varying risk premium /
# state-sigma cross terms).
#
# Reference: Mutschler (2022) perturbation_solver_nonsymmetric_order3.m,
# lines 101-156.  Layout conventions used here:
#   * compound vector z in dynhr's total_cols ordering (one entry per
#     (variable name, lead_lag) pair in dyn_col_map);
#   * ghxx/ghuu/ghxu/ghxuu/ghuuu column layouts: leftmost dim slowest,
#     rightmost dim fastest in col-major flatten (same convention as the
#     existing dynhr order-2/order-3 solvers);
#   * SIGMA2 = vec(Sigma_e) col-major = (k1 slow, k2 fast).
# =====================================================================


#' Bilinear contract of the dense Hessian with two compound transfer
#' matrices: returns n_eq x (ncol(T_a) * ncol(T_b)) with column layout
#' (a slow, b fast) -- i.e., for each row e it is as.vector(t(T_a) %*%
#' H2[e,,] %*% T_b) in col-major (rows-of-product fastest).
#'
#' @noRd
.bilinear_h2 <- function(H2, T_a, T_b) {
  n_eq <- dim(H2)[1]
  out  <- matrix(0, n_eq, ncol(T_a) * ncol(T_b))
  for (e in seq_len(n_eq)) {
    out[e, ] <- as.vector(crossprod(T_a, H2[e, , ] %*% T_b))
  }
  out
}


#' Symmetrize the (zup, zXup) bilinear over the two u'-orderings.
#'
#' Input `raw` has rows indexed by equation e and cols flat-(n_u, n_X, n_u)
#' with the (d1, d2, X) layout flattened col-major (d1 fastest, then d2,
#' then X).  Returns an n_eq x (n_X * n_u * n_u) matrix with cols laid out
#' (X slow, k1 mid, k2 fast) -- matching Mutschler's F-target convention
#' -- equal to raw[d1=k1, d2=k2, X] + raw[d1=k2, d2=k1, X].
#'
#' @noRd
.sym_up_pair <- function(raw, n_X, n_u) {
  n_eq <- nrow(raw)
  out  <- matrix(0, n_eq, n_X * n_u * n_u)
  for (e in seq_len(n_eq)) {
    A   <- array(raw[e, ], dim = c(n_u, n_u, n_X))         # [d1, d2, X]
    sym <- A + aperm(A, c(2L, 1L, 3L))                      # [d1, d2, X] symmetric
    out[e, ] <- as.vector(aperm(sym, c(2L, 1L, 3L)))        # → [k2, k1, X] flat
  }
  out
}


#' Pass-through for .contract_h3(T_X, T_up, T_up) Fxupup output.
#'
#' The actual output layout of .contract_h3 (due to the outer product order
#' `tc_m %o% outer_ab` with tb_j %o% ta_i inside) is:
#'   column = (X-1)*n_u^2 + (up1-1)*n_u + up2
#' i.e. (up2 fastest, up1 middle, X slowest) — which IS already the
#' Mutschler kron(zX, kron(zup, zup)) convention.  No permutation is needed.
#' The prior implementation incorrectly assumed X was fastest and applied an
#' array+aperm that was a no-op for n_u=1 but corrupted the result for n_u>1.
#'
#' @noRd
.h3_xupup_to_mutschler <- function(M, n_X, n_u) {
  M   # contract_h3 already outputs in (up2 fast, up1 mid, X slow) = Mutschler
}


#' Solve for the order-2 stochastic-SS correction `ghss` (Dynare `ghs2`).
#'
#' Canonical solve following Mutschler (2022)
#' `perturbation_solver_nonsymmetric_order3.m` eqs. 95-99; matches
#' Dynare's `oo_.dr.ghs2`:
#'
#'   ghs2 = -(A_L + f_+)^{-1}
#'           [ f_+ · ghuu · vec(Σ)  +  H2(T_up, T_up) · vec(Σ) ]
#'
#' where `T_up = dz/du'` is the jumper block of the compound-variable
#' partial.  Called by `solve_perturbation_order2()` to populate
#' `dr2$ghss`, and reused by `solve_sigma_cross()` via that field.
#'
#' Lives in this file (rather than `solve-perturbation-order2.R`)
#' because it depends on `.bilinear_h2`, which is part of the
#' sigma-cross tensor toolkit.
#'
#' @noRd
#' Two-sided (row + column) equilibrated linear solve.
#'
#' \code{solve()} on a badly-scaled but nonsingular matrix is singular to
#' machine precision (e.g. Caldara_et_al_2012: a value-function variable with
#' SS ~ 2.27e6 makes \code{A_L + fp} span ~1e12, rcond ~5e-19). Scaling both
#' rows and columns to unit max-norm recovers several orders of conditioning
#' (LAPACK's DGEEQU idea). The plain solve is tried first and returned
#' unchanged when it succeeds, so well-conditioned models are byte-identical;
#' equilibration only engages when \code{solve()} fails or returns non-finite.
#' @noRd
.solve_equilibrated <- function(A, b) {
  out <- tryCatch(solve(A, b), error = function(e) NULL)
  if (!is.null(out) && all(is.finite(out))) return(out)
  r <- apply(abs(A), 1L, max); r[!is.finite(r) | r == 0] <- 1   # row scales
  Ar <- A / r
  cc <- apply(abs(Ar), 2L, max); cc[!is.finite(cc) | cc == 0] <- 1  # col scales
  Ac <- sweep(Ar, 2L, cc, `/`)
  br <- b / r
  y  <- tryCatch(solve(Ac, br), error = function(e) qr.solve(Ac, br))
  y / cc
}

.solve_ghss <- function(A_L, fp, ghuu, T_up, H2, Sigma_e) {
  n_eq <- dim(H2)[1]
  n_u  <- ncol(T_up)
  vS   <- as.numeric(Sigma_e)
  term_ghuu <- fp %*% ghuu %*% vS                         # n_eq vector
  term_h2   <- .bilinear_h2(H2, T_up, T_up) %*% vS        # n_eq vector
  AB <- A_L + fp
  as.numeric(-.solve_equilibrated(AB, term_ghuu + term_h2))
}


#' Recompute ghxss and ghuss (Levintal / Mutschler sigma-cross terms)
#'
#' Translates Mutschler (2022) perturbation_solver_nonsymmetric_order3.m
#' lines 101-156.  Requires `dr3` to already carry `ghss` (= ghs2) and
#' the deterministic third-order rules (ghxxx/ghxxu/ghxuu/ghuuu).
#'
#' Uses the model's declared shock covariance `compiled$model$Sigma_e`
#' if available, else the identity.
#'
#' @param dr3       DecisionRules3 from `solve_perturbation_order3()`
#' @param compiled  dynhr_compiled
#' @param ss        Named numeric SS vector
#' @param params    Named numeric parameter vector
#' @param Sigma_e   Optional n_u x n_u shock covariance (default: from
#'                   `compiled$model$Sigma_e` or identity).
#' @return The `dr3` object with `$ghxss` (n_endo x n_state) and `$ghuss`
#'   (n_endo x n_u) added; `$sigma_correction` flag updated.
#' @export
solve_sigma_cross <- function(dr3, compiled, ss, params, Sigma_e = NULL) {
  if (!inherits(dr3, "DecisionRules3")) {
    stop("dr3 must be a DecisionRules3 from solve_perturbation_order3().")
  }
  compiled_order <- compiled$dynamic$max_order
  if (!is.null(compiled_order) && compiled_order < 2L) {
    stop(sprintf(
      paste0("solve_sigma_cross() needs 2nd/3rd-order symbolic derivatives ",
             "but the model was compiled with max_order = %d. Re-run ",
             "compile_model(model, max_order = 2L)."), compiled_order))
  }
  endo_names <- dr3$endo_names
  exo_names  <- dr3$exo_names
  state_idx  <- dr3$state_idx
  n          <- length(endo_names)
  n_s        <- length(state_idx)
  n_u        <- length(exo_names)

  # Sigma_e contract: must match the covariance baked into dr3$ghss /
  # dr3$ghuu (i.e. the one passed to solve_perturbation_order2). If the
  # caller supplies a different Sigma_e they need a fresh order-2 solve
  # first, otherwise ghs2 (from dr3) and the sigma-cross RHS terms (built
  # below from Sigma_e) would use inconsistent covariances.
  if (is.null(Sigma_e)) {
    Sigma_e <- dr3$Sigma_e
    if (is.null(Sigma_e)) {
      stderr  <- .get_shock_stderr(compiled$model, exo_names, params)
      Sigma_e <- diag(stderr^2, n_u, n_u)
    }
  } else if (!is.null(dr3$Sigma_e) &&
             !isTRUE(all.equal(unname(as.matrix(Sigma_e)),
                               unname(as.matrix(dr3$Sigma_e)),
                               tolerance = 1e-12))) {
    stop("Sigma_e passed to solve_sigma_cross() differs from dr3$Sigma_e ",
         "(the covariance baked into dr3$ghss). Re-run ",
         "solve_perturbation_order2() with the new Sigma_e first.")
  }
  if (!is.matrix(Sigma_e) || nrow(Sigma_e) != n_u || ncol(Sigma_e) != n_u) {
    stop(sprintf("Sigma_e must be an n_u x n_u matrix (n_u = %d).", n_u))
  }
  SIGMA2 <- as.numeric(Sigma_e)   # vec col-major: (k1 slow, k2 fast)

  dyn      <- compiled$dynamic
  sys      <- extract_system_matrices(compiled, ss, params)
  has_lead <- sys$is_fwd | sys$is_mixed

  if (!any(has_lead) || n_s == 0L || n_u == 0L) {
    dr3$ghxss <- matrix(0, n, n_s)
    dr3$ghuss <- matrix(0, n, n_u)
    rownames(dr3$ghxss) <- endo_names
    rownames(dr3$ghuss) <- endo_names
    dr3$sigma_correction <- "trivial_no_jumpers_or_no_shocks"
    return(dr3)
  }

  # ----- canonical matrices ------------------------------------------
  ghx  <- dr3$ghx;   ghu  <- dr3$ghu
  ghxx <- dr3$ghxx;  ghxu <- dr3$ghxu; ghuu <- dr3$ghuu
  # dr3$ghxuu is stored in standard Kron layout (state SLOW, u1 MID, u2 FAST),
  # which is exactly the Mutschler/Dynare kron(hx, I_{u^2}) column order — no
  # conversion needed. (Earlier versions of solve_perturbation_order3 emitted
  # an internal (state FAST) layout, which is what the Mutschler converter was
  # designed to flip. After commit "Order-3 perturbation: fix internal/standard
  # col convention", the public field is already in standard layout.)
  ghxuu_m <- dr3$ghxuu
  hx   <- ghx[state_idx, , drop = FALSE]
  hu   <- ghu[state_idx, , drop = FALSE]

  S <- matrix(0, n, n_s)
  for (s in seq_along(state_idx)) S[state_idx[s], s] <- 1
  A_L <- sys$f_zero + sys$f_plus %*% ghx %*% t(S)
  fp  <- sys$f_plus

  # ----- compound-derivative builders --------------------------------
  T_x   <- {
    tm <- .build_transfer_matrices(dyn, ghx, ghu, state_idx, hx, hu,
                                    endo_names, exo_names)
    tm$T_x
  }
  T_u   <- {
    tm <- .build_transfer_matrices(dyn, ghx, ghu, state_idx, hx, hu,
                                    endo_names, exo_names)
    tm$T_u
  }
  T_up   <- .build_T_up  (dyn, ghu,  endo_names, exo_names, has_lead)
  W_xup  <- .build_W_xup (dyn, ghxu, hx, endo_names, exo_names, has_lead)
  W_uup  <- .build_W_uup (dyn, ghxu, hu, endo_names, exo_names, has_lead)
  W_upup <- .build_W_upup(dyn, ghuu, endo_names, exo_names, has_lead)

  # ----- dense Hessian + sparse H3 -----------------------------------
  dy_ss <- .build_dy_ss_o2(compiled, ss)
  H2    <- .compute_model_hessian_symbolic(compiled, dy_ss, params, ss)
  if (is.null(H2)) stop("Symbolic Hessian2 unavailable; required for ghxss/ghuss.")
  H3    <- .compute_model_hessian3_symbolic(compiled, dy_ss, params, ss)
  H_ord <- .reorder_derivatives_to_decl(H2, H3, compiled, n,
                                        label = "sigma-cross derivative",
                                        eq_to_decl = sys$eq_to_decl)
  H2    <- H_ord$H2
  H3    <- H_ord$H3
  n_eq  <- dyn$n_eq

  # ----- reuse dr3$ghss (already canonical Mutschler/Dynare) ---------
  # solve_perturbation_order2() populates dr3$ghss via .solve_ghss with
  # the same Sigma_e (contract enforced above); no need to recompute.
  if (is.null(dr3$ghss)) {
    stop("dr3$ghss is missing; solve_sigma_cross() requires the order-2 ",
         "stochastic-SS correction (run solve_perturbation_order2() first).")
  }
  ghs2        <- dr3$ghss
  names(ghs2) <- endo_names
  ghs2_state  <- ghs2[state_idx]

  T_ss <- .build_T_ss(dyn, ghx, ghs2, state_idx, endo_names, exo_names,
                       has_lead)

  # ==================================================================
  # Fxupup = fp · ghxuu · kron(hx, I_{u^2})
  #        + d2f · kron(zx, zupup)
  #        + d2f · kron(zup, zxup) · (P_up1_x1up2 + P_up2_x1up1)
  #        + d3f · kron(zx, zup, zup)
  # Result layout: n_eq × (n_s × n_u^2), cols (a1 slow, k1 mid, k2 fast).
  # ==================================================================
  K_hx_Iu2 <- hx %x% diag(n_u * n_u)         # (n_s · n_u²) × (n_s · n_u²)
  Fx_term1 <- fp %*% ghxuu_m %*% K_hx_Iu2    # n_eq x (n_s * n_u^2)

  # bilinear_h2(H2, T_a, T_b) gives cols (T_a FAST, T_b SLOW).
  # Target: (a1 state SLOW, k12 exo-pair FAST) → must put W_upup first, T_x second.
  Fx_term2a <- .bilinear_h2(H2, W_upup, T_x)              # cols (u² fast, a1 slow)
  Fx_raw    <- .bilinear_h2(H2, T_up,   W_xup)            # cols (d1 fast, (a1,d2) slow)
  Fx_term2b <- .sym_up_pair(Fx_raw, n_X = n_s, n_u = n_u)

  if (!is.null(H3) && length(H3$triplets) > 0L) {
    # .contract_h3(H3, Ta, Tb, Tc) gives cols (ka SLOW, kb MID, kc FAST).
    # Target (a1 slow, k1 mid, k2 fast): Ta=T_x (a1 slow), Tb=T_up (k1 mid),
    # Tc=T_up (k2 fast).
    Fx_term3 <- .h3_xupup_to_mutschler(
      .contract_h3(H3, T_x, T_up, T_up, n_eq), n_s, n_u)
  } else {
    Fx_term3 <- matrix(0, n_eq, n_s * n_u * n_u)
  }

  Fxupup <- Fx_term1 + Fx_term2a + Fx_term2b + Fx_term3

  # ==================================================================
  # ghxss RHS
  # ==================================================================
  # term A: fp · ghxx · kron(hx, ghs2_state)
  ghs2_state_col <- matrix(ghs2_state, n_s, 1L)
  RHS_xss_A      <- fp %*% ghxx %*% (hx %x% ghs2_state_col)   # n_eq × n_s

  # term B: d2f · kron(zx, zss).  zss is (total_cols × 1), so cols = n_s.
  RHS_xss_B <- .bilinear_h2(H2, T_x, T_ss)   # cols (a1 slow, 1) ⇒ n_s

  # term C: Fxupup · kron(I_s, SIGMA2)
  # Per equation e, slab = matrix(Fxupup[e, ], n_u^2, n_s) col-major
  # ⇒ slab[k_idx, a1] · SIGMA2[k_idx] summed over k_idx = vec result of length n_s.
  RHS_xss_C <- matrix(0, n_eq, n_s)
  for (e in seq_len(n_eq)) {
    slab <- matrix(Fxupup[e, ], n_u * n_u, n_s)
    RHS_xss_C[e, ] <- as.numeric(SIGMA2 %*% slab)
  }

  RHS_xss <- RHS_xss_A + RHS_xss_B + RHS_xss_C

  # ==================================================================
  # Sylvester solve: A_L · X + fp · X · hx = -RHS_xss  (k = 1)
  # Use the compact Schur-based solver instead of forming the dense
  # n·n_s × n·n_s Kronecker matrix, which is singular to machine precision
  # under extreme scaling (Caldara: rcond ~9e-20). Machine-precision
  # identical on well-conditioned models (dense-Schur safety net).
  # ==================================================================
  ghxss <- matrix(.solve_kron_compact(A_L, fp, hx, k = 1L, RHS = -RHS_xss),
                  n, n_s)

  # ==================================================================
  # Fuupup = fp · ghxuu · kron(hu, I_{u^2})
  #        + d2f · kron(zu, zupup)
  #        + d2f · kron(zup, zuup) · (P_up1_u1up2 + P_up2_u1up1)
  #        + d3f · kron(zu, zup, zup)
  # Cols (l slow, k1 mid, k2 fast), total n_u³.
  # ==================================================================
  K_hu_Iu2 <- hu %x% diag(n_u * n_u)
  Fu_term1 <- fp %*% ghxuu_m %*% K_hu_Iu2     # n_eq x n_u^3

  # Target (l slow, k1 mid, k2 fast): put W_upup first (k12 fast), T_u second (l slow).
  Fu_term2a <- .bilinear_h2(H2, W_upup, T_u)              # cols (u² fast, l slow)
  Fu_raw    <- .bilinear_h2(H2, T_up,   W_uup)            # cols (d1 fast, (l,d2) slow)
  Fu_term2b <- .sym_up_pair(Fu_raw, n_X = n_u, n_u = n_u)

  if (!is.null(H3) && length(H3$triplets) > 0L) {
    # Target (l slow, k1 mid, k2 fast): Ta=T_u (l slow), Tb=T_up (k1 mid),
    # Tc=T_up (k2 fast).
    Fu_term3 <- .h3_xupup_to_mutschler(
      .contract_h3(H3, T_u, T_up, T_up, n_eq), n_u, n_u)
  } else {
    Fu_term3 <- matrix(0, n_eq, n_u * n_u * n_u)
  }
  Fuupup <- Fu_term1 + Fu_term2a + Fu_term2b + Fu_term3

  # ==================================================================
  # ghuss RHS
  # ==================================================================
  # term A: fp · (ghxx · kron(hu, ghs2_state) + ghxss · hu)
  RHS_uss_A <- fp %*% (
    ghxx %*% (hu %x% ghs2_state_col) +
      ghxss %*% hu
  )

  # term B: d2f · kron(zu, zss)
  RHS_uss_B <- .bilinear_h2(H2, T_u, T_ss)

  # term C: Fuupup · kron(I_u, SIGMA2)
  RHS_uss_C <- matrix(0, n_eq, n_u)
  for (e in seq_len(n_eq)) {
    slab <- matrix(Fuupup[e, ], n_u * n_u, n_u)
    RHS_uss_C[e, ] <- as.numeric(SIGMA2 %*% slab)
  }

  RHS_uss <- RHS_uss_A + RHS_uss_B + RHS_uss_C

  # Direct solve A_L · ghuss = -RHS_uss (equilibrated for badly-scaled models)
  ghuss <- .solve_equilibrated(A_L, -RHS_uss)
  if (!is.matrix(ghuss)) ghuss <- matrix(ghuss, n, n_u)

  rownames(ghxss) <- endo_names
  rownames(ghuss) <- endo_names
  colnames(ghxss) <- endo_names[state_idx]
  colnames(ghuss) <- exo_names

  dr3$ghxss <- ghxss
  dr3$ghuss <- ghuss
  dr3$Sigma_e_used <- Sigma_e
  dr3$sigma_correction <- if (!is.null(dr3$ghs3) && any(dr3$ghs3 != 0))
                          "full" else "ghxss_ghuss_only"
  dr3
}
