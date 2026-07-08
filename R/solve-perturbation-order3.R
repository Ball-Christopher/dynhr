## R/solve-perturbation-order3.R
## --------------------------------------------------------------------------
## Third-order perturbation solver for DSGE models.
##
## Implements the DETERMINISTIC (sigma^0) third-order policy terms using the
## Binning (2013) Faà di Bruno matrix chain rule:
##
##   y_t = ... + (1/6) g_xxx (x \otimes x \otimes x)
##             + (1/2) g_xxu (x \otimes x \otimes u)
##             + (1/2) g_xuu (x \otimes u \otimes u)
##             + (1/6) g_uuu (u \otimes u \otimes u)
##
## The sigma-correction terms (ghxss, ghuss, ghsss -- time-varying risk
## premia) are NOT yet implemented; they require the augmented-state /
## Levintal (2017) recursion and are slated for a follow-up phase.
##
## Algorithm:
##   1. Reuse order-1/order-2 building blocks: ghx, ghu, ghxx, ghxu, ghuu,
##      hx, hu, hxx, hxu, huu, A_L, fp
##   2. Build transfer matrices T_x, T_u (from order-2 solver)
##   3. Compute the symbolic third-order Hessian tensor F_www (sparse) at SS
##   4. Compute the Faà di Bruno forcing term Phi_xxx via:
##        F_www [T_x, T_x, T_x]                       (direct cubic term)
##      + F_ww  [W_xx, T_x]  summed over 3 pair-singleton permutations
##      + f_+ contributions from KNOWN parts of d^3 y_{t+1}/dx^3
##        (i.e. ghxx [(hxx \otimes hx) sym-permutations])
##   5. Solve Kronecker system  A_L ghxxx + f_+ ghxxx (hx \otimes hx \otimes hx) = -Phi_xxx
##   6. Direct solves (given ghxxx) for ghxxu, ghxuu, ghuuu
##
## References:
##   Binning, A. (2013). Solving Second and Third-Order Approximations to
##     DSGE Models: A Recursive Sylvester Equation Solution. Norges Bank
##     WP 2013/18. -- main algorithmic reference, esp. eq (35)-(43).
##   Andreasen, Fernandez-Villaverde, Rubio-Ramirez (2018). The Pruned
##     State-Space System for Non-Linear DSGE Models. ReStud 85(1).
##   Levintal, O. (2017). Fifth-Order Perturbation Solution to DSGE
##     Models. JEDC 80, 1-16.
## --------------------------------------------------------------------------


# =====================================================================
# Sparse third-order Hessian helpers
# =====================================================================

#' Build the symbolic 4-tensor F_www at the SS from sparse triplets.
#'
#' Returns a list with:
#'   triplets: list of (eq, col1, col2, col3, val) -- canonical ordering
#'             c1 <= c2 <= c3 (Schwarz symmetry; caller expands the orbit
#'             of up to 6 permutations).
#'   n_eq, total_cols: dimensions of the implied dense tensor
#'
#' Lazily evaluated -- never builds the full n_eq * total_cols^3 array.
#'
#' @noRd
.compute_model_hessian3_symbolic <- function(compiled, dy_ss, params, ss) {
  dyn <- compiled$dynamic
  if (is.null(dyn$hessian3_fn) || is.null(dyn$hess3_triplets)) return(NULL)
  if (dyn$n_hess3 == 0L) {
    return(list(triplets = list(), values = numeric(0),
                n_eq = dyn$n_eq, total_cols = dyn$total_cols))
  }
  ## Prefer the C++ value-vector tape over the interpreted hessian3_fn closure;
  ## fall back to the closure when the tape is absent or errors. Bit-parity to
  ## machine precision is asserted by test-hess-tape-parity.R.
  values <- .eval_triplet_tape(dyn, dyn$hess3_tape, dy_ss, params, ss)
  if (is.null(values))
    values <- tryCatch(dyn$hessian3_fn(dy_ss, params, ss),
                      error = function(e) {
                        warning(sprintf(
                          "Symbolic Hessian3 evaluation failed: %s.",
                          conditionMessage(e)))
                        NULL
                      })
  if (is.null(values)) return(NULL)
  list(triplets = dyn$hess3_triplets, values = values,
       n_eq = dyn$n_eq, total_cols = dyn$total_cols)
}


#' Enumerate the symmetry orbit of a (c1, c2, c3) triplet.
#'
#' Returns a list of unique (c1, c2, c3) orderings: 1 for (a,a,a),
#' 3 for (a,a,b) or (a,b,b), 6 for (a,b,c) all distinct.
#'
#' @noRd
.orbit_3 <- function(c1, c2, c3) {
  if (c1 == c2 && c2 == c3) {
    list(c(c1, c2, c3))
  } else if (c1 == c2) {
    list(c(c1, c1, c3), c(c1, c3, c1), c(c3, c1, c1))
  } else if (c2 == c3) {
    list(c(c1, c2, c2), c(c2, c1, c2), c(c2, c2, c1))
  } else {
    # All distinct: 6 permutations
    list(c(c1, c2, c3), c(c1, c3, c2), c(c2, c1, c3),
         c(c2, c3, c1), c(c3, c1, c2), c(c3, c2, c1))
  }
}


# =====================================================================
# Faà di Bruno: third-order forcing terms
# =====================================================================

#' Contract F_www with three rank-1 transfer matrices.
#'
#' For each equation e, compute the n_a * n_b * n_c vector
#'   Phi[e, (a,b,c)] = sum_{i,j,k} F_www[e,i,j,k] * Ta[i,a] * Tb[j,b] * Tc[k,c]
#'
#' Iterates over the canonical sparse triplets and expands the Schwarz
#' orbit explicitly so the resulting Phi is correctly summed over all
#' permutations of (i,j,k).
#'
#' Indexing convention (Kronecker: a slowest, c fastest):
#'   col = (a-1)*n_b*n_c + (b-1)*n_c + c
#'
#' This follows from the outer-product accumulation:
#'   outer_ab = as.vector(Tb[j,] %o% Ta[i,])  -- b fastest, a slowest
#'   contrib  = as.vector(Tc[m,] %o% outer_ab) -- c fastest, then b, then a
#'
#' @noRd
.contract_h3 <- function(h3, Ta, Tb, Tc, n_eq) {
  n_a <- ncol(Ta); n_b <- ncol(Tb); n_c <- ncol(Tc)
  Phi <- matrix(0, n_eq, n_a * n_b * n_c)

  if (length(h3$triplets) == 0L) return(Phi)

  for (k in seq_along(h3$triplets)) {
    t   <- h3$triplets[[k]]
    val <- h3$values[k]
    if (val == 0) next
    e   <- t$eq

    # Expand the FULL Schwarz symmetry orbit of (col1, col2, col3).
    # H3 is stored with canonical (c1 <= c2 <= c3) but is fully symmetric
    # in all three positions.  For any Ta, Tb, Tc (identical or not), we
    # must sum over all distinct permutations of (c1,c2,c3) because the
    # matrices assigned to each slot can differ — e.g. for (T_x, T_up, T_up)
    # the permutation (c3,c1,c1) places the "singleton" index into the T_x
    # slot, contributing T_x[c3]*T_up[c1]^2, which was previously missing
    # for the identical(Tb,Tc) branch.  Using .orbit_3 unconditionally is
    # always correct and fixes the incomplete-orbit bug.
    c1 <- t$col1; c2 <- t$col2; c3 <- t$col3
    orbit <- .orbit_3(c1, c2, c3)
    for (perm in orbit) {
      i <- perm[1]; j <- perm[2]; m <- perm[3]
      # outer product accumulation: Phi[e, (a,b,c)] += val * Ta[i,a] * Tb[j,b] * Tc[m,c]
      ta_i <- Ta[i, ]   # length n_a
      tb_j <- Tb[j, ]   # length n_b
      tc_m <- Tc[m, ]   # length n_c
      # Kronecker outer product: tb_j %o% ta_i is (n_b × n_a), as.vector = b fastest
      outer_ab <- as.vector(tb_j %o% ta_i)               # (a-1)*n_b + b  (b fast)
      # tc_m %o% outer_ab is (n_c × n_a*n_b); as.vector = c fastest
      contrib  <- as.vector(tc_m %o% outer_ab)            # (a-1)*n_b*n_c+(b-1)*n_c+c
      Phi[e, ] <- Phi[e, ] + val * contrib
    }
  }
  Phi
}


#' Build 2nd-order transfer matrices W_xx, W_xu, W_uu at the SS.
#'
#' Each row of W_** is a 2nd-derivative of the compound vector w with
#' respect to (state, state), (state, shock), or (shock, shock).  Used
#' as the "pair" block in the Faà di Bruno chain rule for 3rd-order.
#'
#' For block (variable nm, lead/lag ll):
#'   ll = -1 (lag, state only): W_xx, W_xu, W_uu are all 0
#'      (w_c = x_{s_c} is linear in x_{t-1}; no 2nd derivative)
#'   ll =  0 (current): W_xx = ghxx[j,], W_xu = ghxu[j,], W_uu = ghuu[j,]
#'   ll = +1 (lead):
#'      W_xx = ghxx[j,] (hx \otimes hx) + ghx[j,] hxx     -- size n_s^2
#'      W_xu = ghxx[j,] (hx \otimes hu) + ghx[j,] hxu     -- size n_s*n_u
#'      W_uu = ghxx[j,] (hu \otimes hu) + ghx[j,] huu     -- size n_u^2
#'      (the +0.5 ghss * dsigma^2 piece is omitted for deterministic 3rd order)
#'   ll = shock-block (u_t): all 0 (shock is independent of x_{t-1})
#'
#' @noRd
.build_W2_matrices <- function(dyn, ghx, ghu, ghxx, ghxu, ghuu,
                                hx, hu, hxx, hxu, huu,
                                state_idx, endo_names, exo_names) {
  total_cols <- dyn$total_cols
  n_s        <- length(state_idx)
  n_u        <- length(exo_names)
  dcm        <- dyn$dyn_col_map

  W_xx <- matrix(0, total_cols, n_s * n_s)
  W_xu <- matrix(0, total_cols, n_s * n_u)
  W_uu <- matrix(0, total_cols, n_u * n_u)

  # Precompute Kronecker products for lead-block
  hx_kron_hx <- hx %x% hx   # n_s^2 x n_s^2
  # W_xu cols must match ghxu convention: (state FAST, exo SLOW).
  # (hu %x% hx) gives cols (exo SLOW, state FAST) ← correct.
  # (hx %x% hu) would give (state SLOW, exo FAST) ← wrong for n_u > 1.
  hu_kron_hx <- hu %x% hx   # n_s^2 x (n_u * n_s)
  hu_kron_hu <- hu %x% hu   # n_u^2 x n_u^2

  # Build column index reverse-map
  for (kc in seq_len(nrow(dcm))) {
    c  <- dcm$col[kc]
    nm <- dcm$name[kc]
    ll <- dcm$lead_lag[kc]
    is_exo <- nm %in% exo_names

    if (is_exo || ll == -1L) {
      next   # all zero
    }

    j <- which(endo_names == nm)
    if (length(j) != 1L) next

    if (ll == 0L) {
      if (n_s > 0L)             W_xx[c, ] <- ghxx[j, ]
      if (n_s > 0L && n_u > 0L) W_xu[c, ] <- ghxu[j, ]
      if (n_u > 0L)             W_uu[c, ] <- ghuu[j, ]
    } else if (ll == 1L) {
      # d^2 y_{t+1}/dx^2 = ghxx (hx \otimes hx) + ghx hxx
      if (n_s > 0L) {
        W_xx[c, ] <- as.numeric(ghxx[j, , drop = FALSE] %*% hx_kron_hx) +
                     as.numeric(ghx [j, , drop = FALSE] %*% hxx)
      }
      if (n_s > 0L && n_u > 0L) {
        # d^2 y_{t+1}/dx du = ghxx (hu \otimes hx) + ghx hxu
        W_xu[c, ] <- as.numeric(ghxx[j, , drop = FALSE] %*% hu_kron_hx) +
                     as.numeric(ghx [j, , drop = FALSE] %*% hxu)
      }
      if (n_u > 0L) {
        # d^2 y_{t+1}/du^2 = ghxx (hu \otimes hu) + ghx huu
        W_uu[c, ] <- as.numeric(ghxx[j, , drop = FALSE] %*% hu_kron_hu) +
                     as.numeric(ghx [j, , drop = FALSE] %*% huu)
      }
    }
  }

  list(W_xx = W_xx, W_xu = W_xu, W_uu = W_uu)
}


#' Permute the "abc"-column index of an (n_eq x n_a*n_b*n_c) matrix.
#'
#' Each row of M is a flattened 3-tensor T[a,b,c] in column-major order
#' (a fastest, c slowest).  Returns a new matrix M' whose rows are
#' aperm(T, perm) -- i.e., for perm = c(1,3,2) the new row at (a,b,c)
#' equals the original row at (a, c, b).
#'
#' All three modes must have the same length n (works for n_s^3 case).
#'
#' @noRd
.permute_cube_cols <- function(M, n, perm) {
  if (n == 0L) return(M)
  n_eq <- nrow(M)
  out  <- matrix(0, n_eq, n^3)
  for (e in seq_len(n_eq)) {
    A <- array(M[e, ], dim = c(n, n, n))
    B <- aperm(A, perm)
    out[e, ] <- as.numeric(B)
  }
  out
}


#' Symmetrize an n_eq x n^3 tensor over all three modes.
#' @noRd
.symmetrize_cube_cols <- function(M, n) {
  if (n <= 1L || ncol(M) == 0L) return(M)
  perms <- list(c(1L, 2L, 3L), c(1L, 3L, 2L), c(2L, 1L, 3L),
                c(2L, 3L, 1L), c(3L, 1L, 2L), c(3L, 2L, 1L))
  out <- matrix(0, nrow(M), n^3)
  for (e in seq_len(nrow(M))) {
    A <- array(M[e, ], dim = c(n, n, n))
    S <- array(0, dim = c(n, n, n))
    for (perm in perms) S <- S + aperm(A, perm)
    out[e, ] <- as.numeric(S / length(perms))
  }
  out
}


#' Symmetrize ghxxu over the two state dimensions.
#' @noRd
.symmetrize_xxu_cols <- function(M, n_s, n_u) {
  if (n_s <= 1L || n_u == 0L || ncol(M) == 0L) return(M)
  out <- matrix(0, nrow(M), n_s * n_s * n_u)
  for (e in seq_len(nrow(M))) {
    A <- array(M[e, ], dim = c(n_s, n_s, n_u))
    out[e, ] <- as.numeric((A + aperm(A, c(2L, 1L, 3L))) / 2)
  }
  out
}


#' Symmetrize ghxuu over the two shock dimensions.
#' @noRd
.symmetrize_xuu_cols <- function(M, n_s, n_u) {
  if (n_u <= 1L || n_s == 0L || ncol(M) == 0L) return(M)
  out <- matrix(0, nrow(M), n_s * n_u * n_u)
  for (e in seq_len(nrow(M))) {
    A <- array(M[e, ], dim = c(n_s, n_u, n_u))
    out[e, ] <- as.numeric((A + aperm(A, c(1L, 3L, 2L))) / 2)
  }
  out
}


#' Symmetrize ghuuu over all three shock dimensions.
#' @noRd
.symmetrize_uuu_cols <- function(M, n_u) {
  .symmetrize_cube_cols(M, n_u)
}


# =====================================================================
# Main third-order solver
# =====================================================================



#' Solve the deterministic third-order perturbation of a DSGE model
#'
#' Given the first-/second-order decision rules in `dr2`, computes the
#' third-order terms ghxxx, ghxxu, ghxuu, ghuuu.  Does NOT compute the
#' sigma-correction terms (ghxss, ghuss, ghsss) -- those terms (which
#' capture time-varying risk premia and the third-cumulant correction)
#' are deferred to a follow-up phase that implements the Levintal (2017)
#' augmented-state recursion.
#'
#' @param model    dynhr_mod object from parse_mod()
#' @param compiled dynhr_compiled from compile_model()
#' @param ss       Named numeric steady state vector
#' @param params   Named numeric parameter vector
#' @param dr2      Second-order DecisionRules2 object
#' @param verbose  Print progress
#' @param solver_method Deprecated and ignored (retained for back-compat). The
#'   order-3 Kronecker system now uses the shared `.solve_kron_compact` solver
#'   (the same one orders 4/5 use), which self-selects an eigenbasis fast path
#'   or a robust real-Schur dense fallback; there is no longer a user choice.
#' @param backend Computation backend: `"auto"` selects `"fd"` for
#'   tractable model sizes (when eligible) and `"symbolic"` otherwise;
#'   `"symbolic"` uses the Faà di Bruno chain-rule expansion; `"fd"` uses
#'   a finite-difference oracle (requires single-period leads/lags and no
#'   AUX variables).
#' @param sparse Logical or `NULL` (default). Controls the ghxxx Kronecker
#'   solve. `FALSE` uses the default dense-eligible `.solve_kron_compact`
#'   solver (eigenbasis fast path + dense real-Schur fallback) — unchanged
#'   numerical output. `TRUE` uses the memory-light complex-Schur Kronecker
#'   Bartels–Stewart solver `.solve_kron_compact_sparse`, which never
#'   materialises the ns^3 × ns^3 Kronecker matrix and so breaks the dense
#'   order-3 wall on high-dimensional-but-sparse state blocks (e.g. the
#'   emitted finite HANK, n_s ≈ 32). `NULL` (auto) turns the sparse route on
#'   automatically once `n_state^3` exceeds `sparse_threshold`. The sparse
#'   route is bit-parity to the dense path on well-conditioned models and
#'   residual-verified (with dense fallback) otherwise.
#' @param sparse_threshold Integer. When `sparse = NULL`, the sparse ghxxx
#'   solve is used iff `n_state^3 >= sparse_threshold` (default 8000, i.e.
#'   n_state >= 20). Ignored when `sparse` is `TRUE`/`FALSE`.
#' @return A DecisionRules3 object extending DecisionRules2 with fields
#'   ghxxx, ghxxu, ghxuu, ghuuu (and all lower-order fields preserved)
#'
#' @references
#'   Binning, A. (2013). Underidentified SVAR models: A framework for combining
#'     short-run and long-run restrictions. \emph{Norges Bank Working Paper}.
#'   Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez, J. F. (2018).
#'     The pruned state-space system for non-linear DSGE models.
#'     \emph{Review of Economic Studies}, 85(1), 1-49.
#' @export
solve_perturbation_order3 <- function(model, compiled, ss, params, dr2,
                                       verbose = FALSE,
                                       solver_method = c("auto", "direct", "qz"),
                                       backend = c("auto", "symbolic", "fd"),
                                       sparse = NULL,
                                       sparse_threshold = 8000L) {
  if (!inherits(dr2, "DecisionRules2")) {
    stop("dr2 must be a DecisionRules2 object from solve_perturbation_order2().")
  }
  if (!isTRUE(dr2$bk_satisfied)) {
    stop("First-order Blanchard-Kahn not satisfied; cannot solve to third order.")
  }
  compiled_order <- compiled$dynamic$max_order
  if (!is.null(compiled_order) && compiled_order < 2L) {
    stop(sprintf(
      paste0("solve_perturbation_order3() needs 2nd/3rd-order symbolic ",
             "derivatives but the model was compiled with max_order = %d. ",
             "Re-run compile_model(model, max_order = 2L)."), compiled_order))
  }
  solver_method <- match.arg(solver_method)
  backend <- match.arg(backend)

  # FD backend: use the finite-difference oracle solver for small models
  # (more accurate than the symbolic Faà di Bruno chain rule, which has
  # row-permutation bugs in the Phi_xxx forcing term — see
  # ORDER3_INVESTIGATION.md). Auto-selects FD for tractable model size.
  n_endo <- length(dr2$endo_names)
  n_state <- length(dr2$state_idx)
  n_exo <- length(dr2$exo_names)
  dyn <- compiled$dynamic
  has_aux <- any(grepl("^AUX_", dr2$endo_names))
  has_multi_lead <- any(dyn$dyn_col_map$lead_lag > 1L) ||
                    any(dyn$dyn_col_map$lead_lag < -1L)
  fd_eligible <- !has_aux && !has_multi_lead
  fd_cost <- n_endo * (n_state + n_exo)^3
  use_fd <- (backend == "fd")
  if (use_fd) {
    if (!fd_eligible) stop("FD backend requires single-period leads/lags and no AUX variables.")
    if (verbose) cat("Using FD backend for order-3 (fd_cost=", fd_cost, ").\n")
    return(.solve_perturbation_order3_fd(model, compiled, ss, params, dr2,
                                          verbose = verbose))
  }
  if (backend == "fd") {
    stop("FD backend requires single-period leads/lags and no AUX variables.")
  }

  dyn        <- compiled$dynamic
  n          <- length(dr2$endo_names)
  n_s        <- length(dr2$state_idx)
  n_u        <- length(dr2$exo_names)
  state_idx  <- dr2$state_idx
  endo_names <- dr2$endo_names
  exo_names  <- dr2$exo_names

  ghx  <- dr2$ghx;  ghu  <- dr2$ghu
  ghxx <- dr2$ghxx; ghxu <- dr2$ghxu; ghuu <- dr2$ghuu

  hx  <- ghx[state_idx, , drop = FALSE]
  hu  <- ghu[state_idx, , drop = FALSE]
  hxx <- ghxx[state_idx, , drop = FALSE]
  hxu <- ghxu[state_idx, , drop = FALSE]
  huu <- ghuu[state_idx, , drop = FALSE]

  if (n_s == 0L) {
    if (verbose) message("No state variables; third-order x-terms are zero.")
    return(.trivial_dr3(dr2))
  }

  if (isTRUE(model$model_options$linear)) {
    if (verbose) message("Linear model: all third-order terms are exactly zero; skipping Kronecker solve.")
    return(.linear_dr3(dr2))
  }

  # First-order system matrices (for A_L, fp)
  sys <- extract_system_matrices(compiled, ss, params)
  f0  <- sys$f_zero
  fp  <- sys$f_plus

  S   <- matrix(0, n, n_s)
  for (s in seq_len(n_s)) S[state_idx[s], s] <- 1
  A_L <- f0 + fp %*% ghx %*% t(S)

  if (verbose) {
    cat("Third-order perturbation (deterministic):\n")
    cat(sprintf("  n_endo=%d  n_state=%d  n_exo=%d  n_hess3=%d\n",
                n, n_s, n_u, dyn$n_hess3 %||% 0L))
    cat(sprintf("  Kronecker system size: %d x %d\n",
                n * n_s^3, n * n_s^3))
  }

  # Build dy at SS and transfer matrices
  dy_ss <- .build_dy_ss_o2(compiled, ss)
  tm    <- .build_transfer_matrices(dyn, ghx, ghu, state_idx, hx, hu,
                                     endo_names, exo_names)
  T_x   <- tm$T_x
  T_u   <- tm$T_u

  # Second-order Hessian (sparse triplets, also needed)
  H2 <- .compute_model_hessian_symbolic(compiled, dy_ss, params, ss)
  if (is.null(H2)) stop("Symbolic Hessian2 unavailable; required for order 3.")

  # Third-order Hessian (sparse)
  H3 <- .compute_model_hessian3_symbolic(compiled, dy_ss, params, ss)
  if (is.null(H3)) {
    if (verbose) cat("  Hessian3 unavailable (log-linear model). Phi_direct = 0.\n")
  }

  # ----------------------------------------------------------------
  # Reorder Hessian2 rows from compiled-equation order to
  # declaration-variable order, matching A_L/fp row ordering.
  # Hessian3 triplets use the SAME equation numbering convention
  # as Hessian2, so they must also be remapped for consistency.
  # Without this remapping, Phi_direct (H3) and Phi_pair (H2)
  # contribute to different equations, producing incorrect results.
  #
  # CRITICAL: reuse the SAME mapping that extract_system_matrices()
  # used to reorder A_L/fp (sys$eq_to_decl), NOT a freshly recomputed
  # .build_eq_to_decl(model).  The latter can differ (it omits the
  # f_zero Jacobian-matching step) and silently misaligns the Hessian
  # forcing rows against A_L -- producing spurious order-3 forcing in
  # equations that are actually linear (e.g. BP2014's `r = rbar +
  # eps_tb + eps_r`).  This mirrors the order-2 solver, which already
  # reuses sys$eq_to_decl.  See [[eq-to-decl-consistency-invariant]].
  # ----------------------------------------------------------------
  if (!is.null(compiled$model$equations)) {
    eq_to_decl <- sys$eq_to_decl %||% .build_eq_to_decl(compiled$model)
    if (all(eq_to_decl > 0L) && !identical(eq_to_decl, seq_len(n))) {
      perm <- order(eq_to_decl)
      H2_perm <- array(0, dim = dim(H2))
      for (k in seq_len(n)) H2_perm[k, , ] <- H2[perm[k], , ]
      H2 <- H2_perm
      # Also remap H3 equation indices to match
      if (!is.null(H3) && length(H3$triplets) > 0L) {
        inv_perm <- order(perm)
        for (k in seq_along(H3$triplets)) {
          H3$triplets[[k]]$eq <- inv_perm[H3$triplets[[k]]$eq]
        }
      }
      if (verbose) cat("  Hessian rows reordered: compiled -> declaration order.\n")
    }
  }

  # Build second-order transfer matrices
  W2 <- .build_W2_matrices(dyn, ghx, ghu, ghxx, ghxu, ghuu,
                            hx, hu, hxx, hxu, huu,
                            state_idx, endo_names, exo_names)
  W_xx <- W2$W_xx
  W_xu <- W2$W_xu
  W_uu <- W2$W_uu

  # ----------------------------------------------------------------
  # Phi_xxx: forcing for ghxxx (state-state-state)
  # = F_www[T_x, T_x, T_x]                              (direct cubic)
  # + F_ww[W_xx, T_x] summed over 3 pair-singleton permutations
  # ----------------------------------------------------------------
  if (verbose) cat("  Computing Phi_xxx forcing term...\n")

  ns3 <- n_s^3

  # Direct cubic contribution from F_www
  if (!is.null(H3) && length(H3$triplets) > 0L) {
    Phi_direct <- .contract_h3(H3, T_x, T_x, T_x, n)
  } else {
    Phi_direct <- matrix(0, n, ns3)
  }

  # Mixed (F_ww with W_xx + T_x) contributions, all 3 permutations,
  # all converted to Kronecker storage (slot 1 slowest, slot 3 fastest)
  # so they can be summed with the .contract_h3 direct term (which is
  # already in Kronecker convention -- verified numerically).
  #
  # The raw matrix M = t(W_xx) %*% H_e %*% T_x has shape n_s^2 x n_s.
  # vec(M) col-major: rows fast (n_s^2), cols slow (n_s).  Rows are the
  # W_xx column index = ghxx pair index (a-fast, b-slow).  Cols are
  # T_x column index = singleton state.
  #
  # Mapping for each pair labeling (which two slots go into W_xx, which
  # is the singleton), the array dimensions and conversion aperms are:
  #   pair12: (s1, s2) into W_xx, s3 singleton -> raw array dim (s1, s2, s3)
  #           -> Kronecker (s3, s2, s1) via aperm(_, c(3, 2, 1))
  #   pair13: (s1, s3) into W_xx, s2 singleton -> raw array dim (s1, s3, s2)
  #           -> aperm(_, c(2, 3, 1))
  #   pair23: (s2, s3) into W_xx, s1 singleton -> raw array dim (s2, s3, s1)
  #           -> aperm(_, c(2, 1, 3))
  # All three use the SAME numerical M matrix; only the index re-labeling
  # differs.
  Phi_pair12 <- matrix(0, n, ns3)
  Phi_pair13 <- matrix(0, n, ns3)
  Phi_pair23 <- matrix(0, n, ns3)
  H2_dense   <- H2

  for (e in seq_len(n)) {
    He <- H2_dense[e, , ]
    M  <- t(W_xx) %*% He %*% T_x        # n_s^2 x n_s, raw
    A  <- array(as.vector(M), dim = c(n_s, n_s, n_s))  # raw 3-tensor
    Phi_pair12[e, ] <- as.vector(aperm(A, c(3L, 2L, 1L)))
    Phi_pair13[e, ] <- as.vector(aperm(A, c(2L, 3L, 1L)))
    Phi_pair23[e, ] <- as.vector(aperm(A, c(2L, 1L, 3L)))
  }

  Phi_xxx <- Phi_direct + Phi_pair12 + Phi_pair13 + Phi_pair23

  # ----------------------------------------------------------------
  # Chain-rule cross-terms from y_{t+1} = g(h(x, 0, 0), 0, 0):
  #   d^3(g o h)/dx^3 contains ghxx [(hxx \otimes hx) summed over 3
  #   pair-singleton arrangements], in addition to the ghxxx and h_xxx
  #   pieces (the latter being the unknowns moved to the LHS).
  #
  # Raw `ghxx %*% (hxx %x% hx)` has shape n x n_s^3 with storage
  # (s_C fast, s_A mid, s_B slow) where (s_A, s_B) are the hxx pair
  # (in hxx column convention s_A fast, s_B slow within the pair) and
  # s_C is the hx col.
  #   pair12: (s_A, s_B, s_C) = (s1, s2, s3) -> array dim (s3, s1, s2)
  #           -> Kronecker (s3, s2, s1) via aperm(_, c(1, 3, 2))
  #   pair13: (s_A, s_B, s_C) = (s1, s3, s2) -> array dim (s2, s1, s3)
  #           -> aperm(_, c(3, 1, 2))
  #   pair23: (s_A, s_B, s_C) = (s2, s3, s1) -> array dim (s1, s2, s3)
  #           -> aperm(_, c(3, 2, 1))
  # ----------------------------------------------------------------
  if (verbose) cat("  Adding chain-rule cross-terms from y_{t+1}...\n")

  raw_chain    <- ghxx %*% (hxx %x% hx)      # n x n_s^3 raw
  cross_pair12 <- matrix(0, n, ns3)
  cross_pair13 <- matrix(0, n, ns3)
  cross_pair23 <- matrix(0, n, ns3)
  for (e in seq_len(n)) {
    A <- array(raw_chain[e, ], dim = c(n_s, n_s, n_s))
    cross_pair12[e, ] <- as.vector(aperm(A, c(1L, 3L, 2L)))
    cross_pair13[e, ] <- as.vector(aperm(A, c(3L, 1L, 2L)))
    cross_pair23[e, ] <- as.vector(aperm(A, c(3L, 2L, 1L)))
  }
  cross_known <- cross_pair12 + cross_pair13 + cross_pair23

  Phi_xxx <- Phi_xxx + fp %*% cross_known

  # ----------------------------------------------------------------
  # Solve Kronecker system for ghxxx:
  #   (I_{n_s^3} \otimes A_L + (hx' \otimes hx' \otimes hx') \otimes fp) vec(ghxxx)
  #     = -vec(Phi_xxx)
  # ----------------------------------------------------------------
  if (verbose) cat("  Solving Kronecker system for ghxxx...\n")

  # Shared generalized-Sylvester solver (same one orders 4/5 use):
  # eigenbasis fast path with a kappa(V)^k conditioning guard that falls back
  # to a real-Schur (Bartels-Stewart) dense solve. The Schur fallback never
  # inverts hx's eigenvector matrix, so it is robust to a DEFECTIVE / non-
  # diagonalizable hx (e.g. predetermined-variable lag chains with repeated
  # zero eigenvalues -- Born_Pfeifer_2014). Solves the identical equation as
  # the retired legacy order-3 solver: A_L * X + fp * X * hx^{otimes 3} = RHS
  # (verified bit-parity on rbc/medium_nonlinear_test). `solver_method` is
  # retained for back-compat but is now a no-op (the solver self-selects).
  #
  # OPT-IN sparse route: for high-dimensional-but-sparse state blocks (finite
  # HANK, n_s ~ 32 -> ns^3 ~ 3.3e4 RHS cols) the dense fallback of
  # .solve_kron_compact forms an 8.6 GB ns^3 x ns^3 Kronecker matrix -- a
  # documented >4h / 27 GB wall. .solve_kron_compact_sparse never materialises
  # it (complex-Schur Kronecker Bartels-Stewart). It is bit-parity on
  # well-conditioned models and residual-verified (dense fallback) otherwise.
  use_sparse <- if (is.null(sparse)) (n_s^3 >= sparse_threshold) else isTRUE(sparse)
  if (use_sparse) {
    if (verbose) cat(sprintf(
      "  ghxxx: sparse Kronecker solve (n_s=%d, ns^3=%d cols).\n", n_s, n_s^3))
    ghxxx <- .solve_kron_compact_sparse(A_L, fp, hx, 3L, -Phi_xxx,
                                        verbose = verbose)
  } else {
    ghxxx <- .solve_kron_compact(A_L, fp, hx, 3L, -Phi_xxx, verbose = verbose)
  }

  ghxxx <- .symmetrize_cube_cols(ghxxx, n_s)

  # ----------------------------------------------------------------
  # Solve for ghxxu, ghxuu, ghuuu (direct solves given ghxxx).
  #
  # The same Faà di Bruno chain rule with (T_x, T_x, T_u), (T_x, T_u, T_u),
  # (T_u, T_u, T_u) replacing (T_x, T_x, T_x).  The LHS coefficient
  # collapses to A_L for the "current" block (no Sylvester structure
  # since u_{t+1} is independent of x_{t-1}, so the (h^3) propagation
  # vanishes for u-mixed terms in the lead block of W).
  #
  # NOTE: ghxxu, ghxuu, ghuuu still receive known chain-rule contributions
  # from y_{t+1} = g(h(x,u), ...) involving the propagated hxx, hxu, huu.
  # ----------------------------------------------------------------
  if (verbose) cat("  Solving for ghxxu, ghxuu, ghuuu...\n")

  # Helper: contract H3 once with (Ta, Tb, Tc) summing the appropriate orbit.
  #
  # contract_h3 output layout: (Tc-dim FAST, Tb-dim MID, Ta-dim SLOW)
  #   i.e. col = (a-1)*nb*nc + (b-1)*nc + c   (a=Ta-col SLOW, c=Tc-col FAST)
  #
  # Phi_xxu/Phi_xuu use the INTERNAL convention (s1 FAST, u SLOW), so
  # H3_xxu and H3_xuu must be permuted from contract_h3's layout to match.
  # H3_uuu and H3_xxx already use Kronecker (first SLOW, last FAST) = same
  # as contract_h3, so they are unchanged.
  #
  # The required permutation is always aperm(A, c(3L, 2L, 1L)) on the
  # reshaped 3-tensor: (nc, nb, na) → (na, nb, nc) = INTERNAL order.
  # For n_u = 1 this permutation is the identity on values (due to the full
  # symmetry of H3, which makes (s1,s2) and (s2,s1) entries equal), so
  # single-shock results are unchanged.
  .h3_to_internal <- function(M, n_a, n_b, n_c) {
    if (is.null(M) || nrow(M) == 0L || ncol(M) == 0L) return(M)
    out <- matrix(0, nrow(M), n_a * n_b * n_c)
    for (e in seq_len(nrow(M))) {
      A        <- array(M[e, ], dim = c(n_c, n_b, n_a))   # [c, b, a] = contract_h3 layout
      out[e, ] <- as.vector(aperm(A, c(3L, 2L, 1L)))      # [a, b, c] = internal layout
    }
    out
  }
  H3_xxu <- if (!is.null(H3))
              .h3_to_internal(.contract_h3(H3, T_x, T_x, T_u, n), n_s, n_s, n_u)
            else matrix(0, n, n_s^2 * n_u)
  H3_xuu <- if (!is.null(H3))
              .h3_to_internal(.contract_h3(H3, T_x, T_u, T_u, n), n_s, n_u, n_u)
            else matrix(0, n, n_s * n_u^2)
  H3_uuu <- if (!is.null(H3)) .contract_h3(H3, T_u, T_u, T_u, n) else
            matrix(0, n, n_u^3)

  # Mixed F_ww contributions for x-x-u: 3 pair-singleton permutations,
  # all converted to Kronecker storage of ghxxu's cols (s1, s2, u) with
  # s1 slowest, u fastest -> array dim (u, s2, s1).
  #
  # pair (1,2) sing 3:  M = t(W_xx) He T_u  (n_s^2 x n_u),
  #   raw array dim (s1, s2, u) -> aperm c(3, 2, 1)
  # pair (1,3) sing 2:  M = t(W_xu) He T_x  (n_s*n_u x n_s),
  #   (s_state, k_exo) into W_xu = (s1, u), col = s2.
  #   raw array dim (s1, u, s2) -> aperm c(2, 3, 1)
  # pair (2,3) sing 1:  M = t(W_xu) He T_x  same matrix,
  #   (s_state, k_exo) into W_xu = (s2, u), col = s1.
  #   raw array dim (s2, u, s1) -> aperm c(2, 1, 3)
  Phi_xxu <- H3_xxu
  for (e in seq_len(n)) {
    He <- H2_dense[e, , ]
    M_xx <- t(W_xx) %*% He %*% T_u
    A12  <- array(as.vector(M_xx), dim = c(n_s, n_s, n_u))
    M_xu <- t(W_xu) %*% He %*% T_x
    A_xu  <- array(as.vector(M_xu), dim = c(n_s, n_u, n_s))
    # Explicit indexing: A12[s1,s2,u] = pair(s1,s2) in W_xx, u in T_u
    # A_xu[s_state, u_exo, s_sing] = pair(s_state,u_exo) in W_xu, s_sing in T_x
    for (s1 in seq_len(n_s)) for (s2 in seq_len(n_s)) for (u in seq_len(n_u)) {
      col <- s1 + (s2 - 1L)*n_s + (u - 1L)*n_s^2
      Phi_xxu[e, col] <- Phi_xxu[e, col] +
        A12[s1, s2, u] +          # pair12: W_xx(s1,s2) + T_u(u)
        A_xu[s1, u, s2] +         # pair13: W_xu(s1,u) + T_x(s2)
        A_xu[s2, u, s1]           # pair23: W_xu(s2,u) + T_x(s1)
    }
  }
  # Chain-rule cross-terms from y_{t+1}:
  #   pair12: ghxx %*% (hxx %x% hu) -> raw array dim (u, s1, s2) -> aperm c(1, 3, 2)
  #   pair13: ghxx %*% (hxu %x% hx) -> raw array dim (s2, s1, u) -> aperm c(3, 1, 2)
  #   pair23: ghxx %*% (hxu %x% hx) -> raw array dim (s1, s2, u) -> aperm c(3, 2, 1)
  raw_chain_xxu_a <- ghxx %*% (hxx %x% hu)
  raw_chain_xxu_b <- ghxx %*% (hxu %x% hx)
  cross_xxu_12 <- matrix(0, n, n_s^2 * n_u)
  cross_xxu_13 <- matrix(0, n, n_s^2 * n_u)
  cross_xxu_23 <- matrix(0, n, n_s^2 * n_u)
  for (e in seq_len(n)) {
    A_a <- array(raw_chain_xxu_a[e, ], dim = c(n_u, n_s, n_s))
    A_b <- array(raw_chain_xxu_b[e, ], dim = c(n_s, n_s, n_u))
    for (s1 in seq_len(n_s)) for (s2 in seq_len(n_s)) for (u in seq_len(n_u)) {
      col <- s1 + (s2 - 1L)*n_s + (u - 1L)*n_s^2
      # pair12: hxx(s1,s2) + hu(u)    -> raw col = u+(s1-1)*n_u+(s2-1)*n_s*n_u
      cross_xxu_12[e, col] <- A_a[u, s1, s2]
      # pair13: hxu(s1,u) + hx(s2)    -> raw col = s2+(s1-1)*n_s+(u-1)*n_s^2
      cross_xxu_13[e, col] <- A_b[s2, s1, u]
      # pair23: hxu(s2,u) + hx(s1)    -> raw col = s1+(s2-1)*n_s+(u-1)*n_s^2
      cross_xxu_23[e, col] <- A_b[s1, s2, u]
    }
  }
  Phi_xxu <- Phi_xxu + fp %*% (cross_xxu_12 + cross_xxu_13 + cross_xxu_23)

  # Also: ghxxx propagated via h: f_+ * ghxxx * (hx \otimes hx \otimes hu) is
  # the contribution of the chain-rule UNKNOWN ghxxx to the xxu equation.
  # Since ghxxx is now known, move to RHS.
  # Propagation: ghxxx * (hx\otimes hx \otimes hu) in Kronecker columns
  # (k, sB, sA). aperm(c(3,2,1)) converts to target (sA,sB,k).
  prop_xxu <- ghxxx %*% (hx %x% hx %x% hu)
  prop_xxu_p <- matrix(0, n, n_s^2 * n_u)
  for (e in seq_len(n)) {
    A <- array(prop_xxu[e, ], dim = c(n_u, n_s, n_s))
    prop_xxu_p[e, ] <- as.vector(aperm(A, c(3L, 2L, 1L)))
  }
  rhs_xxu <- -(Phi_xxu + fp %*% prop_xxu_p)
  ghxxu   <- tryCatch(solve(A_L, rhs_xxu),
                      error = function(e) qr.solve(A_L, rhs_xxu))
  ghxxu <- .symmetrize_xxu_cols(ghxxu, n_s, n_u)

  # x-u-u: cols indexed (s1, k1, k2) in Kronecker (s1 slow, k2 fast),
  # array dim (k2, k1, s1).
  # pair (1,2) sing 3: (s1, k1) into W_xu, k2 in T_u.
  #   M = t(W_xu) He T_u (n_s*n_u x n_u); raw array (s1, k1, k2) -> aperm (3,2,1)
  # pair (1,3) sing 2: (s1, k2) into W_xu, k1 in T_u (same M).
  #   raw array (s1, k2, k1) -> aperm (2,3,1)
  # pair (2,3) sing 1: (k1, k2) into W_uu, s1 in T_x.
  #   M = t(W_uu) He T_x (n_u^2 x n_s); raw array (k1, k2, s1) -> aperm (2,1,3)
  Phi_xuu <- H3_xuu
  for (e in seq_len(n)) {
    He <- H2_dense[e, , ]
    M_xu <- t(W_xu) %*% He %*% T_u
    A_xu  <- array(as.vector(M_xu), dim = c(n_s, n_u, n_u))
    M_uu <- t(W_uu) %*% He %*% T_x
    A_uu  <- array(as.vector(M_uu), dim = c(n_u, n_u, n_s))
    for (s1 in seq_len(n_s)) for (u1 in seq_len(n_u)) for (u2 in seq_len(n_u)) {
      col <- s1 + (u1 - 1L)*n_s + (u2 - 1L)*n_s*n_u
      Phi_xuu[e, col] <- Phi_xuu[e, col] +
        A_xu[s1, u1, u2] +         # pair12: W_xu(s1,u1) + T_u(u2)
        A_xu[s1, u2, u1] +         # pair13: W_xu(s1,u2) + T_u(u1)
        A_uu[u1, u2, s1]           # pair23: W_uu(u1,u2) + T_x(s1)
    }
  }
  # Chain-rule cross:
  #   pair12: ghxx %*% (hxu %x% hu) -> raw (k2, s1, k1) -> aperm (1, 3, 2)
  #   pair13: ghxx %*% (hxu %x% hu) -> raw (k1, s1, k2) -> aperm (3, 1, 2)
  #   pair23: ghxx %*% (huu %x% hx) -> raw (s1, k1, k2) -> aperm (3, 2, 1)
  raw_chain_xuu_a <- ghxx %*% (hxu %x% hu)
  raw_chain_xuu_b <- ghxx %*% (huu %x% hx)
  cross_xuu_12 <- matrix(0, n, n_s * n_u^2)
  cross_xuu_13 <- matrix(0, n, n_s * n_u^2)
  cross_xuu_23 <- matrix(0, n, n_s * n_u^2)
  for (e in seq_len(n)) {
    A_a <- array(raw_chain_xuu_a[e, ], dim = c(n_u, n_s, n_u))
    A_b <- array(raw_chain_xuu_b[e, ], dim = c(n_s, n_u, n_u))
    for (s1 in seq_len(n_s)) for (u1 in seq_len(n_u)) for (u2 in seq_len(n_u)) {
      col <- s1 + (u1 - 1L)*n_s + (u2 - 1L)*n_s*n_u
      # pair12: hxu(s1,u1) + hu(u2)  -> raw col = u2+(s1-1)*n_u+(u1-1)*n_s*n_u
      cross_xuu_12[e, col] <- A_a[u2, s1, u1]
      # pair13: hxu(s1,u2) + hu(u1)  -> raw col = u1+(s1-1)*n_u+(u2-1)*n_s*n_u
      cross_xuu_13[e, col] <- A_a[u1, s1, u2]
      # pair23: huu(u1,u2) + hx(s1)  -> raw col = s1+(u1-1)*n_s+(u2-1)*n_s*n_u
      cross_xuu_23[e, col] <- A_b[s1, u1, u2]
    }
  }
  Phi_xuu <- Phi_xuu + fp %*% (cross_xuu_12 + cross_xuu_13 + cross_xuu_23)

  # Propagation: ghxxx * (hx\otimes hu \otimes hu) in Kronecker columns
  # (k2, k1, s). aperm(c(3,2,1)) converts to target (s,k1,k2).
  prop_xuu <- ghxxx %*% (hx %x% hu %x% hu)
  prop_xuu_p <- matrix(0, n, n_s * n_u^2)
  for (e in seq_len(n)) {
    A <- array(prop_xuu[e, ], dim = c(n_u, n_u, n_s))
    prop_xuu_p[e, ] <- as.vector(aperm(A, c(3L, 2L, 1L)))
  }
  rhs_xuu <- -(Phi_xuu + fp %*% prop_xuu_p)
  ghxuu   <- tryCatch(solve(A_L, rhs_xuu),
                      error = function(e) qr.solve(A_L, rhs_xuu))
  ghxuu <- .symmetrize_xuu_cols(ghxuu, n_s, n_u)

  # u-u-u: 3 pair-singleton arrangements with (W_uu, T_u).
  # Phi_uuu cols: (k1 slow, k2 mid, k3 fast) = (k1-1)*n_u^2+(k2-1)*n_u+k3.
  # M = t(W_uu) He T_u (n_u^2 x n_u); W_uu cols in ghuu convention
  # (exo_fast, exo_slow) => A[exo_fast=k1, exo_slow=k2, exo_sing=k3].
  # The 3 pair-singleton arrangements contribute A[k1,k2,k3]+A[k1,k3,k2]+A[k2,k3,k1]
  # at target slot (k1,k2,k3). Each is a different aperm of the same A:
  #   pair12 (fast=k1,slow=k2) sing k3: A[k1,k2,k3] → aperm(A, c(3,2,1))
  #   pair13 (fast=k1,slow=k3) sing k2: A[k1,k3,k2] → aperm(A, c(2,3,1))
  #   pair23 (fast=k2,slow=k3) sing k1: A[k2,k3,k1] → aperm(A, c(2,1,3))
  Phi_uuu <- H3_uuu
  for (e in seq_len(n)) {
    He <- H2_dense[e, , ]
    M_uuu <- t(W_uu) %*% He %*% T_u
    A     <- array(as.vector(M_uuu), dim = c(n_u, n_u, n_u))
    Phi_uuu[e, ] <- Phi_uuu[e, ] +
      as.vector(aperm(A, c(3L, 2L, 1L))) +   # pair12
      as.vector(aperm(A, c(2L, 3L, 1L))) +   # pair13
      as.vector(aperm(A, c(2L, 1L, 3L)))     # pair23
  }
  # Chain-rule cross: ghxx %*% (huu %x% hu).
  # (huu %x% hu) cols: (huu_col SLOW, hu_col FAST) where huu cols = ghuu convention
  # (exo_fast, exo_slow) => raw col encodes (exo_slow SLOWEST, exo_fast MID, exo_sing FAST).
  # A[exo_sing=k3, exo_fast=k1, exo_slow=k2]; 3 arrangements at (k1,k2,k3):
  #   pair12 (fast=k1,slow=k2) sing k3: A[k3,k1,k2] → aperm(A, c(1,3,2))
  #   pair13 (fast=k1,slow=k3) sing k2: A[k2,k1,k3] → aperm(A, c(3,1,2))
  #   pair23 (fast=k2,slow=k3) sing k1: A[k1,k2,k3] → aperm(A, c(3,2,1))
  raw_chain_uuu <- ghxx %*% (huu %x% hu)
  cross_uuu     <- matrix(0, n, n_u^3)
  for (e in seq_len(n)) {
    A <- array(raw_chain_uuu[e, ], dim = c(n_u, n_u, n_u))
    cross_uuu[e, ] <-
      as.vector(aperm(A, c(1L, 3L, 2L))) +   # pair12
      as.vector(aperm(A, c(3L, 1L, 2L))) +   # pair13
      as.vector(aperm(A, c(3L, 2L, 1L)))     # pair23
  }
  Phi_uuu <- Phi_uuu + fp %*% cross_uuu

  # Propagation: ghxxx * (hu\otimes hu \otimes hu) in Kronecker columns
  # (k3, k2, k1). aperm(c(3,2,1)) converts to target (k1,k2,k3).
  prop_uuu <- ghxxx %*% (hu %x% hu %x% hu)
  prop_uuu_p <- matrix(0, n, n_u^3)
  for (e in seq_len(n)) {
    A <- array(prop_uuu[e, ], dim = c(n_u, n_u, n_u))
    prop_uuu_p[e, ] <- as.vector(aperm(A, c(3L, 2L, 1L)))
  }
  rhs_uuu <- -(Phi_uuu + fp %*% prop_uuu_p)
  ghuuu   <- tryCatch(solve(A_L, rhs_uuu),
                      error = function(e) qr.solve(A_L, rhs_uuu))
  ghuuu <- .symmetrize_uuu_cols(ghuuu, n_u)

  # ----------------------------------------------------------------
  # Permute ghxxu and ghxuu from internal convention (state index FAST,
  # exo index SLOW) to standard Kronecker convention (state SLOW, exo FAST)
  # so that  ghxxu %*% (x %x% x %x% u) etc. produces correct Taylor
  # coefficients. The convention bug was invisible for n_u == 1 (e.g. RBC)
  # because the permutation is identity there.
  # ----------------------------------------------------------------
  if (n_s > 0L && n_u > 0L) {
    perm_xxu <- integer(n_s * n_s * n_u)
    for (s1 in seq_len(n_s)) for (s2 in seq_len(n_s)) for (u in seq_len(n_u)) {
      int_col <- s1 + (s2-1L)*n_s + (u-1L)*n_s^2
      std_col <- (s1-1L)*n_s*n_u + (s2-1L)*n_u + u
      perm_xxu[std_col] <- int_col
    }
    ghxxu <- ghxxu[, perm_xxu, drop = FALSE]

    perm_xuu <- integer(n_s * n_u * n_u)
    for (s1 in seq_len(n_s)) for (u1 in seq_len(n_u)) for (u2 in seq_len(n_u)) {
      int_col <- s1 + (u1-1L)*n_s + (u2-1L)*n_s*n_u
      std_col <- (s1-1L)*n_u^2 + (u1-1L)*n_u + u2
      perm_xuu[std_col] <- int_col
    }
    ghxuu <- ghxuu[, perm_xuu, drop = FALSE]
  }

  # ----------------------------------------------------------------
  # Permute dr2$ghxu and dr2$ghxu's downstream copies (inherited into dr3)
  # to standard Kron convention so users get consistent storage across
  # ghxu, ghxxu, ghxuu, ghuuu.  dr2's own ghxu is left untouched (its
  # internal convention is required by other internal solvers).
  # ----------------------------------------------------------------
  ghxu_out <- ghxu
  if (n_s > 0L && n_u > 0L) {
    perm_xu <- integer(n_s * n_u)
    for (s in seq_len(n_s)) for (u in seq_len(n_u)) {
      int_col <- (u-1L)*n_s + s                # internal: state FAST, exo SLOW
      std_col <- (s-1L)*n_u + u                # standard: state SLOW, exo FAST
      perm_xu[std_col] <- int_col
    }
    ghxu_out <- ghxu[, perm_xu, drop = FALSE]
  }

  # ----------------------------------------------------------------
  # Name and assemble the output
  # ----------------------------------------------------------------
  state_vars <- endo_names[state_idx]

  rownames(ghxxx) <- endo_names
  rownames(ghxxu) <- endo_names
  rownames(ghxuu) <- endo_names
  rownames(ghuuu) <- endo_names

  # Column names: outer product of (state, state, state) etc.
  triple_names <- function(a, b, c) {
    out <- character(length(a) * length(b) * length(c))
    idx <- 1L
    for (k in seq_along(c)) for (j in seq_along(b)) for (i in seq_along(a)) {
      out[idx] <- paste(a[i], b[j], c[k], sep = "__x__")
      idx <- idx + 1L
    }
    out
  }
  colnames(ghxxx) <- triple_names(state_vars, state_vars, state_vars)
  colnames(ghxxu) <- triple_names(state_vars, state_vars, exo_names)
  colnames(ghxuu) <- triple_names(state_vars, exo_names, exo_names)
  colnames(ghuuu) <- triple_names(exo_names,  exo_names,  exo_names)

  dr3       <- unclass(dr2)
  dr3$ghxu  <- ghxu_out   # permuted to standard Kron convention for user
  dr3$ghxxx <- ghxxx
  dr3$ghxxu <- ghxxu
  dr3$ghxuu <- ghxuu
  dr3$ghuuu <- ghuuu
  dr3$order <- 3L                                # override dr2$order = 2L
  dr3$third_order_method <- "deterministic_binning_2013"
  dr3$sigma_correction   <- "not_implemented"
  class(dr3) <- c("DecisionRules3", "DecisionRules2", "DecisionRules")

  if (verbose) {
    cat("Third-order solution complete.\n")
    cat(sprintf("  ghxxx: %d x %d  max|.| = %.3g\n",
                nrow(ghxxx), ncol(ghxxx), max(abs(ghxxx))))
    cat(sprintf("  ghxxu: %d x %d  max|.| = %.3g\n",
                nrow(ghxxu), ncol(ghxxu), max(abs(ghxxu))))
    cat(sprintf("  ghxuu: %d x %d  max|.| = %.3g\n",
                nrow(ghxuu), ncol(ghxuu), max(abs(ghxuu))))
    cat(sprintf("  ghuuu: %d x %d  max|.| = %.3g\n",
                nrow(ghuuu), ncol(ghuuu), max(abs(ghuuu))))
  }

  dr3
}


#' Trivial third-order solution (no state variables).
#' @noRd
.trivial_dr3 <- function(dr2) {
  n   <- length(dr2$endo_names)
  n_u <- length(dr2$exo_names)
  dr3       <- unclass(dr2)
  dr3$ghxxx <- matrix(0, n, 0)
  dr3$ghxxu <- matrix(0, n, 0)
  dr3$ghxuu <- matrix(0, n, 0)
  dr3$ghuuu <- matrix(0, n, n_u^3)
  dr3$order <- 3L
  dr3$third_order_method <- "trivial_no_states"
  dr3$sigma_correction   <- "not_implemented"
  class(dr3) <- c("DecisionRules3", "DecisionRules2", "DecisionRules")
  dr3
}


.linear_dr3 <- function(dr2) {
  n          <- length(dr2$endo_names)
  n_s        <- length(dr2$state_idx)
  n_u        <- length(dr2$exo_names)
  endo_names <- dr2$endo_names
  state_vars <- endo_names[dr2$state_idx]
  exo_names  <- dr2$exo_names

  make_zero <- function(nc, nms) {
    m <- matrix(0, n, nc)
    rownames(m) <- endo_names
    if (!is.null(nms)) colnames(m) <- nms
    m
  }

  xxx_nms <- as.vector(outer(outer(state_vars, state_vars, paste, sep = "__x__"),
                              state_vars, paste, sep = "__x__"))
  xxu_nms <- as.vector(outer(outer(state_vars, state_vars, paste, sep = "__x__"),
                              exo_names,  paste, sep = "__x__"))
  xuu_nms <- as.vector(outer(state_vars,
                              as.vector(outer(exo_names, exo_names, paste, sep = "__x__")),
                              paste, sep = "__x__"))
  uuu_nms <- as.vector(outer(outer(exo_names, exo_names, paste, sep = "__x__"),
                              exo_names,  paste, sep = "__x__"))

  dr3       <- unclass(dr2)
  dr3$ghxxx <- make_zero(n_s^3,      xxx_nms)
  dr3$ghxxu <- make_zero(n_s^2*n_u,  xxu_nms)
  dr3$ghxuu <- make_zero(n_s*n_u^2,  xuu_nms)
  dr3$ghuuu <- make_zero(n_u^3,      uuu_nms)
  dr3$order <- 3L
  dr3$third_order_method <- "linear_shortcircuit"
  dr3$sigma_correction   <- "not_implemented"
  class(dr3) <- c("DecisionRules3", "DecisionRules2", "DecisionRules")
  dr3
}


# =====================================================================
# S3 print method
# =====================================================================

#' @export
print.DecisionRules3 <- function(x, ...) {
  cat("Third-order Decision Rules (DecisionRules3)\n")
  cat("  Endogenous variables:", length(x$endo_names), "\n")
  cat("  State variables:     ", x$n_state, "\n")
  cat("  Shocks:              ", x$n_exo, "\n")
  cat("  Perturbation order:  3 (deterministic)\n")
  cat("  Hessian method:     ", x$hessian_method %||% "unknown", "\n")
  cat("  sigma-correction:   ", x$sigma_correction, "\n")
  cat("\nFirst-order:\n")
  cat(sprintf("  ghx: %d x %d   ghu: %d x %d\n",
              nrow(x$ghx), ncol(x$ghx), nrow(x$ghu), ncol(x$ghu)))
  cat("\nSecond-order:\n")
  cat(sprintf("  ghxx: %d x %d max=%.3g   ghxu: %d x %d max=%.3g\n",
              nrow(x$ghxx), ncol(x$ghxx), max(abs(x$ghxx)),
              nrow(x$ghxu), ncol(x$ghxu), max(abs(x$ghxu))))
  cat(sprintf("  ghuu: %d x %d max=%.3g   ghss: len %d max=%.3g\n",
              nrow(x$ghuu), ncol(x$ghuu), max(abs(x$ghuu)),
              length(x$ghss), max(abs(x$ghss))))
  cat("\nThird-order:\n")
  if (!is.null(x$ghxxx) && prod(dim(x$ghxxx)) > 0) {
    cat(sprintf("  ghxxx: %d x %d max=%.3g\n",
                nrow(x$ghxxx), ncol(x$ghxxx), max(abs(x$ghxxx))))
    cat(sprintf("  ghxxu: %d x %d max=%.3g\n",
                nrow(x$ghxxu), ncol(x$ghxxu), max(abs(x$ghxxu))))
    cat(sprintf("  ghxuu: %d x %d max=%.3g\n",
                nrow(x$ghxuu), ncol(x$ghxuu), max(abs(x$ghxuu))))
    cat(sprintf("  ghuuu: %d x %d max=%.3g\n",
                nrow(x$ghuuu), ncol(x$ghuuu), max(abs(x$ghuuu))))
  } else {
    cat("  (no state variables; zero)\n")
  }
  invisible(x)
}


# =====================================================================
# Pruned third-order IRF
# =====================================================================

#' Compute impulse response functions at third order (pruned state space)
#'
#' Implements the Andreasen, Fernandez-Villaverde & Rubio-Ramirez (2018)
#' pruning scheme truncated at third order.  IRFs are deviations from the
#' stochastic steady state: the constant \code{ghss}/\code{ghs3} terms are
#' absorbed into the baseline and cancel, so only the shock-driven deviations
#' are returned.
#'
#' The pruned state tracks three additive layers \eqn{x^f, x^s, x^{rd}}:
#' \itemize{
#'   \item \eqn{x^f} -- first-order (linear) state response
#'   \item \eqn{x^s} -- second-order correction driven by \eqn{(x^f)^2}
#'   \item \eqn{x^{rd}} -- third-order correction driven by
#'     \eqn{(x^f)^3} and \eqn{x^f \otimes x^s} cross terms
#' }
#'
#' Output at each horizon: \eqn{y = y^{(1)} + y^{(2)} + y^{(3)}}.
#'
#' @param dr3       \code{DecisionRules3} object from
#'   \code{\link{solve_perturbation}(order = 3)}.
#' @param model     \code{dynhr_mod} from \code{\link{parse_mod}} (for shock
#'   standard deviations).
#' @param n_periods Number of IRF horizons (default 40).
#' @param shock_size Shock size in standard-deviation units (default 1).
#' @param params    Named numeric parameter vector.  \code{NULL} uses
#'   \code{model$param_values}.
#' @return An \code{IRFCollection}: a named list of \code{n_periods x n_endo}
#'   matrices (one per shock), with \code{attr(., "order") = 3L}.
#' @references Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez,
#'   J. F. (2018). The pruned state-space system for non-linear DSGE models.
#'   \emph{Review of Economic Studies}, 85(1), 1--49.
#' @export
compute_irfs_order3 <- function(dr3, model, n_periods = 40L,
                                 shock_size = 1, params = NULL) {
  if (!inherits(dr3, "DecisionRules3")) {
    stop("dr3 must be a DecisionRules3 object from solve_perturbation(order = 3).")
  }

  ghx   <- dr3$ghx;   ghu   <- dr3$ghu
  ghxx  <- dr3$ghxx;  ghuu  <- dr3$ghuu
  ghxxx <- dr3$ghxxx; ghuuu <- dr3$ghuuu

  # Sigma^2 cross-term matrices (present when solve_perturbation(order=3) was
  # used; NULL when only solve_perturbation_order3() was called directly).
  # ghuss: n_endo x n_exo  — shock x sigma^2 cross term
  # ghxss: n_endo x n_state — state x sigma^2 cross term
  has_uss <- !is.null(dr3$ghuss)
  has_xss <- !is.null(dr3$ghxss)

  endo      <- dr3$endo_names
  exo       <- dr3$exo_names
  state_idx <- dr3$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)
  n_s       <- length(state_idx)

  if (is.null(params)) params <- model$param_values
  shock_stderr <- .get_shock_stderr(model, exo, params)

  # State-row submatrices
  hx   <- ghx  [state_idx, , drop = FALSE]   # n_s x n_s
  hu   <- ghu  [state_idx, , drop = FALSE]   # n_s x n_u
  hxx  <- ghxx [state_idx, , drop = FALSE]   # n_s x n_s^2
  huu  <- ghuu [state_idx, , drop = FALSE]   # n_s x n_u^2
  huuu <- ghuuu[state_idx, , drop = FALSE]   # n_s x n_u^3
  hxxx <- ghxxx[state_idx, , drop = FALSE]   # n_s x n_s^3
  huss <- if (has_uss) dr3$ghuss[state_idx, , drop = FALSE] else NULL  # n_s x n_u
  hxss <- if (has_xss) dr3$ghxss[state_idx, , drop = FALSE] else NULL  # n_s x n_s

  irfs <- list()
  for (k in seq_along(exo)) {
    shock_name <- exo[k]
    irf_mat    <- matrix(0, n_periods, n_endo)
    colnames(irf_mat) <- endo
    rownames(irf_mat) <- paste0("t", seq_len(n_periods))

    eps    <- numeric(n_exo)
    eps[k] <- shock_stderr[shock_name] * shock_size

    # --- Period 1: impact (all previous states are zero) -------------------
    # ghss and ghs3 (constant sigma^2 terms) are absorbed in the baseline and
    # cancel in deviation IRFs.  But ghuss (shock x sigma^2) and ghxss
    # (state x sigma^2) do NOT cancel: ghuss*eps is nonzero at t=1 in the
    # shocked path but zero in the baseline.  See Andreasen, Fernandez-
    # Villaverde & Rubio-Ramirez (2018, RES) pruning Appendix.

    eps2 <- eps %x% eps           # length n_u^2
    eps3 <- eps %x% eps %x% eps   # length n_u^3

    x1 <- as.numeric(hu   %*% eps)
    x2 <- as.numeric(0.5 * huu  %*% eps2)
    # x3 at t=1: (1/6) huuu eps^3 + 0.5 huss eps  (x1_prev=0 so hxss term=0)
    x3 <- as.numeric((1/6) * huuu %*% eps3)
    if (has_uss) x3 <- x3 + 0.5 * as.numeric(huss %*% eps)

    y1 <- as.numeric(ghu  %*% eps)
    y2 <- as.numeric(0.5 * ghuu  %*% eps2)
    # y3 at t=1: (1/6) ghuuu eps^3 + 0.5 ghuss eps
    y3 <- as.numeric((1/6) * ghuuu %*% eps3)
    if (has_uss) y3 <- y3 + 0.5 * as.numeric(dr3$ghuss %*% eps)

    irf_mat[1, ] <- y1 + y2 + y3

    # --- Periods 2..n_periods: no further shocks ---------------------------
    if (n_periods >= 2L) for (t in 2:n_periods) {
      x1p <- x1; x2p <- x2; x3p <- x3

      x1 <- as.numeric(hx %*% x1p)
      x2 <- as.numeric(hx %*% x2p + 0.5 * hxx %*% (x1p %x% x1p))
      # (1/2)(hxx(x1⊗x2 + x2⊗x1)) = hxx(x1⊗x2) by Schwarz symmetry of hxx
      # x3 propagation: add 0.5 hxss x1p (state x sigma^2 cross term)
      x3 <- as.numeric(
        hx  %*% x3p +
        hxx %*% (x1p %x% x2p) +
        (1/6) * hxxx %*% (x1p %x% x1p %x% x1p)
      )
      if (has_xss) x3 <- x3 + 0.5 * as.numeric(hxss %*% x1p)

      y1 <- as.numeric(ghx  %*% x1p)
      y2 <- as.numeric(ghx  %*% x2p + 0.5 * ghxx %*% (x1p %x% x1p))
      # y3 propagation: add 0.5 ghxss x1p (state x sigma^2 cross term)
      y3 <- as.numeric(
        ghx  %*% x3p +
        ghxx %*% (x1p %x% x2p) +
        (1/6) * ghxxx %*% (x1p %x% x1p %x% x1p)
      )
      if (has_xss) y3 <- y3 + 0.5 * as.numeric(dr3$ghxss %*% x1p)

      irf_mat[t, ] <- y1 + y2 + y3
    }

    irfs[[shock_name]] <- irf_mat
  }

  class(irfs) <- "IRFCollection"
  attr(irfs, "n_periods")  <- n_periods
  attr(irfs, "endo_names") <- endo
  attr(irfs, "exo_names")  <- exo
  attr(irfs, "order")      <- 3L
  irfs
}


# =====================================================================
# Pruned third-order stochastic simulator
# =====================================================================

#' Pruned third-order stochastic simulator for DecisionRules3 objects
#'
#' Simulates the model using the pruned third-order perturbation approximation
#' of Andreasen, Fernandez-Villaverde and Rubio-Ramirez (2018).  The pruned
#' state decomposes into three additive layers: x^f (first order), x^s (second
#' order), and x^rd (third order).  The total deviation from steady state is
#' x^f + x^s + x^rd, and observable output is reconstructed from y^(1) + y^(2)
#' + y^(3).
#'
#' @param dr3 A \code{DecisionRules3} object from
#'   \code{solve_perturbation(order = 3)}.
#' @param n_periods Integer.  Number of periods to return (after burn-in).
#' @param shocks Optional \code{(n_periods + burn_in) x n_exo} matrix of
#'   pre-drawn shocks (in model units, i.e. already scaled by stderr).  If
#'   \code{NULL}, shocks are drawn from \eqn{N(0, \sigma_e^2)}.
#' @param model Optional model object used to retrieve shock standard
#'   deviations when \code{shocks = NULL}.
#' @param burn_in Integer.  Number of initial periods to discard.
#' @param pruning Logical.  If \code{TRUE} (default), run the full pruned
#'   three-layer recursion.  If \code{FALSE}, set x^s = x^rd = 0 every
#'   period, reducing to a first-order simulation.
#' @param init_state Optional named numeric vector of initial state deviations
#'   loaded into the first-order pruned component x^f; pair with
#'   \code{burn_in = 0}.  \code{NULL} starts at the steady state.
#' @return Matrix \code{n_periods x n_endo} of deviations from steady state,
#'   with attribute \code{"levels"} giving levels.
#' @references Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez,
#'   J. F. (2018). The pruned state-space system for non-linear DSGE models:
#'   Theory and empirical applications. \emph{Review of Economic Studies},
#'   85(1), 1--49.
#' @export
simulate_model_order3 <- function(dr3, n_periods = 200L, shocks = NULL,
                                   model = NULL, burn_in = 100L,
                                   pruning = TRUE, init_state = NULL) {
  if (!inherits(dr3, "DecisionRules3")) {
    stop("dr3 must be a DecisionRules3 object.")
  }

  # ---- Extract DR2 matrices (needed for all three layers) ----------------
  ghx  <- dr3$ghx;  ghu  <- dr3$ghu
  ghxx <- dr3$ghxx; ghxu <- dr3$ghxu
  ghuu <- dr3$ghuu; ghss <- dr3$ghss

  # ---- Extract DR3 matrices (third-order layer) --------------------------
  ghxxx <- dr3$ghxxx; ghxxu <- dr3$ghxxu
  ghxuu <- dr3$ghxuu; ghuuu <- dr3$ghuuu

  # ---- Optional sigma-correction fields ----------------------------------
  # ghxss / ghuss are present when solve_sigma_cross() has been called
  # (always the case for dr3 from solve_perturbation()).
  # ghs3 is present only when sigma3 != NULL.
  has_xss <- !is.null(dr3$ghxss)
  has_uss <- !is.null(dr3$ghuss)
  has_s3  <- !is.null(dr3$ghs3)

  endo      <- dr3$endo_names
  exo       <- dr3$exo_names
  state_idx <- dr3$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)
  n_s       <- length(state_idx)
  params    <- if (!is.null(model)) model$param_values else NULL

  # ---- State-row submatrices (n_s x *) -----------------------------------
  hx   <- ghx [state_idx, , drop = FALSE]   # n_s x n_s
  hu   <- ghu [state_idx, , drop = FALSE]   # n_s x n_u
  hxx  <- ghxx[state_idx, , drop = FALSE]   # n_s x n_s^2
  hxu  <- ghxu[state_idx, , drop = FALSE]   # n_s x n_s*n_u
  huu  <- ghuu[state_idx, , drop = FALSE]   # n_s x n_u^2
  hss  <- ghss[state_idx]                   # n_s vector

  hxxx <- ghxxx[state_idx, , drop = FALSE]  # n_s x n_s^3
  hxxu <- ghxxu[state_idx, , drop = FALSE]  # n_s x n_s^2*n_u
  hxuu <- ghxuu[state_idx, , drop = FALSE]  # n_s x n_s*n_u^2
  huuu <- ghuuu[state_idx, , drop = FALSE]  # n_s x n_u^3

  hxss <- if (has_xss) dr3$ghxss[state_idx, , drop = FALSE] else NULL  # n_s x n_s
  huss <- if (has_uss) dr3$ghuss[state_idx, , drop = FALSE] else NULL  # n_s x n_u
  hs3  <- if (has_s3)  dr3$ghs3[state_idx]                  else NULL  # n_s vector

  # ---- Shocks ------------------------------------------------------------
  total_periods <- n_periods + burn_in

  shock_stderr <- .get_shock_stderr(model, exo, params)
  if (is.null(shocks)) {
    shocks <- matrix(rnorm(total_periods * n_exo), ncol = n_exo)
    for (k in seq_along(exo)) shocks[, k] <- shocks[, k] * shock_stderr[exo[k]]
  }

  sim <- matrix(0, total_periods, n_endo)
  colnames(sim) <- endo

  # ---- Pruned state components -------------------------------------------
  x1 <- numeric(n_s)   # first-order state component
  x2 <- numeric(n_s)   # second-order state component
  x3 <- numeric(n_s)   # third-order state component (x^rd)

  if (!is.null(init_state)) {
    idx <- match(endo[state_idx], names(init_state))
    ok  <- !is.na(idx)
    if (any(ok)) x1[ok] <- as.numeric(init_state[idx[ok]])
  }

  # ---- Simulation loop ---------------------------------------------------
  # Kronecker conventions (inherited from solve-perturbation-order2.R):
  #   ghxu cols are outer(state_vars, exo_names): state FAST, exo SLOW
  #   => matching kron vector is (e %x% x): exo SLOW outer, state FAST inner
  #
  # ghxxu cols are triple_names(state_vars, state_vars, exo_names):
  #   state1 FASTEST, state2 MIDDLE, exo SLOWEST
  #   => matching kron vector is (e %x% x1 %x% x1): exo SLOW, state FAST
  #
  # ghxuu cols are triple_names(state_vars, exo_names, exo_names):
  #   state FASTEST, exo1 MIDDLE, exo2 SLOWEST
  #   => matching kron vector is (e %x% e %x% x): exo2 SLOW, state FAST

  for (t in seq_len(total_periods)) {
    e       <- shocks[t, ]
    x1_prev <- x1
    x2_prev <- x2
    x3_prev <- x3

    # (A) First-order component
    x1 <- as.numeric(hx %*% x1_prev + hu %*% e)

    if (pruning) {
      # (B) Second-order component
      x2 <- as.numeric(
        hx  %*% x2_prev +
        0.5 * hxx %*% (x1_prev %x% x1_prev) +
        hxu %*% (e %x% x1_prev) +
        0.5 * huu %*% (e %x% e) +
        0.5 * hss
      )

      # (C) Third-order component
      # hxx*(x1⊗x2): factor 1 (not 1/2) — both (x1⊗x2) and (x2⊗x1) terms
      # collapse to this single term by symmetry of hxx (Schwarz symmetry).
      x3_new <- as.numeric(
        hx  %*% x3_prev +
        hxx %*% (x1_prev %x% x2_prev) +
        hxu %*% (e %x% x2_prev) +
        0.5 * hxxu %*% (e %x% x1_prev %x% x1_prev) +
        0.5 * hxuu %*% (e %x% e %x% x1_prev) +
        (1/6) * hxxx %*% (x1_prev %x% x1_prev %x% x1_prev) +
        (1/6) * huuu %*% (e %x% e %x% e)
      )
      if (!is.null(hxss)) x3_new <- x3_new + 0.5 * as.numeric(hxss %*% x1_prev)
      if (!is.null(huss)) x3_new <- x3_new + 0.5 * as.numeric(huss %*% e)
      if (!is.null(hs3))  x3_new <- x3_new + (1/6) * hs3
      x3 <- x3_new
    } else {
      x2 <- numeric(n_s)
      x3 <- numeric(n_s)
    }

    # ---- Observable reconstruction ----------------------------------------
    # y^(1): linear terms using x^f_{t-1} and e_t
    y1 <- as.numeric(ghx %*% x1_prev + ghu %*% e)

    # y^(2): second-order correction
    y2 <- as.numeric(
      ghx  %*% x2_prev +
      0.5 * ghxx %*% (x1_prev %x% x1_prev) +
      ghxu %*% (e %x% x1_prev) +
      0.5 * ghuu %*% (e %x% e) +
      0.5 * ghss
    )

    # y^(3): third-order correction
    y3_val <- as.numeric(
      ghx  %*% x3_prev +
      ghxx %*% (x1_prev %x% x2_prev) +
      ghxu %*% (e %x% x2_prev) +
      0.5 * ghxxu %*% (e %x% x1_prev %x% x1_prev) +
      0.5 * ghxuu %*% (e %x% e %x% x1_prev) +
      (1/6) * ghxxx %*% (x1_prev %x% x1_prev %x% x1_prev) +
      (1/6) * ghuuu %*% (e %x% e %x% e)
    )
    if (has_xss) y3_val <- y3_val + 0.5 * as.numeric(dr3$ghxss %*% x1_prev)
    if (has_uss) y3_val <- y3_val + 0.5 * as.numeric(dr3$ghuss %*% e)
    if (has_s3)  y3_val <- y3_val + (1/6) * dr3$ghs3

    sim[t, ] <- y1 + y2 + y3_val
  }

  sim <- sim[(burn_in + 1L):total_periods, , drop = FALSE]

  # Add SS levels
  sim_levels <- sim
  for (j in seq_along(endo)) {
    ss_val <- dr3$ys[endo[j]]
    if (!is.na(ss_val)) sim_levels[, j] <- sim[, j] + ss_val
  }
  attr(sim, "levels") <- sim_levels
  sim
}
