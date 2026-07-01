## R/solve-perturbation-order2.R
## --------------------------------------------------------------------------
## Second-order perturbation solver for DSGE models.
##
## Implements the Schmitt-Grohé & Uribe (2004) algorithm using the
## Binning (2013) matrix chain rule formulation:
##   1. Obtain model Hessian H[eq, col1, col2] at steady state
##        PRIMARY: symbolic (from compiled$dynamic$hessian2_fn + hess2_triplets)
##        FALLBACK: numerical finite-differences of the first-order Jacobian
##   2. Build "transfer matrices" mapping (state, shock) perturbations to
##      the compound variable vector (y_{t-1}, y_t, y_{t+1}, u_t)
##   3. Contract Hessian with transfer matrices (Faà di Bruno chain rule)
##      to get forcing terms Phi_xx, Phi_xu, Phi_uu
##   4. Solve Kronecker linear system for ghxx (Sylvester equation)
##   5. Direct solves for ghxu, ghuu, ghss (uncertainty correction)
##
## References:
##   Schmitt-Grohé, S. & Uribe, M. (2004). Solving Dynamic General
##     Equilibrium Models Using a Second-Order Approximation to the
##     Policy Function. J. Economic Dynamics and Control 28(4): 755-775.
##   Binning, A. (2013). Solving Second and Third-Order Approximations to
##     DSGE Models: A Recursive Sylvester Equation Solution.
##     Norges Bank WP 2013/18.
##   Andreasen, M.M., Fernandez-Villaverde, J. & Rubio-Ramirez, J.F.
##     (2018). The Pruned State-Space System for Non-Linear DSGE Models.
##     Review of Economic Studies 85(1): 1-49.
## --------------------------------------------------------------------------


# =====================================================================
# Internal helpers
# =====================================================================

#' Build the compound variable vector dy at the steady state
#'
#' Mirrors the logic in extract_system_matrices() but returns only the
#' named dy vector (not the full system matrices), suitable for
#' numerical differentiation.
#'
#' @param compiled dynhr_compiled
#' @param ss       Named numeric steady state
#' @return Named numeric vector for jacobian_fn()
#' @noRd
.build_dy_ss_o2 <- function(compiled, ss) {
  dyn  <- compiled$dynamic
  exo  <- compiled$model$varexo_names
  dy   <- numeric(nrow(dyn$dyn_col_map))
  keys <- character(nrow(dyn$dyn_col_map))

  for (k in seq_len(nrow(dyn$dyn_col_map))) {
    nm  <- dyn$dyn_col_map$name[k]
    ll  <- dyn$dyn_col_map$lead_lag[k]
    sfx <- if (ll == 0L) "__0" else if (ll > 0L) paste0("__p", ll) else paste0("__m", abs(ll))
    keys[k] <- paste0(nm, sfx)
    dy[k]   <- if (nm %in% names(ss)) ss[[nm]] else 0
  }
  names(dy) <- keys
  dy
}


#' Build equation-to-declaration mapping from model AST
#'
#' Extracts the LHS variable name from each equation's AST and maps it to
#' the declaration-order variable index.  This mirrors the same mapping
#' used by extract_system_matrices() to reorder Jacobian rows, ensuring
#' that Hessian rows are permuted consistently.
#'
#' @param model dynhr_mod from parse_mod()
#' @return Integer vector of length n_eq: eq_to_decl[i] gives the
#'   declaration index of the LHS variable of equation i.
#' @noRd
.build_eq_to_decl <- function(model, f_zero = NULL) {
  endo <- model$var_names
  n_eq <- length(model$equations)
  ## Delegate to the shared package-level LHS-variable helper (.lhs_endo_var in
  ## solve-extract-system.R) -- fixes the dead "uniop" branch and eliminates
  ## the duplicate copy of this logic.
  used <- logical(length(endo)); names(used) <- endo
  eq_to_decl <- integer(n_eq)
  name_mapped <- logical(n_eq)
  for (i in seq_len(n_eq)) {
    v <- .lhs_endo_var(model$equations[[i]]$lhs)
    if (!is.null(v) && nzchar(v) && v %in% endo && !used[[v]]) {
      used[[v]] <- TRUE
      eq_to_decl[i] <- match(v, endo)
      name_mapped[i] <- TRUE
    }
  }
  ## Positional fill (base behaviour): remaining equations -> remaining
  ## variables in declaration order.  PREFERRED because it preserves the
  ## historically-validated mapping for every well-formed model -- the order-2
  ## solver, the order-3/4/5 solvers, and the analytic solution-derivative code
  ## all permute the model Hessian by this mapping and must agree.
  unassigned <- which(!(seq_along(endo) %in% eq_to_decl[eq_to_decl > 0L]))
  unmapped <- which(eq_to_decl == 0L)
  for (k in seq_along(unmapped)) {
    if (k <= length(unassigned)) eq_to_decl[unmapped[k]] <- unassigned[k]
  }
  ## H5 refinement (only with a numeric static Jacobian whose rows are in
  ## ORIGINAL equation order): a POSITIONALLY-assigned equation whose assigned
  ## variable carries ~zero contemporaneous weight is a compound-LHS
  ## mis-assignment (e.g. `a*c = ...` positionally grabbed `b`).  Free ONLY
  ## those provably-bad assignments and re-match them by max-|entry|.  Well-
  ## formed models are untouched -- their positional assignments DO load on the
  ## assigned variable -- so no validated mapping changes, while genuine
  ## compound-LHS models get the correct owner.
  if (!is.null(f_zero)) {
    bad <- vapply(seq_len(n_eq), function(i) {
      !name_mapped[i] && eq_to_decl[i] >= 1L &&
        abs(f_zero[i, eq_to_decl[i]]) <= 1e-9
    }, logical(1))
    if (any(bad)) {
      eq_to_decl[bad] <- 0L
      eq_to_decl <- .jacobian_match_unmapped(eq_to_decl, f_zero, endo)
      unassigned <- which(!(seq_along(endo) %in% eq_to_decl[eq_to_decl > 0L]))
      unmapped <- which(eq_to_decl == 0L)
      for (k in seq_along(unmapped))
        if (k <= length(unassigned)) eq_to_decl[unmapped[k]] <- unassigned[k]
    }
  }
  eq_to_decl
}


#' Compute numerical second-order Hessian of model residuals
#'
#' Finite-differences each column of the Jacobian to get d^2F/(dw_i dw_j)
#' for all pairs of compound-variable columns.
#'
#' @param dyn    compiled$dynamic (list with jacobian_fn, dyn_col_map, n_eq, total_cols)
#' @param dy_ss  Named numeric vector from .build_dy_ss_o2()
#' @param params Named numeric parameter vector
#' @param ss     Named numeric steady state
#' @param h      Step size for central finite differences (default 1e-4)
#' @return Array of dimension c(n_eq, total_cols, total_cols)
#' @noRd
.compute_model_hessian <- function(dyn, dy_ss, params, ss, h = 1e-4) {
  n_eq       <- dyn$n_eq
  total_cols <- dyn$total_cols
  n_dy       <- length(dy_ss)  # may be < total_cols if exo appended separately

  # The jacobian_fn takes a named dy vector of length n_dy (endo + exo),
  # indexed by name. total_cols == n_dy here.
  H <- array(0, dim = c(n_eq, total_cols, total_cols))

  for (c2 in seq_len(total_cols)) {
    dy_p <- dy_ss
    dy_m <- dy_ss
    # Only perturb columns that exist in the named dy vector
    if (c2 <= n_dy && c2 >= 1L) {
      step    <- max(h, abs(dy_ss[c2]) * h)
      dy_p[c2] <- dy_ss[c2] + step
      dy_m[c2] <- dy_ss[c2] - step

      Jp <- tryCatch(dyn$jacobian_fn(dy_p, params, ss),
                     error = function(e) NULL)
      Jm <- tryCatch(dyn$jacobian_fn(dy_m, params, ss),
                     error = function(e) NULL)

      if (!is.null(Jp) && !is.null(Jm)) {
        H[, , c2] <- (Jp - Jm) / (2 * step)
      }
    }
  }

  # Symmetrize to reduce numerical noise: H[e,i,j] and H[e,j,i] should agree
  for (e in seq_len(n_eq)) {
    He <- H[e, , ]
    H[e, , ] <- (He + t(He)) / 2
  }

  H
}


#' Compute symbolic second-order Hessian using precompiled expressions
#'
#' Uses the hessian2_fn and hess2_triplets stored in the compiled model to
#' evaluate the exact symbolic Hessian at the steady state.  Much faster than
#' numerical finite-differences for large models; exact for any step size.
#'
#' Falls back silently to NULL if the compiled model was built before symbolic
#' Hessian support was added (hessian2_fn absent).
#'
#' @param compiled dynhr_compiled with dynamic$hessian2_fn
#' @param dy_ss    Named numeric SS compound vector from .build_dy_ss_o2()
#' @param params   Named numeric parameter vector
#' @param ss       Named numeric steady state
#' @return Array c(n_eq, total_cols, total_cols), or NULL if unavailable
#' @noRd
.compute_model_hessian_symbolic <- function(compiled, dy_ss, params, ss) {
  dyn <- compiled$dynamic

  if (is.null(dyn$hessian2_fn) || is.null(dyn$hess2_triplets)) return(NULL)

  n_eq       <- dyn$n_eq
  total_cols <- dyn$total_cols
  n_hess     <- dyn$n_hess %||% length(dyn$hess2_triplets)

  H <- array(0, dim = c(n_eq, total_cols, total_cols))

  if (n_hess == 0L) return(H)  # log-linear model: all-zero Hessian

  ## Prefer the C++ value-vector tape (positional indexing) over the interpreted
  ## hessian2_fn closure (named-vector string lookups); fall back to the closure
  ## when the tape is absent (unsupported construct / use_rcpp off) or errors.
  ## Bit-parity to machine precision is asserted by test-hess-tape-parity.R.
  values <- .eval_triplet_tape(dyn, dyn$hess2_tape, dy_ss, params, ss)
  if (is.null(values))
    values <- tryCatch(
      dyn$hessian2_fn(dy_ss, params, ss),
      error = function(e) {
        warning(sprintf("Symbolic Hessian evaluation failed: %s. Falling back to numerical.",
                        conditionMessage(e)))
        NULL
      }
    )
  if (is.null(values)) return(NULL)

  for (k in seq_len(n_hess)) {
    t  <- dyn$hess2_triplets[[k]]
    v  <- values[k]
    H[t$eq, t$col1, t$col2] <- v
    if (t$col1 != t$col2) H[t$eq, t$col2, t$col1] <- v  # symmetry
  }

  H
}


#' Build transfer matrices from (state x, shock u) to compound-variable columns
#'
#' T_x[c, s] = d(w_c) / d(x_s):  how column c of dy responds to state var s
#' T_u[c, k] = d(w_c) / d(u_k):  how column c of dy responds to shock k
#'
#' For the three timing groups:
#'   - lag (ll == -1): y_{t-1,j} = x_{t-1,j} for state vars j  =>  T_x[c, s] = 1 if j==state_idx[s]
#'   - current (ll == 0): y_{t,j} responds via first-order policy ghx[j,], ghu[j,]
#'   - lead (ll == +1):   y_{t+1,j} responds via ghx[j,]*hx (state) and ghx[j,]*hu (shock)
#'
#' @param dyn        compiled$dynamic
#' @param ghx        Decision rules matrix (n_endo x n_state) from first-order solution
#' @param ghu        Decision rules matrix (n_endo x n_exo) from first-order solution
#' @param state_idx  Integer vector: indices of state variables in endo list
#' @param hx         n_state x n_state state transition  (= ghx[state_idx, ])
#' @param hu         n_state x n_exo shock-to-state      (= ghu[state_idx, ])
#' @param endo_names Character vector of endogenous variable names
#' @param exo_names  Character vector of exogenous variable names
#' @return List with T_x (total_cols x n_state) and T_u (total_cols x n_exo)
#' @noRd
.build_transfer_matrices <- function(dyn, ghx, ghu, state_idx, hx, hu,
                                     endo_names, exo_names) {
  n_s        <- length(state_idx)
  n_u        <- length(exo_names)
  total_cols <- dyn$total_cols
  n_dyn      <- dyn$n_dyn_cols
  dcm        <- dyn$dyn_col_map  # data.frame: name, lead_lag, col

  T_x <- matrix(0, total_cols, n_s)
  T_u <- matrix(0, total_cols, n_u)

  ghx_hx <- ghx %*% hx  # n_endo x n_state: y_{t+1} response to state via transition
  ghx_hu <- ghx %*% hu  # n_endo x n_u:     y_{t+1} response to shock via transition

  # Build reverse map: col index -> (variable name, lead_lag, is_exo, endo_j, exo_k)
  # dyn_col_map rows are ordered by col, so row k has col = dcm$col[k]
  col_name <- character(total_cols)
  col_ll   <- integer(total_cols)
  col_exo  <- logical(total_cols)

  for (k in seq_len(nrow(dcm))) {
    c  <- dcm$col[k]
    nm <- dcm$name[k]
    ll <- dcm$lead_lag[k]
    col_name[c] <- nm
    col_ll[c]   <- ll
    col_exo[c]  <- nm %in% exo_names
  }

  for (c in seq_len(total_cols)) {
    nm <- col_name[c]
    ll <- col_ll[c]

    if (col_exo[c]) {
      # Exogenous variable at time t: direct shock impact
      k_exo <- which(exo_names == nm)
      if (length(k_exo) == 1L) T_u[c, k_exo] <- 1
      # T_x stays zero (shocks are independent of lagged state)

    } else if (nchar(nm) > 0L) {
      # Endogenous variable: look up index in endo list
      j <- which(endo_names == nm)
      if (length(j) != 1L) next

      if (ll == -1L) {
        # y_{t-1}: only state variables have non-zero derivative wrt x_{t-1}
        s <- which(state_idx == j)
        if (length(s) == 1L) T_x[c, s] <- 1
        # T_u stays zero (lagged value independent of current shock)

      } else if (ll == 0L) {
        # y_t: first-order policy response
        T_x[c, ] <- ghx[j, ]
        T_u[c, ] <- ghu[j, ]

      } else if (ll == 1L) {
        # y_{t+1}: response through state transition (ghx*hx for x, ghx*hu for u)
        T_x[c, ] <- ghx_hx[j, ]
        T_u[c, ] <- ghx_hu[j, ]

      } else {
        # Higher-order lags/leads: aux-expansion should prevent this, but
        # handle gracefully by leaving T_x and T_u as zero and warning once
        # (caller issues warning if needed)
      }
    }
  }

  list(T_x = T_x, T_u = T_u)
}


#' Contract Hessian with transfer matrices to get forcing terms
#'
#' For each equation e:
#'   Phi_xx[e, ] = vec( T_x' H[e,,] T_x )  -- n_state^2 column
#'   Phi_xu[e, ] = vec( T_x' H[e,,] T_u )  -- n_state * n_exo
#'   Phi_uu[e, ] = vec( T_u' H[e,,] T_u )  -- n_exo^2
#'
#' @param H    Array c(n_eq, total_cols, total_cols) from .compute_model_hessian
#' @param T_x  Matrix total_cols x n_state
#' @param T_u  Matrix total_cols x n_exo
#' @param n_eq Number of equations
#' @return List with Phi_xx, Phi_xu, Phi_uu (each n_eq x appropriate size)
#' @noRd
.compute_phi_matrices <- function(H, T_x, T_u, n_eq) {
  n_s  <- ncol(T_x)
  n_u  <- ncol(T_u)

  Phi_xx <- matrix(0, n_eq, n_s * n_s)
  Phi_xu <- matrix(0, n_eq, n_s * n_u)
  Phi_uu <- matrix(0, n_eq, n_u * n_u)

  Tx_t <- t(T_x)  # n_state x total_cols
  Tu_t <- t(T_u)  # n_exo   x total_cols

  for (e in seq_len(n_eq)) {
    He <- H[e, , ]  # total_cols x total_cols
    if (n_s > 0L) {
      Phi_xx[e, ] <- as.vector(Tx_t %*% He %*% T_x)
    }
    if (n_s > 0L && n_u > 0L) {
      Phi_xu[e, ] <- as.vector(Tx_t %*% He %*% T_u)
    }
    if (n_u > 0L) {
      Phi_uu[e, ] <- as.vector(Tu_t %*% He %*% T_u)
    }
  }

  list(Phi_xx = Phi_xx, Phi_xu = Phi_xu, Phi_uu = Phi_uu)
}


# =====================================================================
# Main second-order solver
# =====================================================================

#' Solve the second-order perturbation of a DSGE model
#'
#' Given the first-order decision rules (ghx, ghu), computes the
#' second-order terms ghxx, ghxu, ghuu and the uncertainty correction
#' ghss using the Schmitt-Grohé & Uribe (2004) method.
#'
#' The algorithm:
#'   1. Compute numerical Hessian of model residuals at the steady state
#'   2. Build transfer matrices T_x, T_u mapping (state, shock) perturbations
#'      to the compound variable vector used by the Jacobian
#'   3. Contract Hessian with T_x, T_u to get forcing terms Phi_xx, Phi_xu, Phi_uu
#'   4. Solve the Kronecker linear system for ghxx (Sylvester equation)
#'   5. Direct linear solves for ghxu, ghuu (given ghxx)
#'   6. Direct solve for ghss using shock covariance
#'
#' @param model    dynhr_mod object from parse_mod()
#' @param compiled dynhr_compiled from compile_model()
#' @param ss       Named numeric steady state vector
#' @param params   Named numeric parameter vector
#' @param dr1      First-order DecisionRules from solve_perturbation()
#' @param Sigma_e  n_exo x n_exo shock covariance matrix. NULL = diagonal from
#'                 model shocks block (standard deviation values).
#' @param h        Step size for numerical Hessian (default 1e-4)
#' @param verbose  Print progress
#' @return A DecisionRules2 object extending DecisionRules with fields
#'   ghxx, ghxu, ghuu, ghss (and all first-order fields preserved)
#'
#' @references
#'   Schmitt-Grohé, S., & Uribe, M. (2004). Solving dynamic general equilibrium
#'     models using a second-order approximation to the policy function.
#'     \emph{Journal of Economic Dynamics and Control}, 28(4), 755-775.
#'   Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez, J. F. (2018).
#'     The pruned state-space system for non-linear DSGE models.
#'     \emph{Review of Economic Studies}, 85(1), 1-49.
#' @export
solve_perturbation_order2 <- function(model, compiled, ss, params,
                                       dr1, Sigma_e = NULL,
                                       h = 1e-4, verbose = FALSE) {
  if (!inherits(dr1, "DecisionRules")) {
    stop("dr1 must be a DecisionRules object from solve_perturbation()")
  }
  if (!isTRUE(dr1$bk_satisfied)) {
    stop("First-order Blanchard-Kahn conditions not satisfied; cannot proceed to second order.")
  }

  dyn        <- compiled$dynamic
  n          <- dr1$n_state + length(which(!dr1$endo_names %in% dr1$state_vars))
  n          <- length(dr1$endo_names)   # n_endo
  n_s        <- length(dr1$state_idx)    # n_state (= n_minus)
  n_u        <- length(dr1$exo_names)    # n_exo
  state_idx  <- dr1$state_idx
  endo_names <- dr1$endo_names
  exo_names  <- dr1$exo_names

  ghx <- dr1$ghx   # n x n_s
  ghu <- dr1$ghu   # n x n_u
  hx  <- ghx[state_idx, , drop = FALSE]  # n_s x n_s  (state transition)
  hu  <- ghu[state_idx, , drop = FALSE]  # n_s x n_u  (shock-to-state)

  if (n_s == 0L) {
    if (verbose) message("No state variables; second-order terms are all zero.")
    return(.trivial_dr2(dr1))
  }

  if (isTRUE(model$model_options$linear)) {
    if (verbose) message("Linear model: all second-order terms are exactly zero; skipping Kronecker solve.")
    return(.linear_dr2(dr1, n_s, n_u, endo_names, state_idx, exo_names, model, params))
  }

  # ----------------------------------------------------------------
  # First-order system matrices (needed for A_L)
  # ----------------------------------------------------------------
  sys <- extract_system_matrices(compiled, ss, params)
  f0  <- sys$f_zero   # n x n
  fp  <- sys$f_plus   # n x n

  # ----------------------------------------------------------------
  # Effective feedback matrix A_L = f0 + fp * ghx * S'
  # where S[j,s] = 1 iff variable j is state variable s
  # ----------------------------------------------------------------
  S   <- matrix(0, n, n_s)
  for (s in seq_len(n_s)) S[state_idx[s], s] <- 1
  A_L <- f0 + fp %*% ghx %*% t(S)   # n x n

  if (verbose) {
    cat("Second-order perturbation:\n")
    cat("  n_endo =", n, "  n_state =", n_s, "  n_exo =", n_u, "\n")
    cat("  Kronecker system size:", n * n_s^2, "x", n * n_s^2, "\n")
    n_hess_sym <- compiled$dynamic$n_hess %||% 0L
    cat("  Hessian: symbolic triplets =", n_hess_sym, "\n")
  }

  # ----------------------------------------------------------------
  # Build dy at SS and compute numerical Hessian
  # ----------------------------------------------------------------
  dy_ss <- .build_dy_ss_o2(compiled, ss)

  # Warn if model has higher-order lags/leads (aux expansion should prevent)
  all_ll <- dyn$dyn_col_map$lead_lag
  if (any(abs(all_ll) > 1L)) {
    warning(paste(
      "Model has leads/lags beyond +/-1 in dynamic Jacobian columns.",
      "Transfer matrices for |lead_lag| > 1 are set to zero.",
      "Ensure aux-expansion was applied to the model."))
  }

  # Try symbolic Hessian first (compiled at parse time via hessian2_fn);
  # fall back to numerical finite-differences if unavailable or broken.
  H <- .compute_model_hessian_symbolic(compiled, dy_ss, params, ss)
  hessian_method <- "symbolic"
  if (is.null(H)) {
    if (verbose) cat("  Symbolic Hessian unavailable; using numerical (h=", h, ")...\n")
    H <- .compute_model_hessian(dyn, dy_ss, params, ss, h)
    hessian_method <- "numerical"
  }
  if (verbose) cat("  Hessian method:", hessian_method, "\n")

  # ----------------------------------------------------------------
  # Reorder Hessian rows from compiled-equation order to
  # declaration-variable order.  The Jacobian in extract_system_matrices
  # undergoes the same reordering so that A_L, fp, fm have rows aligned
  # with declaration order.  Without this step, the Hessian rows would
  # reference equations in compiled order while A_L/fp reference them
  # in declaration order, breaking the Phi_xx forcing term.
  #
  # This applies to BOTH the symbolic (hessian2_fn) and the numerical
  # (.compute_model_hessian -> jacobian_fn FD) paths: both produce rows in
  # compiled-equation order.  Previously this was gated on the symbolic
  # path only, so the numerical fallback returned a row-misaligned ghxx.
  # ----------------------------------------------------------------
  if (!is.null(compiled$model$equations)) {
    ## Reuse the mapping that extract_system_matrices() already used to reorder
    ## the Jacobian, so the Hessian is permuted consistently.  Recomputing it
    ## here from sys$f_zero would be WRONG for compound-LHS models: that f_zero
    ## is already row-reordered, so the Jacobian fallback would see misaligned
    ## rows and return the identity, leaving the Hessian unpermuted relative to
    ## the reordered Jacobian (swapped second-order terms -- H5).
    eq_to_decl <- sys$eq_to_decl %||% .build_eq_to_decl(compiled$model)
    if (all(eq_to_decl > 0L) && !identical(eq_to_decl, seq_len(n))) {
      perm <- order(eq_to_decl)
      H_perm <- array(0, dim = dim(H))
      for (k in seq_len(n)) H_perm[k, , ] <- H[perm[k], , ]
      H <- H_perm
      if (verbose) cat("  Hessian rows reordered: compiled -> declaration order.\n")
    }
  }

  if (verbose) cat("  Building transfer matrices...\n")

  # ----------------------------------------------------------------
  # Transfer matrices
  # ----------------------------------------------------------------
  tm  <- .build_transfer_matrices(dyn, ghx, ghu, state_idx, hx, hu,
                                   endo_names, exo_names)
  T_x <- tm$T_x  # total_cols x n_s
  T_u <- tm$T_u  # total_cols x n_u

  # ----------------------------------------------------------------
  # Forcing terms Phi_xx, Phi_xu, Phi_uu
  # ----------------------------------------------------------------
  if (verbose) cat("  Computing forcing terms (Phi matrices)...\n")
  phi <- .compute_phi_matrices(H, T_x, T_u, n)

  Phi_xx <- phi$Phi_xx  # n x n_s^2
  Phi_xu <- phi$Phi_xu  # n x n_s*n_u
  Phi_uu <- phi$Phi_uu  # n x n_u^2

  # ----------------------------------------------------------------
  # Solve for ghxx: Kronecker system
  #   (I_{n_s^2} ⊗ A_L + (hx'⊗hx') ⊗ fp) * vec(ghxx) = -vec(Phi_xx)
  # ----------------------------------------------------------------
  if (verbose) cat("  Solving Kronecker system for ghxx (compact Sylvester)...\n")

  ns2 <- n_s * n_s

  # Solve the generalized-Sylvester system
  #   A_L·ghxx + fp·ghxx·(hx ⊗ hx) = -Phi_xx
  # via the compact Schur-based solver instead of forming the dense
  # n·ns² × n·ns² Kronecker matrix and LU-solving it. Forming K_xx squares the
  # model's dynamic range and is singular to machine precision on badly-scaled
  # models (e.g. Caldara_et_al_2012: gamma=40, value-function variable with
  # SS ~ 2.27e6 → rcond ~8e-21). .solve_kron_compact never forms K_xx and falls
  # back to an unconditionally-stable dense Schur solve, so it recovers these
  # cases while staying machine-precision identical on well-conditioned models.
  ghxx <- tryCatch(
    .solve_kron_compact(A_L, fp, hx, k = 2L, RHS = -Phi_xx, verbose = verbose),
    error = function(e) {
      stop(sprintf(
        "Failed to solve Kronecker system for ghxx (size %dx%d). %s",
        n * ns2, n * ns2, conditionMessage(e)))
    })
  ghxx <- matrix(ghxx, n, ns2)  # ensure n x n_s^2 layout

  # ----------------------------------------------------------------
  # Solve for ghxu (direct solve given ghxx)
  #   A_L * ghxu = -(Phi_xu + fp * ghxx * (hu ⊗ hx))
  #
  # Phi_xu cols: (state FAST, exo SLOW) = col (exo-1)*n_s + state.
  # (hu %x% hx) cols: (hu-col=exo SLOW, hx-col=state FAST) → same. ✓
  # The naive (hx %x% hx) would give (state SLOW, exo FAST) — wrong for n_u > 1.
  # ----------------------------------------------------------------
  if (verbose) cat("  Solving for ghxu...\n")
  rhs_xu <- -(Phi_xu + fp %*% ghxx %*% (hu %x% hx))  # n x n_s*n_u
  ghxu   <- tryCatch(
    solve(A_L, rhs_xu),
    error = function(e) {
      warning("solve(A_L, rhs_xu) failed; using least-squares fallback.")
      qr.solve(A_L, rhs_xu)
    })

  # ----------------------------------------------------------------
  # Solve for ghuu (direct solve given ghxx)
  #   A_L * ghuu = -(Phi_uu + fp * ghxx * (hu ⊗ hu))
  # ----------------------------------------------------------------
  if (verbose) cat("  Solving for ghuu...\n")
  rhs_uu <- -(Phi_uu + fp %*% ghxx %*% (hu %x% hu))  # n x n_u^2
  ghuu   <- tryCatch(
    solve(A_L, rhs_uu),
    error = function(e) {
      warning("solve(A_L, rhs_uu) failed; using least-squares fallback.")
      qr.solve(A_L, rhs_uu)
    })

  # ----------------------------------------------------------------
  # Shock covariance matrix
  # ----------------------------------------------------------------
  if (is.null(Sigma_e)) {
    # Full covariance from model shocks block (includes off-diagonals from
    # corr/var-pair entries; diagonal when no off-diagonals are declared).
    Sigma_e <- .get_shock_cov(model, exo_names, params)
  }

  # ----------------------------------------------------------------
  # Solve for ghss (uncertainty correction).
  #   (A_L + fp) * ghss = -(fp * ghuu * vec(Sigma_e)
  #                       +  d2f · kron(T_up, T_up) * vec(Sigma_e))
  # Mutschler (2022) eqs. 95-99; matches Dynare oo_.dr.ghs2. See
  # .solve_ghss() in R/solve-perturbation-order3-sigma.R.
  # ----------------------------------------------------------------
  if (verbose) cat("  Solving for ghss (uncertainty correction)...\n")

  has_lead <- sys$is_fwd | sys$is_mixed
  T_up <- .build_T_up(dyn, ghu, endo_names, exo_names, has_lead)
  ghss <- tryCatch(
    .solve_ghss(A_L, fp, ghuu, T_up, H, Sigma_e),
    error = function(e) {
      warning(sprintf(".solve_ghss failed (%s); ghss set to zero.",
                      conditionMessage(e)))
      numeric(n)
    })

  # ----------------------------------------------------------------
  # Name the output matrices
  # ----------------------------------------------------------------
  rownames(ghxx) <- endo_names
  rownames(ghxu) <- endo_names
  rownames(ghuu) <- endo_names
  names(ghss)    <- endo_names

  state_vars <- endo_names[state_idx]

  # Kronecker column naming convention: "FAST__x__SLOW"
  # outer(A, B) has rows = A (vary fast in as.vector), cols = B (vary slow).
  # Col index col = (j-1)*n_A + i  →  colname = A[i]__x__B[j].
  # Because ghxx/ghuu/ghxxx/ghuuu are fully symmetric in same-type indices,
  # ghxx[e, "k__x__a"] == ghxx[e, "a__x__k"] — users can look up any ordering.
  # For mixed-type matrices (ghxu, ghxuu etc.) the variable types (state vs
  # shock) make the combination unambiguous regardless of fast/slow ordering.
  sv_pairs <- outer(state_vars, state_vars, paste, sep = "__x__")
  colnames(ghxx) <- as.vector(sv_pairs)

  sx_pairs <- outer(state_vars, exo_names, paste, sep = "__x__")
  colnames(ghxu) <- as.vector(sx_pairs)

  uu_pairs <- outer(exo_names, exo_names, paste, sep = "__x__")
  colnames(ghuu) <- as.vector(uu_pairs)

  # ----------------------------------------------------------------
  # Assemble DecisionRules2 object (inherits from DecisionRules)
  # ----------------------------------------------------------------
  dr2 <- c(
    unclass(dr1),          # all first-order fields
    list(
      ghxx              = ghxx,
      ghxu              = ghxu,
      ghuu              = ghuu,
      ghss              = ghss,
      Sigma_e           = Sigma_e,
      order             = 2L,
      hessian_method    = hessian_method,
      hessian_step      = h,
      second_order_ok   = TRUE
    )
  )
  class(dr2) <- c("DecisionRules2", "DecisionRules")

  if (verbose) {
    cat("Second-order solution complete.\n")
    cat("  ghxx:", nrow(ghxx), "x", ncol(ghxx), "\n")
    cat("  ghxu:", nrow(ghxu), "x", ncol(ghxu), "\n")
    cat("  ghuu:", nrow(ghuu), "x", ncol(ghuu), "\n")
    cat("  max|ghss|:", round(max(abs(ghss)), 6), "\n")
  }

  dr2
}


#' Second-order solution for linear models (all higher-order terms exactly zero)
#' @noRd
.linear_dr2 <- function(dr1, n_s, n_u, endo_names, state_idx, exo_names,
                        model, params) {
  n   <- length(endo_names)
  state_vars <- endo_names[state_idx]

  ghxx <- matrix(0, n, n_s^2)
  ghxu <- matrix(0, n, n_s * n_u)
  ghuu <- matrix(0, n, n_u^2)

  colnames(ghxx) <- as.vector(outer(state_vars, state_vars, paste, sep = "__x__"))
  colnames(ghxu) <- as.vector(outer(state_vars, exo_names,  paste, sep = "__x__"))
  colnames(ghuu) <- as.vector(outer(exo_names,  exo_names,  paste, sep = "__x__"))
  rownames(ghxx) <- rownames(ghxu) <- rownames(ghuu) <- endo_names

  stderr  <- .get_shock_stderr(model, exo_names, params)
  Sigma_e <- diag(stderr^2, n_u, n_u)
  if (!is.null(exo_names)) colnames(Sigma_e) <- rownames(Sigma_e) <- exo_names

  dr2 <- c(unclass(dr1), list(
    ghxx           = ghxx,
    ghxu           = ghxu,
    ghuu           = ghuu,
    ghss           = setNames(numeric(n), endo_names),
    Sigma_e        = Sigma_e,
    order          = 2L,
    hessian_method = "linear_shortcircuit",
    hessian_step   = NA_real_,
    second_order_ok = TRUE
  ))
  class(dr2) <- c("DecisionRules2", "DecisionRules")
  dr2
}


#' Trivial second-order solution (no state variables)
#' @noRd
.trivial_dr2 <- function(dr1) {
  n   <- length(dr1$endo_names)
  n_u <- length(dr1$exo_names)
  dr2 <- c(unclass(dr1), list(
    ghxx = matrix(0, n, 0), ghxu = matrix(0, n, 0),
    ghuu = matrix(0, n, n_u^2), ghss = numeric(n),
    Sigma_e = diag(n_u), order = 2L, hessian_step = NA_real_,
    second_order_ok = TRUE
  ))
  class(dr2) <- c("DecisionRules2", "DecisionRules")
  dr2
}


# =====================================================================
# S3 methods for DecisionRules2
# =====================================================================

#' Print method for second-order decision rules
#' @param x   A \code{DecisionRules2} object.
#' @param ... Unused; included for S3 compatibility.
#' @export
print.DecisionRules2 <- function(x, ...) {
  cat("Second-order Decision Rules (DecisionRules2)\n")
  cat("  Endogenous variables:", length(x$endo_names), "\n")
  cat("  State variables:     ", x$n_state, "\n")
  cat("  Shocks:              ", x$n_exo, "\n")
  cat("  BK satisfied:        ", x$bk_satisfied, "\n")
  cat("  Perturbation order:  2\n")
  cat("  Hessian method:     ", x$hessian_method %||% "unknown", "\n")
  cat("\nFirst-order decision rules:\n")
  cat("  ghx:", nrow(x$ghx), "x", ncol(x$ghx), "\n")
  cat("  ghu:", nrow(x$ghu), "x", ncol(x$ghu), "\n")
  cat("\nSecond-order terms:\n")
  if (!is.null(x$ghxx) && prod(dim(x$ghxx)) > 0) {
    cat("  ghxx:", nrow(x$ghxx), "x", ncol(x$ghxx),
        " max|.| =", format(max(abs(x$ghxx)), digits = 4), "\n")
    cat("  ghxu:", nrow(x$ghxu), "x", ncol(x$ghxu),
        " max|.| =", format(max(abs(x$ghxu)), digits = 4), "\n")
    cat("  ghuu:", nrow(x$ghuu), "x", ncol(x$ghuu),
        " max|.| =", format(max(abs(x$ghuu)), digits = 4), "\n")
    cat("  ghss: length", length(x$ghss),
        " max|.| =", format(max(abs(x$ghss)), digits = 4), "\n")
  } else {
    cat("  (no state variables; all zero)\n")
  }
  invisible(x)
}


# =====================================================================
# Second-order simulation and IRFs (pruned state space)
# =====================================================================

#' Impulse response functions at second order (pruned approximation)
#'
#' Computes IRFs using the pruned state-space system of Andreasen et al.
#' (2018). At each step:
#'   \eqn{x1_t = hx \cdot x1_{t-1} + hu \cdot \varepsilon_t}  (first-order state)
#'   \eqn{x2_t = hx \cdot x2_{t-1} + (1/2) hxx \cdot (x1 \otimes x1) + hxu \cdot (x1 \otimes \varepsilon) + (1/2) huu \cdot (\varepsilon \otimes \varepsilon) + (1/2) hss}
#'   \eqn{y_t  = y1_t + y2_t}
#' where \eqn{y1_t = ghx \cdot x1 + ghu \cdot \varepsilon} and y2_t is the second-order correction.
#'
#' @param dr2        DecisionRules2 object
#' @param model      dynhr_mod (for shock variances)
#' @param n_periods  Number of IRF periods
#' @param shock_size Shock size in std dev units
#' @param params     Named parameter vector
#' @param pruning    Use pruned state space (TRUE, recommended for order=2)
#' @return List of matrices, one per shock. Each n_periods x n_endo.
#'
#' @details
#' **Steady-state baseline (dynhr vs. Dynare).**
#' dynhr reports IRFs as deviations from the **deterministic** steady state
#' (\eqn{y^*}, the solution to the non-stochastic model).  Dynare's
#' \code{stoch_simul} order-2 IRFs are instead reported as deviations from the
#' **stochastic** steady state \eqn{y^* + 0.5 \cdot \texttt{ghss}}, which
#' includes the constant second-order correction for uncertainty.
#'
#' Consequently, a direct period-by-period comparison between dynhr and Dynare
#' order-2 IRFs will show a systematic offset of approximately
#' \eqn{0.5 \cdot \texttt{ghss}} (applied to the relevant endogenous variables
#' at every horizon).  The discrepancy is largest where \code{ghss} / \code{ghs2}
#' is largest, and is typically on the order of 1e-4 to 1e-3 for standard
#' calibrations.
#'
#' To reconcile with Dynare, either:
#' \itemize{
#'   \item Subtract \eqn{0.5 \cdot \texttt{ghss}} from each dynhr IRF column
#'         (shifting to the stochastic-SS baseline), or
#'   \item Compare the **dynamics** (period-over-period changes) rather than
#'         levels, since the offset is constant across horizons.
#' }
#' @export
compute_irfs_order2 <- function(dr2, model, n_periods = 40L,
                                 shock_size = 1,
                                 params = NULL,
                                 pruning = TRUE) {
  if (!inherits(dr2, "DecisionRules2")) {
    stop("dr2 must be a DecisionRules2 object.")
  }

  ghx  <- dr2$ghx; ghu <- dr2$ghu
  ghxx <- dr2$ghxx; ghxu <- dr2$ghxu
  ghuu <- dr2$ghuu; ghss <- dr2$ghss

  endo      <- dr2$endo_names
  exo       <- dr2$exo_names
  state_idx <- dr2$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)
  n_s       <- length(state_idx)

  if (is.null(params)) params <- model$param_values
  shock_stderr <- .get_shock_stderr(model, exo, params)

  hx  <- ghx[state_idx, , drop = FALSE]
  hu  <- ghu[state_idx, , drop = FALSE]

  if (!is.null(ghxx) && nrow(ghxx) == n_endo) {
    hxx <- ghxx[state_idx, , drop = FALSE]  # n_s x n_s^2
    hxu <- ghxu[state_idx, , drop = FALSE]  # n_s x n_s*n_exo
    huu <- ghuu[state_idx, , drop = FALSE]  # n_s x n_exo^2
    hss <- ghss[state_idx]                  # n_s
  } else {
    hxx <- matrix(0, n_s, n_s^2)
    hxu <- matrix(0, n_s, n_s * n_exo)
    huu <- matrix(0, n_s, n_exo^2)
    hss <- numeric(n_s)
  }

  irfs <- list()
  for (k in seq_along(exo)) {
    shock_name <- exo[k]
    irf_mat    <- matrix(0, n_periods, n_endo)
    colnames(irf_mat) <- endo
    rownames(irf_mat) <- paste0("t", seq_len(n_periods))

    eps    <- numeric(n_exo)
    eps[k] <- shock_stderr[shock_name] * shock_size

    # Pruned state-space: compute DEVIATIONS from the stochastic steady state.
    # ghss/hss are constant mean corrections already embedded in the baseline,
    # so they cancel out and are NOT included in deviation IRFs.
    x1 <- numeric(n_s)   # first-order state deviation
    x2 <- numeric(n_s)   # second-order state correction (deviation from mean)

    # Period 1: impact
    eps_vec <- as.numeric(eps)

    # First-order state: x1 = hu * eps (from x1_prev = 0)
    x1 <- as.numeric(hu %*% eps_vec)

    # Second-order state deviation (x1_prev=0 at SS, so hxx and hxu vanish)
    x2 <- as.numeric(0.5 * huu %*% (eps_vec %x% eps_vec))

    # First-order output deviation
    y1 <- as.numeric(ghu %*% eps_vec)
    # Second-order output deviation (no ghss: it's the SS-level correction)
    y2 <- as.numeric(0.5 * ghuu %*% (eps_vec %x% eps_vec))

    irf_mat[1, ] <- y1 + y2

    # Periods 2..n_periods: no further shocks (pure propagation)
    for (t in 2:n_periods) {
      x1_prev <- x1
      x2_prev <- x2

      x1 <- as.numeric(hx %*% x1_prev)
      x2 <- as.numeric(
        hx %*% x2_prev +
        0.5 * hxx %*% (x1_prev %x% x1_prev)
        # hxu and huu terms are zero (no shock after period 1)
        # hss is excluded: it's the constant deviation of ergodic mean from SS
      )

      y1 <- as.numeric(ghx %*% x1_prev)
      y2 <- as.numeric(
        ghx %*% x2_prev +
        0.5 * ghxx %*% (x1_prev %x% x1_prev)
        # ghss excluded for same reason
      )

      irf_mat[t, ] <- y1 + y2
    }

    irfs[[shock_name]] <- irf_mat
  }

  class(irfs) <- "IRFCollection"
  attr(irfs, "n_periods")  <- n_periods
  attr(irfs, "endo_names") <- endo
  attr(irfs, "exo_names")  <- exo
  attr(irfs, "order")      <- 2L
  irfs
}


#' Simulate the model forward at second order (pruned state space)
#'
#' @param dr2         DecisionRules2 object
#' @param n_periods   Number of simulation periods
#' @param shocks      Matrix (n_periods + burn_in) x n_exo. If NULL, draws random.
#' @param model       dynhr_mod (for shock variances)
#' @param burn_in     Burn-in periods to discard
#' @param pruning     Use pruned state-space (default TRUE)
#' @param init_state  Optional named numeric vector of initial state deviations
#'   loaded into the pruned first-order component; pair with \code{burn_in = 0}.
#'   \code{NULL} starts at the steady state.
#' @return Matrix n_periods x n_endo (deviations from SS)
#' @export
simulate_model_order2 <- function(dr2, n_periods = 200L, shocks = NULL,
                                   model = NULL, burn_in = 100L,
                                   pruning = TRUE, init_state = NULL) {
  if (!inherits(dr2, "DecisionRules2")) {
    stop("dr2 must be a DecisionRules2 object.")
  }

  ghx  <- dr2$ghx; ghu  <- dr2$ghu
  ghxx <- dr2$ghxx; ghxu <- dr2$ghxu
  ghuu <- dr2$ghuu; ghss <- dr2$ghss

  endo      <- dr2$endo_names
  exo       <- dr2$exo_names
  state_idx <- dr2$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)
  n_s       <- length(state_idx)
  params    <- if (!is.null(model)) model$param_values else NULL

  hx  <- ghx[state_idx, , drop = FALSE]
  hu  <- ghu[state_idx, , drop = FALSE]
  hxx <- ghxx[state_idx, , drop = FALSE]
  hxu <- ghxu[state_idx, , drop = FALSE]
  huu <- ghuu[state_idx, , drop = FALSE]
  hss <- ghss[state_idx]

  ## Near-unit-root guard: pruned order-2 simulation mixes slowly when the
  ## dominant eigenvalue of hx is close to 1.  The ergodic mean may require
  ## tens of thousands of periods to be well estimated; use
  ## compute_moments_order2()$mean for the analytically exact ergodic mean.
  max_eig <- max(Mod(eigen(hx, only.values = TRUE)$values))
  ## Threshold 0.99: a genuine near-unit-root cutoff.  A lower value (e.g.
  ## 0.95) false-alarms on standard calibrations like RBC (rho = 0.95), where
  ## the default burn-in + n_periods is adequate.
  if (isTRUE(getOption("dynhr.warn_near_unit_root", TRUE)) &&
      max_eig > 0.99) {
    warning(sprintf(
      paste0("simulate_model_order2(): near-unit-root state ",
             "(max |eigenvalue(hx)| = %.4f > 0.95). ",
             "Pruned simulation may require many more than %d periods to reach ",
             "the ergodic distribution. ",
             "Use compute_moments_order2()$mean for the analytical ergodic mean."),
      max_eig, n_periods),
      call. = FALSE)
  }

  total_periods <- n_periods + burn_in

  shock_stderr <- .get_shock_stderr(model, exo, params)
  if (is.null(shocks)) {
    shocks <- matrix(rnorm(total_periods * n_exo), ncol = n_exo)
    for (k in seq_along(exo)) shocks[, k] <- shocks[, k] * shock_stderr[exo[k]]
  }

  sim <- matrix(0, total_periods, n_endo)
  colnames(sim) <- endo

  x1 <- numeric(n_s)
  x2 <- numeric(n_s)

  ## Condition on a custom starting state: load the initial deviation into the
  ## first-order pruned component x1 (ordered as endo[state_idx]); the second-
  ## order correction x2 starts at its steady value 0.  Pair with burn_in = 0.
  if (!is.null(init_state)) {
    idx <- match(endo[state_idx], names(init_state))
    ok  <- !is.na(idx)
    if (any(ok)) x1[ok] <- as.numeric(init_state[idx[ok]])
  }

  for (t in seq_len(total_periods)) {
    e  <- shocks[t, ]
    x1_prev <- x1
    x2_prev <- x2

    # Update states.  ghxu / hxu cols are (state FAST, exo SLOW) from
    # outer(state_vars, exo_names), so the matching Kronecker vector is
    # (e %x% x), giving rows (exo SLOW, state FAST).  (x %x% e) would mix
    # the indices and silently produce wrong values for n_s>1, n_u>1.
    x1 <- as.numeric(hx %*% x1_prev + hu %*% e)
    if (pruning) {
      x2 <- as.numeric(
        hx  %*% x2_prev +
        0.5 * hxx %*% (x1_prev %x% x1_prev) +
        hxu %*% (e %x% x1_prev) +
        0.5 * huu %*% (e %x% e) +
        0.5 * hss
      )
    } else {
      x2 <- numeric(n_s)
    }

    # Output: first-order + second-order correction
    y1 <- as.numeric(ghx %*% x1_prev + ghu %*% e)
    y2 <- as.numeric(
      ghx  %*% x2_prev +
      0.5 * ghxx %*% (x1_prev %x% x1_prev) +
      ghxu %*% (e %x% x1_prev) +
      0.5 * ghuu %*% (e %x% e) +
      0.5 * ghss
    )

    sim[t, ] <- y1 + y2
  }

  sim <- sim[(burn_in + 1L):total_periods, , drop = FALSE]

  # Add SS levels
  sim_levels <- sim
  for (j in seq_along(endo)) {
    ss_val <- dr2$ys[endo[j]]
    if (!is.na(ss_val)) sim_levels[, j] <- sim[, j] + ss_val
  }
  attr(sim, "levels") <- sim_levels
  sim
}
