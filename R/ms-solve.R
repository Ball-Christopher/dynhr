## R/ms-solve.R
## --------------------------------------------------------------------------
## Structural Markov-switching DSGE perturbation solver.
##
## Implements the Maih (2015) functional iteration for the coupled first-order
## decision rules {ghx_s, ghu_s} in a 2+ regime MS-DSGE model where structural
## parameters (and hence the Jacobian blocks) differ across regimes.
##
## Reference: Maih, J. (2015). Efficient perturbation methods for solving
##   regime-switching DSGE models. Norges Bank Working Paper 1/2015.
##
## COUPLED SYSTEM (for state-block G_s = ghx_s[state_idx, ]):
##
##   For each regime s, simultaneously:
##     (A^s_0 + A^s_+ * Phi_s) * G_s + A^s_- = 0
##
##   where Phi_s = sum_{s'} P[s,s'] * G_{s'} (transition-weighted coupling).
##
## FUNCTIONAL ITERATION (Eq. C in the brief):
##   Initialise G_s^{(0)} from single-regime QZ per regime.
##   At each iteration k+1:
##     Phi_s^{(k)} = sum_{s'} P[s,s'] * G_{s'}^{(k)}
##     Solve: (A^s_0 + A^s_+ * Phi_s^{(k)}) * G_s^{(k+1)} = -A^s_-
##   Until max_s max|G_s^{(k+1)} - G_s^{(k)}| < tol.
##
## REDUCTION ORACLE: when all params_by_regime are identical, every G_s must
##   converge to the single-regime QZ solution to 1e-10, for any P.
## --------------------------------------------------------------------------


#' Solve structural Markov-switching DSGE perturbation (coupled decision rules)
#'
#' Computes regime-specific first-order decision rules (\code{ghx_s},
#' \code{ghu_s}) for an MS-DSGE model where structural parameters differ across regimes,
#' using the Maih (2015) functional iteration.
#'
#' @param model      dynhr_mod object (same model for all regimes; params
#'                   differ per regime via \code{params_by_regime}).
#' @param compiled   dynhr_compiled (same Jacobian functions used for all
#'                   regimes; only the evaluated param/ss point changes).
#' @param ss_by_regime  List of length h; each element is a named numeric
#'                   steady-state vector (output of \code{solve_steady}).
#' @param params_by_regime  List of length h; each element is a named numeric
#'                   parameter vector for that regime.
#' @param P          h x h row-stochastic transition matrix; entry
#'                   \code{P[s, sp]} is the probability of regime sp at t+1
#'                   given regime s at t.
#' @param max_iter   Maximum functional-iteration steps (default 500).
#' @param tol        Convergence tolerance: max over regimes of max absolute
#'                   change in G_s (default 1e-10).
#' @param verbose    Logical; print iteration progress (default FALSE).
#'
#' @return An object of class \code{"MsDecisionRules"} with fields:
#'   \describe{
#'     \item{\code{dr}}{Length-h list of \code{DecisionRules} objects, one per
#'       regime. Each has \code{ghx}, \code{ghu}, \code{ys}, \code{state_idx},
#'       \code{bk_satisfied}, plus \code{c_const} (constant drift term).}
#'     \item{\code{P}}{The transition matrix.}
#'     \item{\code{pi0}}{Ergodic distribution of P.}
#'     \item{\code{ss_list}}{List of per-regime steady states.}
#'     \item{\code{converged}}{Logical: TRUE if iteration converged.}
#'     \item{\code{n_iter}}{Number of iterations until convergence (or
#'       \code{max_iter} if not converged).}
#'     \item{\code{mss_ok}}{Logical: TRUE if the mean-square stability
#'       Kronecker condition holds (spectral radius < 1).}
#'   }
#'
#' @details
#' \strong{NON-CONVERGENCE:} If the iteration does not converge (max_iter
#' hit, divergence, or non-finite iterates), the function \code{stop()}s with
#' a diagnostic message. It never silently returns an unconverged solution.
#'
#' \strong{Reduction oracle:} when all \code{params_by_regime} elements are
#' identical, the converged \code{ghx_s} must equal the single-regime
#' \code{solve_perturbation} output to within 1e-10.
#'
#' @export
solve_ms_perturbation <- function(model, compiled,
                                   ss_by_regime, params_by_regime,
                                   P,
                                   max_iter = 500L,
                                   tol      = 1e-10,
                                   verbose  = FALSE) {

  ## ---- input validation ---------------------------------------------------
  h <- nrow(P)
  if (!is.matrix(P) || nrow(P) != ncol(P))
    stop("solve_ms_perturbation: P must be a square matrix.", call. = FALSE)
  if (h < 2L)
    stop("solve_ms_perturbation: P must be at least 2x2.", call. = FALSE)
  rs <- rowSums(P)
  if (any(abs(rs - 1) > 1e-10))
    stop(sprintf(
      "solve_ms_perturbation: P rows must sum to 1; max deviation %.2e.",
      max(abs(rs - 1))), call. = FALSE)
  if (any(P < 0))
    stop("solve_ms_perturbation: P must be non-negative.", call. = FALSE)

  if (!is.list(ss_by_regime) || length(ss_by_regime) != h)
    stop(sprintf(
      "solve_ms_perturbation: ss_by_regime must be a list of length %d.", h),
      call. = FALSE)
  if (!is.list(params_by_regime) || length(params_by_regime) != h)
    stop(sprintf(
      "solve_ms_perturbation: params_by_regime must be a list of length %d.", h),
      call. = FALSE)

  ## ---- cache system structure (shared across all regimes and iterations) ---
  sys_cache <- cache_system_structure(compiled)
  n_endo    <- sys_cache$n_endo
  n_exo     <- sys_cache$n_exo
  endo      <- sys_cache$endo
  exo       <- sys_cache$exo

  ## ---- Step 1: per-regime system matrices + starting values from QZ --------
  ## For each regime s: extract (f_minus^s, f_zero^s, f_plus^s, f_exo^s) and
  ## run the single-regime QZ solver. The QZ solution G_s^{(0)} is the best
  ## starting value for the functional iteration.

  sys_list  <- vector("list", h)   # raw system matrices per regime
  dr0_list  <- vector("list", h)   # initial (single-regime) decision rules

  for (s in seq_len(h)) {
    ss_s     <- ss_by_regime[[s]]
    params_s <- params_by_regime[[s]]

    ## Validate SS
    if (!is.numeric(ss_s) || is.null(names(ss_s)))
      stop(sprintf(
        "solve_ms_perturbation: ss_by_regime[[%d]] must be a named numeric vector.", s),
        call. = FALSE)
    if (any(!is.finite(ss_s)))
      stop(sprintf(
        "solve_ms_perturbation: ss_by_regime[[%d]] has non-finite values; steady-state solve may have failed.", s),
        call. = FALSE)

    sys_s <- extract_system_matrices_fast(sys_cache, ss_s, params_s)
    sys_list[[s]] <- sys_s

    ## Single-regime QZ solution as starting value
    dr0_s <- tryCatch(
      .solve_from_system(sys_s, model, compiled, ss_s, params_s, FALSE),
      error = function(e) {
        stop(sprintf(
          "solve_ms_perturbation: single-regime QZ solve failed for regime %d: %s",
          s, conditionMessage(e)), call. = FALSE)
      }
    )

    if (!isTRUE(dr0_s$bk_satisfied))
      warning(sprintf(
        "solve_ms_perturbation: BK condition violated in regime %d starting values; iteration may not converge.", s),
        call. = FALSE)

    dr0_list[[s]] <- dr0_s
    if (verbose)
      cat(sprintf("  Regime %d: QZ init done, BK=%s, n_state=%d\n",
                  s, dr0_s$bk_satisfied, dr0_s$n_state))
  }

  ## Recover shared dimension info from regime 1 (all regimes share the same
  ## model structure; only params/SS differ).
  state_idx <- dr0_list[[1L]]$state_idx
  n_state   <- length(state_idx)

  if (n_state == 0L) {
    ## Degenerate case: no state variables. Return QZ solutions directly.
    warning("solve_ms_perturbation: no state variables; returning single-regime QZ solutions.", call. = FALSE)
    dr_list <- dr0_list
    for (s in seq_len(h)) {
      dr_list[[s]]$c_const     <- numeric(n_endo)
      dr_list[[s]]$mss_satisfied <- TRUE
    }
    return(.build_ms_dr(dr_list, P, ss_by_regime, converged = TRUE,
                         n_iter = 0L, mss_ok = TRUE))
  }

  ## ---- Step 2: extract the "reduced" system blocks needed for iteration ----
  ## For the functional iteration, we work with the dynamic block (post-QR
  ## static elimination) in the reduced n_d x n_minus form. The brief shows
  ## that after static elimination, the relevant matrices are:
  ##   A^s_+ (dynamic rows of Qf_plus, columns for forward vars)  -- n_d x n_plus
  ##   A^s_0 (dynamic rows of Qf_zero, all dynamic cols)          -- n_d x n_d
  ##   A^s_- (dynamic rows of Qf_minus, columns for state vars)   -- n_d x n_minus
  ##
  ## HOWEVER: because we're doing a direct linear solve
  ##   (A^s_0 + A^s_+ * Phi_s) * G_s = -A^s_-
  ## we need these in the companion-pencil format. Instead of re-implementing
  ## all of .solve_from_system's QR machinery, we use a smarter approach:
  ## the QZ starting value G_s^{(0)} already lives in the correct space.
  ##
  ## KEY INSIGHT: In the single-regime QZ, the pencil is:
  ##   f_plus^s * ghx_s * ghx_s + f_zero^s * ghx_s + f_minus^s = 0  (state block)
  ## In the MS iteration, we replace ghx_s * ghx_s with Phi_s * ghx_s:
  ##   f_plus^s * Phi_s * ghx_s + f_zero^s * ghx_s + f_minus^s = 0
  ## i.e. (f_zero^s + f_plus^s * Phi_s) * ghx_s = -f_minus^s  (state cols only)
  ##
  ## We work in the FULL n_endo space using the original (un-QR'd) matrices,
  ## solving the state column of ghx directly. This matches the brief's Eq. C
  ## and avoids reconstructing the static elimination machinery.
  ##
  ## More precisely: ghx is n_endo x n_state. The system for the state block is:
  ##   (f_zero^s + f_plus^s * Phi_s_full) * ghx_s = -f_minus_state^s
  ## where Phi_s_full[i,j] = ghx_s[j, ] for state variables (full n_endo x n_endo
  ## matrix whose (i,j)-entry is ghx_{s'}[state_j -> endo_i]).
  ##
  ## We form Phi_s as an n_endo x n_endo matrix (zeros for non-state columns)
  ## and solve the full n_endo x n_state system.

  ## Build Phi_s_full: n_endo x n_endo. Only the state_idx columns are nonzero.
  ## Phi_s_full[, state_idx] = sum_{s'} P[s,s'] * ghx_{s'}[, state_idx via state rows]
  ## Actually Phi_s in the brief is n_state x n_state (state-to-state block).
  ## The full-system analogue for the ghu solve is:
  ##   (f_zero^s + f_plus^s * Phi_s_embed) * ghu_s = -f_exo^s
  ## where Phi_s_embed is n_endo x n_endo with Phi_s_embed[, c] = 0 for non-state c,
  ## and Phi_s_embed[:, state_idx] = (sum P[s,s'] ghx_{s'})[:,state_idx] -- but
  ## that's just the full ghx P-average.

  ## Initialise G_s = current ghx full (n_endo x n_state), from QZ starting values.
  G_list <- lapply(dr0_list, function(dr) dr$ghx)  # each n_endo x n_state

  ## Precompute: f_minus[:, state_idx] is the "right-hand side forcing" for the
  ## state block. This is the state-column sub-matrix of f_minus.
  ## For full-system solve, RHS = -f_minus_state where f_minus_state = f_minus[, state_idx].
  ## But state_idx refers to endo indices. f_minus is n_eq x n_endo.
  ## The system (f_zero + f_plus * Phi_embed) * X = -f_minus is n_eq x n_state,
  ## solving for X = ghx (n_endo x n_state) -- BUT the system is rectangular
  ## (n_eq rows = n_endo rows since n_eq == n_endo after reconciliation).

  ## Actually: from solve_perturbation line 1042-1050, after ghx is known:
  ##   P_mat[, state_idx] <- ghx   (n_endo x n_endo, zeros elsewhere)
  ##   M = f_zero + f_plus %*% P_mat
  ##   ghx_cols = solve(M, -f_minus[, state_idx])  -- but this is NOT how ghx was found
  ## ghx comes from QZ, not from a linear solve in the full space.
  ##
  ## For the MS iteration, the iteration step at regime s is:
  ##   Build Phi_embed = matrix(0, n_endo, n_endo); Phi_embed[, state_idx] = t(Phi_s @ t(ghx))
  ##   No -- Phi_s (n_state x n_state) acts on the state columns of ghx.
  ##   The MS Sylvester (Eq. C in brief, after static elim) is a n_d x n_minus system.
  ##
  ## CORRECT FORMULATION for full-space direct solve:
  ## Following the brief §3 step 2: "Solve (A^s_0 + A^s_+ * Phi_s^{(k)}) * G_s = -A^s_-"
  ## These are the reduced (post-QR) matrices. We need to replicate the static elim.
  ##
  ## SIMPLER ALTERNATIVE: use solve_perturbation_fast() with a MODIFIED compiled
  ## model where f_plus is replaced by f_plus * Phi_s (i.e. treat Phi_s as a
  ## "forced" one-step-ahead slope). We pass a custom sys to .solve_from_system.
  ##
  ## This is the cleanest path: build a modified sys_s where f_plus_modified =
  ## f_plus^s * Phi_embed and call .solve_from_system with that modified system.
  ## The QZ then solves: (f_plus_mod * ghx + f_zero) * ghx + f_minus = 0
  ## which is wrong -- QZ expects a quadratic, not our linear system.
  ##
  ## CORRECT APPROACH (direct linear solve on reduced system):
  ## We re-implement the static elimination inline for each iteration to get
  ## (A^s_0 + A^s_+ * Phi_s) and -A^s_-, then solve directly.
  ## See: .ms_reduced_system() helper below.

  if (verbose) cat(sprintf("  MS iteration: h=%d, n_state=%d, tol=%.1e, max_iter=%d\n",
                            h, n_state, tol, max_iter))

  ## ---- Step 3: Functional iteration ----------------------------------------
  iter_result <- .ms_fixed_point_iteration(
    sys_list    = sys_list,
    G_init      = G_list,
    P           = P,
    state_idx   = state_idx,
    n_endo      = n_endo,
    h           = h,
    max_iter    = max_iter,
    tol         = tol,
    verbose     = verbose
  )

  if (!iter_result$converged) {
    stop(sprintf(paste0(
      "solve_ms_perturbation: functional iteration did NOT converge after %d iterations.\n",
      "  Final max change = %.4e (tol = %.1e).\n",
      "  Regimes may be structurally dissimilar or near the MS-determinacy boundary.\n",
      "  Try: (1) increase max_iter, (2) different params, or (3) check BK per regime."),
      max_iter, iter_result$final_change, tol), call. = FALSE)
  }

  G_conv <- iter_result$G_list   # converged ghx, each n_endo x n_state

  if (verbose)
    cat(sprintf("  Converged in %d iterations, final change = %.2e\n",
                iter_result$n_iter, iter_result$final_change))

  ## ---- Step 4: Recover ghu_s for each regime (Eq. B) ----------------------
  ## ghu_s = -solve(f_zero^s + f_plus^s * Phi_s, f_exo^s)
  ## where Phi_s is the converged coupling (n_endo x n_endo embed matrix).

  dr_list <- vector("list", h)

  for (s in seq_len(h)) {
    ## Build converged Phi_s embed (n_endo x n_endo)
    Phi_s_embed <- .ms_build_Phi_embed(G_conv, P, s, state_idx, n_endo)

    sys_s    <- sys_list[[s]]
    f_zero_s <- sys_s$f_zero
    f_plus_s <- sys_s$f_plus
    f_exo_s  <- sys_s$f_exo

    M_s  <- f_zero_s + f_plus_s %*% Phi_s_embed
    ghu_s <- tryCatch(
      solve(M_s, -f_exo_s),
      error = function(e) .safe_inv(M_s) %*% (-f_exo_s)
    )

    rownames(ghu_s) <- endo
    if (n_exo > 0) colnames(ghu_s) <- exo

    ## ghx_s from converged G_s (already named from QZ path)
    ghx_s <- G_conv[[s]]
    if (!is.null(rownames(dr0_list[[s]]$ghx))) rownames(ghx_s) <- endo
    if (!is.null(colnames(dr0_list[[s]]$ghx))) colnames(ghx_s) <- colnames(dr0_list[[s]]$ghx)

    ## Constant drift term c_s = ghx_s * sum_{s'} P[s,s'] * (ys_{s'} - ys_s)
    ## (FRWZ 2016 Proposition 2 — regime-specific mean shift)
    ys_s <- ss_by_regime[[s]]
    c_s_arg <- numeric(n_state)
    for (sp in seq_len(h)) {
      diff_ss <- ss_by_regime[[sp]][endo[state_idx]] - ys_s[endo[state_idx]]
      diff_ss[is.na(diff_ss)] <- 0
      c_s_arg <- c_s_arg + P[s, sp] * diff_ss
    }
    c_const_s <- drop(ghx_s[state_idx, , drop = FALSE] %*% c_s_arg)

    ## Per-regime BK check (spectral radius of converged ghx state block)
    state_block_s <- ghx_s[state_idx, , drop = FALSE]
    sr_s <- tryCatch(
      max(Mod(eigen(state_block_s, only.values = TRUE)$values)),
      error = function(e) NA_real_)
    bk_s <- is.finite(sr_s) && sr_s <= 1 + 1e-6

    dr_s <- list(
      ghx          = ghx_s,
      ghu          = ghu_s,
      ys           = ys_s,
      endo_names   = endo,
      exo_names    = exo,
      state_vars   = dr0_list[[1L]]$state_vars,
      state_idx    = state_idx,
      n_state      = n_state,
      n_exo        = n_exo,
      bk_satisfied = bk_s,
      c_const      = c_const_s
    )
    class(dr_s) <- "DecisionRules"
    dr_list[[s]] <- dr_s
  }

  ## ---- Step 5: Mean-square stability check ---------------------------------
  ## M_Kron = sum_{s,s'} P[s,s'] * (ghx_{s'} kron ghx_{s'}) -- spectral radius < 1
  mss_ok <- .ms_check_mss(G_conv, P, state_idx, h)
  if (verbose)
    cat(sprintf("  Mean-square stability: %s\n", if (mss_ok) "OK" else "VIOLATED"))

  for (s in seq_len(h))
    dr_list[[s]]$mss_satisfied <- mss_ok

  .build_ms_dr(dr_list, P, ss_by_regime, converged = TRUE,
                n_iter = iter_result$n_iter, mss_ok = mss_ok)
}


## ============================================================================
## Internal: fixed-point iteration loop
## ============================================================================

## Run the Maih (2015) functional iteration.
##
## At each step k:
##   For each regime s:
##     1. Phi_s^{(k)} = sum_{s'} P[s,s'] * G_{s'}^{(k)}  (n_endo x n_state embed)
##     2. Solve (f_zero^s + f_plus^s * Phi_embed^{(k)}) * G_s^{(k+1)} = -f_minus_state^s
##
## Step 2 uses the FULL space solve (n_endo x n_state RHS).
## This avoids re-implementing static QR elimination per iteration.
## The full-space system is square (n_endo x n_endo LHS, n_endo x n_state RHS).
##
## @return list(G_list, converged, n_iter, final_change)
## @noRd
.ms_fixed_point_iteration <- function(sys_list, G_init, P, state_idx,
                                       n_endo, h, max_iter, tol, verbose) {
  G_list <- G_init

  ## Precompute RHS: -f_minus[, state_idx] per regime (fixed across iterations)
  rhs_list <- lapply(sys_list, function(sys_s)
    -sys_s$f_minus[, state_idx, drop = FALSE])

  n_state <- length(state_idx)
  final_change <- Inf

  for (k in seq_len(max_iter)) {
    G_new_list <- vector("list", h)

    for (s in seq_len(h)) {
      ## Build Phi_s embed (n_endo x n_endo): columns at non-state positions are 0;
      ## at state_idx positions: Phi_s_embed[, state_idx] = sum_{s'} P[s,s'] * G_{s'}
      ## where G_{s'} is n_endo x n_state.
      ## So Phi_s_embed has n_endo rows and n_endo cols, nonzero only in state cols.
      Phi_s_embed <- .ms_build_Phi_embed(G_list, P, s, state_idx, n_endo)

      sys_s    <- sys_list[[s]]
      f_zero_s <- sys_s$f_zero
      f_plus_s <- sys_s$f_plus

      ## LHS = f_zero^s + f_plus^s * Phi_s_embed  (n_endo x n_endo)
      M_s <- f_zero_s + f_plus_s %*% Phi_s_embed

      ## Solve: M_s * G_s_new = -f_minus[, state_idx]  (n_endo x n_state)
      G_s_new <- tryCatch(
        solve(M_s, rhs_list[[s]]),
        error = function(e) NULL
      )

      if (is.null(G_s_new) || !all(is.finite(G_s_new))) {
        ## Non-finite or singular: iteration diverged
        return(list(G_list = G_list, converged = FALSE,
                    n_iter = k, final_change = Inf))
      }

      G_new_list[[s]] <- G_s_new
    }

    ## Check convergence: max over regimes of max|G_new - G_old|
    changes <- vapply(seq_len(h), function(s)
      max(abs(G_new_list[[s]] - G_list[[s]])), numeric(1))
    final_change <- max(changes)

    G_list <- G_new_list

    if (verbose && (k %% 50L == 0L || k <= 5L))
      cat(sprintf("    iter %d: max_change = %.3e\n", k, final_change))

    if (final_change < tol) {
      return(list(G_list = G_list, converged = TRUE,
                  n_iter = k, final_change = final_change))
    }

    ## Divergence guard: if change grows beyond 1e8, bail early
    if (final_change > 1e8) {
      return(list(G_list = G_list, converged = FALSE,
                  n_iter = k, final_change = final_change))
    }
  }

  ## max_iter exhausted without convergence
  list(G_list = G_list, converged = FALSE,
       n_iter = max_iter, final_change = final_change)
}


## ============================================================================
## Internal helpers
## ============================================================================

## Build the Phi_s embed matrix (n_endo x n_endo).
## Only the columns corresponding to state_idx are non-zero.
## Phi_s_embed[, state_idx] = sum_{s'} P[s,s'] * G_{s'}
## where G_{s'} has shape n_endo x n_state (ghx for regime s').
## @noRd
.ms_build_Phi_embed <- function(G_list, P, s, state_idx, n_endo) {
  h <- length(G_list)
  n_state <- length(state_idx)
  Phi_state <- matrix(0, n_endo, n_state)
  for (sp in seq_len(h)) {
    Phi_state <- Phi_state + P[s, sp] * G_list[[sp]]
  }
  ## Embed into n_endo x n_endo: column j of Phi_embed = Phi_state[, j] at state_idx[j]
  Phi_embed <- matrix(0, n_endo, n_endo)
  Phi_embed[, state_idx] <- Phi_state
  Phi_embed
}

## Check mean-square stability — returns a list with both the per-regime BK
## sufficient condition and the true Maih (2015) Proposition 1 Kronecker
## companion spectral radius.
##
## IMPORTANT: the CORRECT Maih (2015) MSS criterion is:
##   Build the (h * n_state^2) x (h * n_state^2) block matrix M where
##   block M[s, s'] = P[s,s'] * kron(Gs, Gs)   (n_state^2 x n_state^2)
##   MSS holds iff spectral_radius(M) < 1.
##
## Note on the "simple" formula sum P[s,s'] (G_{s'} kron G_{s'}):
##   That formula sums over the FROM regime s' for a fixed TO regime — it gives
##   the row-sum of the Kronecker block matrix. The spectral radius of the full
##   block matrix need NOT equal the max row-norm (only for non-negative matrices
##   via Perron-Frobenius). The full eigen is required for correctness.
##
## Per-regime BK (max_s sr(G_s) < 1) is sufficient but not necessary for MSS.
## A model where per-regime BK fails can still be MS-stable if the weighted
## Kronecker companion sr < 1 (the transition dynamics dampen the instability).
## Conversely, a model where all per-regime BK pass can fail MSS if the
## Kronecker coupling amplifies shocks across regimes.
##
## @noRd
.ms_check_mss <- function(G_list, P, state_idx, h) {
  n_state <- length(state_idx)
  if (n_state == 0L) return(TRUE)

  ## Per-regime spectral radius (per-regime BK sufficient condition)
  sr_per_regime <- vapply(seq_len(h), function(s) {
    Gs <- G_list[[s]][state_idx, , drop = FALSE]
    tryCatch(
      max(Mod(eigen(Gs, only.values = TRUE)$values)),
      error = function(e) NA_real_)
  }, numeric(1))

  ## True Kronecker companion check (Maih 2015, Prop. 1)
  ## Block matrix M of size (h * n_state^2) x (h * n_state^2)
  ## M[(s-1)*n2 + (1:n2), (sp-1)*n2 + (1:n2)] = P[s,sp] * kron(Gs, Gs)
  n2 <- n_state^2L
  dim_M <- h * n2
  M_kron <- matrix(0, dim_M, dim_M)
  for (s in seq_len(h)) {
    Gs <- G_list[[s]][state_idx, , drop = FALSE]
    KG <- kronecker(Gs, Gs)   # n_state^2 x n_state^2
    row_idx <- (s - 1L) * n2 + seq_len(n2)
    for (sp in seq_len(h)) {
      col_idx <- (sp - 1L) * n2 + seq_len(n2)
      M_kron[row_idx, col_idx] <- P[s, sp] * KG
    }
  }
  sr_kron <- tryCatch(
    max(Mod(eigen(M_kron, only.values = TRUE)$values)),
    error = function(e) NA_real_)

  ## mss_ok: TRUE if either the Kronecker sr < 1 (authoritative) or all
  ## per-regime BK pass (sufficient). Use Kronecker as authoritative gate.
  kron_ok <- is.finite(sr_kron) && sr_kron < 1 + 1e-6
  bk_ok   <- all(is.finite(sr_per_regime)) && max(sr_per_regime) < 1 + 1e-6

  ## Return the per-regime check (backward-compatible: TRUE/FALSE scalar) but
  ## attach diagnostics as attributes for tests and ms_foc_residual.
  ## The authoritative criterion is the Kronecker one; we fall back to BK if
  ## the Kronecker computation fails (non-finite eigenvalues).
  result <- if (is.finite(sr_kron)) kron_ok else bk_ok
  attr(result, "sr_kron")        <- sr_kron
  attr(result, "sr_per_regime")  <- sr_per_regime
  attr(result, "bk_ok")          <- bk_ok
  result
}

#' Mean-square-stability diagnostics for a Markov-switching solution
#'
#' Extracts the mean-square-stability (MSS) diagnostics computed during
#' \code{\link{solve_ms_perturbation}} from an \code{MsDecisionRules} object.
#' The authoritative criterion is the spectral radius of the Maih (2015)
#' Kronecker companion block; the per-regime Blanchard-Kahn radii are a
#' weaker sufficient condition.
#'
#' @param ms_dr An \code{MsDecisionRules} object (from
#'   \code{\link{solve_ms_perturbation}}).
#' @return A named list with elements \code{sr_kron} (spectral radius of the
#'   Maih 2015 Kronecker companion block), \code{sr_per_regime} (per-regime
#'   spectral radii, length \code{h}), \code{mss_ok} (\code{TRUE} if
#'   \code{sr_kron < 1}, authoritative), and \code{bk_ok} (\code{TRUE} if all
#'   per-regime radii are below 1, sufficient condition).
#' @seealso \code{\link{solve_ms_perturbation}}
#' @export
ms_mss_diagnostics <- function(ms_dr) {
  if (!inherits(ms_dr, "MsDecisionRules"))
    stop("ms_mss_diagnostics: ms_dr must be an MsDecisionRules object.", call. = FALSE)
  h       <- length(ms_dr$dr)
  state_idx <- ms_dr$dr[[1L]]$state_idx
  G_list  <- lapply(ms_dr$dr, function(dr) dr$ghx)
  flag    <- .ms_check_mss(G_list, ms_dr$P, state_idx, h)
  list(
    sr_kron        = attr(flag, "sr_kron"),
    sr_per_regime  = attr(flag, "sr_per_regime"),
    mss_ok         = isTRUE(flag),
    bk_ok          = isTRUE(attr(flag, "bk_ok"))
  )
}

## Compute FOC residual for each regime (for diagnostics / oracle check).
## Returns max-regime Frobenius-norm of the FOC residual matrix.
## FOC residual for regime s: f_plus^s * Phi_s_embed * G_s + f_zero^s * G_s + f_minus_state^s
## @noRd
.ms_foc_residual <- function(dr_list, sys_list, P, state_idx, h) {
  n_endo    <- nrow(sys_list[[1L]]$f_zero)
  G_list    <- lapply(dr_list, function(dr) dr$ghx)   # each n_endo x n_state
  max_resid <- 0
  for (s in seq_len(h)) {
    Phi_s_embed <- .ms_build_Phi_embed(G_list, P, s, state_idx, n_endo)
    sys_s     <- sys_list[[s]]
    G_s       <- G_list[[s]]
    rhs_state <- sys_s$f_minus[, state_idx, drop = FALSE]
    resid_s   <- sys_s$f_plus %*% Phi_s_embed %*% G_s +
                 sys_s$f_zero %*% G_s +
                 rhs_state
    max_resid <- max(max_resid, norm(resid_s, "F"))
  }
  max_resid
}

## Build the MsDecisionRules return object.
## @noRd
.build_ms_dr <- function(dr_list, P, ss_list, converged, n_iter, mss_ok) {
  structure(
    list(
      dr        = dr_list,
      P         = P,
      pi0       = .ms_ergodic_dist(P),
      ss_list   = ss_list,
      converged = converged,
      n_iter    = n_iter,
      mss_ok    = mss_ok
    ),
    class = c("MsDecisionRules", "list")
  )
}


#' @export
#' @noRd
print.MsDecisionRules <- function(x, ...) {
  cat(sprintf("<MsDecisionRules>  %d regimes, converged=%s, n_iter=%d, mss_ok=%s\n",
              length(x$dr), x$converged, x$n_iter, x$mss_ok))
  cat("  Ergodic dist:", paste(round(x$pi0, 4), collapse = ", "), "\n")
  for (s in seq_along(x$dr)) {
    dr_s <- x$dr[[s]]
    cat(sprintf("  Regime %d: BK=%s, mss=%s, n_state=%d\n",
                s, isTRUE(dr_s$bk_satisfied), isTRUE(dr_s$mss_satisfied),
                dr_s$n_state))
  }
  invisible(x)
}


## ============================================================================
## FOC residual diagnostic (exported for testing)
## ============================================================================

#' Compute per-regime FOC residual for a converged MsDecisionRules object
#'
#' Returns the maximum Frobenius norm of the FOC residual across all regimes.
#' At convergence this should be near zero (< 1e-8 for well-conditioned problems).
#'
#' @param ms_dr   Output of \code{\link{solve_ms_perturbation}}.
#' @param compiled  dynhr_compiled (same as used in solve_ms_perturbation).
#' @param ss_by_regime  List of per-regime steady states.
#' @param params_by_regime  List of per-regime parameter vectors.
#' @return Scalar: max over regimes of Frobenius norm of the FOC residual matrix.
#' @export
ms_foc_residual <- function(ms_dr, compiled, ss_by_regime, params_by_regime) {
  if (!inherits(ms_dr, "MsDecisionRules"))
    stop("ms_foc_residual: ms_dr must be an MsDecisionRules object.", call. = FALSE)
  h         <- length(ms_dr$dr)
  sys_cache <- cache_system_structure(compiled)
  sys_list  <- lapply(seq_len(h), function(s)
    extract_system_matrices_fast(sys_cache, ss_by_regime[[s]], params_by_regime[[s]]))
  state_idx <- ms_dr$dr[[1L]]$state_idx
  n_endo    <- sys_cache$n_endo
  .ms_foc_residual(ms_dr$dr, sys_list, ms_dr$P, state_idx, h)
}
