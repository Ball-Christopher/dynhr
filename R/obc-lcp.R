## R/obc-lcp.R
## --------------------------------------------------------------------------
## LCP (Linear Complementarity Problem) OBC simulation solver.
##
## Provides:
##   lcp_comp_residual()  -- complementarity mismatch residual of a regime path
##                           on its constrained path (drives the Newton method)
##   lcp_jacobian_col()   -- finite-difference column of the Newton Jacobian J
##   lcp_build_jacobian() -- assemble the T*n_spec × T*n_spec Jacobian
##   lcp_build_system()   -- the stacked LCP(q, M) of the piecewise-linear path
##   lcp_lemke()          -- Lemke's complementarity pivoting
##   solve_obc_lcp()      -- Lemke (default) or Newton regime search
##
## THE STACKED LCP (Holden 2016 / DynareOBC; fixed 2026-09-25, W48)
##
##   Keep every tagged equation and add a multiplier z_jt >= 0 to it,
##   F_j(t) = sign_j * z_jt, known from period 1 (a "news" shock).  The model
##   stays linear, so the constrained variables are linear in z:
##
##     w = sign * (x - bound) = q + M z,   z >= 0,  w >= 0,  z'w = 0,
##
##   q from the all-slack path and M[:, (j,t)] the response to z_jt = 1.
##   M is NOT lower-triangular: a multiplier in period t moves every earlier
##   period through expectations.  (Before W48 M was built by simulating
##   one-period binding policies that assumed the next period slack, which
##   is block-lower-triangular but not the model: a spell of 2+ periods was
##   solved with the wrong expectation.)  Its solution is the exact
##   piecewise-linear OccBin path; the binding set is {z > 0}.
##
## RESIDUAL OF THE NEWTON METHOD
##
##   For each (j, t) on the constrained path of the current regime path:
##     SLACK period:   R_jt = max(0, -sign_j (x_j(t) - bound_j))
##                     (> 0 when the bound is violated -> should bind)
##     BINDING period: R_jt = max(0, -sign_j F_j(t))
##                     (> 0 when the multiplier has the wrong sign -> release)
##   R = 0 iff the regime path is an equilibrium.  (The old residual used the
##   SLACK-policy prediction, so a bound violated only because another
##   constraint binds was never flagged.)
##
## Shock sequences are surprises: each re-solve (at period 1 and at every
## period with a non-zero shock) is one stacked LCP, as in
## boehl_solve_regime_path().
##
## References: Holden T. (2016), "Computation of solutions to dynamic models
##   with occasionally binding constraints"; Holden T. & Paetz M. (2012),
##   Norges Bank WP 2013/18.
## --------------------------------------------------------------------------


# =============================================================================
# Complementarity mismatch residual
# =============================================================================

#' Complementarity mismatch residual for an OBC regime path
#'
#' Evaluated on the constrained path \code{sim} of \code{regime_path}:
#'   - Slack period: positive when bound j is violated on the path.
#'   - Binding period: positive when the multiplier of the tagged equation
#'     has the wrong sign (the constraint should be released).
#'
#' @param sim          List from boehl_simulate() under regime_path
#' @param shock_seq    n_exo x T shock matrix
#' @param regime_cache Policy cache seeded by obc_ensure_policy()
#' @param regime_path  Integer vector (length T): current bitfield regime
#' @param state_init   Initial state (NULL = zero)
#' @return n_spec x T numeric matrix (>= 0); zero iff the regime path is an
#'         equilibrium
#' @noRd
lcp_comp_residual <- function(sim, shock_seq, regime_cache, regime_path,
                              state_init = NULL) {
  .lcp_comp_eval(sim, shock_seq, regime_cache, regime_path, state_init)$res
}

#' Complementarity mismatch and the scale-free equilibrium test
#'
#' $res is lcp_comp_residual()'s matrix (in the units of the model).
#' $consistent is TRUE when the complementarity check of .obc_pwl_check()
#' reproduces regime_path, i.e. every slack period is inside its bound and
#' every binding period has a multiplier of the right sign up to the
#' round-off band rtol RELATIVE to |x| + |b| and to the tagged equation's
#' term magnitudes (.obc_gap_binds() / .obc_mult_keeps()).  This replaces an
#' absolute test on the norm of $res (W79): a mismatch norm below 1e-8 in
#' level units passed at every size of the violation relative to the path.
#' @noRd
.lcp_comp_eval <- function(sim, shock_seq, regime_cache, regime_path,
                           state_init = NULL, rtol = 1e-8) {
  ctx   <- .obc_pwl_cache_context(regime_cache)
  rules <- .obc_pwl_rules(ctx, regime_path)
  chk   <- .obc_pwl_check(ctx, rules, sim, shock_seq, regime_path, state_init,
                          rtol)
  B     <- t(.obc_bitfield_to_binding(regime_path, ctx$n_spec))
  res   <- ifelse(B, pmax(0, -chk$mult), pmax(0, -chk$gap))
  list(res = matrix(res, ctx$n_spec, ncol(shock_seq)),
       consistent = identical(chk$regime_path, as.integer(regime_path)))
}


# =============================================================================
# Jacobian column: finite-difference sensitivity along one regime flip
# =============================================================================

#' Finite-difference Jacobian column for the complementarity mismatch residual
#'
#' Flips spec j_src at period t_src, re-simulates, and returns the change in
#' the stacked T*n_spec complementarity mismatch vector.  This is column
#' (j_src, t_src) of the Newton Jacobian J.  A flip in period t_src changes
#' the rules of every period up to t_src (expectations), so the column is
#' generally non-zero in earlier rows too.
#'
#' Column ordering: col = (t - 1)*n_spec + j (period-major, consistent with
#' as.vector(matrix[n_spec, T])).
#'
#' @param j_src        Integer: spec index to flip (1-based)
#' @param t_src        Integer: period index to flip (1-based)
#' @param res_vec      Numeric vector (length T*n_spec): current complementarity
#'                     mismatch from as.vector(lcp_comp_residual(...))
#' @param shock_seq    n_exo x T shock matrix
#' @param dr_slack     Slack-regime DecisionRules
#' @param regime_cache Policy cache environment (seeded by obc_ensure_policy)
#' @param sys          System matrices (unused; the cache carries them)
#' @param regime_path  Integer vector (length T): current regime
#' @param specs        OBC spec list
#' @param obs_idx      Integer vector of observable positions (unused)
#' @param state_init   Numeric state vector or NULL (zero initial state)
#' @return Numeric vector (length T*n_spec): R_flipped − R_current
#' @noRd
lcp_jacobian_col <- function(j_src, t_src, res_vec,
                              shock_seq, dr_slack, regime_cache,
                              sys, regime_path, specs, obs_idx, state_init) {
  n_spec <- length(specs)

  flags <- obc_regime_flags(regime_path[t_src], n_spec)
  flags[j_src] <- !flags[j_src]
  rp_mod        <- regime_path
  rp_mod[t_src] <- obc_regime_idx(flags)

  sim_mod <- boehl_simulate(shock_seq, dr_slack, regime_cache, rp_mod, state_init)
  res_mod <- lcp_comp_residual(sim_mod, shock_seq, regime_cache, rp_mod,
                               state_init)

  as.vector(res_mod) - res_vec
}


# =============================================================================
# Jacobian assembly
# =============================================================================

#' Build the T*n_spec × T*n_spec Newton Jacobian for the stacked system
#'
#' Assembles J by calling lcp_jacobian_col() for every (j, t) pair.
#' Column index: col = (t - 1)*n_spec + j.
#'
#' @param res_vec      Numeric vector (length T*n_spec): current residuals
#' @param shock_seq    n_exo x T shock matrix
#' @param dr_slack     Slack-regime DecisionRules
#' @param regime_cache Policy cache environment
#' @param sys          System matrices
#' @param regime_path  Integer vector (length T)
#' @param specs        OBC spec list
#' @param obs_idx      Integer vector of observable positions
#' @param state_init   Numeric state vector or NULL
#' @return Numeric matrix (T*n_spec × T*n_spec) — the Newton Jacobian J
#' @noRd
lcp_build_jacobian <- function(res_vec, shock_seq, dr_slack, regime_cache,
                                sys, regime_path, specs, obs_idx, state_init) {
  n_T    <- ncol(shock_seq)
  n_spec <- length(specs)
  n_lcp  <- n_T * n_spec

  J <- matrix(0, n_lcp, n_lcp)

  for (t in seq_len(n_T)) {
    for (j in seq_len(n_spec)) {
      col <- (t - 1L) * n_spec + j
      J[, col] <- lcp_jacobian_col(j, t, res_vec,
                                    shock_seq, dr_slack, regime_cache,
                                    sys, regime_path, specs, obs_idx, state_init)
    }
  }
  J
}


# =============================================================================
# Stacked LCP system construction (for Lemke)
# =============================================================================

#' Responses of the constrained variables to anticipated multipliers
#'
#' Column (t-1)*n_spec + j is the response of sign_j' x_j'(t') (rows, same
#' ordering) to a unit multiplier z_jt in the tagged equation of spec j,
#' F_j(t) = sign_j z_jt, known from period 1, under the slack rule (zero
#' initial state).  With A = F0 + Fp ghx_slack E_s the multiplier enters the
#' constant of the rules c_t = A^{-1}(sign_j e_j) and c_u = -A^{-1} Fp c_{u+1}
#' for u < t (the backward recursion of .obc_pwl_rules() with no regime
#' change).
#' @noRd
.lcp_multiplier_matrix <- function(ctx, n_T) {
  n_spec <- ctx$n_spec
  n_lcp  <- n_T * n_spec
  M <- matrix(0, n_lcp, n_lcp)
  if (n_spec == 0L || n_T == 0L) return(M)

  A <- ctx$F0
  A[, ctx$si] <- A[, ctx$si] + ctx$Fp %*% ctx$ghx
  B <- -solve(A, ctx$Fp)
  G <- ctx$ghx

  for (j in seq_len(n_spec)) {
    e <- numeric(ctx$n_endo)
    e[ctx$eq[j]] <- ctx$sgn[j]
    ## d[, m + 1] = B^m A^{-1} sign_j e_j: the constant m periods before z
    d <- matrix(0, ctx$n_endo, n_T)
    d[, 1L] <- solve(A, e)
    for (m in seq_len(n_T - 1L)) d[, m + 1L] <- drop(B %*% d[, m])

    for (t_src in seq_len(n_T)) {
      col <- (t_src - 1L) * n_spec + j
      s <- numeric(ctx$n_state)
      for (tau in seq_len(n_T)) {
        y <- drop(G %*% s)
        if (tau <= t_src) y <- y + d[, t_src - tau + 1L]
        M[(tau - 1L) * n_spec + seq_len(n_spec), col] <- ctx$sgn * y[ctx$var]
        s <- y[ctx$si]
      }
    }
  }
  M
}

#' Build the stacked OccBin LCP system (q, M) from a shock sequence
#'
#' Constructs the T*n_spec dimensional LCP(q, M)
#'
#'   find z >= 0, w >= 0 with w = M z + q and z'w = 0
#'
#' whose solution is the piecewise-linear OBC path of the shock sequence
#' (surprise shocks, regime path known from period 1; exact for a single
#' shock in period 1):
#'   z_k > 0  ↔  constraint k = (spec j, period t) binds (z = its multiplier)
#'   w_k > 0  ↔  constraint k is slack (w = sign_j (x_j(t) - bound_j))
#'
#' q is the all-slack path: q_k = sign_j * (x_j^slack(t) - bound_j).
#' M is the response of w to the multipliers (.lcp_multiplier_matrix()); it
#' is not lower-triangular in t (expectations).
#'
#' @param shock_seq    n_exo x T shock matrix
#' @param dr_slack     Slack-regime DecisionRules
#' @param sys          System matrices
#' @param specs        OBC spec list
#' @param obs_idx      Integer vector of observable positions (unused)
#' @param state_init   Numeric state vector or NULL (zero initial state)
#' @return Named list with:
#'   $q    -- numeric vector (length T*n_spec): all-slack complementarity gaps
#'   $M    -- numeric matrix (T*n_spec × T*n_spec): multiplier responses
#'   $sim0 -- all-slack simulation (list with $paths and $states)
#' @noRd
lcp_build_system <- function(shock_seq, dr_slack, sys, specs, obs_idx,
                              state_init = NULL) {
  n_T <- ncol(shock_seq)
  ctx <- .obc_pwl_context(sys, dr_slack, specs)

  sim0 <- .obc_pwl_forward(ctx, vector("list", n_T), shock_seq, state_init)
  q <- as.vector(ctx$sgn * (sim0$paths[ctx$var, , drop = FALSE] - ctx$bnd))

  list(q = q, M = .lcp_multiplier_matrix(ctx, n_T), sim0 = sim0)
}


# =============================================================================
# Lemke's complementarity pivoting algorithm
# =============================================================================

#' Lemke's principal pivoting algorithm for LCP(q, M)
#'
#' Solves the linear complementarity problem: find z ≥ 0, w ≥ 0 such that
#'   w = M z + q    and    z' w = 0
#' using Lemke's complementarity pivoting method with an artificial variable z_0.
#'
#' Algorithm outline:
#'   1. If q ≥ 0: trivial solution z = 0, w = q.
#'   2. Introduce z_0 (covering variable) with column d (default d = 1_n).
#'      Drive z_0 into the basis at the row with the most negative q entry.
#'   3. Pivot: at each step bring the complement of the leaving variable into
#'      the basis (complementary pivoting).  Stop when z_0 leaves the basis.
#'
#' Lemke terminates with a solution when M is a P-matrix (every principal
#' minor positive), which is Holden's (2016) condition for the OccBin
#' stacked LCP to have a unique solution for every q.  M is not triangular
#' (a multiplier moves earlier periods through expectations), so the pivot
#' count is not bounded by n; a ray termination is reported as
#' converged = FALSE.
#'
#' Scale-free (W79, 2026-09): the pivots run on the equivalent LCP
#' (q / s_q, M / s_M) with s_q, s_M the powers of two nearest max|q| and
#' max|M| (z = z_n s_q / s_M, w = w_n s_q), so every tolerance below is
#' RELATIVE to the problem's scale.  Scaling by a power of two is exact in
#' floating point and commutes with the Gauss-Jordan pivots: away from a
#' tolerance decision the solution is bit-identical to pivoting on (q, M)
#' itself.  (The tolerances used to be absolute: with the shocks and the
#' bound in units 1e-6 times smaller every q entry passed q >= -1e-8 and
#' the binding LCP came back trivial.)
#'
#' @param q        Numeric vector (length n): constant term of the LCP
#' @param M        Numeric matrix (n × n): LCP coefficient matrix
#' @param d        Covering vector (length n, all positive); NULL → ones
#' @param max_iter Maximum pivot steps; NULL → 50*n
#' @param tol      Relative numerical tolerance (default 1e-10): q is
#'   trivial when q >= -tol max|q| (to the nearest power of two), and the
#'   pivot and ratio tests use tol on the normalised tableau
#' @return Named list:
#'   $z         -- numeric vector (length n): primal solution (z ≥ 0)
#'   $w         -- numeric vector (length n): dual solution (w = Mz + q ≥ 0)
#'   $converged -- logical
#'   $n_iter    -- integer: pivot steps taken
#'   $message   -- character: termination reason
#'   $z_scale   -- the natural scale of z, s_q / s_M (a power of two): a
#'                 multiplier is positive when z > tol * z_scale
#' @noRd
lcp_lemke <- function(q, M, d = NULL, max_iter = NULL, tol = 1e-10) {
  n <- length(q)
  if (is.null(d))        d        <- rep(1.0, n)
  if (is.null(max_iter)) max_iter <- 50L * n

  stopifnot(nrow(M) == n, ncol(M) == n, length(d) == n, all(d > 0))

  # Normalise by powers of two (exact): w_n = q_n + M_n z_n + d z0_n with
  # q_n = q / s_q, M_n = M / s_M, z_n = z s_M / s_q, z0_n = z0 / s_q.
  s_q <- .lcp_pow2(max(abs(q)))
  s_M <- .lcp_pow2(max(abs(M)))
  z_scale <- s_q / s_M
  q_in <- q
  M_in <- M
  q <- q / s_q
  M <- M / s_M

  # Trivial solution: q >= 0 already (relative to max|q|).
  if (all(q >= -tol))
    return(list(z = numeric(n), w = pmax(q_in, 0),
                converged = TRUE, n_iter = 0L, message = "trivial",
                z_scale = z_scale))

  # --- Variable encoding (column indices in the tableau) ---
  # Cols  1..n   : z_1 .. z_n     (LCP decision variables)
  # Col   n+1    : z_0            (artificial covering variable)
  # Cols  n+2..2n+1: w_1 .. w_n  (LCP slack variables)
  # Col   2n+2   : RHS
  #
  # Complementary pairs: (z_k, w_k) where z_k = col k, w_k = col n+1+k.
  # comp(z_k)  = col n+1+k    (k = 1..n)
  # comp(w_k)  = col k        (w_k is at col n+1+k)
  # z_0 has no complement.

  comp_col <- function(col) {
    if (col >= 1L && col <= n)           return(n + 1L + col)   # z_k -> w_k
    if (col >= n + 2L && col <= 2L*n+1L) return(col - n - 1L)  # w_k -> z_k
    stop("z_0 (col n+1) has no complement")
  }

  # Initial full tableau: n rows × (2n+2) columns.
  # Convention A: basic_i = RHS_i - sum_j tab[i,j]*nonbasic_j.
  # The system is: w - Mz - dz_0 = q, so nonbasic columns store -M and -d.
  # Row i (basic = w_i): w_i = q_i - (-M[i,:])*z - (-d_i)*z_0
  #   => stored coefficients are -M (for z cols) and -d (for z_0 col).
  # After pivoting z_0 into row t_0 (piv_elem = -d[t_0] < 0):
  #   RHS_new[t_0] = q[t_0] / (-d[t_0]) = -q[t_0]/d[t_0] > 0 ✓
  tab <- cbind(-M, -d, diag(n), q)   # n × (2n+2)
  bas <- (n + 2L):(2L*n + 1L)      # initial basic variables: w_1..w_n

  # Helper: Gauss-Jordan pivot on (row r, column c).
  pivot <- function(tab, r, c) {
    piv <- tab[r, c]
    if (abs(piv) < tol * 1e-4)
      stop(sprintf("degenerate pivot element %.2e at row %d col %d", piv, r, c))
    tab[r, ]  <- tab[r, ] / piv
    for (i in seq_len(nrow(tab)))
      if (i != r) tab[i, ] <- tab[i, ] - tab[i, c] * tab[r, ]
    tab
  }

  # Step 1: drive z_0 into the basis at the most-negative-q row.
  t0 <- which.min(q)                  # row index
  tab <- pivot(tab, t0, n + 1L)       # z_0 (col n+1) enters at row t0
  leaving_var  <- bas[t0]             # w_{t0} leaves
  bas[t0]      <- n + 1L              # z_0 is now basic at row t0
  enter_col    <- comp_col(leaving_var) # complement of w_{t0} = z_{t0} (col t0)

  for (iter in seq_len(max_iter)) {
    # Minimum-ratio test for the entering column.
    col_vals <- tab[, enter_col]
    rhs_vals <- tab[, 2L * n + 2L]

    pos_rows <- which(col_vals > tol)
    if (length(pos_rows) == 0L)
      return(list(z = numeric(n), w = numeric(n),
                  converged = FALSE, n_iter = iter,
                  message = "secondary ray: problem is unbounded",
                  z_scale = z_scale))

    ratios    <- rhs_vals[pos_rows] / col_vals[pos_rows]
    r         <- pos_rows[which.min(ratios)]
    leaving_var <- bas[r]

    tab    <- pivot(tab, r, enter_col)
    bas[r] <- enter_col

    # z_0 left the basis → solution found.
    if (leaving_var == n + 1L) {
      z <- numeric(n)
      w <- numeric(n)
      for (i in seq_len(n)) {
        bv  <- bas[i]
        val <- tab[i, 2L * n + 2L]
        if (bv >= 1L  && bv <= n)          z[bv]           <- max(0, val)
        if (bv >= n+2L && bv <= 2L*n+1L)  w[bv - n - 1L]  <- max(0, val)
      }
      # Back to the caller's units; recompute w for numerical accuracy.
      z <- z * z_scale
      w <- pmax(drop(M_in %*% z) + q_in, 0)
      return(list(z = z, w = w,
                  converged = TRUE, n_iter = iter,
                  message = "z_0 left basis", z_scale = z_scale))
    }

    enter_col <- comp_col(leaving_var)
  }

  # Extract best approximation at max_iter.
  z <- numeric(n)
  for (i in seq_len(n)) {
    bv <- bas[i]
    if (bv >= 1L && bv <= n) z[bv] <- max(0, tab[i, 2L*n+2L])
  }
  z <- z * z_scale
  list(z = z, w = pmax(drop(M_in %*% z) + q_in, 0),
       converged = FALSE, n_iter = max_iter, message = "max iterations reached",
       z_scale = z_scale)
}

#' Power of two nearest a positive scale (1 for zero / non-finite)
#'
#' Dividing by a power of two is exact in floating point, so normalising a
#' problem by .lcp_pow2(scale) changes no digit of the arithmetic -- only
#' which side of a tolerance a value falls.
#' @noRd
.lcp_pow2 <- function(s) {
  if (!is.finite(s) || s <= 0) return(1)
  2^round(log2(s))
}




# =============================================================================
# Lemke-based OBC regime solver
# =============================================================================

#' Extract an OBC regime path from a Lemke LCP solution
#'
#' Given the z solution of LCP(q, M) built by lcp_build_system(), maps z back
#' to a binary integer regime path: z_k > tol * max|z| means constraint k
#' binds (relative, so the classification does not depend on the units of
#' the model; an absolute threshold before W79).
#' Then forward-simulates the model under that regime path.
#'
#' @param z        Numeric vector (length T*n_spec): Lemke primal solution
#' @param shock_seq n_exo x T shock matrix
#' @param dr_slack Slack-regime DecisionRules
#' @param sys      System matrices
#' @param specs    OBC spec list
#' @param obs_idx  Integer vector of observable positions
#' @param state_init Numeric state vector or NULL
#' @param tol      Relative threshold for treating z_k as positive
#'   (default 1e-6, times max|z|)
#' @return Named list with $regime_path, $paths, $states
#' @noRd
lcp_regime_from_z <- function(z, shock_seq, dr_slack, sys, specs, obs_idx,
                               state_init = NULL, tol = 1e-6) {
  n_T    <- ncol(shock_seq)
  n_spec <- length(specs)

  rc <- new.env(parent = emptyenv(), hash = TRUE)
  obc_ensure_policy(0L, rc, sys, dr_slack, specs, obs_idx)

  Z <- matrix(z > tol * max(abs(z)), n_spec, n_T)
  regime_path <- .obc_binding_to_bitfield(t(Z))

  sim <- boehl_simulate(shock_seq, dr_slack, rc, regime_path, state_init)
  list(regime_path = regime_path, paths = sim$paths, states = sim$states)
}


# =============================================================================
# Main solver
# =============================================================================

#' LCP OBC simulation solver
#'
#' Finds the binding/slack regime path and the piecewise-linear path for a
#' deterministic shock sequence using one of two methods selected by the
#' \code{method} argument.  Both return the same equilibrium as
#' \code{\link{boehl_solve_regime_path}}: time-varying OccBin rules along the
#' regime path, every slack period inside every bound and every binding
#' period with a multiplier of the right sign.  Shocks are surprises: at
#' period 1 and at every period with a non-zero shock the expected regime
#' path is re-solved from the current state.
#'
#' \strong{method = "lemke"} (default): Lemke's complementarity pivoting on
#'   the stacked LCP(q, M) of \code{lcp_build_system()} -- every tagged
#'   equation carries a non-negative multiplier known from the start
#'   (Holden 2016), which keeps the model linear in the multipliers.
#'
#' \strong{method = "newton"}: Newton iteration on the stacked complementarity
#'   mismatch residual.  At each step the T*n_spec × T*n_spec finite-difference
#'   Jacobian of regime flips is assembled (T*n_spec re-simulations), and a QR
#'   Newton step gives a regime update.
#'
#' @param shock_seq        n_exo x T numeric matrix of structural shocks
#' @param dr_slack         Slack-regime DecisionRules (from solve_perturbation)
#' @param sys              System matrices (from extract_system_matrices_fast)
#' @param specs            OBC spec list (from obc_parse_tags / obc_collect_specs)
#' @param obs_idx          Integer vector: observable positions in endo vector;
#'                         NULL → all endo indices
#' @param method           Character: \code{"lemke"} (default) or \code{"newton"}
#' @param max_iter         Maximum pivot steps (lemke) or Newton iterations
#'                         (newton) per re-solve; NULL → solver-specific
#'                         defaults
#' @param tol              Relative complementarity tolerance (default 1e-8):
#'                         the Lemke pivots and multiplier classification
#'                         are relative to max|q| and max|M| of the stacked
#'                         LCP, and a regime path is an equilibrium when every
#'                         bound holds up to \code{tol} times \eqn{|x| + |b|}
#'                         and every multiplier has the right sign up to
#'                         \code{tol} times its equation's term magnitudes
#'                         (\code{max(tol, 1e-6)} for the \code{converged}
#'                         flag), so the solution does not depend on the
#'                         units of the model (absolute before W79)
#' @param regime_path_init Integer vector (length T) for warm-starting the
#'                         Newton solver; ignored for method = "lemke"
#' @param state_init       Numeric vector (length n_state) initial model state;
#'                         NULL → zero initial state (standard IRF convention)
#' @return Named list with:
#'   $regime_path -- integer vector (length T): final bitfield regime per period
#'   $paths       -- n_endo x T matrix of simulated endogenous variable paths
#'   $fb_norm     -- scalar: complementarity mismatch norm at the returned path
#'   $converged   -- logical
#'   $n_iter      -- integer: pivot steps (lemke) or Newton iterations (newton),
#'                   the largest over the re-solves
#'
#' @references
#'   Boehl, G. (2020). Efficient solution and computation of models with
#'     occasionally binding constraints. \emph{Deutsche Bundesbank Discussion
#'     Paper}, 38/2020.
#'   Holden, T. D. (2016). Computation of solutions to dynamic models with
#'     occasionally binding constraints. \emph{EconStor Preprints}, 130142.
#'   Cottle, R. W., Pang, J.-S., & Stone, R. E. (2009).
#'     \emph{The Linear Complementarity Problem}. SIAM.
#' @export
solve_obc_lcp <- function(shock_seq, dr_slack, sys, specs,
                           obs_idx          = NULL,
                           method           = c("lemke", "newton"),
                           max_iter         = NULL,
                           tol              = 1e-8,
                           regime_path_init = NULL,
                           state_init       = NULL) {
  method <- match.arg(method)
  n_T    <- ncol(shock_seq)

  if (is.null(obs_idx)) obs_idx <- seq_len(nrow(dr_slack$ghx))

  regime_cache <- new.env(parent = emptyenv(), hash = TRUE)
  obc_ensure_policy(0L, regime_cache, sys, dr_slack, specs, obs_idx)
  ctx <- .obc_pwl_cache_context(regime_cache)

  ## One re-solve: a shock in the first column of `e` only, from state `s`.
  ## Returns the list shape of .obc_pwl_solve() plus the mismatch matrix.
  ## Every test is RELATIVE to the scale of the problem (W79): the
  ## equilibrium check is .lcp_comp_eval()'s round-off band (|x| + |b| and the
  ## tagged equation's term magnitudes), a Lemke multiplier is positive above
  ## tol times its natural scale.  (Absolute before: a mismatch norm below
  ## max(tol, 1e-6) and z > tol in level units.)
  finish <- function(e, s, regime_path, n_iter, ok) {
    sim  <- boehl_simulate(e, dr_slack, regime_cache, regime_path, s)
    comp <- .lcp_comp_eval(sim, e, regime_cache, regime_path, s,
                           rtol = max(tol, 1e-6))
    list(regime_path = regime_path, paths = sim$paths, states = sim$states,
         phi = comp$res, converged = ok && comp$consistent,
         n_iter = n_iter)
  }

  inner_lemke <- function(ctx, e, s, ini) {
    lcp_sys <- lcp_build_system(e, dr_slack, sys, specs, obs_idx, s)
    lemke   <- lcp_lemke(lcp_sys$q, lcp_sys$M,
                         max_iter = if (!is.null(max_iter)) as.integer(max_iter),
                         tol = tol)
    Z  <- matrix(lemke$z > tol * lemke$z_scale, length(specs), n_T)
    finish(e, s, .obc_binding_to_bitfield(t(Z)), lemke$n_iter,
           lemke$converged)
  }

  inner_newton <- function(ctx, e, s, ini) {
    newton_max  <- if (!is.null(max_iter)) as.integer(max_iter) else 20L
    regime_path <- if (!is.null(ini)) as.integer(ini) else integer(n_T)
    for (iter in seq_len(newton_max)) {
      sim      <- boehl_simulate(e, dr_slack, regime_cache, regime_path, s)
      comp     <- .lcp_comp_eval(sim, e, regime_cache, regime_path, s,
                                 rtol = tol)
      comp_vec <- as.vector(comp$res)
      if (comp$consistent)
        return(finish(e, s, regime_path, iter - 1L, TRUE))

      # Finite-difference Jacobian of regime flips, one re-simulation per
      # column; least-squares QR Newton step (a flip that changes nothing
      # gives a zero column: its aliased coefficient is set to 0, no flip).
      J <- lcp_build_jacobian(comp_vec, e, dr_slack, regime_cache, sys,
                              regime_path, specs, obs_idx, s)
      delta <- qr.coef(qr(J, tol = 1e-12), -comp_vec)
      delta[is.na(delta)] <- 0

      # delta_{jt} > 0.5 → flip constraint j in period t
      flip <- matrix(delta > 0.5, length(specs), n_T)
      B    <- t(.obc_bitfield_to_binding(regime_path, length(specs)))
      new_regime <- .obc_binding_to_bitfield(t(xor(B, flip)))
      if (identical(new_regime, regime_path))
        return(finish(e, s, regime_path, iter, FALSE))
      regime_path <- new_regime
    }
    finish(e, s, regime_path, newton_max, FALSE)
  }

  res <- .obc_pwl_solve_surprise(
    ctx, shock_seq, state_init,
    if (method == "newton") regime_path_init,
    if (method == "lemke") inner_lemke else inner_newton)

  list(
    regime_path = res$regime_path,
    paths       = res$paths,
    fb_norm     = sqrt(sum(res$phi^2)),
    converged   = res$converged,
    n_iter      = res$n_iter
  )
}
