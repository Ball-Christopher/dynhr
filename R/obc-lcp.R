## R/obc-lcp.R
## --------------------------------------------------------------------------
## LCP (Linear Complementarity Problem) OBC simulation solver.
##
## Provides:
##   lcp_comp_residual()  -- complementarity mismatch residual (non-zero when
##                           regime is inconsistent with the slack-policy
##                           predictions; drives the Newton iteration)
##   lcp_jacobian_col()   -- finite-difference column of the Newton Jacobian J
##   lcp_build_jacobian() -- assemble the T*n_spec × T*n_spec Jacobian
##   solve_obc_lcp()      -- Newton iteration on the stacked complementarity system
##
## RESIDUAL DEFINITION
##
##   For each (j, t) the complementarity mismatch is:
##
##     For a SLACK period (regime bit j = 0 at t):
##       R_{jt} = max(0, sign_j * (bound_j − x_j^slack(t)))
##       > 0 when the slack policy would violate bound j → regime should be binding.
##
##     For a BINDING period (regime bit j = 1 at t):
##       R_{jt} = max(0, sign_j * (x_j^slack(t) − bound_j))
##       > 0 when the slack policy satisfies bound j → binding is unnecessary.
##
##   R = 0 everywhere iff the regime path is the unique complementarity solution:
##   every slack period has a non-violated slack-policy prediction, and every
##   binding period genuinely needs to bind.
##
##   NOTE: The Boehl FB residual (a+b−sqrt(a²+b²)) is always identically zero
##   for an all-slack simulation — because b_jt = max(0, x_cur−bound) = 0 when
##   x_cur = x_slack < bound — and is therefore not a suitable convergence
##   criterion or Newton objective.  The complementarity mismatch R does not
##   have this degeneracy.
##
## ALGORITHM
##
##   Newton iteration on the stacked T*n_spec complementarity system:
##
##     1. Initialise regime_path (warm-start or all-slack).
##     2. Simulate; compute R and ||R||.
##     3. If ||R|| < tol: CONVERGED.
##     4. Build Jacobian J (T*n_spec × T*n_spec):
##          J[:, (j,t)] = R(regime_flipped_j_at_t) − R(current_regime)
##        Block-lower-triangular in t (causal: flipping period t only affects t'≥t).
##     5. Newton step: solve J * delta = −R via QR (SVD regularised fallback).
##     6. Update regime: for each (j, t), do the flip if delta_{jt} > 0.5.
##     7. Ensure new policies are cached; goto 2.
##
## RELATION TO BOEHL ITERATIVE
##   Boehl iterative is a fixed-point iteration that checks each period's
##   slack-policy prediction independently.  The LCP Newton uses the full
##   T-period Jacobian and simultaneously accounts for all inter-period
##   state propagation, preventing the flip-flop cycling that can occur in
##   multi-spell scenarios.
##
## Reference: Holden T. & Paetz M. (2012), Norges Bank WP 2013/18.
## --------------------------------------------------------------------------


# =============================================================================
# Complementarity mismatch residual
# =============================================================================

#' Complementarity mismatch residual for an OBC regime path
#'
#' For each (j, t) returns a non-negative scalar that is positive iff the
#' current regime assignment is inconsistent with the slack-policy prediction:
#'   - Slack period: positive when the slack policy would violate bound j.
#'   - Binding period: positive when binding is unnecessary (slack would satisfy).
#'
#' This is the proper convergence metric and Newton objective for the LCP
#' solver.  Unlike the Boehl FB residual, it is non-zero at incorrect regimes.
#'
#' @param states       n_state x T matrix: post-period state from boehl_simulate
#' @param shock_seq    n_exo x T shock matrix
#' @param specs        OBC spec list
#' @param dr_slack     Slack-regime DecisionRules
#' @param regime_path  Integer vector (length T): current bitfield regime
#' @return n_spec x T numeric matrix (≥ 0); zero iff regime is complementarity-
#'         consistent with the current state trajectory
#' @noRd
lcp_comp_residual <- function(states, shock_seq, specs, dr_slack, regime_path) {
  n_T     <- ncol(shock_seq)
  n_spec  <- length(specs)
  n_state <- nrow(states)
  res     <- matrix(0, n_spec, n_T)
  s_prev  <- numeric(n_state)

  for (t in seq_len(n_T)) {
    eps   <- shock_seq[, t]
    flags <- obc_regime_flags(regime_path[t], n_spec)

    for (j in seq_len(n_spec)) {
      s <- specs[[j]]

      # Slack-policy prediction of x_j at period t from state s_{t-1}
      x_slack <- sum(dr_slack$ghx[s$var_idx, ] * s_prev) +
                 sum(dr_slack$ghu[s$var_idx, ] * eps)

      # gap > 0  means slack policy VIOLATES the bound (should bind)
      # gap <= 0 means slack policy SATISFIES the bound (can be slack)
      gap <- if (s$op == ">") s$bound - x_slack else x_slack - s$bound

      if (flags[j]) {
        # Currently binding: residual is positive only when binding is unneeded
        res[j, t] <- max(0, -gap)   # positive when gap < 0 (slack policy OK)
      } else {
        # Currently slack: residual is positive only when bound is violated
        res[j, t] <- max(0,  gap)   # positive when gap > 0 (slack policy violates)
      }
    }
    s_prev <- states[, t]
  }
  res
}


# =============================================================================
# Jacobian column: finite-difference sensitivity along one regime flip
# =============================================================================

#' Finite-difference Jacobian column for the complementarity mismatch residual
#'
#' Flips spec j_src at period t_src, re-simulates, and returns the change in
#' the stacked T*n_spec complementarity mismatch vector.  This is column
#' (j_src, t_src) of the Newton Jacobian J.
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
#' @param regime_cache Policy cache environment
#' @param sys          System matrices
#' @param regime_path  Integer vector (length T): current regime
#' @param specs        OBC spec list
#' @param obs_idx      Integer vector of observable positions
#' @param state_init   Numeric state vector or NULL (zero initial state)
#' @return Numeric vector (length T*n_spec): R_flipped − R_current
#' @noRd
lcp_jacobian_col <- function(j_src, t_src, res_vec,
                              shock_seq, dr_slack, regime_cache,
                              sys, regime_path, specs, obs_idx, state_init) {
  n_spec <- length(specs)

  flags_cur <- obc_regime_flags(regime_path[t_src], n_spec)
  flags_flp <- flags_cur
  flags_flp[j_src] <- !flags_cur[j_src]
  r_flipped <- obc_regime_idx(flags_flp)

  obc_ensure_policy(r_flipped, regime_cache, sys, dr_slack, specs, obs_idx)

  rp_mod        <- regime_path
  rp_mod[t_src] <- r_flipped

  sim_mod <- boehl_simulate(shock_seq, dr_slack, regime_cache, rp_mod, state_init)
  res_mod <- lcp_comp_residual(sim_mod$states, shock_seq, specs, dr_slack, rp_mod)

  as.vector(res_mod) - res_vec
}


# =============================================================================
# Jacobian assembly
# =============================================================================

#' Build the T*n_spec × T*n_spec Newton Jacobian for the stacked system
#'
#' Assembles J by calling lcp_jacobian_col() for every (j, t) pair.
#' Column index: col = (t - 1)*n_spec + j.
#' J is block-lower-triangular in t: flipping period t_src only affects
#' the residuals at periods t' ≥ t_src (causal model structure).
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

#' Build the stacked OccBin LCP system (q, M) from a shock sequence
#'
#' Constructs the T*n_spec dimensional LCP(q, M) that characterises the OBC
#' regime path for a given deterministic shock sequence.
#'
#' LCP: find z ≥ 0, w ≥ 0 such that  w = M z + q  and  z' w = 0.
#'
#' Interpretation:
#'   z_k > 0  ↔  constraint k = (spec j, period t) is BINDING (shadow price > 0)
#'   w_k > 0  ↔  constraint k is SLACK (positive gap)
#'
#' q is computed from the all-slack simulation:
#'   q_k = sign_j * (x_j^slack(t) - bound_j)
#'   Positive when the slack policy satisfies the bound; negative when it
#'   would be violated (these periods will be binding in the solution).
#'
#' M[:,col_k] is the finite-difference sensitivity of z to adding binding at k:
#'   M[:,col_k] = z(only_k_binding) - z(all_slack)
#' M is block-lower-triangular in t (causal structure), so Lemke's algorithm
#' for this system reduces to a single-pass forward substitution and terminates
#' in exactly T*n_spec pivots.
#'
#' @param shock_seq    n_exo x T shock matrix
#' @param dr_slack     Slack-regime DecisionRules
#' @param sys          System matrices
#' @param specs        OBC spec list
#' @param obs_idx      Integer vector of observable positions
#' @param state_init   Numeric state vector or NULL (zero initial state)
#' @return Named list with:
#'   $q    -- numeric vector (length T*n_spec): all-slack complementarity gaps
#'   $M    -- numeric matrix (T*n_spec × T*n_spec): sensitivity (lower triangular)
#'   $sim0 -- all-slack simulation (list with $paths and $states)
#' @noRd
lcp_build_system <- function(shock_seq, dr_slack, sys, specs, obs_idx,
                              state_init = NULL) {
  n_T    <- ncol(shock_seq)
  n_spec <- length(specs)
  n_lcp  <- n_T * n_spec
  n_state <- length(dr_slack$state_idx)

  # Seed a local cache with the slack policy and one individual-binding policy
  # per spec (bitfield = 2^(j-1) for spec j).
  rc <- new.env(parent = emptyenv(), hash = TRUE)
  obc_ensure_policy(0L, rc, sys, dr_slack, specs, obs_idx)
  for (j in seq_len(n_spec))
    obc_ensure_policy(2L^(j - 1L), rc, sys, dr_slack, specs, obs_idx)

  # All-slack simulation
  sim0 <- boehl_simulate(shock_seq, dr_slack, rc, integer(n_T), state_init)

  # q vector: complementarity gaps under the all-slack simulation.
  # q_k = sign_j * (x_j^slack(t) - bound_j).
  q       <- numeric(n_lcp)
  s_prev  <- if (!is.null(state_init) && length(state_init) == n_state)
               as.numeric(state_init) else numeric(n_state)

  for (t in seq_len(n_T)) {
    eps <- shock_seq[, t]
    for (j in seq_len(n_spec)) {
      s   <- specs[[j]]
      x_s <- sum(dr_slack$ghx[s$var_idx, ] * s_prev) +
             sum(dr_slack$ghu[s$var_idx, ] * eps)
      sign_j        <- if (s$op == ">") 1 else -1
      q[(t - 1L) * n_spec + j] <- sign_j * (x_s - s$bound)
    }
    s_prev <- sim0$states[, t]
  }

  # M matrix: finite-difference sensitivity at the all-slack baseline.
  # Column col_k = (t_src-1)*n_spec + j_src: simulate with only spec j_src
  # binding at period t_src; M[:,col_k] = z_mod - z_0 (unclipped gaps).
  M <- matrix(0, n_lcp, n_lcp)

  for (t_src in seq_len(n_T)) {
    for (j_src in seq_len(n_spec)) {
      col   <- (t_src - 1L) * n_spec + j_src
      r_bind <- 2L^(j_src - 1L)

      rp_mod         <- integer(n_T)
      rp_mod[t_src]  <- r_bind

      sim_mod <- boehl_simulate(shock_seq, dr_slack, rc, rp_mod, state_init)

      # z_0[k'] = sign_{j'} * (paths_0[var_{j'}, t'] - bound_{j'})
      # z_mod[k'] = sign_{j'} * (paths_mod[var_{j'}, t'] - bound_{j'})
      # M[k', col] = z_mod[k'] - z_0[k']
      for (t_prime in seq_len(n_T)) {
        for (j_prime in seq_len(n_spec)) {
          s_prime    <- specs[[j_prime]]
          sign_prime <- if (s_prime$op == ">") 1 else -1
          row        <- (t_prime - 1L) * n_spec + j_prime
          M[row, col] <- sign_prime * (sim_mod$paths[s_prime$var_idx, t_prime] -
                                        sim0$paths[s_prime$var_idx, t_prime])
        }
      }
    }
  }

  list(q = q, M = M, sim0 = sim0)
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
#' For the OccBin stacked LCP, M is block-lower-triangular in t (causal model
#' structure) with positive diagonal at binding periods.  This makes M a
#' P-matrix on the binding sub-block, guaranteeing monotone termination in
#' at most T*n_spec pivots (single forward substitution pass).
#'
#' @param q        Numeric vector (length n): constant term of the LCP
#' @param M        Numeric matrix (n × n): LCP coefficient matrix
#' @param d        Covering vector (length n, all positive); NULL → ones
#' @param max_iter Maximum pivot steps; NULL → 50*n
#' @param tol      Numerical tolerance (default 1e-10)
#' @return Named list:
#'   $z         -- numeric vector (length n): primal solution (z ≥ 0)
#'   $w         -- numeric vector (length n): dual solution (w = Mz + q ≥ 0)
#'   $converged -- logical
#'   $n_iter    -- integer: pivot steps taken
#'   $message   -- character: termination reason
#' @noRd
lcp_lemke <- function(q, M, d = NULL, max_iter = NULL, tol = 1e-10) {
  n <- length(q)
  if (is.null(d))        d        <- rep(1.0, n)
  if (is.null(max_iter)) max_iter <- 50L * n

  stopifnot(nrow(M) == n, ncol(M) == n, length(d) == n, all(d > 0))

  # Trivial solution: q >= 0 already.
  if (all(q >= -tol))
    return(list(z = numeric(n), w = pmax(q, 0),
                converged = TRUE, n_iter = 0L, message = "trivial"))

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
                  message = "secondary ray: problem is unbounded"))

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
      # Recompute w for numerical accuracy.
      w <- pmax(drop(M %*% z) + q, 0)
      return(list(z = z, w = w,
                  converged = TRUE, n_iter = iter,
                  message = "z_0 left basis"))
    }

    enter_col <- comp_col(leaving_var)
  }

  # Extract best approximation at max_iter.
  z <- numeric(n)
  for (i in seq_len(n)) {
    bv <- bas[i]
    if (bv >= 1L && bv <= n) z[bv] <- max(0, tab[i, 2L*n+2L])
  }
  list(z = z, w = pmax(drop(M %*% z) + q, 0),
       converged = FALSE, n_iter = max_iter, message = "max iterations reached")
}


# =============================================================================
# Lemke-based OBC regime solver
# =============================================================================

#' Extract an OBC regime path from a Lemke LCP solution
#'
#' Given the z solution of LCP(q, M) built by lcp_build_system(), maps z back
#' to a binary integer regime path: z_k > tol means constraint k binds.
#' Then forward-simulates the model under that regime path.
#'
#' @param z        Numeric vector (length T*n_spec): Lemke primal solution
#' @param shock_seq n_exo x T shock matrix
#' @param dr_slack Slack-regime DecisionRules
#' @param sys      System matrices
#' @param specs    OBC spec list
#' @param obs_idx  Integer vector of observable positions
#' @param state_init Numeric state vector or NULL
#' @param tol      Threshold for treating z_k as positive (default 1e-6)
#' @return Named list with $regime_path, $paths, $states
#' @noRd
lcp_regime_from_z <- function(z, shock_seq, dr_slack, sys, specs, obs_idx,
                               state_init = NULL, tol = 1e-6) {
  n_T    <- ncol(shock_seq)
  n_spec <- length(specs)

  rc <- new.env(parent = emptyenv(), hash = TRUE)
  obc_ensure_policy(0L, rc, sys, dr_slack, specs, obs_idx)

  regime_path <- integer(n_T)
  for (t in seq_len(n_T)) {
    flags <- logical(n_spec)
    for (j in seq_len(n_spec)) {
      k <- (t - 1L) * n_spec + j
      flags[j] <- z[k] > tol
    }
    r <- obc_regime_idx(flags)
    if (r != 0L) obc_ensure_policy(r, rc, sys, dr_slack, specs, obs_idx)
    regime_path[t] <- r
  }

  sim <- boehl_simulate(shock_seq, dr_slack, rc, regime_path, state_init)
  list(regime_path = regime_path, paths = sim$paths, states = sim$states)
}


# =============================================================================
# Main solver
# =============================================================================

#' LCP OBC simulation solver
#'
#' Finds the binding/slack regime path for a given deterministic shock sequence
#' using one of two methods selected by the \code{method} argument.
#'
#' \strong{method = "lemke"} (default): Lemke's complementarity pivoting.
#'   Builds the stacked T*n_spec LCP(q, M) from an all-slack simulation baseline
#'   and solves it with the Lemke principal-pivot algorithm.  For the OccBin
#'   piecewise-linear model M is block-lower-triangular (causal structure), so
#'   Lemke terminates in exactly T*n_spec pivots (one forward substitution pass)
#'   with provable monotone guarantees.
#'
#' \strong{method = "newton"}: Newton iteration on the stacked complementarity
#'   mismatch residual.  At each step the T*n_spec × T*n_spec finite-difference
#'   Jacobian is assembled (T*n_spec re-simulations), and a QR Newton step
#'   gives a regime update.  Useful as a fallback when Lemke does not converge.
#'
#' @param shock_seq        n_exo x T numeric matrix of structural shocks
#' @param dr_slack         Slack-regime DecisionRules (from solve_perturbation)
#' @param sys              System matrices (from extract_system_matrices_fast)
#' @param specs            OBC spec list (from obc_parse_tags / obc_collect_specs)
#' @param obs_idx          Integer vector: observable positions in endo vector;
#'                         NULL → all endo indices
#' @param method           Character: \code{"lemke"} (default) or \code{"newton"}
#' @param max_iter         Maximum pivot steps (lemke) or Newton iterations
#'                         (newton); NULL → solver-specific defaults
#' @param tol              Complementarity tolerance (default 1e-8)
#' @param regime_path_init Integer vector (length T) for warm-starting the
#'                         Newton solver; ignored for method = "lemke"
#' @param state_init       Numeric vector (length n_state) initial model state;
#'                         NULL → zero initial state (standard IRF convention)
#' @return Named list with:
#'   $regime_path -- integer vector (length T): final bitfield regime per period
#'   $paths       -- n_endo x T matrix of simulated endogenous variable paths
#'   $fb_norm     -- scalar: complementarity mismatch norm at the returned path
#'   $converged   -- logical
#'   $n_iter      -- integer: pivot steps (lemke) or Newton iterations (newton)
#'
#' @references
#'   Boehl, G. (2020). Efficient solution and computation of models with
#'     occasionally binding constraints. \emph{Deutsche Bundesbank Discussion
#'     Paper}, 38/2020.
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
  n_spec <- length(specs)

  if (is.null(obs_idx)) obs_idx <- seq_len(nrow(dr_slack$ghx))

  # ------------------------------------------------------------------
  # LEMKE PATH
  # ------------------------------------------------------------------
  if (method == "lemke") {
    lemke_max <- if (!is.null(max_iter)) as.integer(max_iter) else NULL

    lcp_sys <- lcp_build_system(shock_seq, dr_slack, sys, specs, obs_idx,
                                 state_init)
    lemke_res <- lcp_lemke(lcp_sys$q, lcp_sys$M,
                            max_iter = lemke_max, tol = tol)

    # Map Lemke z solution to binary regime path and simulate.
    sol <- lcp_regime_from_z(lemke_res$z, shock_seq, dr_slack, sys, specs,
                              obs_idx, state_init, tol = tol)

    # Compute complementarity mismatch norm of the resulting regime.
    rc_tmp <- new.env(parent = emptyenv(), hash = TRUE)
    obc_ensure_policy(0L, rc_tmp, sys, dr_slack, specs, obs_idx)
    for (r in unique(sol$regime_path[sol$regime_path != 0L]))
      obc_ensure_policy(r, rc_tmp, sys, dr_slack, specs, obs_idx)

    sim_final  <- boehl_simulate(shock_seq, dr_slack, rc_tmp,
                                  sol$regime_path, state_init)
    comp_mat   <- lcp_comp_residual(sim_final$states, shock_seq, specs,
                                     dr_slack, sol$regime_path)
    comp_norm  <- sqrt(sum(comp_mat^2))

    return(list(
      regime_path = sol$regime_path,
      paths       = sim_final$paths,
      fb_norm     = comp_norm,
      converged   = lemke_res$converged && comp_norm < max(tol, 1e-6),
      n_iter      = lemke_res$n_iter
    ))
  }

  # ------------------------------------------------------------------
  # NEWTON PATH
  # ------------------------------------------------------------------
  newton_max <- if (!is.null(max_iter)) as.integer(max_iter) else 20L

  regime_cache <- new.env(parent = emptyenv(), hash = TRUE)
  obc_ensure_policy(0L, regime_cache, sys, dr_slack, specs, obs_idx)

  regime_path <- if (!is.null(regime_path_init) && length(regime_path_init) == n_T)
                   as.integer(regime_path_init)
                 else integer(n_T)

  for (r in unique(regime_path[regime_path != 0L]))
    obc_ensure_policy(r, regime_cache, sys, dr_slack, specs, obs_idx)

  sim <- NULL

  for (iter in seq_len(newton_max)) {
    sim      <- boehl_simulate(shock_seq, dr_slack, regime_cache, regime_path,
                               state_init)
    comp_mat <- lcp_comp_residual(sim$states, shock_seq, specs, dr_slack,
                                  regime_path)
    comp_vec <- as.vector(comp_mat)
    comp_norm <- sqrt(sum(comp_vec^2))

    if (comp_norm < tol) {
      return(list(
        regime_path = regime_path,
        paths       = sim$paths,
        fb_norm     = comp_norm,
        converged   = TRUE,
        n_iter      = iter - 1L
      ))
    }

    # Build the T*n_spec × T*n_spec Jacobian by finite differencing each
    # binary regime flip direction (one re-simulation per column).
    J <- lcp_build_jacobian(comp_vec, shock_seq, dr_slack, regime_cache,
                             sys, regime_path, specs, obs_idx, state_init)

    # Newton step: solve J * delta = −comp_vec.
    # SVD regularised pseudoinverse when J is rank-deficient.
    delta <- qr.solve(J, -comp_vec, tol = 1e-12)

    # Map delta to binary regime updates: delta_{jt} > 0.5 → do the flip.
    new_regime <- integer(n_T)
    for (t in seq_len(n_T)) {
      flags <- obc_regime_flags(regime_path[t], n_spec)
      for (j in seq_len(n_spec)) {
        col <- (t - 1L) * n_spec + j
        if (delta[col] > 0.5) flags[j] <- !flags[j]
      }
      r_new <- obc_regime_idx(flags)
      if (r_new != 0L)
        obc_ensure_policy(r_new, regime_cache, sys, dr_slack, specs, obs_idx)
      new_regime[t] <- r_new
    }

    if (identical(new_regime, regime_path)) {
      # Newton step did not change the regime — fixed point reached.
      return(list(
        regime_path = regime_path,
        paths       = sim$paths,
        fb_norm     = comp_norm,
        converged   = comp_norm < tol,
        n_iter      = iter
      ))
    }

    regime_path <- new_regime
  }

  # Maximum iterations reached.
  sim      <- boehl_simulate(shock_seq, dr_slack, regime_cache, regime_path,
                             state_init)
  comp_mat <- lcp_comp_residual(sim$states, shock_seq, specs, dr_slack, regime_path)
  comp_norm <- sqrt(sum(comp_mat^2))

  list(
    regime_path = regime_path,
    paths       = sim$paths,
    fb_norm     = comp_norm,
    converged   = FALSE,
    n_iter      = newton_max
  )
}
