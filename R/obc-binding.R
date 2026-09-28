## R/obc-binding.R
## --------------------------------------------------------------------------
## OBC binding-regime algebra.
##
## Provides:
##   obc_build_binding_sys()   -- substitute constraint equations into Jacobians
##   obc_solve_binding()       -- OccBin terminal-substitution policy solve
##   .obc_spec_rows()          -- system row of each spec's tagged equation
##   .obc_pwl_*()              -- piecewise-linear OccBin engine: time-varying
##                                rules along a regime path, forward
##                                simulation, complementarity check,
##                                guess-and-verify, surprise-shock re-solves
##
## The OccBin approach (Guerrieri-Iacoviello 2015) constructs the binding-regime
## policy by replacing each tagged equation with its constraint equality and
## substituting the SLACK policy as the next-period terminal condition.  This
## gives a determinate linear system even when the Taylor principle is lost
## (e.g. at the ZLB).  It is the rule of the LAST period of a spell only; a
## spell of 2+ periods needs the backward recursion of .obc_pwl_rules().
##
## Consumers of the engine below: the deterministic path solvers
## (boehl_solve_regime_path, solve_obc_lcp, obc_simulate, compute_irfs_obc,
## ramsey_obc_pwlinear, occbin_solve_path(method = "pwlinear"), W48) and,
## since W49 (2026-09-25), the piecewise-linear Kalman filter
## kalman_filter_obc_pkf() (Dynare's OccBin PKF) with everything built on it
## (obc_guess_verify, kalman_filter_obc, the smoother, the historical
## decomposition, the PKF posterior and conditional forecasts), the
## inversion filter (R/obc-inversion-filter.R, per-period rules along the
## regime path) and, since W50, the particle filters (PPF / COPF,
## R/obc-ppf.R: each particle follows the regime sequence solved from its own
## state).  No filter uses the one-period policies of obc_ensure_policy()
## (next period slack) any more.
##
## Design note: a near-zero lag coefficient (1e-10) is retained on the
## constrained variable so its classification as a state variable is preserved
## and the Kalman filter state dimension is unchanged across regimes.
## --------------------------------------------------------------------------


# =============================================================================
# Binding system construction
# =============================================================================

#' Build the binding-regime system matrices
#'
#' For each OBC spec, replaces row eq_idx of the system Jacobians with the
#' constraint equation: var_idx_t = bound. To preserve the state-space
#' dimension (critical for the Kalman filter), a near-zero lag coefficient
#' (1e-10) is retained for the constrained variable so the classifier keeps
#' it as a state variable. The bound itself is captured in const_b.
#'
#' @param sys   System matrices list (from extract_system_matrices_fast)
#' @param specs List of OBC specs from obc_parse_tags
#' @return List:
#'   $sys_b   -- modified system matrices (same structure as sys)
#'   $const_b -- numeric vector (length n_endo), nonzero at bound equations
#' @noRd
obc_build_binding_sys <- function(sys, specs) {
  f_minus_b <- sys$f_minus
  f_zero_b  <- sys$f_zero
  f_plus_b  <- sys$f_plus
  f_exo_b   <- sys$f_exo
  const_b   <- numeric(sys$n_endo)
  rows      <- .obc_spec_rows(sys, specs)

  for (k in seq_along(specs)) {
    s <- specs[[k]]
    i <- rows[[k]]
    j <- s$var_idx
    b <- s$bound

    # Replace row i with the binding constraint: var_j_t = bound
    # f_zero[i, j] = 1 (contemporaneous var_j)
    # const_b[i]   = b (the bound, which becomes the RHS forcing)
    #
    # Near-zero lag on j (1e-10 << any economic quantity) retains var_j in
    # the state-variable classification, keeping T_b and T_s the same size.
    f_minus_b[i, ] <- 0
    f_zero_b[i, ]  <- 0
    f_plus_b[i, ]  <- 0
    f_exo_b[i, ]   <- 0

    f_zero_b[i, j]  <- 1
    f_minus_b[i, j] <- 1e-10   # keeps var_j classified as state variable
    const_b[i]      <- b
  }

  sys_b         <- sys
  sys_b$f_minus <- f_minus_b
  sys_b$f_zero  <- f_zero_b
  sys_b$f_plus  <- f_plus_b
  sys_b$f_exo   <- f_exo_b

  list(sys_b = sys_b, const_b = const_b)
}


#' Row of the system matrices that holds each spec's tagged equation
#'
#' \code{spec$eq_idx} is the equation's position in \code{model$equations}.
#' \code{extract_system_matrices()} permutes the Jacobian rows into
#' declaration-variable order (\code{J[order(eq_to_decl), ]}) and returns
#' \code{eq_to_decl}; \code{extract_system_matrices_fast()} keeps the model
#' order and returns no \code{eq_to_decl}.  Replacing row \code{eq_idx} of a
#' permuted system replaced ANOTHER equation whenever the equations were not
#' written in declaration order (fixed 2026-09-25, W48).  A spec that already
#' carries its system row (\code{row_idx}, e.g. the OccBin block specs of
#' .occbin_build_pwlinear_drs()) is used as is.
#'
#' @param sys   System matrices
#' @param specs OBC spec list
#' @return Integer vector: system row per spec
#' @noRd
.obc_spec_rows <- function(sys, specs) {
  n   <- nrow(sys$f_zero)
  e   <- sys$eq_to_decl
  pos <- seq_len(n)
  if (!is.null(e) && length(e) == n && all(e > 0L))
    pos[order(e)] <- seq_len(n)
  vapply(specs, function(s) {
    if (!is.null(s$row_idx)) as.integer(s$row_idx) else pos[[s$eq_idx]]
  }, integer(1))
}


# =============================================================================
# Binding policy function
# =============================================================================

#' Compute the binding-regime policy function via the OccBin terminal substitution
#'
#' A permanently pegged rate (or any binding OBC) makes the linearised model
#' indeterminate: the Taylor principle is lost and there is no nominal anchor.
#' The OccBin solution (Guerrieri-Iacoviello 2015) avoids this by substituting
#' the SLACK policy function as the terminal condition for the next period:
#'
#'   E\[y_{t+1}\] ~= dr_slack$ghx * s_t    (next period is expected to be slack)
#'
#' Substituting into the binding-regime system:
#'   (f_zero_b + f_plus_b * ghx_slack * E_s) * y_t = const_b - f_minus_b\[:,state_idx\] * s_{t-1}
#'
#' where E_s is the n_state x n_endo selection matrix that picks state variables.
#' This linear system is generically non-singular (slack policy provides the anchor)
#' and gives a well-defined, determinate binding-regime policy function.
#'
#' Returns a pseudo-DecisionRules object with the same interface as dr_slack
#' (ghx, ghu, state_idx, etc.) so kalman_filter_obc can use it directly.
#'
#' OBC bounds in DEVIATION-FROM-STEADY-STATE units
#'
#' THE package convention (fixed 0.9.4, ledger A6): an \code{[mcp = 'v > b']}
#' bound is a LEVEL, in the same units as the model variable \code{v}.  Every
#' internal comparison happens in deviations (\code{ghx}/\code{ghu} produce
#' deviations, and \code{boehl_simulate()}'s \code{paths} are deviations), so
#' the bound must be shifted by the steady state:
#'
#'   v_t >= b   <=>   (v_t - v_ss) >= b - v_ss
#'
#' \code{obc_solve_binding()} already pinned at \code{b - v_ss}, while
#' \code{R/obc-boehl.R} compared raw deviations against \code{b}.  With a
#' non-zero steady state the two disagreed and the Boehl slackness test
#' declared the constraint binding in EVERY period.  Both now call this.
#'
#' @param specs     OBC spec list from obc_parse_tags
#' @param dr_slack  Slack-regime DecisionRules (supplies \code{$ys})
#' @return Numeric vector, one deviation-form bound per spec.
#' @noRd
.obc_bound_dev <- function(specs, dr_slack) {
  ys <- if (inherits(dr_slack$ys, "dynhr_steady")) dr_slack$ys$values else dr_slack$ys
  vapply(specs, function(s) {
    lvl <- if (!is.null(ys) && length(ys) >= s$var_idx) as.numeric(ys[[s$var_idx]]) else 0
    if (!is.finite(lvl)) lvl <- 0
    s$bound - lvl
  }, numeric(1))
}


#' @param sys       System matrices (from extract_system_matrices_fast)
#' @param dr_slack  Slack-regime DecisionRules (provides terminal condition)
#' @param specs     OBC spec list from obc_parse_tags
#' @param obs_idx   Integer vector: observable positions in endo vector
#' @return List($dr, $c_state, $c_obs, $c_full), or NULL if linear solve fails
#' @noRd
obc_solve_binding <- function(sys, dr_slack, specs, obs_idx) {
  built   <- obc_build_binding_sys(sys, specs)
  sys_b   <- built$sys_b
  const_b <- built$const_b

  # Correct const_b to deviation form (see .obc_bound_dev below).
  bound_dev <- .obc_bound_dev(specs, dr_slack)
  rows      <- .obc_spec_rows(sys, specs)
  for (j in seq_along(specs)) {
    const_b[rows[[j]]] <- bound_dev[[j]]
  }

  n_endo    <- sys$n_endo
  n_exo     <- sys$n_exo
  state_idx <- dr_slack$state_idx
  n_state   <- length(state_idx)

  # Selection matrix E_s (n_state x n_endo): picks state variables from y_t
  # so that y_t[state_idx] = E_s %*% y_t
  E_s <- matrix(0, n_state, n_endo)
  for (i in seq_len(n_state)) E_s[i, state_idx[i]] <- 1

  # Binding-regime coefficient matrix with slack terminal condition substituted:
  #   A_b = f_zero_b + f_plus_b * ghx_slack * E_s   (n_endo x n_endo)
  A_b <- sys_b$f_zero + sys_b$f_plus %*% dr_slack$ghx %*% E_s

  # Solve for the three components of the binding policy function:
  #   y_t = ghx_b * s_{t-1} + ghu_b * eps_t + c_b
  # where:
  #   ghx_b (n_endo x n_state) = -A_b^{-1} * f_minus_b[:, state_idx]
  #   ghu_b (n_endo x n_exo)   = -A_b^{-1} * f_exo_b
  #   c_b   (n_endo)           =  A_b^{-1} * const_b
  rhs_ghx <- -sys_b$f_minus[, state_idx, drop = FALSE]
  rhs_ghu <- -sys_b$f_exo
  rhs_c   <- const_b

  # Check condition number before solve; fall back to SVD pseudoinverse if singular
  rcond_Ab <- rcond(A_b)
  if (is.finite(rcond_Ab) && rcond_Ab > .Machine$double.eps) {
    soln <- list(
      ghx = solve(A_b, rhs_ghx),
      ghu = solve(A_b, rhs_ghu),
      c   = solve(A_b, rhs_c)
    )
  } else {
    # Fallback: pseudoinverse via SVD
    sv    <- svd(A_b)
    d_inv <- ifelse(sv$d > 1e-10, 1 / sv$d, 0)
    A_inv <- sv$v %*% diag(d_inv, nrow = length(d_inv)) %*% t(sv$u)
    soln <- list(ghx = A_inv %*% rhs_ghx, ghu = A_inv %*% rhs_ghu, c = A_inv %*% rhs_c)
  }

  rownames(soln$ghx) <- dr_slack$endo_names
  colnames(soln$ghx) <- dr_slack$endo_names[state_idx]
  rownames(soln$ghu) <- dr_slack$endo_names
  if (length(dr_slack$exo_names) > 0) colnames(soln$ghu) <- dr_slack$exo_names

  # Construct a pseudo-DecisionRules object for the binding regime.
  # Uses the same state_idx as dr_slack so kalman_filter_obc dimensions match.
  dr_b <- structure(list(
    ghx          = soln$ghx,
    ghu          = soln$ghu,
    ys           = dr_slack$ys,
    endo_names   = dr_slack$endo_names,
    exo_names    = dr_slack$exo_names,
    state_vars   = dr_slack$state_vars,
    state_idx    = state_idx,
    n_state      = n_state,
    n_stable     = n_state,
    n_exo        = n_exo,
    eigenvalues  = NULL,
    n_unstable   = 0L,
    bk_satisfied = TRUE
  ), class = "DecisionRules")

  c_full  <- soln$c
  c_state <- c_full[state_idx]
  c_obs   <- c_full[obs_idx]

  list(dr = dr_b, c_state = c_state, c_obs = c_obs, c_full = c_full)
}


# =============================================================================
# Piecewise-linear OccBin engine: time-varying rules along a regime path
# =============================================================================
#
# obc_solve_binding() above gives the rule of a period that binds while the
# NEXT period is slack.  That is the right rule only for the LAST period of a
# spell.  For a regime path r_1..r_H (slack after H) the rational-expectations
# piecewise-linear solution (Guerrieri-Iacoviello 2015, Dynare OccBin) has a
# rule per period, obtained by a backward recursion from the last binding
# period L:
#
#   E_t y_{t+1} = G_{t+1} s_t + c_{t+1}           (G_{L+1} = ghx_slack, c = 0)
#   (F0_r + Fp_r G_{t+1} E_s) y_t
#        = k_r - Fp_r c_{t+1} - Fm_r s_{t-1} - Fe_r eps_t
#
# where (Fm_r, F0_r, Fp_r, Fe_r, k_r) is the linear system of regime r_t (the
# bound row x_j = b_j replacing the tagged equation of every binding spec).
# The recursion runs through SLACK periods before L too: a period that is
# slack but anticipates a later binding period does not follow ghx_slack.
#
# The regime path is verified against the CONSTRAINED path it produces
# (OccBin guess-and-verify): a slack period binds when x_j violates its
# bound; a binding period is released when the multiplier of its relaxed
# equation has the wrong sign (the residual F = lhs - rhs of the tagged
# equation, evaluated on the path with E_t y_{t+1}; Dynare's MCP / lmmcp
# convention F >= 0 at a lower bound, F <= 0 at an upper bound -- the same
# test as pf_newton_solve()).  A constraint that is violated only because
# another one binds is therefore imposed.
# =============================================================================

#' Linear-system context for the piecewise-linear OccBin engine
#'
#' @param sys      System matrices (f_minus, f_zero, f_plus, f_exo; n_endo rows)
#' @param dr_slack Slack-regime DecisionRules
#' @param specs    OBC spec list (eq_idx, var_idx, bound, and op)
#' @return List of the matrices, spec indices, signs and deviation bounds
#' @noRd
.obc_pwl_context <- function(sys, dr_slack, specs) {
  si <- dr_slack$state_idx
  list(
    Fm      = sys$f_minus[, si, drop = FALSE],
    F0      = sys$f_zero,
    Fp      = sys$f_plus,
    Fe      = sys$f_exo,
    ghx     = dr_slack$ghx,
    ghu     = dr_slack$ghu,
    si      = si,
    n_endo  = nrow(dr_slack$ghx),
    n_state = length(si),
    n_exo   = ncol(dr_slack$ghu),
    n_spec  = length(specs),
    eq      = .obc_spec_rows(sys, specs),
    var     = vapply(specs, function(s) as.integer(s$var_idx), integer(1)),
    sgn     = vapply(specs, function(s) if (identical(s$op, "<")) -1 else 1,
                     numeric(1)),
    bnd     = .obc_bound_dev(specs, dr_slack),
    endo_names = dr_slack$endo_names
  )
}

#' Engine context stored in a regime cache by obc_ensure_policy()
#' @noRd
.obc_pwl_cache_context <- function(regime_cache) {
  if (exists(".pwl_ctx", envir = regime_cache, inherits = FALSE))
    return(get(".pwl_ctx", envir = regime_cache, inherits = FALSE))
  if (!exists(".pwl_src", envir = regime_cache, inherits = FALSE))
    .dynhr_abort("the OBC regime cache was not seeded by obc_ensure_policy(), ",
                 "so the system matrices of the piecewise-linear solution ",
                 "are unknown.", class = "dynhr_error_obc_cache")
  src <- get(".pwl_src", envir = regime_cache, inherits = FALSE)
  ctx <- .obc_pwl_context(src$sys, src$dr_slack, src$specs)
  assign(".pwl_ctx", ctx, envir = regime_cache)
  ctx
}

#' Linear system of one regime (bound rows replacing the tagged equations)
#' @noRd
.obc_pwl_regime_sys <- function(ctx, regime_idx) {
  Fm <- ctx$Fm; F0 <- ctx$F0; Fp <- ctx$Fp; Fe <- ctx$Fe
  k  <- numeric(ctx$n_endo)
  if (regime_idx != 0L) {
    for (j in which(obc_regime_flags(regime_idx, ctx$n_spec))) {
      i <- ctx$eq[j]
      Fm[i, ] <- 0; F0[i, ] <- 0; Fp[i, ] <- 0; Fe[i, ] <- 0
      F0[i, ctx$var[j]] <- 1
      k[i] <- ctx$bnd[j]
    }
  }
  list(Fm = Fm, F0 = F0, Fp = Fp, Fe = Fe, k = k)
}

#' Time-varying decision rules along a regime path (backward recursion)
#'
#' @param ctx         Context from .obc_pwl_context()
#' @param regime_path Integer bitfield per period (slack after the path)
#' @param strict      TRUE: abort on a singular regime system; FALSE: return
#'   NULL instead (the Kalman filter treats such a regime guess as a failed
#'   candidate, as Dynare's OccBin filter does)
#' @return List (length of regime_path): element t is NULL (slack rule
#'   ghx_slack/ghu_slack, no constant) or list(ghx, ghu, c) with
#'   y_t = ghx s_{t-1} + ghu eps_t + c (deviations)
#' @noRd
.obc_pwl_rules <- function(ctx, regime_path, strict = TRUE) {
  regime_path <- as.integer(regime_path)
  rules <- vector("list", length(regime_path))
  bind  <- which(regime_path != 0L)
  if (length(bind) == 0L) return(rules)

  ns <- ctx$n_state
  ne <- ctx$n_exo
  G  <- ctx$ghx                     # rule of period t+1 (slack after L)
  cc <- numeric(ctx$n_endo)
  pieces <- list()
  for (t in seq.int(max(bind), 1L)) {
    key <- as.character(regime_path[t])
    if (is.null(pieces[[key]]))
      pieces[[key]] <- .obc_pwl_regime_sys(ctx, regime_path[t])
    m <- pieces[[key]]
    A <- m$F0
    A[, ctx$si] <- A[, ctx$si] + m$Fp %*% G
    rc <- rcond(A)
    if (!strict && (!is.finite(rc) || rc <= .Machine$double.eps))
      return(NULL)
    if (!is.finite(rc) || rc <= .Machine$double.eps)
      .dynhr_abort(sprintf(paste0(
        "the piecewise-linear OBC system of regime %d in period %d is ",
        "singular (rcond %.2e): this regime path has no determinate ",
        "solution."), regime_path[t], t, rc),
        class = "dynhr_error_obc_singular_regime")
    sol <- solve(A, cbind(-m$Fm, -m$Fe, m$k - drop(m$Fp %*% cc)))
    G  <- sol[, seq_len(ns), drop = FALSE]
    cc <- sol[, ns + ne + 1L]
    rules[[t]] <- list(ghx = G, ghu = sol[, ns + seq_len(ne), drop = FALSE],
                       c = cc)
  }
  rules
}

#' Forward simulation under time-varying rules
#'
#' Shocks are surprises: eps_t enters through the period-t rule only.
#' @return list(paths = n_endo x T, states = n_state x T)
#' @noRd
.obc_pwl_forward <- function(ctx, rules, shock_seq, state_init = NULL) {
  n_T    <- ncol(shock_seq)
  paths  <- matrix(NA_real_, ctx$n_endo, n_T)
  states <- matrix(NA_real_, ctx$n_state, n_T)
  rownames(paths) <- ctx$endo_names
  s <- if (!is.null(state_init) && length(state_init) == ctx$n_state)
         as.numeric(state_init) else numeric(ctx$n_state)
  for (t in seq_len(n_T)) {
    ru  <- if (t <= length(rules)) rules[[t]] else NULL
    eps <- shock_seq[, t]
    y <- if (is.null(ru)) drop(ctx$ghx %*% s) + drop(ctx$ghu %*% eps)
         else drop(ru$ghx %*% s) + drop(ru$ghu %*% eps) + ru$c
    paths[, t]  <- y
    s           <- y[ctx$si]
    states[, t] <- s
  }
  list(paths = paths, states = states)
}

#' Scale-free complementarity decisions of the OccBin guess-and-verify
#'
#' Dynare 7.1's OccBin decides with the PLAIN inequalities of the
#' occbin_constraints block (the preprocessor-generated occbin_difference.m:
#' binding = x < b, relax = the relax condition, e.g. x > b or a shadow value
#' above the bound; no tolerance), and solve_one_constraint.m updates
#' B <- (B | binding) & !(B & relax).  dynhr takes the same decisions up to a
#' round-off band RELATIVE to the magnitude of the quantities compared, so the
#' regimes do not depend on the units of the model.  (Before W78 (2026-09) the
#' band was an absolute 1e-8 in level units: with every shock std, the data
#' and the bound scaled by 1e-4 a violation of 1 percent of the bound was
#' ignored.)  A slack period binds when gap = sgn (x - b) < -rtol (|x| + |b|);
#' a binding period stays binding while mult = sgn F >= -rtol M, with M the
#' sum of the absolute terms of the tagged equation's residual F.
#' @param gap,x  sign-adjusted slack and the constrained variable's value
#' @param bnd    deviation-form bound, recycled down the rows (one entry per
#'   constraint row)
#' @param rtol   relative round-off band
#' @noRd
.obc_gap_binds <- function(gap, x, bnd, rtol) {
  gap < -rtol * (abs(x) + abs(bnd))
}

#' Binding period kept binding (see .obc_gap_binds())
#' @param mult,mag sign-adjusted residual of the tagged equation and the sum
#'   of its absolute terms
#' @noRd
.obc_mult_keeps <- function(mult, mag, rtol) {
  mult >= -rtol * mag
}

#' Complementarity check of a regime path against its constrained path
#'
#' @param tol relative round-off band of the decisions (.obc_gap_binds())
#' @return list(regime_path = the verified (updated) regime path,
#'   gap = n_spec x T sign-adjusted slack sgn*(x - b),
#'   mult = n_spec x T sign-adjusted multiplier sgn*F of the tagged equation,
#'   phi = Fischer-Burmeister residual of (gap, mult), zero at a solution)
#' @noRd
.obc_pwl_check <- function(ctx, rules, sim, shock_seq, regime_path,
                           state_init = NULL, tol = 1e-8) {
  n_T    <- ncol(shock_seq)
  n_spec <- ctx$n_spec
  gap  <- matrix(0, n_spec, n_T)
  mult <- matrix(0, n_spec, n_T)
  if (n_spec == 0L)
    return(list(regime_path = integer(n_T), gap = gap, mult = mult,
                phi = gap))
  eq  <- ctx$eq
  Fm  <- ctx$Fm[eq, , drop = FALSE]
  F0  <- ctx$F0[eq, , drop = FALSE]
  Fp  <- ctx$Fp[eq, , drop = FALSE]
  Fe  <- ctx$Fe[eq, , drop = FALSE]
  s_prev <- if (!is.null(state_init) && length(state_init) == ctx$n_state)
              as.numeric(state_init) else numeric(ctx$n_state)
  mag <- matrix(0, n_spec, n_T)
  for (t in seq_len(n_T)) {
    y_t <- sim$paths[, t]
    s_t <- sim$states[, t]
    ## E_t y_{t+1}: the period-(t+1) rule at s_t with no (surprise) shock
    ru  <- if (t < length(rules)) rules[[t + 1L]] else NULL
    Ey1 <- if (is.null(ru)) drop(ctx$ghx %*% s_t)
           else drop(ru$ghx %*% s_t) + ru$c
    Fr  <- drop(Fm %*% s_prev) + drop(F0 %*% y_t) + drop(Fp %*% Ey1) +
           drop(Fe %*% shock_seq[, t])
    gap[, t]  <- ctx$sgn * (y_t[ctx$var] - ctx$bnd)
    mult[, t] <- ctx$sgn * Fr
    mag[, t]  <- drop(abs(Fm) %*% abs(s_prev)) + drop(abs(F0) %*% abs(y_t)) +
                 drop(abs(Fp) %*% abs(Ey1)) +
                 drop(abs(Fe) %*% abs(shock_seq[, t]))
    s_prev <- s_t
  }
  B     <- t(.obc_bitfield_to_binding(regime_path, n_spec))   # n_spec x T
  X     <- sim$paths[ctx$var, , drop = FALSE]
  B_new <- ifelse(B, .obc_mult_keeps(mult, mag, tol),
                  .obc_gap_binds(gap, X, ctx$bnd, tol))
  w     <- 2L^(seq_len(n_spec) - 1L)
  list(regime_path = as.integer(colSums(B_new * w)),
       gap = gap, mult = mult,
       phi = gap + mult - sqrt(gap^2 + mult^2))
}

#' Piecewise-linear OccBin solution for a known shock path (guess-and-verify)
#'
#' Iterates regime path -> time-varying rules -> constrained path ->
#' complementarity check until the regime path reproduces itself.  The shocks
#' of shock_seq are surprises; with a shock in the first period only this is
#' the exact perfect-foresight piecewise-linear solution (slack after the
#' horizon).
#' @param strict passed to .obc_pwl_rules(): FALSE returns converged = FALSE
#'   (paths NULL) when a guessed regime path has a singular system
#' @return list(regime_path, paths, states, phi, converged, n_iter)
#' @noRd
.obc_pwl_solve <- function(ctx, shock_seq, state_init = NULL, init = NULL,
                           max_iter = 50L, tol = 1e-8, strict = TRUE) {
  n_T    <- ncol(shock_seq)
  regime <- if (!is.null(init) && length(init) == n_T) as.integer(init)
            else integer(n_T)
  seen <- character(0)
  out  <- NULL
  for (iter in seq_len(max_iter)) {
    rules <- .obc_pwl_rules(ctx, regime, strict = strict)
    if (is.null(rules) && any(regime != 0L))      # singular guess (strict = FALSE)
      return(list(regime_path = regime, paths = NULL, states = NULL,
                  phi = NULL, converged = FALSE, n_iter = iter))
    sim   <- .obc_pwl_forward(ctx, rules, shock_seq, state_init)
    chk   <- .obc_pwl_check(ctx, rules, sim, shock_seq, regime, state_init,
                            tol)
    out   <- list(regime_path = regime, paths = sim$paths,
                  states = sim$states, phi = chk$phi, converged = FALSE,
                  n_iter = iter)
    if (identical(chk$regime_path, regime)) {
      out$converged <- TRUE
      return(out)
    }
    seen   <- c(seen, paste(regime, collapse = ","))
    regime <- chk$regime_path
    if (paste(regime, collapse = ",") %in% seen) break      # cycle
  }
  out
}

#' OccBin simulation of a sequence of SURPRISE shocks
#'
#' At period 1 and at every period with a non-zero shock the regime path
#' expected from then on (no further shocks, horizon ncol(shock_seq)) is
#' solved by inner(); its path is followed until the next shock arrives.
#' This is Dynare's occbin_solver treatment of a shock sequence.
#'
#' @param inner function(ctx, shock_seq, state_init, init) returning a list
#'   like .obc_pwl_solve()
#' @return list(regime_path, paths, states, phi, converged, n_iter); n_iter
#'   is the largest iteration count of any re-solve
#' @noRd
.obc_pwl_solve_surprise <- function(ctx, shock_seq, state_init = NULL,
                                    init = NULL, inner) {
  n_T    <- ncol(shock_seq)
  hit    <- which(colSums(abs(shock_seq)) > 0)
  starts <- sort(unique(c(1L, hit)))
  ends   <- c(starts[-1L] - 1L, n_T)

  paths  <- matrix(NA_real_, ctx$n_endo, n_T)
  states <- matrix(NA_real_, ctx$n_state, n_T)
  rownames(paths) <- ctx$endo_names
  regime <- integer(n_T)
  phi    <- matrix(0, ctx$n_spec, n_T)
  conv   <- TRUE
  n_iter <- 0L
  s <- if (!is.null(state_init) && length(state_init) == ctx$n_state)
         as.numeric(state_init) else numeric(ctx$n_state)
  if (!is.null(init) && length(init) != n_T) init <- NULL

  expected <- NULL          # regime path expected by the previous re-solve
  for (k in seq_along(starts)) {
    a <- starts[k]
    b <- ends[k]
    e <- matrix(0, ctx$n_exo, n_T)
    e[, 1L] <- shock_seq[, a]
    ## Initial guess: the previous re-solve's expectation from period a on;
    ## a supplied init overrides the period-a regime (its realised one).
    ini <- if (is.null(expected)) integer(n_T)
           else c(expected[(a - starts[k - 1L] + 1L):n_T],
                  integer(a - starts[k - 1L]))
    if (!is.null(init)) ini[1L] <- init[a]
    if (k == 1L && !is.null(init)) ini <- init
    res <- inner(ctx, e, s, ini)
    expected <- res$regime_path
    idx <- seq_len(b - a + 1L)
    paths[, a:b]  <- res$paths[, idx]
    states[, a:b] <- res$states[, idx]
    regime[a:b]   <- res$regime_path[idx]
    phi[, a:b]    <- res$phi[, idx]
    conv   <- conv && isTRUE(res$converged)
    n_iter <- max(n_iter, as.integer(res$n_iter))
    s <- states[, b]
  }
  list(regime_path = regime, paths = paths, states = states, phi = phi,
       converged = conv, n_iter = n_iter)
}
