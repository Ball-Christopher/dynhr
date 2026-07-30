## R/hank-egm2d.R
## --------------------------------------------------------------------------
## Two-asset household with DISCRETE ADJUSTMENT: each period the household
## either pays a fixed cost F_adj (in liquid units) to rebalance the illiquid
## account -- the ADJUST branch, which is exactly the smooth two-asset problem
## with cash-on-hand shifted by -F_adj -- or leaves it to roll over passively
## at its return (a' = (1+ra) a) -- the NO-ADJUST branch, a one-dimensional
## liquid savings problem per illiquid slice. The choice is smoothed by iid
## extreme-value taste shocks with scale sigma_taste (Iskhakov-Jorgensen-Rust-
## Schjerning), so policies and aggregates stay smooth in prices, which is
## what the fake-news Jacobian pipeline requires. sigma_taste -> 0 is the
## deterministic-inaction limit.
##
## WHY THIS EXISTS (briefs/19 section 9, F21/F22): the smooth convex cost
## structurally cannot generate empirical wealthy hand-to-mouth shares --
## stiffening it REDUCES inaction, because a convex cost leaves small
## continuous rebalancing cheap. Inaction requires a non-convexity at zero
## adjustment; the fixed cost supplies it, and this discrete-choice
## formulation is the standard discrete-time implementation (a literal
## kinked-FOC solver would break the crossing search's monotone marginal
## cost and is numerically dominated).
##
## RECURSION. State (e, b, a); V denotes the value LEVEL, (Vb, Va) the
## marginals. The level is the genuinely new object relative to the smooth
## solver -- the discrete choice needs levels to compare branches:
##   V^A(e,b,a): the smooth two-stage EGM at income y - F_adj (policies AND
##               level u(c^A) + beta E V'(e', b'^A, a'^A));
##   V^N(e,b,a): a' = (1+ra) a fixed; one-asset EGM in b per a-slice with
##               continuations interpolated along a at (1+ra) a;
##   V = sigma * logsumexp(V^A / sigma, V^N / sigma),
##   P = P(adjust) = logit((V^A - V^N) / sigma).
## Under EV taste shocks the envelope across the discrete choice is exact:
##   dV/dstate = P dV^A/dstate + (1-P) dV^N/dstate  (no dP term),
## so Vb = P Vb^A + (1-P) Vb^N and likewise Va, with
##   Vb^N = (1+rb) u'(c^N)  and  Va^N = (1+ra) * Wa(e, b'^N, (1+ra)a)
## (the no-adjust illiquid envelope compounds through the passive rollover).
##
## REUSE: the adjust branch calls .hank_egm2_step (the validated smooth
## kernel) verbatim; the no-adjust branch calls .hank_egm_step (the validated
## one-asset kernel) per slice. Both reductions are therefore also the
## oracles: F_adj = 0, sigma -> 0 must reproduce hank_egm2_solve, and
## F_adj -> large with ra = 0 must reproduce hank_egm_solve per slice.
##
## The tiny consumption floor on the ADJUST branch can bind at cells where
## paying F_adj is infeasible (poorest cells at large F_adj). Unlike the
## smooth solver -- where a binding floor at convergence means the returned
## policy is wrong (F10) -- here it is BENIGN BY DESIGN: those cells' V^A is
## astronomically negative, the logit sends P to 0, and the no-adjust branch
## (always feasible) carries the cell. hank_egm2d_solve checks instead that
## no cell with MATERIAL adjust probability sits on the floor.
##
## R reference implementation only in this wave; C++ kernels follow the
## established pattern once the R oracles are green.
## --------------------------------------------------------------------------


#' CRRA utility level
#'
#' \eqn{u(c) = c^{1-1/eis}/(1-1/eis)}, with the \eqn{eis = 1} log limit.
#' Levels are needed (not just marginals) because the discrete choice
#' compares branch VALUES.
#' @keywords internal
.hank_crra_u <- function(c, eis) {
  if (abs(eis - 1) < 1e-12) log(c)
  else c^(1 - 1 / eis) / (1 - 1 / eis)
}


#' Interpolation weights along the illiquid grid at one off-grid point
#'
#' Bracketing index and lower weight for the scalar query \code{aq} in
#' \code{a_grid}, extrapolating linearly above the top (unclamped weight,
#' clamped index -- the \code{.hank_interp1} convention). Used to evaluate
#' continuation objects at the passive rollover point \code{(1+ra) a}.
#' @keywords internal
.hank2d_a_coord <- function(a_grid, aq) {
  n_a <- length(a_grid)
  i <- min(max(findInterval(aq, a_grid), 1L), n_a - 1L)
  list(i = i, p = (a_grid[i + 1L] - aq) / (a_grid[i + 1L] - a_grid[i]))
}


#' One backward step of the discrete-adjustment two-asset household
#'
#' See the file header for the recursion. Inputs are next-period LEVEL
#' \code{V_p} and marginals \code{Vb_p, Va_p} (each \code{n_e x n_b x n_a});
#' returns the updated triple plus both branches' policies and the adjust
#' probability \code{P}.
#'
#' @keywords internal
.hank_egm2d_step <- function(V_p, Vb_p, Va_p, b_grid, a_grid, k_grid, y,
                             rb, ra, beta, eis, chi0, chi1, chi2,
                             F_adj, sigma_taste, Pi, Psi1_grid,
                             phi_contrib = 0) {
  n_e <- length(y); n_b <- length(b_grid); n_a <- length(a_grid)

  ## Discounted expectations shared by both branches (contract Pi once).
  EV <- beta * (Pi %*% matrix(V_p,  n_e, n_b * n_a))
  Wa <- beta * (Pi %*% matrix(Va_p, n_e, n_b * n_a))
  dim(EV) <- dim(Wa) <- c(n_e, n_b, n_a)

  ## ---- ADJUST branch: the validated smooth kernel at y - F_adj ------------
  adj <- .hank_egm2_step(Vb_p, Va_p, b_grid, a_grid, k_grid, y - F_adj,
                         rb, ra, beta, eis, chi0, chi1, chi2, Pi,
                         Psi1_grid = Psi1_grid)
  ## level: u(c^A) + beta E V'(e, b'^A, a'^A)
  V_A <- array(0, c(n_e, n_b, n_a))
  for (e in seq_len(n_e))
    V_A[e, , ] <- .hank_crra_u(adj$c[e, , ], eis) +
      .hank2_bilin(matrix(EV[e, , ], n_b, n_a), b_grid, a_grid,
                   as.numeric(adj$b[e, , ]), as.numeric(adj$a[e, , ]))

  ## ---- NO-ADJUST branch: one-asset EGM in b per illiquid slice ------------
  ## NOT .hank_egm_step: the mixed discrete-choice envelope is a smoothed max
  ## of two concave functions, hence NOT concave, so the branch's endogenous
  ## liquid grid can be non-monotone in transients -- .hank_egm_step's
  ## findInterval ERRORS there (observed at n = 32), while the F8 forward
  ## sweep (.hank2_interp_coord_rows) degrades gracefully exactly as the
  ## smooth two-asset kernel does. Same EGM algebra, tolerant inversion.
  ##
  ## KIWISAVER CONTRIBUTION (briefs/19 sections 9.5-9.7): a fraction
  ## phi_contrib of labor income flows into the locked account regardless of
  ## the adjustment choice -- the F26 participation-trap fix, and the NZ
  ## default-enrolment institution. The no-adjust rollover becomes
  ## a' = (1+ra) a + phi*y_e (now e-DEPENDENT) and the liquid budget keeps
  ## (1-phi) y_e. The ADJUST branch is deliberately unchanged: with the
  ## account open, a mandated contribution is fungible with the household's
  ## own portfolio choice, so it would have no effect there.
  tiny <- tiny_floor()
  y_liq <- (1 - phi_contrib) * y
  b_N <- array(0, c(n_e, n_b, n_a)); c_N <- b_N
  V_N <- b_N; Vb_N <- b_N; Va_N <- b_N; a_N <- b_N
  for (k in seq_len(n_a)) {
    for (e in seq_len(n_e)) {
      ap_ek <- (1 + ra) * a_grid[k] + phi_contrib * y[e]   # rollover, off-grid
      co <- .hank2d_a_coord(a_grid, ap_ek)
      ## continuation objects on the b-grid at a' = ap_ek. The rollover point
      ## is the DECIDER's (today's e enters through the contribution), but the
      ## continuation is still averaged over TOMORROW's e' by Pi.
      Vb_sl <- co$p * Vb_p[, , co$i] + (1 - co$p) * Vb_p[, , co$i + 1L]
      EV_sl <- co$p * EV[e, , co$i] + (1 - co$p) * EV[e, , co$i + 1L]
      Wa_sl <- co$p * Wa[e, , co$i] + (1 - co$p) * Wa[e, , co$i + 1L]
      Wb_ek <- pmax(beta * as.numeric(Pi[e, ] %*% Vb_sl), tiny)   # over b'
      c_endo <- Wb_ek^(-eis)
      b_endo <- (c_endo + b_grid - y_liq[e]) / (1 + rb)
      cb <- .hank2_interp_coord_rows(matrix(b_endo, 1L), b_grid)
      bp <- pmax(as.numeric(.hank2_apply_coord_vec(cb$i, cb$p, b_grid)),
                 b_grid[1L])
      cN <- pmax((1 + rb) * b_grid + y_liq[e] - bp, tiny)   # exact budget
      b_N[e, , k]  <- bp
      c_N[e, , k]  <- cN
      a_N[e, , k]  <- ap_ek
      Vb_N[e, , k] <- (1 + rb) * cN^(-1 / eis)
      ib <- pmin(pmax(findInterval(bp, b_grid), 1L), n_b - 1L)
      pb <- (b_grid[ib + 1L] - bp) / (b_grid[ib + 1L] - b_grid[ib])
      V_N[e, , k]  <- .hank_crra_u(cN, eis) +
        pb * EV_sl[ib] + (1 - pb) * EV_sl[ib + 1L]
      Va_N[e, , k] <- (1 + ra) * (pb * Wa_sl[ib] + (1 - pb) * Wa_sl[ib + 1L])
    }
  }

  ## ---- discrete choice with EV taste shocks -------------------------------
  ## Overflow-safe logit/logsumexp: dV can be ~1e12 where the adjust branch
  ## is infeasible (floored c), which is exactly the P -> 0 limit.
  dV <- (V_A - V_N) / sigma_taste
  P  <- stats::plogis(dV)
  m  <- pmax(V_A, V_N)
  V  <- m + sigma_taste * log(exp((V_A - m) / sigma_taste) +
                                exp((V_N - m) / sigma_taste))
  Vb <- P * adj$Vb + (1 - P) * Vb_N
  Va <- P * adj$Va + (1 - P) * Va_N

  list(V = V, Vb = Vb, Va = Va, P = P,
       b_A = adj$b, a_A = adj$a, c_A = adj$c, chi_A = adj$chi,
       b_N = b_N, a_N = a_N, c_N = c_N)
}


#' Solve the discrete-adjustment two-asset household to a stationary policy
#'
#' Validation record (briefs/19 sections 9.4-9.8): per-branch budget
#' identities at machine precision; the F -> large / ra = 0 reduction onto
#' \code{hank_egm_solve} per slice (8e-11); weak dominance over a brute-force
#' node-restricted VFI at every cell; the comparative Bellman gate (no worse
#' than the fully-validated smooth solver on the identical instrument -- the
#' original F24 "upper envelope" blocker was RETRACTED after the smooth
#' control failed that instrument identically, F25); the F26 participation
#' trap measured and closed by \code{phi_contrib}; and the fake-news Jacobian
#' of the block built on this solver validated against brute-force ND
#' (\code{test-hank-het2d-block.R}).  C++ kernels are a documented follow-up;
#' this is the R reference implementation.
#'
#' The fixed point of \code{\link{.hank_egm2d_step}}: a household that pays a
#' fixed cost \code{F_adj} to rebalance its illiquid account, with the choice
#' smoothed by extreme-value taste shocks of scale \code{sigma_taste}.  See
#' the file header of \code{R/hank-egm2d.R} for the recursion and
#' \code{briefs/19-twoasset-hank-scope.md} section 9 for why this variant
#' exists: the smooth convex cost cannot generate empirical wealthy
#' hand-to-mouth shares (stiffening it REDUCES inaction), while the fixed
#' cost's non-convexity at zero adjustment can.
#'
#' Nests the smooth household: at \code{F_adj = 0} the adjust branch weakly
#' dominates (its choice set contains the no-adjust point), so as
#' \code{sigma_taste} shrinks the solution converges to
#' \code{\link{hank_egm2_solve}} -- one of the two reduction oracles in
#' \code{test-hank-egm2d.R}.  At large \code{F_adj} adjustment never happens
#' and each illiquid slice is a one-asset problem
#' (\code{\link{hank_egm_solve}}), the other reduction.
#'
#' @param b_grid,a_grid,y,rb,ra,beta,eis,chi0,chi1,chi2,Pi,k_grid,n_k,k_max
#'   As in \code{\link{hank_egm2_solve}} (the convex cost is retained on the
#'   adjust branch for interior optimality).
#' @param F_adj Fixed adjustment cost, in liquid (cash-on-hand) units;
#'   \code{>= 0}.
#' @param sigma_taste Extreme-value taste-shock scale \code{> 0}. Smaller is
#'   closer to deterministic inaction but kinks the aggregates; the fake-news
#'   Jacobian needs it comfortably positive (reference calibrations use
#'   1e-3..1e-1 in consumption units).
#' @param phi_contrib Fraction of labour income paid DIRECTLY into the locked
#'   (illiquid) account, bypassing cash-on-hand: liquid income is
#'   \code{(1 - phi_contrib) * y} and the remainder arrives as an illiquid
#'   contribution regardless of the adjust/no-adjust decision. \code{0} (the
#'   default) restores the contribution-free household exactly.
#' @param tol,maxit Convergence tolerance on policies and the adjust
#'   probability, and the iteration cap.
#'
#' @return A list with the mixed envelopes \code{V}, \code{Vb}, \code{Va},
#'   the adjust probability \code{P}, per-branch policies (\code{b_A},
#'   \code{a_A}, \code{c_A}, \code{chi_A}, \code{b_N}, \code{a_N},
#'   \code{c_N}), \code{iterations}, \code{converged}, and the echoed
#'   calibration.
#' @seealso \code{\link{hank_egm2_solve}} (smooth), \code{\link{hank_egm_solve}}
#'   (one-asset), \code{\link{hank_het2d_block}}
#' @export
hank_egm2d_solve <- function(b_grid, a_grid, y, rb, ra, beta, eis,
                             chi0, chi1, chi2, F_adj, sigma_taste, Pi,
                             k_grid = NULL, n_k = 50L, k_max = 1,
                             tol = 1e-9, maxit = 5000L, phi_contrib = 0) {
  ## Reuse the smooth wrapper's full validation by SOLVING the F = 0 smooth
  ## problem first: it validates every shared argument identically for both
  ## solvers AND supplies consistent (Vb, Va) starting values. The extra cost
  ## is one smooth solve, which is cheap next to the discrete iteration.
  smooth <- hank_egm2_solve(b_grid, a_grid, y = y, rb = rb, ra = ra,
                            beta = beta, eis = eis, chi0 = chi0, chi1 = chi1,
                            chi2 = chi2, Pi = Pi, k_grid = k_grid, n_k = n_k,
                            k_max = k_max)
  if (!is.numeric(F_adj) || length(F_adj) != 1L || !is.finite(F_adj) ||
      F_adj < 0)
    stop("hank_egm2d_solve: 'F_adj' must be a finite scalar >= 0.")
  if (!is.numeric(phi_contrib) || length(phi_contrib) != 1L ||
      !is.finite(phi_contrib) || phi_contrib < 0 || phi_contrib >= 1)
    stop("hank_egm2d_solve: 'phi_contrib' must be a scalar in [0, 1) (the ",
         "mandatory contribution rate out of labor income; 0 disables it).")
  if (!is.numeric(sigma_taste) || length(sigma_taste) != 1L ||
      !is.finite(sigma_taste) || sigma_taste <= 0)
    stop("hank_egm2d_solve: 'sigma_taste' must be a finite scalar > 0 (the ",
         "sigma_taste -> 0 deterministic limit kinks the aggregates; use a ",
         "small positive value).")
  n_e <- length(y); n_b <- length(b_grid); n_a <- length(a_grid)
  k_grid <- smooth$k_grid
  Psi1_grid <- .hank_psi(matrix(a_grid, n_a, n_a),
                         matrix(a_grid, n_a, n_a, byrow = TRUE),
                         ra, chi0, chi1, chi2)$Psi1

  ## Levels initialised from the smooth solution's consumption at the
  ## stationary-consumption bound u(c)/(1-beta); marginals from the smooth
  ## fixed point.
  V  <- .hank_crra_u(pmax(smooth$c, tiny_floor()), eis) / (1 - beta)
  Vb <- smooth$Vb; Va <- smooth$Va

  converged <- FALSE; it <- 0L; prev <- NULL; step <- NULL
  for (it in seq_len(as.integer(maxit))) {
    step <- .hank_egm2d_step(V, Vb, Va, b_grid, a_grid, k_grid, y, rb, ra,
                             beta, eis, chi0, chi1, chi2, F_adj, sigma_taste,
                             Pi, Psi1_grid, phi_contrib = phi_contrib)
    V <- step$V; Vb <- step$Vb; Va <- step$Va
    cur <- c(step$b_A, step$a_A, step$b_N, step$P)
    if (!is.null(prev)) {
      gap <- max(abs(cur - prev))
      if (!is.finite(gap)) break
      if (gap < tol) { converged <- TRUE; break }
    }
    prev <- cur
  }

  if (!all(is.finite(step$P)) || !all(is.finite(step$b_N)) ||
      !all(is.finite(step$c_N)))
    stop("hank_egm2d_solve: non-finite iterate after ", it, " iterations.")
  ## The adjust branch's tiny consumption floor is benign ONLY where the
  ## choice kills the branch (see the file header). A floored cell that
  ## retains material adjust probability means the discretization is genuinely
  ## inadequate -- same failure class as the smooth solver's F10.
  floored_active <- step$c_A <= tiny_floor() & step$P > 1e-6
  if (any(floored_active))
    stop("hank_egm2d_solve: ", sum(floored_active), " cell(s) have floored ",
         "adjust-branch consumption but non-negligible adjust probability. ",
         "The calibration/grids cannot represent the adjust decision there; ",
         "reduce F_adj, refine the grids, or raise sigma_taste.")

  c(step, list(iterations = it, converged = converged,
               b_grid = b_grid, a_grid = a_grid, k_grid = k_grid, y = y,
               rb = rb, ra = ra, beta = beta, eis = eis,
               chi0 = chi0, chi1 = chi1, chi2 = chi2,
               F_adj = F_adj, sigma_taste = sigma_taste,
               phi_contrib = phi_contrib, Pi = Pi))
}
