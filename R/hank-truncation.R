## R/hank-truncation.R
## --------------------------------------------------------------------------
## Finite-state, Kalman-filterable HANK via a SCALAR-ANCHORED COARSE asset
## grid (Tier 18 flagship, milestone m5; scratch derivation + verification in
## .claude/orchestration/truncation/, sessions 2026-07-05..08).
##
## Construction (m4b, orchestrator-verified):
##   1. Solve the model on a FINE grid (the truth): `hank_ks_steady`.
##   2. Re-solve the household block on a COARSE grid (n_a ~ 16-40) at the
##      SAME prices, with ONE scalar discount wedge beta_eff = beta*exp(lx)
##      chosen so the coarse aggregate assets equal the fine-model target
##      (`hank_coarse_anchored`). This makes the coarse steady state EXACT in
##      the aggregate (A, and C via the aggregate budget) and — the
##      quantitative finding — also collapses the coarse grid's DYNAMICS
##      errors by 1-2 orders of magnitude (persistence eigenvalue, Jacobian
##      impact), because the stationary asset level and the forward
##      operator's second eigenvalue are both controlled by effective
##      patience beta_eff*(1+r).
##   3. Assemble the finite-dimensional Reiter (2009) system around that
##      anchored steady state and solve it by Klein (2000) QZ into a linear
##      state space x_{t+1} = T x_t + R eps, y_t = Z x_t
##      (`hank_reiter_statespace`), on which dynhr's Kalman filter runs
##      directly (`hank_reiter_kalman_loglik`).
##
## Numerical facts to keep in mind (from the m4b verification):
##   - The structural lead matrix A of the Reiter system is SINGULAR:
##     borrowing-constrained nodes have exactly-zero G_V rows (their dVa is a
##     static function of prices), and end-of-period nodes that no interior
##     policy interpolates through give zero G_V columns (dangling jumps).
##     Klein/QZ absorbs both as infinite generalized eigenvalues; a plain
##     eigen(solve(A, B)) Blanchard-Kahn solve FAILS. Do not "simplify" back.
##   - Mass conservation gives a unit generalized eigenvalue whose eigenvector
##     has sum(dD) != 0; the dynamics live on the invariant sum(dD) = 0
##     subspace. `drop_dist_coord = TRUE` (default) removes the redundant
##     distribution coordinate so T is strictly stable and the KF can use the
##     stationary (Lyapunov) initial covariance.
##   - All vectorized household objects use the package distribution order
##     .hank_mat_to_vec(X) = as.numeric(t(X)) (asset index fastest).
## --------------------------------------------------------------------------


#' Scalar-anchored coarse-grid household block
#'
#' Solves the household block on \code{a_grid} at fixed prices \code{(r, w)}
#' with a single scalar discount wedge \code{beta_eff = beta * exp(lx)},
#' root-finding \code{lx} so the block's aggregate end-of-period assets equal
#' \code{A_target} (typically the FINE-grid model's aggregate at the same
#' prices). See the file header for why this one scalar also repairs the
#' coarse grid's aggregate dynamics.
#'
#' @param a_grid Coarse asset grid (see \code{\link{hank_asset_grid}}).
#' @param Pi,e Income transition matrix and levels.
#' @param beta,eis Baseline discount factor and EIS.
#' @param r,w Prices at which to anchor (the fine model's steady state).
#' @param A_target Aggregate assets the anchored block must reproduce.
#' @param lx_bracket Initial search bracket for \code{lx} (expanded downward
#'   automatically; the upper end is capped so \code{beta_eff * (1 + r) < 1}).
#' @param tol Root-finding tolerance on \code{lx}.
#' @param Pi_fn,Pi_inputs Optional transition-matrix constructor and its
#'   steady-state inputs (see \code{\link{hank_het_block}}), carried onto the
#'   anchored block so the HANK+SAM endogenous-Pi machinery
#'   (\code{\link{hank_het_jacobian}} transition-input columns,
#'   \code{\link{hank_sam_reiter_statespace}}) works on the coarse block.
#'
#' @return A list with the anchored \code{block} (a
#'   \code{\link{hank_het_block}} solved at \code{beta_eff}), \code{beta_eff},
#'   \code{lx}, and \code{A_err = block$A - A_target}.
#' @export
hank_coarse_anchored <- function(a_grid, Pi, e, beta, eis, r, w, A_target,
                                 lx_bracket = c(-0.05, 0.005), tol = 1e-12,
                                 Pi_fn = NULL, Pi_inputs = NULL) {
  mk <- function(lx) suppressWarnings(
    hank_het_block(a_grid, Pi, e, beta = beta * exp(lx), eis = eis,
                   r = r, w = w, Pi_fn = Pi_fn, Pi_inputs = Pi_inputs))
  f <- function(lx) mk(lx)$A - A_target

  ## cap the upper end so beta_eff * (1 + r) < 1 (EGM diverges otherwise)
  lx_hi <- min(lx_bracket[2],
               if (r > -1) log(1 / (beta * (1 + r))) - 1e-6 else lx_bracket[2])
  lx_lo <- min(lx_bracket[1], lx_hi - 1e-3)
  f_hi <- f(lx_hi)
  if (f_hi < 0)
    stop(sprintf(paste0("hank_coarse_anchored(): A_target %.6g unreachable -- ",
                        "even beta_eff*(1+r) -> 1 gives A = %.6g. ",
                        "Coarse grid amax too small?"), A_target, f_hi + A_target))
  ## expand the lower end until the bracket straddles the root
  f_lo <- f(lx_lo)
  n_expand <- 0L
  while (f_lo > 0 && n_expand < 12L) {
    lx_lo <- lx_lo * 2 - 0.01
    f_lo <- f(lx_lo)
    n_expand <- n_expand + 1L
  }
  if (f_lo > 0)
    stop("hank_coarse_anchored(): could not bracket the anchoring root.")

  lx <- stats::uniroot(f, c(lx_lo, lx_hi), f.lower = f_lo, f.upper = f_hi,
                       tol = tol)$root
  blk <- mk(lx)
  list(block = blk, beta_eff = beta * exp(lx), lx = lx,
       A_err = blk$A - A_target)
}


#' Anchored coarse-grid version of a Krusell-Smith steady state
#'
#' Given a (fine-grid) \code{\link{hank_ks_steady}} solution, builds the
#' scalar-anchored coarse-grid economy at the SAME prices and capital stock.
#' Because the anchor reproduces \code{A = K} exactly, the coarse economy
#' clears the capital market at the fine model's \code{r}, and the returned
#' object is a fully valid \code{hank_ks}: all downstream GE machinery
#' (\code{\link{hank_ks_ge_jacobian}}, \code{\link{hank_ks_linear_irf}},
#' \code{\link{hank_ks_nonlinear_irf}}) applies unchanged.
#'
#' @param ks A \code{hank_ks} object (the fine-grid truth).
#' @param n_a Number of coarse asset gridpoints (16-40 is the useful range;
#'   ignored if \code{a_grid} is given).
#' @param a_grid Optional explicit coarse grid; default
#'   \code{hank_asset_grid(amax = max(fine grid), n = n_a, amin = min(fine grid))}.
#'
#' @return A \code{hank_ks} object for the anchored coarse economy, with extra
#'   fields \code{beta_base} (the unanchored beta), \code{anchor_lx}, and
#'   \code{n_a_coarse}. Its \code{beta} field IS \code{beta_eff}.
#'
#' @section Per-draw use in estimation (frozen-lx recipe):
#' When a posterior draw moves HOUSEHOLD parameters (not just \code{rho_z},
#' \code{sigma_z}), do NOT re-solve the fine model per draw. Calibrate
#' \code{lx} ONCE here at a reference \eqn{\theta}, then per draw solve the
#' coarse GE directly at the wedged discount factor -- i.e.
#' \code{hank_ks_steady(a_grid_coarse, ..., beta = beta_draw * exp(lx_ref))}
#' -- and feed that to \code{\link{hank_reiter_linearize}}. Verified
#' 2026-07-08 (\code{.claude/orchestration/truncation/m5_reanchor.R}, 3x3
#' grid, beta +-0.005 x EIS +-25\%, n_a = 24 vs fine n_a = 500): \code{lx} is
#' essentially INVARIANT in beta (drift < 1e-6 over +-0.005, so the frozen
#' wedge is SS-exact in the beta direction) and mildly linear in EIS
#' (d lx / d EIS ~ 8.7e-4); the frozen wedge beats the unanchored coarse GE
#' ~4x on every steady-state metric and matches the per-draw ideal re-anchor
#' on the persistence eigenvalue (whose residual error is the irreducible
#' coarse-grid error, ~1.5e-4 at n_a = 24).
#' @export
hank_ks_coarse_anchored <- function(ks, n_a = 24L, a_grid = NULL) {
  if (!inherits(ks, "hank_ks"))
    stop("hank_ks_coarse_anchored(): `ks` must be a hank_ks object.")
  fine_grid <- ks$block$a_grid
  if (is.null(a_grid))
    a_grid <- hank_asset_grid(amax = max(fine_grid), n = as.integer(n_a),
                              amin = min(fine_grid))
  anch <- hank_coarse_anchored(a_grid, ks$block$Pi, ks$block$e,
                               beta = ks$beta, eis = ks$eis,
                               r = ks$r, w = ks$w, A_target = ks$K,
                               Pi_fn = ks$block$Pi_fn,
                               Pi_inputs = ks$block$Pi_inputs)
  structure(
    list(r = ks$r, w = ks$w, K = ks$K, Z = ks$Z,
         alpha = ks$alpha, delta = ks$delta,
         beta = anch$beta_eff, eis = ks$eis, block = anch$block,
         mkt_residual = anch$block$A - ks$K,
         beta_base = ks$beta, anchor_lx = anch$lx,
         n_a_coarse = length(a_grid)),
    class = "hank_ks")
}


#' Anchored coarse-grid version of a TAXED Krusell-Smith steady state
#'
#' Additive taxed-household counterpart of \code{\link{hank_ks_coarse_anchored}},
#' added to support re-anchoring the finite-HANK capital-tax instrument at a
#' RAMSEY steady state (see R/hank-finite-ramsey-ss.R). Given a
#' \code{hank_ks_taxed} object (\code{hank_ks_taxed_steady()} output -- a
#' household solved under a linear capital tax \code{tau} with balanced-budget
#' lump-sum rebate \code{tr = tau*r*K}), builds the scalar-anchored coarse-grid
#' economy at the SAME (tau, r, w, K): one scalar discount wedge
#' \code{beta_eff = beta*exp(lx)} is root-found so the coarse TAXED household
#' block's aggregate assets reproduce the fine model's \code{K}, mirroring
#' \code{\link{hank_coarse_anchored}} exactly except that the household solve
#' underneath is the taxed EGM problem
#' (net return \code{(1-tau)*r}, income \code{w*e+tr}) rather than the
#' untaxed one. At \code{tau = 0} this reproduces
#' \code{\link{hank_ks_coarse_anchored}} to machine precision, since the
#' taxed and untaxed household blocks coincide there.
#'
#' @param ks_taxed A \code{hank_ks_taxed} object (see
#'   \code{\link{hank_ks_taxed_steady}}).
#' @param n_a Number of coarse asset gridpoints (as
#'   \code{\link{hank_ks_coarse_anchored}}).
#' @param a_grid Optional explicit coarse grid; default as in
#'   \code{\link{hank_ks_coarse_anchored}}.
#' @param lx_bracket,tol Passed to the taxed anchoring root-find.
#'
#' @return A \code{hank_ks_taxed}/\code{hank_ks} object (both classes; see
#'   \code{\link{hank_finite_mod}}'s taxed-initval hook, which reads
#'   \code{$tau}/\code{$tr} when present) for the anchored coarse taxed
#'   economy, with the same extra fields as
#'   \code{\link{hank_ks_coarse_anchored}} (\code{beta_base}, \code{anchor_lx},
#'   \code{n_a_coarse}), plus \code{tau} and \code{tr} carried through from
#'   \code{ks_taxed}.
#' @export
hank_ks_taxed_coarse_anchored <- function(ks_taxed, n_a = 24L, a_grid = NULL,
                                          lx_bracket = c(-0.05, 0.005),
                                          tol = 1e-12) {
  if (!inherits(ks_taxed, "hank_ks_taxed"))
    stop("hank_ks_taxed_coarse_anchored(): `ks_taxed` must be a hank_ks_taxed object.")
  fine_grid <- ks_taxed$block$a_grid
  if (is.null(a_grid))
    a_grid <- hank_asset_grid(amax = max(fine_grid), n = as.integer(n_a),
                              amin = min(fine_grid))
  tau <- ks_taxed$tau; r <- ks_taxed$r; w <- ks_taxed$w; K <- ks_taxed$K
  beta <- ks_taxed$beta; eis <- ks_taxed$eis
  Pi <- ks_taxed$block$Pi; e <- ks_taxed$block$e

  mk <- function(lx) suppressWarnings(
    .hank_taxed_het_block(a_grid, Pi, e, beta = beta * exp(lx), eis = eis,
                          tau = tau, r = r, w = w, K = K))
  f <- function(lx) mk(lx)$A - K

  r_net <- (1 - tau) * r
  lx_hi <- min(lx_bracket[2],
               if (r_net > -1) log(1 / (beta * (1 + r_net))) - 1e-6 else lx_bracket[2])
  lx_lo <- min(lx_bracket[1], lx_hi - 1e-3)
  f_hi <- f(lx_hi)
  if (f_hi < 0)
    stop(sprintf(paste0("hank_ks_taxed_coarse_anchored(): A_target %.6g unreachable -- ",
                        "even beta_eff*(1+r_net) -> 1 gives A = %.6g. ",
                        "Coarse grid amax too small?"), K, f_hi + K))
  f_lo <- f(lx_lo)
  n_expand <- 0L
  while (f_lo > 0 && n_expand < 12L) {
    lx_lo <- lx_lo * 2 - 0.01
    f_lo <- f(lx_lo)
    n_expand <- n_expand + 1L
  }
  if (f_lo > 0)
    stop("hank_ks_taxed_coarse_anchored(): could not bracket the anchoring root.")

  lx <- stats::uniroot(f, c(lx_lo, lx_hi), f.lower = f_lo, f.upper = f_hi,
                       tol = tol)$root
  blk <- mk(lx)

  structure(
    list(r = r, w = w, K = K, tau = tau, tr = blk$tr, Z = ks_taxed$Z,
         alpha = ks_taxed$alpha, delta = ks_taxed$delta,
         beta = beta * exp(lx), eis = eis, block = blk,
         mkt_residual = blk$A - K,
         beta_base = beta, anchor_lx = lx,
         n_a_coarse = length(a_grid)),
    class = c("hank_ks_taxed", "hank_ks"))
}


#' Reiter (2009) linear state space of a Krusell-Smith economy
#'
#' FD-linearizes the household block (one EGM backward step + the Young
#' forward operator) and the firm block around the steady state of \code{ks},
#' and solves the resulting rational-expectations system
#' \eqn{A x_{t+1} = B x_t}, \eqn{x = (dD, dK, dz, dVa)}, by Klein (2000)
#' generalized-Schur decomposition into
#' \deqn{s_{t+1} = T s_t + R \epsilon_{t+1}, \qquad y_t = Z s_t,}
#' with predetermined state \eqn{s = (dD, dK, dz)} and TFP following
#' \eqn{dz_{t+1} = \rho_z dz_t + \epsilon_{t+1}} (log deviation). Observation
#' rows are provided for \code{A}, \code{C} (household aggregates, dated by
#' the beginning-of-period distribution), \code{K} (capital in use at t),
#' \code{r}, \code{w}, \code{z}, and \code{Y} (output).
#'
#' Meant for the ANCHORED coarse economy (\code{\link{hank_ks_coarse_anchored}});
#' it works on any \code{hank_ks}, but the state dimension is
#' \code{n_e * n_a + 1} so keep the grid coarse.
#'
#' @param x A \code{hank_ks} object, or a precomputed
#'   \code{\link{hank_reiter_linearize}} result. Pass the latter when
#'   re-solving many times (e.g. per posterior draw): only the Klein solve
#'   depends on \code{rho_z}; the FD linearization does not.
#' @param rho_z TFP log-AR(1) persistence.
#' @param sigma_z TFP innovation standard deviation (stored for the KF).
#' @param drop_dist_coord Drop the redundant last distribution coordinate
#'   (mass conservation makes \code{sum(dD) = 0} invariant), yielding a
#'   strictly stable \code{T}. Default \code{TRUE}; keep it on for filtering.
#' @param delta_fd Relative finite-difference step for the linearization
#'   (ignored when \code{x} is already a linearization).
#'
#' @return An object of class \code{hank_reiter_ss}: \code{T_mat},
#'   \code{R_mat}, \code{Z_mat} (named rows), \code{state_names},
#'   \code{spectral_radius}, \code{n_stable}, \code{n_infinite},
#'   \code{gen_eig_mod}, \code{Qj} (jump loading \code{dVa_t = Qj s_t} in the
#'   FULL coordinates), \code{rho_z}, \code{sigma_z}, and the calibration.
#' @export
hank_reiter_statespace <- function(x, rho_z, sigma_z = 0.01,
                                   drop_dist_coord = TRUE,
                                   delta_fd = 1e-6) {
  lin <- if (inherits(x, "hank_reiter_lin")) x
         else hank_reiter_linearize(x, delta_fd = delta_fd)
  .hank_reiter_assemble(lin, rho_z = rho_z, sigma_z = sigma_z,
                        drop_dist_coord = drop_dist_coord)
}


#' FD linearization of the household and firm blocks around a hank_ks steady state
#'
#' The expensive, \code{rho_z}-independent half of
#' \code{\link{hank_reiter_statespace}}: finite-difference derivatives of the
#' EGM backward step (w.r.t. next-period marginal value and prices), of the
#' Young forward push (w.r.t. the savings policy), and the analytic firm-block
#' price derivatives. Compute once and pass to
#' \code{hank_reiter_statespace(lin, rho_z = ...)} inside estimation loops.
#'
#' @param ks A \code{hank_ks} object.
#' @param delta_fd Relative finite-difference step.
#' @return An object of class \code{hank_reiter_lin} (all linearization
#'   matrices plus the calibration; see the source for fields).
#' @export
hank_reiter_linearize <- function(ks, delta_fd = 1e-6) {
  if (!inherits(ks, "hank_ks"))
    stop("hank_reiter_linearize(): `ks` must be a hank_ks object.")
  blk <- ks$block
  n_e <- blk$n_e; n_a <- blk$n_a; n <- n_e * n_a
  Pi <- blk$Pi; e <- blk$e
  r <- ks$r; w <- ks$w; K_ss <- ks$K
  alpha <- ks$alpha; delta <- ks$delta; Z_ss <- ks$Z

  ## stationary labour supply + firm-FOC consistency (fail loud on a bad ks)
  s_e <- { v <- Re(eigen(t(Pi))$vectors[, 1]); v / sum(v) }
  L <- sum(s_e * e)
  r_foc <- alpha * Z_ss * K_ss^(alpha - 1) * L^(1 - alpha) - delta
  w_foc <- (1 - alpha) * Z_ss * K_ss^alpha * L^(-alpha)
  if (abs(r_foc - r) > 1e-6 * max(1, abs(r)) || abs(w_foc - w) > 1e-6 * w)
    stop(sprintf(paste0("hank_reiter_linearize(): firm FOCs inconsistent with ",
                        "(r, w, K, Z, L): r_foc = %.8f vs r = %.8f, ",
                        "w_foc = %.6f vs w = %.6f."), r_foc, r, w_foc, w))

  ## firm-block price derivatives (dr_t/dw_t wrt capital IN USE at t and dz_t)
  r_K <- alpha * (alpha - 1) * Z_ss * K_ss^(alpha - 2) * L^(1 - alpha)
  r_Z <- alpha * K_ss^(alpha - 1) * L^(1 - alpha) * Z_ss   # z is a LOG deviation
  w_K <- alpha * (1 - alpha) * Z_ss * K_ss^(alpha - 1) * L^(-alpha)
  w_Z <- (1 - alpha) * K_ss^alpha * L^(-alpha) * Z_ss

  vc   <- .hank_mat_to_vec
  unvc <- function(v) matrix(v, n_e, byrow = TRUE)
  Va_ss <- blk$Va; a_ss <- blk$a; c_ss <- blk$c; D_ss <- blk$D
  a_grid <- blk$a_grid; amin <- a_grid[1L]

  ## ---- FD linearization of the household step (G_V, ga_V, gc_V, price derivs) --
  G_V <- matrix(0, n, n); ga_V <- matrix(0, n, n); gc_V <- matrix(0, n, n)
  vVa <- vc(Va_ss)
  for (j in seq_len(n)) {
    h <- delta_fd * (1 + abs(vVa[j]))
    Vp <- vVa; Vp[j] <- Vp[j] + h
    Vm <- vVa; Vm[j] <- Vm[j] - h
    sp <- .hank_block_step(blk, unvc(Vp), r, w)
    sm <- .hank_block_step(blk, unvc(Vm), r, w)
    G_V[, j]  <- (vc(sp$Va) - vc(sm$Va)) / (2 * h)
    ga_V[, j] <- (vc(sp$a)  - vc(sm$a))  / (2 * h)
    gc_V[, j] <- (vc(sp$c)  - vc(sm$c))  / (2 * h)
  }
  hp <- delta_fd
  sp <- .hank_block_step(blk, Va_ss, r + hp, w)
  sm <- .hank_block_step(blk, Va_ss, r - hp, w)
  gVa_r <- (vc(sp$Va) - vc(sm$Va)) / (2 * hp)
  ga_r  <- (vc(sp$a)  - vc(sm$a))  / (2 * hp)
  gc_r  <- (vc(sp$c)  - vc(sm$c))  / (2 * hp)
  sp <- .hank_block_step(blk, Va_ss, r, w + hp)
  sm <- .hank_block_step(blk, Va_ss, r, w - hp)
  gVa_w <- (vc(sp$Va) - vc(sm$Va)) / (2 * hp)
  ga_w  <- (vc(sp$a)  - vc(sm$a))  / (2 * hp)
  gc_w  <- (vc(sp$c)  - vc(sm$c))  / (2 * hp)

  ## Da = d(Lambda(a)' D_ss)/da_k. One-sided UP step at constrained nodes:
  ## a = amin sits ON a Young-lottery kink and below it is infeasible.
  ## MATRIX-FREE (0.9.0.0026): this loop calls push() up to 2n times, and the
  ## sparse route built a whole n x n Lambda for a single matvec each time --
  ## quadratic in the grid size for a product that is linear in it. Measured at
  ## n_a = 60, n_e = 3 (installed -O2) the loop was 73% of
  ## hank_reiter_linearize; see briefs/21 section 12. .hank_forward_push()
  ## contracts the same identity instead of materializing it, so this is an
  ## algebraic reassociation (round-off-level agreement, asserted in
  ## test-hank-forward-push.R), not a different discretization.
  push <- function(A_mat) .hank_forward_push(A_mat, a_grid, Pi, D_ss)
  p0 <- push(a_ss)
  va <- vc(a_ss)
  Da <- matrix(0, n, n)
  for (k in seq_len(n)) {
    hk <- delta_fd * (1 + abs(va[k]))
    ap <- va; ap[k] <- ap[k] + hk
    if (va[k] <= amin + 1e-12) {
      Da[, k] <- (push(unvc(ap)) - p0) / hk
    } else {
      am <- va; am[k] <- am[k] - hk
      Da[, k] <- (push(unvc(ap)) - push(unvc(am))) / (2 * hk)
    }
  }

  structure(
    list(G_V = G_V, ga_V = ga_V, gc_V = gc_V,
         gVa_r = gVa_r, ga_r = ga_r, gc_r = gc_r,
         gVa_w = gVa_w, ga_w = ga_w, gc_w = gc_w,
         Da = Da, va = va,
         r_K = r_K, r_Z = r_Z, w_K = w_K, w_Z = w_Z,
         L = L, n = n, n_e = n_e, n_a = n_a, delta_fd = delta_fd, ks = ks),
    class = "hank_reiter_lin")
}


#' @keywords internal
.hank_reiter_assemble <- function(lin, rho_z, sigma_z, drop_dist_coord) {
  ks <- lin$ks; blk <- ks$block
  n_e <- lin$n_e; n_a <- lin$n_a; n <- lin$n; L <- lin$L
  G_V <- lin$G_V; ga_V <- lin$ga_V; gc_V <- lin$gc_V
  gVa_r <- lin$gVa_r; ga_r <- lin$ga_r; gc_r <- lin$gc_r
  gVa_w <- lin$gVa_w; ga_w <- lin$ga_w; gc_w <- lin$gc_w
  Da <- lin$Da; va <- lin$va
  r_K <- lin$r_K; r_Z <- lin$r_Z; w_K <- lin$w_K; w_Z <- lin$w_Z
  D_ss <- blk$D; c_ss <- blk$c
  alpha <- ks$alpha; Z_ss <- ks$Z; K_ss <- ks$K
  vc <- .hank_mat_to_vec

  ## ---- structural system A x_{t+1} = B x_t, x = [dD; dK; dz; dVa] --------------
  m  <- 2L * n + 2L
  iD <- 1:n; iK <- n + 1L; iz <- n + 2L; iV <- n + 2L + 1:n
  npred <- n + 2L
  A <- matrix(0, m, m); B <- matrix(0, m, m)
  A[iD, iD] <- diag(n)
  A[iD, iV] <- -Da %*% ga_V
  B[iD, iD] <- as.matrix(Matrix::t(blk$Lambda))
  B[iD, iK] <- Da %*% (ga_r * r_K + ga_w * w_K)
  B[iD, iz] <- Da %*% (ga_r * r_Z + ga_w * w_Z)
  A[iK, iK] <- 1
  A[iK, iV] <- -as.numeric(D_ss %*% ga_V)
  B[iK, iD] <- va
  B[iK, iK] <- sum(D_ss * (ga_r * r_K + ga_w * w_K))
  B[iK, iz] <- sum(D_ss * (ga_r * r_Z + ga_w * w_Z))
  A[iz, iz] <- 1
  B[iz, iz] <- rho_z
  A[iV, iV] <- G_V
  B[iV, iV] <- diag(n)
  B[iV, iK] <- -(gVa_r * r_K + gVa_w * w_K)
  B[iV, iz] <- -(gVa_r * r_Z + gVa_w * w_Z)

  ## ---- Klein (2000): qz.dgges(B, A) => growth eigenvalues mu solve
  ## det(B - mu A) = 0; A's zero rows/columns become mu = Inf (see header) ------
  qzd <- QZ::qz.dgges(B, A)
  mu  <- complex(real = qzd$ALPHAR, imaginary = qzd$ALPHAI) / qzd$BETA
  sel <- is.finite(Mod(mu)) & (Mod(mu) < 1 + 1e-6)
  n_s <- sum(sel)
  if (n_s != npred)
    stop(sprintf(paste0("hank_reiter_statespace(): Blanchard-Kahn failure -- ",
                        "%d stable generalized eigenvalues, need %d ",
                        "(predetermined states)."), n_s, npred))
  o <- QZ::qz.dtgsen(qzd$S, qzd$T, qzd$Q, qzd$Z, sel)
  Z11 <- o$Z[1:npred, 1:npred]
  Z21 <- o$Z[npred + 1:n, 1:npred]
  S11 <- o$S[1:npred, 1:npred]
  T11 <- o$T[1:npred, 1:npred]
  Qj   <- Z21 %*% solve(Z11)                       # dVa_t   = Qj s_t
  Tmat <- Z11 %*% solve(T11, S11) %*% solve(Z11)   # s_{t+1} = Tmat s_t
  Rvec <- c(rep(0, n), 0, 1)

  ## ---- observation rows (functions of s_t = (dD_t, dK_t, dz_t)) ----------------
  selK <- c(rep(0, n), 1, 0); selz <- c(rep(0, n), 0, 1)
  selD <- rbind(diag(n), 0, 0)
  dr_row <- r_K * selK + r_Z * selz
  dw_row <- w_K * selK + w_Z * selz
  da_map <- ga_V %*% (Qj %*% Tmat) + outer(ga_r, dr_row) + outer(ga_w, dw_row)
  dc_map <- gc_V %*% (Qj %*% Tmat) + outer(gc_r, dr_row) + outer(gc_w, dw_row)
  ZA <- as.numeric(va %*% t(selD)) + as.numeric(D_ss %*% da_map)
  ZC <- as.numeric(vc(c_ss) %*% t(selD)) + as.numeric(D_ss %*% dc_map)
  ## output: Y_t = Z_t K_t^alpha L^(1-alpha), K_t = capital in use at t
  ZY <- alpha * Z_ss * K_ss^(alpha - 1) * L^(1 - alpha) * selK +
        Z_ss * K_ss^alpha * L^(1 - alpha) * selz
  Zmat <- rbind(A = ZA, C = ZC, K = selK, z = selz, r = dr_row, w = dw_row,
                Y = ZY)
  state_names <- c(paste0("dD", seq_len(n)), "dK", "dz")

  ## ---- optionally drop the redundant distribution coordinate -------------------
  ## sum(dD) = 0 is invariant; embed/project exactly on that subspace.
  if (drop_dist_coord) {
    keep <- c(seq_len(n - 1L), n + 1L, n + 2L)     # drop dD_n
    E <- matrix(0, npred, npred - 1L)              # embed: dD_n = -sum(others)
    E[keep, ] <- diag(npred - 1L)
    E[n, 1:(n - 1L)] <- -1
    Tmat <- (Tmat %*% E)[keep, , drop = FALSE]
    Rvec <- Rvec[keep]
    Zmat <- Zmat %*% E
    Qj   <- Qj %*% E
    state_names <- state_names[keep]
  }
  sr <- max(Mod(eigen(Tmat, only.values = TRUE)$values))
  if (drop_dist_coord && sr >= 1 - 1e-10)
    warning(sprintf("hank_reiter_statespace(): spectral radius %.8f >= 1 after
  dropping the distribution coordinate -- state space is not stationary.", sr))

  structure(
    list(T_mat = Tmat, R_mat = matrix(Rvec, ncol = 1L), Z_mat = Zmat,
         state_names = state_names, obs_names = rownames(Zmat),
         spectral_radius = sr, n_stable = n_s,
         n_infinite = sum(abs(qzd$BETA) < 1e-12),
         gen_eig_mod = Mod(mu), Qj = Qj,
         rho_z = rho_z, sigma_z = sigma_z,
         drop_dist_coord = drop_dist_coord,
         n_e = n_e, n_a = n_a, ks = ks),
    class = "hank_reiter_ss")
}


#' Impulse responses of a Reiter state space
#'
#' @param rss A \code{\link{hank_reiter_statespace}} object.
#' @param T_h Horizon.
#' @param shock Innovation size (log-TFP units) hitting \code{dz} at t = 1.
#' @return A \code{T_h x n_obs} matrix of observable paths (columns named as
#'   \code{rss$obs_names}).
#' @export
hank_reiter_irf <- function(rss, T_h = 100L, shock = 0.01) {
  s <- as.numeric(rss$R_mat) * shock
  out <- matrix(0, T_h, nrow(rss$Z_mat),
                dimnames = list(NULL, rss$obs_names))
  for (t in seq_len(T_h)) {
    out[t, ] <- as.numeric(rss$Z_mat %*% s)
    s <- as.numeric(rss$T_mat %*% s)
  }
  out
}


#' Kalman log-likelihood of aggregate data on a Reiter state space
#'
#' Runs dynhr's Kalman filter core on the \code{\link{hank_reiter_statespace}}
#' state space with the exact stationary (Lyapunov) initial covariance —
#' the KF-filterable finite-state HANK payoff.
#'
#' @param Y \code{T x n_obs} matrix of DEMEANED observations (deviations from
#'   steady state), columns in the order of \code{observables}.
#' @param rss A \code{hank_reiter_ss} (build with \code{drop_dist_coord = TRUE};
#'   the Lyapunov initialization requires a strictly stable \code{T}).
#' @param observables Which observation rows to use (subset of
#'   \code{rss$obs_names}).
#' @param me_var Measurement-error variance added to each observable (needed
#'   when \code{length(observables) > 1}: one structural shock, so more than
#'   one observable is stochastically singular without it).
#'
#' @return Scalar Gaussian log-likelihood.
#' @export
hank_reiter_kalman_loglik <- function(Y, rss, observables = "A",
                                      me_var = 0) {
  if (!inherits(rss, "hank_reiter_ss"))
    stop("hank_reiter_kalman_loglik(): `rss` must be a hank_reiter_ss object.")
  if (!all(observables %in% rss$obs_names))
    stop("hank_reiter_kalman_loglik(): unknown observable(s): ",
         paste(setdiff(observables, rss$obs_names), collapse = ", "))
  if (rss$spectral_radius >= 1 - 1e-10)
    stop("hank_reiter_kalman_loglik(): T is not strictly stable; build the ",
         "state space with drop_dist_coord = TRUE.")
  Y <- as.matrix(Y)
  if (ncol(Y) != length(observables))
    stop("hank_reiter_kalman_loglik(): ncol(Y) must match observables.")

  TT <- rss$T_mat
  RR <- rss$R_mat
  ZZ <- rss$Z_mat[observables, , drop = FALSE]
  ns <- nrow(TT)

  ## stationary covariance by doubling: P = T P T' + sigma_z^2 R R'
  P  <- rss$sigma_z^2 * (RR %*% t(RR))
  Ak <- TT
  for (i in 1:60) {
    P  <- P + Ak %*% P %*% t(Ak)
    Ak <- Ak %*% Ak
    if (max(abs(Ak)) < 1e-14) break
  }

  out <- .kf_univariate_dispatch(
    Y_minus_d = t(Y),
    ZZ = ZZ, TT = TT, RR = RR,
    DD = matrix(0, length(observables), 1L),
    Sigma_e = matrix(rss$sigma_z^2, 1L, 1L),
    s0 = rep(0, ns),
    P_state = P,
    P_inf_state = NULL,
    me_variance = me_var)
  out$loglik
}


#' PSKF (skewed-shock) log-likelihood on a Reiter state space
#'
#' Runs the pruned-skewed (CSN) Kalman filter on the
#' \code{\link{hank_reiter_statespace}} state space, letting the TFP
#' innovation be SKEW-normal with shape \code{alpha_z} instead of Gaussian:
#' \deqn{s_{t+1} = T s_t + R \epsilon_{t+1}, \qquad
#'       \epsilon \sim \mathrm{CSN}(\cdot;\ \sigma_z^2,\ \alpha_z),}
#' mean-corrected so \eqn{E[\epsilon] = 0} (the anchored steady state is
#' preserved). At \code{alpha_z = 0} this is EXACTLY the Gaussian filter --
#' \code{\link{hank_reiter_kalman_loglik}} -- up to the shared stationary
#' (Lyapunov) initialization, which both use.
#'
#' The scalar-shock lift is a special case of the package's DSGE CSN
#' construction (one skewness dimension, seed covariance \eqn{\sigma_z^2},
#' \eqn{\Gamma_e = \alpha_z / \sigma_z}); because the single shock loads
#' fully into the state, the pseudoinverse \eqn{\Gamma_\eta} lift and the
#' \eqn{\Delta_\eta = I} Schur fallback are both EXACT here (zero residual
#' shock variance off the state range).
#'
#' @param Y \code{T x n_obs} matrix of DEMEANED observations (deviations from
#'   steady state), columns in the order of \code{observables}.
#' @param rss A \code{hank_reiter_ss} (build with \code{drop_dist_coord =
#'   TRUE}; the Lyapunov initialization requires a strictly stable \code{T}).
#' @param alpha_z Skew-normal shape of the TFP innovation (0 = Gaussian;
#'   positive = right-skewed).
#' @param observables Which observation rows to use (subset of
#'   \code{rss$obs_names}).
#' @param me_var Measurement-error variance added to each observable (needed
#'   when \code{length(observables) > 1}, as in the Gaussian filter).
#' @param cut_tol,max_q Skewness-dimension pruning controls of the underlying
#'   PSKF recursion (see the pinned \code{max_q = 5} Miwa-exact cap;
#'   memory: pskf-multishock-pruning-bias). With ONE shock the growth is one
#'   dimension per period, so the defaults are ample.
#'
#' @return Scalar log-likelihood.
#' @seealso \code{\link{hank_reiter_kalman_loglik}}
#' @export
hank_reiter_pskf_loglik <- function(Y, rss, alpha_z = 0, observables = "A",
                                    me_var = 0, cut_tol = 0.01, max_q = 5L) {
  if (!inherits(rss, "hank_reiter_ss"))
    stop("hank_reiter_pskf_loglik(): `rss` must be a hank_reiter_ss object.")
  if (!all(observables %in% rss$obs_names))
    stop("hank_reiter_pskf_loglik(): unknown observable(s): ",
         paste(setdiff(observables, rss$obs_names), collapse = ", "))
  if (rss$spectral_radius >= 1 - 1e-10)
    stop("hank_reiter_pskf_loglik(): T is not strictly stable; build the ",
         "state space with drop_dist_coord = TRUE.")
  if (!is.numeric(alpha_z) || length(alpha_z) != 1L || !is.finite(alpha_z))
    stop("hank_reiter_pskf_loglik(): `alpha_z` must be a finite scalar.")
  Y <- as.matrix(Y)
  if (ncol(Y) != length(observables))
    stop("hank_reiter_pskf_loglik(): ncol(Y) must match observables.")

  ZZ <- rss$Z_mat[observables, , drop = FALSE]
  lift <- .csn_state_noise_lift(
    RR = rss$R_mat,
    DD = matrix(0, length(observables), 1L),
    Sigma_e = matrix(rss$sigma_z^2, 1L, 1L),
    alpha = alpha_z,
    me_variance = me_var)

  .pskf_filter(t(Y), TT = rss$T_mat, ZZ = ZZ,
               mu_eta = lift$mu_eta, Sigma_eta = lift$Sigma_eta,
               Gamma_eta = lift$Gamma_eta, nu_eta = lift$nu_eta,
               Delta_eta = lift$Delta_eta, mu_eps = lift$mu_eps,
               Sigma_eps = lift$Sigma_eps,
               cut_tol = cut_tol, max_q = max_q)
}


#' Build the \code{.kf_loglik_adjoint} state-space list from a Reiter state space
#'
#' Internal mapping shared by \code{\link{hank_reiter_kalman_grad}}: identical
#' to the assembly \code{\link{hank_reiter_kalman_loglik}} feeds to
#' \code{.kf_univariate_dispatch} (same \code{ZZ} row subset, same
#' \code{Sigma_e = sigma_z^2}, \code{DD = 0}, \code{d = 0}).
#' @noRd
.hank_reiter_ss_to_kf <- function(rss, observables, me_var) {
  ZZ <- rss$Z_mat[observables, , drop = FALSE]
  list(TT = rss$T_mat, RR = rss$R_mat, ZZ = ZZ,
       DD = matrix(0, length(observables), 1L),
       d  = rep(0, length(observables)),
       Sigma_e = matrix(rss$sigma_z^2, 1L, 1L))
}


#' Adjoint-KF x central-FD hybrid gradient of the Reiter HANK loglik (P3 #3)
#'
#' Gradient of \code{\link{hank_reiter_kalman_loglik}} with respect to a named
#' vector of model parameters \code{theta}, computed as the PRAGMATIC HYBRID
#' validated by the NZSIM E-wave: an EXACT adjoint (reverse-mode) Kalman-filter
#' gradient with respect to the state-space matrices (\code{T, R, Z, Sigma_e}),
#' contracted against CENTRAL-FD derivatives of those matrices through the
#' Klein solve / linearization for the structural parameters in
#' \code{theta}. The full analytic parameter -> Klein-solve chain is
#' deliberately GATED on Tier-18 A2 (adjoint-of-QZ); until that lands, the
#' FD-of-solve leg is the honest cost: for \code{rho_z}/\code{sigma_z} it is a
#' state-space REBUILD only (\code{hank_reiter_statespace(lin, rho_z, sigma_z)}
#' on an already-linearized economy, ~ms), while for a household parameter
#' like \code{beta} it is a coarse GE re-solve + relinearize (the frozen-lx
#' recipe of \code{\link{hank_ks_coarse_anchored}}, ~10s-100s of ms per
#' evaluation).
#'
#' The adjoint leg (\code{.kf_loglik_adjoint}) requires a NA-free \code{Y}
#' (it has no univariate/missing-data dispatch); use
#' \code{\link{hank_reiter_kalman_loglik}} directly (its
#' \code{.kf_univariate_dispatch} path) when \code{Y} has gaps.
#'
#' @section IMPORTANT -- \code{me_var} convention with \code{length(observables) > 1}:
#' \code{\link{hank_reiter_kalman_loglik}} always routes through
#' \code{.kf_univariate_dispatch}, which (per the documented convention in
#' \code{R/kalman-filter.R}) treats \code{me_var > 0} as TRUE iid diagonal
#' measurement noise (a real noise model). \code{.kf_loglik_adjoint} runs the
#' plain multivariate recursion, which instead adds \code{me_var} to \code{F}
#' ONLY as a positive-definiteness regularizer (never as a noise term in the
#' state/covariance update). The two conventions are IDENTICAL at
#' \code{me_var = 0} (verified here to ~1e-8, both exact evaluations of the
#' same state space) but DIFFER by \code{O(me_var)} otherwise -- confirmed
#' empirically: with 2 observables and \code{me_var = 1e-8} the two logliks
#' differ by ~2 nats, NOT shrinking as \code{me_var -> 0} faster than that
#' \code{O(me_var)} term implies. This function's \code{loglik}/\code{grad}
#' are therefore an EXACT gradient of the adjoint's OWN (regularizer-
#' convention) likelihood, not of \code{hank_reiter_kalman_loglik}'s, whenever
#' \code{me_var > 0} and \code{length(observables) > 1}. Do not mix the two
#' functions' outputs in the same estimation when that condition holds; the
#' single-observable case (\code{length(observables) == 1}, where \code{F} is
#' scalar and stochastic singularity does not arise) is unaffected.
#'
#' @param Y \code{T x n_obs} matrix (or vector when \code{n_obs = 1}) of
#'   DEMEANED observations, columns in the order of \code{observables}. Must
#'   not contain \code{NA}.
#' @param rss_fn \code{function(theta)} returning a \code{hank_reiter_ss}
#'   object (e.g. \code{\link{hank_reiter_statespace}} closed over a fixed
#'   linearization for \code{rho_z}/\code{sigma_z}-only theta, or the
#'   frozen-lx per-draw recipe -- re-solving \code{hank_ks_steady} at
#'   \code{beta_draw * exp(lx_ref)} then \code{hank_reiter_linearize} +
#'   \code{hank_reiter_statespace} -- when \code{theta} includes \code{beta}
#'   or another household parameter).
#' @param theta NAMED numeric vector of parameter values at which to
#'   differentiate.
#' @param observables Which observation rows of the state space to use
#'   (subset of \code{rss_fn(theta)$obs_names}).
#' @param me_var Measurement-error variance added to each observable -- see
#'   the convention warning above when combined with multiple observables.
#' @param fd_step Relative central-FD step: \code{h_j = fd_step *
#'   max(|theta_j|, 1e-2)}. Default \code{1e-5}, i.e. a few times
#'   \code{.Machine$double.eps^(1/3) ~ 6e-6} (the classical central-difference
#'   optimum for smooth functions): the FD leg here differentiates THROUGH the
#'   Klein/QZ solve, which is smooth in \code{theta} away from
#'   Blanchard-Kahn-boundary/eigenvalue-crossing points, so the textbook step
#'   is a safe default; callers doing the beta leg (coarser, noisier GE
#'   re-solve) may want a larger step.
#'
#' @return \code{list(loglik, grad)}: \code{loglik} is the adjoint path's
#'   log-likelihood at \code{theta} (matches
#'   \code{\link{hank_reiter_kalman_loglik}} to ~1e-8 when \code{me_var = 0}
#'   or \code{length(observables) == 1} -- both exact Gaussian-KF evaluations
#'   of the SAME state space; see the \code{me_var} convention section
#'   otherwise); \code{grad} is a numeric vector named as \code{theta}.
#'
#' @section Cost model:
#' \code{2 * length(theta) + 1} evaluations of \code{rss_fn} (central FD per
#' parameter, plus the one at \code{theta} itself) -- this IS the cost model;
#' there is no caching beyond that baseline.
#'
#' @seealso \code{\link{hank_reiter_kalman_loglik}}, \code{\link{hank_ks_coarse_anchored}}
#' @export
hank_reiter_kalman_grad <- function(Y, rss_fn, theta, observables = "A",
                                    me_var = 0, fd_step = NULL) {
  if (is.null(names(theta)) || any(!nzchar(names(theta))))
    stop("hank_reiter_kalman_grad(): `theta` must be a NAMED numeric vector.")
  Y <- as.matrix(Y)
  if (anyNA(Y))
    stop("hank_reiter_kalman_grad(): Y must not contain missing values ",
         "(the adjoint path has no univariate/NA dispatch; use ",
         "hank_reiter_kalman_loglik() directly for data with gaps).")
  if (ncol(Y) != length(observables))
    stop("hank_reiter_kalman_grad(): ncol(Y) must match observables.")
  if (is.null(fd_step)) fd_step <- 1e-5

  par_names <- names(theta)
  n_par <- length(theta)

  rss0 <- rss_fn(theta)
  if (!inherits(rss0, "hank_reiter_ss"))
    stop("hank_reiter_kalman_grad(): rss_fn(theta) must return a ",
         "hank_reiter_ss object.")
  if (!all(observables %in% rss0$obs_names))
    stop("hank_reiter_kalman_grad(): unknown observable(s): ",
         paste(setdiff(observables, rss0$obs_names), collapse = ", "))
  if (rss0$spectral_radius >= 1 - 1e-10)
    stop("hank_reiter_kalman_grad(): T is not strictly stable at theta; ",
         "build the state space with drop_dist_coord = TRUE.")

  ss0 <- .hank_reiter_ss_to_kf(rss0, observables, me_var)

  ## -- central FD through rss_fn for each theta_j: dTT/dRR/dZZ/dSigma_e -----
  ## (dSigma_e = d(sigma_z^2) falls out automatically because rss carries
  ## sigma_z and .hank_reiter_ss_to_kf squares it).
  d_ss_list <- vector("list", n_par)
  for (j in seq_len(n_par)) {
    h <- fd_step * max(abs(theta[[j]]), 1e-2)
    theta_p <- theta; theta_p[[j]] <- theta[[j]] + h
    theta_m <- theta; theta_m[[j]] <- theta[[j]] - h

    rss_p <- rss_fn(theta_p)
    rss_m <- rss_fn(theta_m)
    ss_p <- .hank_reiter_ss_to_kf(rss_p, observables, me_var)
    ss_m <- .hank_reiter_ss_to_kf(rss_m, observables, me_var)

    d_ss_list[[j]] <- list(
      dTT      = (ss_p$TT - ss_m$TT) / (2 * h),
      dRR      = (ss_p$RR - ss_m$RR) / (2 * h),
      dZZ      = (ss_p$ZZ - ss_m$ZZ) / (2 * h),
      dSigma_e = (ss_p$Sigma_e - ss_m$Sigma_e) / (2 * h))
  }

  out <- .kf_loglik_adjoint(t(Y), ss0, d_ss_list, me_variance = me_var)
  names(out$grad) <- par_names
  out
}
