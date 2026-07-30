## R/hank-het-block2d.R
## --------------------------------------------------------------------------
## The DISCRETE-ADJUSTMENT two-asset household as a sequence-space block: the
## P-mixed forward operator, the steady-state block constructor, and the
## nonlinear perfect-foresight transition carrying the (V, Vb, Va) triple.
## Companion to R/hank-egm2d.R (the solver) and R/hank-jacobian2d.R (the
## fake-news Jacobian); design and gate history in briefs/19 sections 9.3-9.8.
##
## The block's transition is a MIXTURE: a household at cell x adjusts with
## probability P(x) (moving to the adjust-branch policy pair) and rolls over
## passively otherwise, so
##     Lambda = diag(P) %*% Lambda_A + diag(1 - P) %*% Lambda_N,
## a row mixture of two joint Young lotteries -- still row-stochastic, so
## hank_stationary_dist() applies unchanged.
##
## AGGREGATE ADJUSTMENT SPENDING: the block's CHI output is
## E[P * (chi_A + F_adj)] -- the convex cost AND the fixed cost are both real
## resource costs, paid only by adjusters. With that definition the household
## budget aggregates to the SAME resource identity as the smooth block,
##     C + CHI = Y + rb * B + ra * A,
## because the KiwiSaver contribution is an internal transfer (income splits
## (1-phi) y to cash and phi y into the locked account; the two budget lines
## sum to the undivided budget on both branches). That identity is this
## block's headline gate, exactly as it was for hank_het2_block.
##
## CLASS: "hank_het2d_block", inheriting NEITHER hank_het_block NOR
## hank_het2_block. The het2 machinery reads fields (b, a, c, Vb, Va) that
## this block does not carry (it has per-branch policies b_A/b_N, ... and the
## triple V/Vb/Va), so silent acceptance anywhere would die on NULLs at best;
## every entry point rejects loudly instead (the .hank_reject_* pattern).
## --------------------------------------------------------------------------


#' P-mixed forward operator for a discrete-adjustment two-asset household
#'
#' The row-stochastic transition over the joint \code{(e, b, a)} cell space
#' when each cell adjusts with probability \code{P} (moving by the
#' adjust-branch joint lottery) and rolls over passively otherwise:
#' \code{diag(P) Lambda_A + diag(1-P) Lambda_N}.
#'
#' @param sol A converged \code{\link{hank_egm2d_solve}} solution (or a
#'   \code{\link{hank_het2d_block}}): fields \code{b_A}, \code{a_A},
#'   \code{b_N}, \code{a_N}, \code{P} are used.
#' @param b_grid,a_grid,Pi As in \code{\link{hank_forward_operator2}}.
#' @return A sparse row-stochastic \code{dgCMatrix} over the joint cell space.
#' @seealso \code{\link{hank_forward_operator2}}, \code{\link{hank_het2d_block}}
#' @export
hank_forward_operator2d <- function(sol, b_grid, a_grid, Pi) {
  need <- c("b_A", "a_A", "b_N", "a_N", "P")
  if (!is.list(sol) || !all(need %in% names(sol)))
    stop("hank_forward_operator2d: 'sol' must carry the per-branch policies ",
         "and P (a hank_egm2d_solve solution or hank_het2d_block).")
  L_A <- hank_forward_operator2(sol$b_A, sol$a_A, b_grid, a_grid, Pi)
  L_N <- hank_forward_operator2(sol$b_N, sol$a_N, b_grid, a_grid, Pi)
  Pv  <- .hank2_arr_to_vec(sol$P)
  Matrix::Diagonal(x = Pv) %*% L_A + Matrix::Diagonal(x = 1 - Pv) %*% L_N
}


#' Construct a discrete-adjustment two-asset household block at steady state
#'
#' Solves the fixed-cost/taste-shock household
#' (\code{\link{hank_egm2d_solve}}) and its stationary joint distribution at
#' fixed aggregate prices, and stores everything the sequence-space Jacobian
#' (\code{\link{hank_het2d_jacobian}}) and nonlinear transition
#' (\code{\link{hank_td2d_nonlinear}}) need.  This is the household for the
#' wealthy hand-to-mouth calibrations the smooth convex cost cannot reach
#' (briefs/19, F21/F22): the fixed cost generates INACTION, and the
#' KiwiSaver-style contribution \code{phi_contrib} keeps the illiquid
#' participation trap closed (F26 -- without an inflow the stationary
#' distribution is initial-condition-dependent; the constructor verifies
#' uniqueness by a two-start check).
#'
#' @param Pi_fn,Pi_inputs Optional transition-probability machinery, exactly
#'   as in \code{\link{hank_het_block}} / \code{\link{hank_het2_block}}:
#'   supplying them makes the named inputs (e.g. the job-finding rate \code{f}
#'   and separation rate \code{s} of \code{\link{hank_employment_income}})
#'   perturbable aggregate inputs of this block, with their own
#'   \code{\link{hank_het2d_jacobian}} columns and
#'   \code{\link{hank_td2d_nonlinear}} paths. The discrete adjust/no-adjust
#'   choice is orthogonal to this: \code{Pi} enters the backward step and the
#'   forward operator, and the taste-shock mixture rides along unchanged.
#' @param b_grid,a_grid,Pi,e,beta,eis,rb,ra,w,chi0,chi1,chi2,n_k,k_max,Tr As in
#'   \code{\link{hank_het2_block}}.
#' @param F_adj,sigma_taste,phi_contrib As in \code{\link{hank_egm2d_solve}}.
#' @param tol,maxit Solver controls.
#' @param dist_tol,dist_maxit Stationary-distribution controls.
#' @param two_start_tol Maximum allowed gap between stationary distributions
#'   computed from two different starting points; a breach means the F26
#'   participation trap is open at this calibration (raise \code{phi_contrib}
#'   or reduce \code{F_adj}).
#'
#' @return An object of class \code{hank_het2d_block}: per-branch policies
#'   (\code{b_A}, \code{a_A}, \code{c_A}, \code{chi_A}, \code{b_N},
#'   \code{a_N}, \code{c_N}), the adjust probability \code{P}, the mixed
#'   envelopes \code{V}, \code{Vb}, \code{Va}, MIXED policies (\code{b},
#'   \code{a}, \code{c}: probability-weighted, for reporting), \code{Lambda},
#'   \code{D}, aggregates \code{B}, \code{A}, \code{C}, \code{CHI} (convex +
#'   fixed cost, adjusters only), \code{adj_freq} (aggregate adjustment
#'   probability), dims, and the calibration.
#' @seealso \code{\link{hank_het2_block}} (smooth),
#'   \code{\link{hank_egm2d_solve}}, \code{\link{hank_het2d_jacobian}}
#' @export
hank_het2d_block <- function(b_grid, a_grid, Pi, e, beta, eis, rb, ra, w,
                             chi0 = 0.25, chi1 = 6.5, chi2 = 2,
                             F_adj, sigma_taste, phi_contrib = 0.06,
                             n_k = 50L, k_max = 1, Tr = 0,
                             Pi_fn = NULL, Pi_inputs = NULL,
                             tol = 1e-9, maxit = 5000L,
                             dist_tol = 1e-13, dist_maxit = 200000L,
                             two_start_tol = 1e-8) {
  if (!is.numeric(e) || !all(is.finite(e)))
    stop("hank_het2d_block: 'e' must be a finite numeric vector.")
  if (!is.numeric(w) || length(w) != 1L || !is.finite(w))
    stop("hank_het2d_block: 'w' must be a finite numeric scalar.")
  if (!is.numeric(Tr) || length(Tr) != 1L || !is.finite(Tr))
    stop("hank_het2d_block: 'Tr' must be a finite numeric scalar.")
  .hank_check_pi_fn(Pi_fn, Pi_inputs, Pi,
                    reserved = c("rb", "ra", "w", "Tr"),
                    caller = "hank_het2d_block")
  hh <- hank_egm2d_solve(b_grid, a_grid, y = w * e + Tr, rb = rb, ra = ra,
                         beta = beta, eis = eis, chi0 = chi0, chi1 = chi1,
                         chi2 = chi2, F_adj = F_adj,
                         sigma_taste = sigma_taste, Pi = Pi, n_k = n_k,
                         k_max = k_max, tol = tol, maxit = maxit,
                         phi_contrib = phi_contrib)
  if (!hh$converged)
    warning("hank_het2d_block: household solve did not converge at steady state")

  Lam <- hank_forward_operator2d(hh, b_grid, a_grid, Pi)
  sd  <- hank_stationary_dist(Lam, tol = dist_tol, maxit = dist_maxit)
  ## F26 guard: uniqueness of the invariant distribution, checked from a
  ## second start (all mass at the richest cell). The participation trap is
  ## exactly an initial-condition dependence here, and it is grid- and
  ## calibration-dependent, so it is CHECKED rather than assumed.
  n_cell <- length(sd$d)
  d2 <- rep(0, n_cell); d2[n_cell] <- 1
  for (i in seq_len(as.integer(dist_maxit))) {
    dn <- as.numeric(Matrix::t(Lam) %*% d2)
    if (max(abs(dn - d2)) < dist_tol) { d2 <- dn; break }
    d2 <- dn
  }
  gap2 <- max(abs(sd$d - d2))
  if (gap2 > two_start_tol)
    warning("hank_het2d_block: stationary distributions from two different ",
            "starts disagree by ", format(gap2), " -- the F26 participation ",
            "trap is open at this calibration (near-multiple invariant ",
            "sets). Raise 'phi_contrib' or reduce 'F_adj'; aggregates below ",
            "are start-dependent.", call. = FALSE)
  D <- sd$d

  Pv <- .hank2_arr_to_vec(hh$P)
  ## mixed (probability-weighted) policies, for reporting and aggregates
  b_mix <- hh$P * hh$b_A + (1 - hh$P) * hh$b_N
  a_mix <- hh$P * hh$a_A + (1 - hh$P) * hh$a_N
  c_mix <- hh$P * hh$c_A + (1 - hh$P) * hh$c_N
  chi_mix <- hh$P * (hh$chi_A + F_adj)          # real resource costs, adjusters
  structure(
    c(hh[c("b_A", "a_A", "c_A", "chi_A", "b_N", "a_N", "c_N", "P",
           "V", "Vb", "Va", "k_grid", "iterations", "converged")],
      list(b = b_mix, a = a_mix, c = c_mix, chi = chi_mix,
           b_grid = b_grid, a_grid = a_grid, Pi = Pi, e = e,
           beta = beta, eis = eis, rb = rb, ra = ra, w = w, Tr = Tr,
           Pi_fn = Pi_fn, Pi_inputs = Pi_inputs,
           chi0 = chi0, chi1 = chi1, chi2 = chi2,
           F_adj = F_adj, sigma_taste = sigma_taste,
           phi_contrib = phi_contrib,
           Lambda = Lam, D = D,
           B = hank_aggregate2(D, b_mix),
           A = hank_aggregate2(D, a_mix),
           C = hank_aggregate2(D, c_mix),
           CHI = hank_aggregate2(D, chi_mix),
           adj_freq = sum(D * Pv),
           n_e = length(e), n_b = length(b_grid), n_a = length(a_grid),
           dist_converged = sd$converged, two_start_gap = gap2)),
    class = "hank_het2d_block")
}


#' One nonlinear discrete-adjustment backward step at given aggregate prices
#'
#' Maps block inputs \code{(rb, ra, w, Tr)} to \code{\link{.hank_egm2d_step}};
#' the primitive the fake-news sweep and nonlinear transition are built from.
#' \code{Psi1_grid} depends on \code{ra} and is rebuilt per call (the F14
#' trap: caching it would silently corrupt every \code{ra} column).
#' @keywords internal
.hank_block_step2d <- function(block, V_p, Vb_p, Va_p, rb, ra, w,
                               Tr = NULL, Pi = NULL) {
  if (is.null(Pi)) Pi <- block$Pi
  if (is.null(Tr)) Tr <- .hank_block_tr(block)
  n_a <- block$n_a
  Psi1 <- .hank_psi(matrix(block$a_grid, n_a, n_a),
                    matrix(block$a_grid, n_a, n_a, byrow = TRUE),
                    ra, block$chi0, block$chi1, block$chi2)$Psi1
  .hank_egm2d_step(V_p, Vb_p, Va_p, block$b_grid, block$a_grid, block$k_grid,
                   y = w * block$e + Tr, rb = rb, ra = ra, beta = block$beta,
                   eis = block$eis, chi0 = block$chi0, chi1 = block$chi1,
                   chi2 = block$chi2, F_adj = block$F_adj,
                   sigma_taste = block$sigma_taste, Pi = Pi,
                   Psi1_grid = Psi1, phi_contrib = block$phi_contrib)
}


#' Nonlinear perfect-foresight transition of a discrete-adjustment block
#'
#' The \code{\link{hank_td2_nonlinear}} analogue for the fixed-cost household:
#' the backward pass carries the \code{(V, Vb, Va)} TRIPLE (the level is what
#' the discrete choice compares), and the forward pass pushes the distribution
#' through the per-period P-mixed operator.
#'
#' @param block A \code{\link{hank_het2d_block}}.
#' @param pi_input_paths Optional named list of length-\code{T_h} LEVEL paths
#'   for (a subset of) the block's transition-probability inputs
#'   (\code{names(block$Pi_inputs)}); missing inputs stay at their
#'   steady-state values. Requires a block built with
#'   \code{Pi_fn}/\code{Pi_inputs}. \code{Pi_t} is the transition between
#'   \code{t} and \code{t+1}, entering the date-\code{t} backward step AND
#'   the date-\code{t} forward push -- the same timing convention as every
#'   other tier. Default \code{NULL} (constant steady-state \code{Pi}).
#' @param rb_path,ra_path,w_path,Tr_path Optional length-\code{T_h} LEVEL
#'   paths; defaults hold each at its steady-state value.
#' @param T_h Integer horizon.
#' @param keep_policies When \code{TRUE}, also return per-period branch
#'   policies, \code{P}, and the sparse operators.
#' @param D0 Optional initial distribution (package two-asset cell order).
#'
#' @return A list with paths \code{B}, \code{A}, \code{C}, \code{CHI},
#'   \code{ADJ} (aggregate adjustment probability), and \code{Dpath}.
#' @seealso \code{\link{hank_het2d_block}}, \code{\link{hank_td2_nonlinear}}
#' @export
hank_td2d_nonlinear <- function(block, rb_path = NULL, ra_path = NULL,
                                w_path = NULL, Tr_path = NULL, T_h = NULL,
                                keep_policies = FALSE, D0 = NULL,
                                pi_input_paths = NULL) {
  if (!inherits(block, "hank_het2d_block"))
    stop("hank_td2d_nonlinear: 'block' must be a hank_het2d_block.")
  if (!is.null(D0)) {
    D0 <- as.numeric(D0)
    if (length(D0) != length(block$D) || any(D0 < -1e-12) ||
        abs(sum(D0) - 1) > 1e-8)
      stop("hank_td2d_nonlinear: D0 must be a nonnegative length-",
           length(block$D), " distribution summing to 1.", call. = FALSE)
  }
  if (is.null(T_h)) {
    lens <- c(length(rb_path), length(ra_path), length(w_path),
              length(Tr_path),
              if (!is.null(pi_input_paths))
                vapply(pi_input_paths, length, integer(1)))
    T_h <- max(lens, if (all(lens == 0L)) 1L else 0L)
  }
  if (is.null(rb_path)) rb_path <- rep(block$rb, T_h)
  if (is.null(ra_path)) ra_path <- rep(block$ra, T_h)
  if (is.null(w_path))  w_path  <- rep(block$w,  T_h)
  if (is.null(Tr_path)) Tr_path <- rep(.hank_block_tr(block), T_h)
  stopifnot(length(rb_path) == T_h, length(ra_path) == T_h,
            length(w_path) == T_h, length(Tr_path) == T_h)

  ## Per-period transition matrices; NULL = steady-state Pi everywhere,
  ## which is the zero-allocation path every pre-existing call takes.
  Pi_path <- .hank_pi_path(block, pi_input_paths, T_h)
  Pi_at   <- function(t) if (is.null(Pi_path)) block$Pi else Pi_path[[t]]

  steps <- vector("list", T_h)
  V <- block$V; Vb <- block$Vb; Va <- block$Va
  for (t in T_h:1L) {
    st <- .hank_block_step2d(block, V, Vb, Va, rb_path[t], ra_path[t],
                             w_path[t], Tr = Tr_path[t], Pi = Pi_at(t))
    steps[[t]] <- st
    V <- st$V; Vb <- st$Vb; Va <- st$Va
  }

  D <- if (is.null(D0)) block$D else D0
  B <- numeric(T_h); A <- numeric(T_h); C <- numeric(T_h)
  CHI <- numeric(T_h); ADJ <- numeric(T_h)
  Dpath <- matrix(0, length(D), T_h)
  Lam_keep <- if (keep_policies) vector("list", T_h) else NULL
  for (t in seq_len(T_h)) {
    st <- steps[[t]]
    Dpath[, t] <- D
    b_mix <- st$P * st$b_A + (1 - st$P) * st$b_N
    a_mix <- st$P * st$a_A + (1 - st$P) * st$a_N
    c_mix <- st$P * st$c_A + (1 - st$P) * st$c_N
    B[t]   <- hank_aggregate2(D, b_mix)
    A[t]   <- hank_aggregate2(D, a_mix)
    C[t]   <- hank_aggregate2(D, c_mix)
    CHI[t] <- hank_aggregate2(D, st$P * (st$chi_A + block$F_adj))
    ADJ[t] <- sum(D * .hank2_arr_to_vec(st$P))
    Lam <- hank_forward_operator2d(st, block$b_grid, block$a_grid, Pi_at(t))
    if (keep_policies) Lam_keep[[t]] <- Lam
    D <- as.numeric(Matrix::t(Lam) %*% D)
  }
  out <- list(B = B, A = A, C = C, CHI = CHI, ADJ = ADJ, Dpath = Dpath)
  if (keep_policies) {
    out$steps <- steps
    out$Lambda <- Lam_keep
  }
  out
}
