## R/hank-egm2.R
## --------------------------------------------------------------------------
## Two-asset (liquid/illiquid) household solver: the SSJ-standard
## Auclert-Bardoczy-Rognlie-Straub (2021) two-asset household, solved by the
## two-stage endogenous-gridpoints method.
##
## Household problem (recursive), with LIQUID b (return rb) and ILLIQUID a
## (return ra, convex adjustment cost Psi):
##
##   V(e, b, a) = max_{b' >= b_grid[1], a'} u(c) + beta E[V(e', b', a') | e]
##   s.t.  c + b' + a' + Psi(a', a) = (1 + rb) b + (1 + ra) a + y(e)
##         u(c) = c^{1-1/eis} / (1 - 1/eis)
##
## The illiquid asset earns more (ra > rb) but moving it costs Psi, so
## households hold liquid buffer stock and a large illiquid position: the
## "wealthy hand-to-mouth" of Kaplan-Moll-Violante.  This is the mechanism a
## one-asset HANK cannot produce (its constrained households are also poor).
##
## PORTED FROM, and validated against, the reference implementation:
## shade-econ/sequence-jacobian, src/sequence_jacobian/hetblocks/hh_twoasset.py
## (pinned 2026-07-15; the transcription, the deviations below, and the
## reference calibration are recorded in briefs/19-twoasset-hank-scope.md
## section 3 -- READ THAT BEFORE CHANGING ANY FORMULA HERE).
##
## THREE DELIBERATE DEVIATIONS from the reference (brief section 3.3):
##   F1 was a reference SUBTLETY, not a deviation: k_grid must be DECREASING
##      (SSJ builds it as agrid(...)[::-1]).  b_endo is decreasing in the
##      multiplier kappa, so a decreasing kappa grid makes b_endo INCREASING
##      along the axis -- which is what the ascending-x interpolation in step 6
##      requires.  hank_egm2_solve() enforces and documents this.
##   F3 Wa/Wb/c are floored at `tiny` (SSJ does not floor).  This is dynhr's
##      one-asset convention (.hank_egm_step) -- a transient-NaN guard that is
##      an exact no-op at any converged solution; test-hank-egm2.R asserts the
##      floors are inactive at convergence.
##   Expectations: SSJ takes the Pi expectation OUTSIDE the backward function
##      (its @het(exogenous='Pi') decorator does it), so its Va_p is already
##      Pi %*% Va_next.  dynhr's convention takes it INSIDE the step, matching
##      .hank_egm_step.  Mathematically identical; a real contract difference.
##
## GRID/ARRAY ORDER: everything is (e, b, a) -- income slowest, ILLIQUID a
## FASTEST -- matching SSJ's numpy (z, b, a) row-major layout, so reference
## dumps transpose trivially.  Policies/marginals are n_e x n_b x n_a arrays.
## Because e and b are the leading axes, the (e,b')-by-a' unfolding used by the
## crossing search is a free reshape: matrix(X, n_e*n_b, n_a).
## --------------------------------------------------------------------------


## The transient floor shared by .hank_egm2_step (Wa/Wb/c) and the compiled
## kernel. Defined once so the step, the kernel and the convergence-time
## feasibility check below cannot drift apart.
tiny_floor <- function() 1e-12


#' Adjustment cost of moving the illiquid asset, and its derivatives
#'
#' The SSJ/ABRS convex adjustment cost and its partials, evaluated
#' element-wise on conformable arrays (or an \eqn{(a', a)} matrix pair):
#' \deqn{\Psi(a', a) = \frac{\chi_1}{\chi_2}\,|a' - (1+r_a)a|
#'       \left(\frac{|a' - (1+r_a)a|}{(1+r_a)a + \chi_0}\right)^{\chi_2 - 1}.}
#' The \eqn{\chi_0} shift in the denominator (rather than a \code{max(a, .)}
#' guard) keeps it finite and positive at \eqn{a = 0}; the cost is zero at the
#' no-adjustment point \eqn{a' = (1+r_a)a}.
#'
#' @param ap Numeric: next-period illiquid holdings \eqn{a'}.
#' @param a Numeric: current illiquid holdings \eqn{a} (conformable with
#'   \code{ap}).
#' @param ra Numeric: illiquid return.
#' @param chi0,chi1,chi2 Adjustment-cost parameters: \code{chi0 > 0} (the
#'   denominator shift), \code{chi1 >= 0} (the scale), \code{chi2 > 1} (the
#'   curvature; \code{chi2 = 2} is the reference calibration and makes
#'   \eqn{\Psi} smooth in \eqn{a'}).
#'
#' @return List with \code{Psi}, \code{Psi1} (\eqn{\partial\Psi/\partial a'})
#'   and \code{Psi2} (\eqn{\partial\Psi/\partial a}), each shaped like the
#'   recycled \code{ap}/\code{a}.
#' @keywords internal
.hank_psi <- function(ap, a, ra, chi0, chi1, chi2) {
  a_with_return <- (1 + ra) * a
  a_change      <- ap - a_with_return
  abs_change    <- abs(a_change)
  adj_denom     <- a_with_return + chi0
  core          <- (abs_change / adj_denom)^(chi2 - 1)
  Psi  <- (chi1 / chi2) * abs_change * core
  Psi1 <- chi1 * sign(a_change) * core
  Psi2 <- -(1 + ra) * (Psi1 + (chi2 - 1) * Psi / adj_denom)
  list(Psi = Psi, Psi1 = Psi1, Psi2 = Psi2)
}


#' Interpolation coordinates against increasing knots (SSJ interpolate_coord)
#'
#' Row-wise: for each row \code{l} of \code{X} (strictly increasing knots),
#' the bracketing index and lower weight representing each query in \code{xq},
#' so that \code{xq = p * X[l, i] + (1-p) * X[l, i+1]}.  Out-of-domain queries
#' EXTRAPOLATE linearly (\code{p} escapes \code{[0,1]}), matching SSJ and
#' \code{\link{.hank_interp1}} -- unlike \code{\link{.hank_lottery}}, which
#' clamps (a lottery must not create negative mass).
#'
#' This is SSJ's FORWARD SWEEP, not a binary search, and the difference is
#' load-bearing.  The sweep carries the bracket index across the (increasing)
#' query points, so it is O(n_x + n_q); on strictly increasing knots it returns
#' exactly what \code{findInterval} would (they disagree on which side of an
#' exact tie to bracket, but both then interpolate to the same value).  Where
#' they differ is on NON-monotone knots: \code{findInterval} errors and a
#' binary search silently returns a wrong bracket, whereas the sweep marches
#' forward and degrades gracefully.
#'
#' That tolerance is required, not incidental.  \code{b_endo} is monotone in
#' \eqn{b'} only once the marginal values are well-behaved, and the FIRST
#' iterate from SSJ's own \code{hh_init} guess is not: at plausible
#' calibrations (e.g. \code{chi2 = 3}) it is non-monotone, and the transient
#' washes out within a few iterations.  A monotonicity requirement here kills
#' those solves outright even though the fixed point is perfectly well-behaved
#' -- the same class of transient the one-asset \code{tiny} floors exist to
#' absorb (see \code{\link{.hank_egm_step}}).  Genuine breakdowns are caught
#' instead by \code{hank_egm2_solve}'s non-finite post-check, which fires for
#' both backends identically.
#'
#' \code{src/hank_egm2.cpp} carries the identical sweep, so the two backends
#' agree on non-monotone input as well as monotone.
#'
#' @param X Numeric \code{n_lead x n_x} matrix: row \code{l} is one knot vector.
#' @param xq Numeric length-\code{n_q}: INCREASING query points (shared across
#'   rows).
#' @return List with \code{i} and \code{p}, each \code{n_lead x n_q}.
#' @keywords internal
.hank2_interp_coord_rows <- function(X, xq) {
  n_lead <- nrow(X); n_x <- ncol(X); n_q <- length(xq)
  iout <- matrix(0L, n_lead, n_q)
  pout <- matrix(0,  n_lead, n_q)
  for (l in seq_len(n_lead)) {
    x <- X[l, ]
    xi <- 1L; x_low <- x[1L]; x_high <- x[2L]
    for (q in seq_len(n_q)) {
      xqc <- xq[q]
      while (xi < n_x - 1L) {
        if (x_high >= xqc) break
        xi <- xi + 1L; x_low <- x_high; x_high <- x[xi + 1L]
      }
      iout[l, q] <- xi
      pout[l, q] <- (x_high - xqc) / (x_high - x_low)
    }
  }
  list(i = iout, p = pout)
}


#' Apply interpolation coordinates to a SHARED knot-value vector
#'
#' \code{yq = p * y[i] + (1-p) * y[i+1]}, preserving the shape of \code{i}.
#' @keywords internal
.hank2_apply_coord_vec <- function(i, p, y) {
  array(p * y[i] + (1 - p) * y[i + 1L], dim = dim(i))
}


#' Apply interpolation coordinates to PER-ROW knot values
#'
#' Row-wise \code{yq[l, ] = p[l, ] * Y[l, i[l, ]] + (1-p[l, ]) * Y[l, i[l, ]+1]}.
#' @keywords internal
.hank2_apply_coord_rows <- function(i, p, Y) {
  out <- matrix(0, nrow(i), ncol(i))
  for (l in seq_len(nrow(i))) {
    ii <- i[l, ]
    out[l, ] <- p[l, ] * Y[l, ii] + (1 - p[l, ]) * Y[l, ii + 1L]
  }
  out
}


#' Locate where a decreasing lhs crosses an increasing rhs (SSJ
#' lhs_equals_rhs_interpolate)
#'
#' For each row \code{l} of \code{lhs} and each column \code{j} of \code{rhs},
#' find \code{i} with \code{lhs[l, i] > rhs[i, j]} and
#' \code{lhs[l, i+1] < rhs[i+1, j]}, and the weight \code{p} with
#' \code{p*(lhs[l,i] - rhs[i,j]) + (1-p)*(lhs[l,i+1] - rhs[i+1,j]) == 0}, i.e.
#' the linearly-interpolated crossing point.
#'
#' Assumes (as SSJ documents) that \code{lhs - rhs} is DECREASING in \code{i}
#' and the solution is increasing in \code{j}.  Here that holds because
#' \code{lhs} is \eqn{W_a/W_b} (decreasing in \eqn{a'} by concavity) and
#' \code{rhs} is \eqn{1 + \Psi_1(a', a)} (increasing in \eqn{a'} by convexity
#' of \eqn{\Psi}).
#'
#' Under that monotonicity, "first \code{i} with \code{lhs - rhs < 0}" equals
#' \code{colSums(D >= 0) + 1}, which is how this vectorizes.  SSJ instead
#' carries \code{i} ACROSS \code{j} as a monotone sweep; the two agree whenever
#' the documented assumption holds, and the per-column form used here is
#' strictly more robust when it does not (brief 19, finding F5).  SSJ also
#' bounds its scan with \code{nj} where \code{ni} is meant (F4) -- invisible
#' there because \eqn{\Psi_1} is square; the correct bound is used here.
#'
#' Corners, matching SSJ exactly:
#' \itemize{
#'   \item no crossing at the bottom (\code{lhs[l,1] < rhs[1,j]}): return
#'     \code{i = 1, p = 1}, i.e. the illiquid floor \eqn{a' = }\code{a_grid[1]};
#'   \item no crossing by the top: return \code{i = n_i - 1} with \code{p}
#'     from the last two points, i.e. linear extrapolation above the grid.
#' }
#'
#' @param lhs Numeric \code{n_lead x n_i} matrix (rows = broadcast leading
#'   states, columns = the search axis).
#' @param rhs Numeric \code{n_i x n_j} matrix.
#' @return List with \code{i} (integer) and \code{p}, each \code{n_lead x n_j}.
#' @keywords internal
## TWO IMPLEMENTATIONS, DISPATCHED ON A MEASURED SIZE THRESHOLD (0.9.0.0028).
##
## The row loop below (`.hank_lhs_eq_rhs_loop`) does eight R-level operations
## on an n_i x n_j matrix per leading row, and n_i = n_j = n_a, so on the grids
## the fake-news sweep actually uses those matrices are TINY and the loop is R
## call overhead rather than arithmetic. Profiled on the installed -O2 build it
## was ~40% of a whole two-asset het Jacobian at n_a = 16 -- the largest single
## item left in that route.
##
## `.hank_lhs_eq_rhs_vec` issues the same comparisons and subtractions as a
## handful of array operations instead of 8*n_lead small ones. Same FLOPs, same
## elementwise order, so it is BIT-IDENTICAL to the loop rather than merely
## close -- asserted with expect_identical in test-hank-egm2-lhs-eq-rhs.R.
##
## BUT IT IS NOT UNIFORMLY FASTER, and shipping it unconditionally would have
## been a large regression on big grids. Measured speedup vs the loop, by the
## per-row matrix size n_i*n_j (briefs/21 section 14):
##
##     n_i*n_j     144    256     625    1024     1600    2500    4096
##     speedup    8.0x   2-7x   2.5-3x   1.5x   1.0-1.2x  0.8x   0.67x
##
## Once the per-row matrix is big enough to amortize the loop's own overhead,
## materializing the n_i x n_j x n_lead intermediate costs more than it saves.
## Hence the threshold: vectorize at or below 1024 (n_a <= 32), loop above it.
## 1024 rather than the ~1600 break-even, deliberately: at 1600 the end-to-end
## Jacobian measured 0.99x, i.e. the isolated win had already been eaten, so
## the threshold is set where the gain is still unambiguous rather than where
## the curves cross. Both paths are kept and both are gated against each other.
.hank_lhs_eq_rhs <- function(lhs, rhs) {
  if (ncol(lhs) * ncol(rhs) <= 1024L) .hank_lhs_eq_rhs_vec(lhs, rhs)
  else .hank_lhs_eq_rhs_loop(lhs, rhs)
}

## The original row loop; the large-n_i*n_j path.
.hank_lhs_eq_rhs_loop <- function(lhs, rhs) {
  n_lead <- nrow(lhs); n_i <- ncol(lhs); n_j <- ncol(rhs)
  iout <- matrix(0L, n_lead, n_j)
  pout <- matrix(0,  n_lead, n_j)
  jj   <- seq_len(n_j)
  for (l in seq_len(n_lead)) {
    D  <- lhs[l, ] - rhs                  # (n_i x n_j), decreasing down rows
    i1 <- pmin(colSums(D >= 0) + 1L, n_i) # first i with D < 0, capped at the top
    at_floor <- i1 == 1L
    lo <- pmax(i1 - 1L, 1L)
    Dl <- D[cbind(lo, jj)]
    Du <- D[cbind(i1, jj)]
    p  <- -Du / (Dl - Du)
    ## Bottom corner: the crossing is below the grid -> sit exactly on a_grid[1].
    p[at_floor]  <- 1
    lo[at_floor] <- 1L
    iout[l, ] <- lo
    pout[l, ] <- p
  }
  list(i = iout, p = pout)
}

## Chunk-vectorized over the leading axis; the small-n_i*n_j path.
.hank_lhs_eq_rhs_vec <- function(lhs, rhs) {
  n_lead <- nrow(lhs); n_i <- ncol(lhs); n_j <- ncol(rhs)
  ## Empty input: match the loop's types exactly (it never assigns, so its
  ## integer `iout` is never promoted to double the way a populated one is).
  if (n_lead == 0L)
    return(list(i = matrix(0L, 0L, n_j), p = matrix(0, 0L, n_j)))
  iout <- matrix(0, n_lead, n_j)
  pout <- matrix(0, n_lead, n_j)
  ## Chunked so the intermediate stays bounded (~1e6 doubles) even though the
  ## threshold above already keeps n_i*n_j small: n_lead is unbounded.
  chunk <- max(1L, min(n_lead, as.integer(1e6 %/% max(n_i * n_j, 1L))))
  jj0 <- (seq_len(n_j) - 1L) * n_i
  for (start in seq.int(1L, n_lead, by = chunk)) {
    idx <- start:min(start + chunk - 1L, n_lead)
    m <- length(idx)
    ## D[i, j, l] = lhs[idx[l], i] - rhs[i, j], flattened to n_i x (n_j*m) with
    ## j fastest -- the layout colSums() and the linear indexing below read.
    D <- aperm(array(lhs[idx, , drop = FALSE], c(m, n_i, n_j)),
               c(2L, 3L, 1L)) - array(rhs, c(n_i, n_j, m))
    dim(D) <- c(n_i, n_j * m)
    i1 <- pmin(colSums(D >= 0) + 1L, n_i)  # first i with D < 0, capped at top
    at_floor <- i1 == 1L
    lo <- pmax(i1 - 1L, 1L)
    off <- rep(jj0, times = m) +
      rep((seq_len(m) - 1L) * (n_i * n_j), each = n_j)
    Du <- D[i1 + off]
    p  <- -Du / (D[lo + off] - Du)
    p[at_floor]  <- 1
    lo[at_floor] <- 1
    ## `lo`/`p` run j-fastest, l-slower; the outputs are n_lead x n_j.
    iout[idx, ] <- matrix(lo, m, n_j, byrow = TRUE)
    pout[idx, ] <- matrix(p,  m, n_j, byrow = TRUE)
  }
  list(i = iout, p = pout)
}


## Broadcast an (e, b, a)-shaped array whose value depends on one axis only.
## Cheap enough to rebuild per step; kept explicit for readability.
.hank2_bcast_e <- function(v, n_e, n_mid, n_a)
  array(rep(v, times = n_mid * n_a), c(n_e, n_mid, n_a))
.hank2_bcast_mid <- function(v, n_e, n_mid, n_a)
  array(rep(rep(v, each = n_e), times = n_a), c(n_e, n_mid, n_a))
.hank2_bcast_a <- function(v, n_e, n_mid, n_a)
  array(rep(v, each = n_e * n_mid), c(n_e, n_mid, n_a))


#' One two-asset EGM backward step
#'
#' The two-stage EGM of SSJ's \code{hh_twoasset.hh} (steps 2-7; see the file
#' header and brief 19 section 3.2).  Given next-period marginal values, returns
#' updated marginal values and this period's liquid/illiquid/consumption
#' policies.
#'
#' @param Vb_p,Va_p Numeric \code{n_e x n_b x n_a} arrays: next-period marginal
#'   value of LIQUID and ILLIQUID assets.  The \code{Pi} expectation is taken
#'   HERE (dynhr convention; SSJ takes it outside -- see the file header).
#' @param b_grid,a_grid Numeric: increasing liquid / illiquid grids.
#' @param k_grid Numeric: the DECREASING multiplier grid for the
#'   liquid-constrained branch (see the file header, F1).
#' @param y Numeric length-\code{n_e}: labour income by income state.
#' @param rb,ra Numeric: liquid and illiquid returns.
#' @param beta,eis Discount factor and elasticity of intertemporal substitution.
#' @param chi0,chi1,chi2 Adjustment-cost parameters (see \code{\link{.hank_psi}}).
#' @param Pi Numeric \code{n_e x n_e} row-stochastic income transition matrix.
#' @param Psi1_grid Optional precomputed \code{n_a x n_a} matrix
#'   \eqn{\Psi_1[i, j] = \Psi_1(a' = }\code{a_grid[i]}\eqn{, a = }\code{a_grid[j]}\eqn{)} (SSJ's
#'   \code{marginal_cost_grid}); rebuilt when \code{NULL}.  It depends only on
#'   \code{(a_grid, ra, chi*)}, so callers iterating the step at fixed prices
#'   pass it once.
#'
#' @param theta_coll Scalar LTV in \eqn{[0, 1)} on end-of-period illiquid
#'   collateral (D1, brief 19 section 9.12).  When positive, the liquid AXIS of
#'   every array is reinterpreted as the GAP COORDINATE \eqn{x = b + \theta a}
#'   (liquid wealth plus collateral capacity), in which the collateral
#'   constraint \eqn{b' \ge b_{grid}[1] - \theta a'} is the CONSTANT floor
#'   \eqn{x' \ge b_{grid}[1]} -- a grid line, so the validated constant-floor
#'   machinery applies verbatim.  The transform only re-weights the budget and
#'   FOC: illiquid outlay costs \eqn{(1-\theta)} liquid units, illiquid arrival
#'   resources carry \eqn{\rho_a = (1+r_a) - \theta(1+r_b)}, and the illiquid
#'   FOC compares \eqn{W_a/W_x} to \eqn{(1-\theta) + \Psi_1}.  \eqn{\Psi} stays
#'   a function of \eqn{(a', a)} unchanged.  At \code{theta_coll = 0} every
#'   expression reduces bit-for-bit to the pre-collateral form
#'   (\eqn{x \equiv b}).  True liquid is \eqn{b = x - \theta a} (recovered by
#'   the callers; this step works in \eqn{x} throughout).
#' @param dtheta_next Scalar \eqn{\theta_{t+1} - \theta_t} for TIME-VARYING
#'   collateral paths.  \eqn{x} is defined with the CONTEMPORANEOUS
#'   \eqn{\theta}, so a household choosing \eqn{x'} today arrives at
#'   \eqn{x_{t+1} = x' + (\theta_{t+1}-\theta_t) a'} tomorrow; the next-period
#'   marginals are pre-shifted along the liquid axis by that amount (one-sided
#'   linear extrapolation at the edges -- tightening pushes legacy borrowers
#'   below tomorrow's floor, the forced-deleveraging case).  Exact no-op at 0.
#'
#' @return List with \code{Vb}, \code{Va}, \code{b}, \code{a}, \code{c} and
#'   \code{chi} (the per-cell adjustment cost), each \code{n_e x n_b x n_a}.
#'   Under collateral, \code{Vb}/\code{b} are the marginal value of and policy
#'   for the gap coordinate \eqn{x}.
#' @keywords internal
.hank_egm2_step <- function(Vb_p, Va_p, b_grid, a_grid, k_grid, y, rb, ra,
                            beta, eis, chi0, chi1, chi2, Pi,
                            Psi1_grid = NULL, theta_coll = 0,
                            dtheta_next = 0) {
  n_e <- dim(Va_p)[1L]; n_b <- dim(Va_p)[2L]; n_a <- dim(Va_p)[3L]
  n_k <- length(k_grid)
  tiny <- 1e-12

  if (is.null(Psi1_grid))
    Psi1_grid <- .hank_psi(matrix(a_grid, n_a, n_a),                # a' down rows
                           matrix(a_grid, n_a, n_a, byrow = TRUE),  # a across cols
                           ra, chi0, chi1, chi2)$Psi1
  rhs <- (1 - theta_coll) + Psi1_grid
  rho_a <- (1 + ra) - theta_coll * (1 + rb)   # illiquid arrival, x-coordinates

  ## --- Step 1b (theta PATH only): arrival-coordinate shift ------------------
  ## x is defined with the contemporaneous theta, so tomorrow's state is
  ## x_next = x' + dtheta_next * a'.  Shift the next-period marginals along
  ## the liquid axis per a-slice, then the rest of the step is the standard
  ## constant-floor algorithm.  Indices clamp to [1, n_b-1] and the weight is
  ## left unclamped: one-sided LINEAR EXTRAPOLATION at both edges.
  if (dtheta_next != 0) {
    for (i in seq_len(n_a)) {
      xq <- b_grid + dtheta_next * a_grid[i]
      j  <- pmin(pmax(findInterval(xq, b_grid), 1L), n_b - 1L)
      pw <- matrix((b_grid[j + 1L] - xq) / (b_grid[j + 1L] - b_grid[j]),
                   n_e, n_b, byrow = TRUE)
      Vb_p[, , i] <- Vb_p[, j, i] * pw + Vb_p[, j + 1L, i] * (1 - pw)
      Va_p[, , i] <- Va_p[, j, i] * pw + Va_p[, j + 1L, i] * (1 - pw)
    }
  }

  ## --- Step 2: discounted expected marginal values, W(e, b', a') ------------
  ## Pi contracts over TODAY's e.  Because e is the leading axis, the unfolding
  ## to n_e x (n_b*n_a) is a free reshape of the same memory.
  Wb <- pmax(beta * (Pi %*% matrix(Vb_p, n_e, n_b * n_a)), tiny)
  Wa <- pmax(beta * (Pi %*% matrix(Va_p, n_e, n_b * n_a)), tiny)
  dim(Wb) <- dim(Wa) <- c(n_e, n_b, n_a)
  W_ratio <- Wa / Wb

  ## Broadcast helpers on the (e, b', a) and (e, kappa, a) shapes.
  Y_e   <- .hank2_bcast_e(y, n_e, n_b, n_a)
  B_mid <- .hank2_bcast_mid(b_grid, n_e, n_b, n_a)
  A_cur <- .hank2_bcast_a(a_grid, n_e, n_b, n_a)
  Y_ek   <- .hank2_bcast_e(y, n_e, n_k, n_a)
  A_cur_k <- .hank2_bcast_a(a_grid, n_e, n_k, n_a)

  ## --- Step 3: a'(e, b', a) for the LIQUID-UNCONSTRAINED branch -------------
  ## Illiquid FOC: Wa/Wb == 1 + Psi1(a', a).
  cr  <- .hank_lhs_eq_rhs(matrix(W_ratio, n_e * n_b, n_a), rhs)
  a_endo_unc <- array(.hank2_apply_coord_vec(cr$i, cr$p, a_grid),
                      c(n_e, n_b, n_a))
  ## Floor the INTERPOLATED Wb before the power, not just Wb itself (F3). At
  ## the top corner the crossing search returns a weight outside [0, 1], i.e.
  ## it EXTRAPOLATES; since Wb falls in a', a far-enough extrapolation drives
  ## the interpolated value negative, and (negative)^(-eis) is NaN. Flooring
  ## the base keeps every intermediate finite and is an exact no-op wherever
  ## the extrapolation stays positive -- i.e. at every converged solution and
  ## every well-posed calibration.
  c_endo_unc <- pmax(array(.hank2_apply_coord_rows(cr$i, cr$p,
                                                   matrix(Wb, n_e * n_b, n_a)),
                           c(n_e, n_b, n_a)), tiny)^(-eis)

  ## --- Step 4: invert the budget to get b'(e, b, a), a'(e, b, a) ------------
  ## x-budget: c + x' + (1-theta) a' + Psi = y + (1+rb) x + rho_a a.
  psi_unc <- .hank_psi(a_endo_unc, A_cur, ra, chi0, chi1, chi2)$Psi
  b_endo  <- (c_endo_unc + (1 - theta_coll) * a_endo_unc + B_mid - Y_e -
                rho_a * A_cur + psi_unc) / (1 + rb)
  ## Interpolate the b_endo -> b' map at the fixed b_grid.  b' is the MIDDLE
  ## axis, so swap it last (SSJ's swapaxes(1, 2)), interpolate, swap back.
  BE <- matrix(aperm(b_endo,     c(1L, 3L, 2L)), n_e * n_a, n_b)
  AE <- matrix(aperm(a_endo_unc, c(1L, 3L, 2L)), n_e * n_a, n_b)
  co <- .hank2_interp_coord_rows(BE, b_grid)
  a_unc <- aperm(array(.hank2_apply_coord_rows(co$i, co$p, AE),
                       c(n_e, n_a, n_b)), c(1L, 3L, 2L))
  b_unc <- aperm(array(.hank2_apply_coord_vec(co$i, co$p, b_grid),
                       c(n_e, n_a, n_b)), c(1L, 3L, 2L))

  ## --- Step 5: a'(e, kappa, a) for the LIQUID-CONSTRAINED branch ------------
  ## At x' = b_grid[1] the liquid Euler holds with multiplier kappa >= 0:
  ## u'(c) = Wb*(1+kappa), and the illiquid FOC becomes
  ## W_ratio/(1+kappa) == rhs (= (1-theta) + Psi1; the x'-floor does not
  ## involve a', so collateral adds NO multiplier shift here -- that was the
  ## point of the gap coordinate).
  Wr1 <- matrix(W_ratio[, 1L, ], n_e, n_a)
  Wb1 <- matrix(Wb[, 1L, ],      n_e, n_a)
  lhs_con <- array(0, c(n_e, n_k, n_a))
  for (kk in seq_len(n_k)) lhs_con[, kk, ] <- Wr1 / (1 + k_grid[kk])
  crk <- .hank_lhs_eq_rhs(matrix(lhs_con, n_e * n_k, n_a), rhs)
  a_endo_con <- array(.hank2_apply_coord_vec(crk$i, crk$p, a_grid),
                      c(n_e, n_k, n_a))
  ## y is the b'=first slice of Wb, broadcast over kappa (rows index (e, kappa)
  ## with e fastest, so replicate the e-rows n_k times).
  Wb1_rep <- Wb1[rep(seq_len(n_e), times = n_k), , drop = FALSE]
  Kfac <- .hank2_bcast_mid(1 + k_grid, n_e, n_k, n_a)
  c_endo_con <- Kfac^(-eis) *
    pmax(array(.hank2_apply_coord_rows(crk$i, crk$p, Wb1_rep),
               c(n_e, n_k, n_a)), tiny)^(-eis)   # floored base, as in step 3

  ## --- Step 6: map kappa -> b to get a'(e, b, a) on the constrained branch --
  psi_con  <- .hank_psi(a_endo_con, A_cur_k, ra, chi0, chi1, chi2)$Psi
  b_endo_k <- (c_endo_con + (1 - theta_coll) * a_endo_con + b_grid[1L] - Y_ek -
                 rho_a * A_cur_k + psi_con) / (1 + rb)
  ## b_endo_k is DECREASING in kappa; k_grid is decreasing, so b_endo_k is
  ## INCREASING along the kappa AXIS -- which the ascending-knot interpolation
  ## below requires (file header, F1).
  BEK <- matrix(aperm(b_endo_k,   c(1L, 3L, 2L)), n_e * n_a, n_k)
  AEK <- matrix(aperm(a_endo_con, c(1L, 3L, 2L)), n_e * n_a, n_k)
  cok <- .hank2_interp_coord_rows(BEK, b_grid)
  a_con <- aperm(array(.hank2_apply_coord_rows(cok$i, cok$p, AEK),
                       c(n_e, n_a, n_b)), c(1L, 3L, 2L))

  ## --- Step 7: combine, then consumption from the budget residual -----------
  ## Take the constrained branch wherever the unconstrained one would violate
  ## the liquid floor.  (SSJ clamps only b; a' is left unclamped and the Young
  ## lottery handles any overshoot downstream -- brief 19, F2.)
  a_pol <- a_unc; b_pol <- b_unc
  con <- b_unc <= b_grid[1L]
  b_pol[con] <- b_grid[1L]
  a_pol[con] <- a_con[con]

  ps    <- .hank_psi(a_pol, A_cur, ra, chi0, chi1, chi2)
  c_pol <- pmax(Y_e + (1 + rb) * B_mid + rho_a * A_cur - ps$Psi -
                  (1 - theta_coll) * a_pol - b_pol, tiny)
  uc <- c_pol^(-1 / eis)
  list(Vb = (1 + rb) * uc,
       Va = (rho_a - ps$Psi2) * uc,
       b = b_pol, a = a_pol, c = c_pol, chi = ps$Psi)
}


#' Unclamped bilinear interpolation of one income-state's grid function
#'
#' Evaluates \code{W[e, , ]} at query points \code{(bq, aq)} by bilinear
#' interpolation with LINEAR EXTRAPOLATION outside the grid (weights are not
#' clamped, only the bracketing indices) -- the \code{\link{.hank_interp1}}
#' convention, deliberately NOT the \code{\link{.hank_lottery}} one.  A lottery
#' must clamp (mass cannot leave the grid), but an oracle that evaluated the
#' marginal value AT THE TOP GRIDPOINT instead of at the policy would charge
#' the solver a spurious residual precisely at the cells where the policy
#' overshoots.
#' @keywords internal
.hank2_bilin <- function(W_e, b_grid, a_grid, bq, aq) {
  n_b <- length(b_grid); n_a <- length(a_grid)
  bi <- pmin(pmax(findInterval(bq, b_grid), 1L), n_b - 1L)
  ai <- pmin(pmax(findInterval(aq, a_grid), 1L), n_a - 1L)
  pb <- (b_grid[bi + 1L] - bq) / (b_grid[bi + 1L] - b_grid[bi])
  pa <- (a_grid[ai + 1L] - aq) / (a_grid[ai + 1L] - a_grid[ai])
  pb * pa * W_e[cbind(bi, ai)] +
    (1 - pb) * pa * W_e[cbind(bi + 1L, ai)] +
    pb * (1 - pa) * W_e[cbind(bi, ai + 1L)] +
    (1 - pb) * (1 - pa) * W_e[cbind(bi + 1L, ai + 1L)]
}


#' First-order-condition residuals of a converged two-asset policy
#'
#' Self-contained correctness oracle for \code{\link{hank_egm2_solve}}, the
#' two-asset counterpart of \code{\link{hank_euler_residual}}.  Two conditions
#' are reported, evaluated at the CHOSEN portfolio \eqn{(b', a')}:
#' \itemize{
#'   \item the LIQUID Euler \eqn{u'(c) = W_b(e, b', a')}, an equality only
#'     where the liquid constraint is slack (at the floor the multiplier
#'     \eqn{\kappa \ge 0} makes it \eqn{u'(c) \ge W_b});
#'   \item the ILLIQUID FOC \eqn{W_a(e, b', a') = u'(c)\,(1 + \Psi_1(a', a))},
#'     an equality wherever \eqn{a'} is interior.
#' }
#' The continuation values are RECONSTRUCTED FROM THE POLICIES rather than
#' read from the stored marginals: next-period consumption \eqn{c'} and
#' illiquid choice \eqn{a''} are interpolated at \eqn{(b', a')} per future
#' income state, and the envelope forms
#' \eqn{V_b' = (1+r_b)\,u'(c')} and
#' \eqn{V_a' = (1 + r_a - \Psi_2(a'', a'))\,u'(c')} are then applied
#' analytically before the \eqn{\Pi}-expectation.  Two reasons.  First,
#' independence: the oracle then checks the POLICIES against the FOCs without
#' trusting the solver's own \eqn{V_b}/\eqn{V_a} bookkeeping.  Second,
#' accuracy: policies are nearly linear in wealth while \eqn{u'(c) = c^{-1/eis}}
#' is violently convex, so interpolating the composite instead of the policy
#' commits far larger error (measured ~60x at the reference grids) and charges
#' it to the solver.  This mirrors \code{\link{hank_euler_residual}}, which
#' interpolates the consumption policy for the same reason.
#'
#' \strong{Why \eqn{u'(c)} appears on BOTH right-hand sides.}  The naive
#' illiquid form \eqn{W_a = W_b (1 + \Psi_1)} holds only where the liquid
#' constraint is slack: at the floor the solver's optimality conditions are
#' \eqn{u'(c) = (1+\kappa) W_b} and \eqn{W_a = (1+\kappa) W_b (1 + \Psi_1)}
#' with \eqn{\kappa > 0}, so comparing \eqn{W_a} to \eqn{W_b(1+\Psi_1)} charges
#' an \eqn{O(\kappa)} "residual" to cells that are exactly optimal.  Written
#' against \eqn{u'(c)} the illiquid condition is branch-universal.  (This is
#' the defect that sank the first, cut implementation of this function -- its
#' illiquid residual GREW under grid refinement because the constrained region
#' resolves better as the grid refines; see briefs/19, F20.)
#'
#' The residuals are not zero by construction: the two-stage EGM interpolates
#' its endogenous grids back onto the fixed ones, so what remains measures the
#' DISCRETIZATION error, and watching it fall under refinement is the intended
#' use (a grid-sizing diagnostic complementing
#' \code{\link{hank_twoasset_grid_check}}).  The solver-correctness oracle
#' proper is the degenerate reduction in \code{test-hank-egm2.R}.
#'
#' @param hh A solved household from \code{\link{hank_egm2_solve}} or a
#'   \code{\link{hank_het2_block}} (both carry the required fields).
#' @param constraint_tol Cells whose liquid policy is within
#'   \code{constraint_tol} of \code{b_grid[1]} are excluded from the liquid
#'   residual; cells whose illiquid policy is within it of \code{a_grid[1]}
#'   are excluded from the illiquid residual.
#'
#' @return A list with \code{max_abs_liquid} / \code{max_abs_illiquid} (max
#'   absolute residuals over unconstrained cells) and the full
#'   \code{n_e x n_b x n_a} arrays \code{liquid} / \code{illiquid} (\code{NA}
#'   where the relevant constraint binds).
#' @seealso \code{\link{hank_egm2_solve}}, \code{\link{hank_euler_residual}}
#'   (one-asset), \code{\link{hank_twoasset_grid_check}}
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' hh  <- hank_egm2_solve(hank_asset_grid(20, 15, 0), hank_asset_grid(30, 15, 0),
#'                        y = inc$e, rb = 0.005, ra = 0.02, beta = 0.95,
#'                        eis = 0.5, chi0 = 0.25, chi1 = 6.5, chi2 = 2,
#'                        Pi = inc$Pi, n_k = 10)
#' hank_euler2_residual(hh)$max_abs_liquid
#' @export
hank_euler2_residual <- function(hh, constraint_tol = 1e-8) {
  ## A het2d block carries every field this check looks for (its MIXED
  ## policies), so the field test alone would silently accept it -- and the
  ## smooth FOCs are the wrong conditions for a discrete-choice household
  ## (the mixture satisfies them nowhere). Reject loudly (the F19 lesson).
  if (inherits(hh, "hank_het2d_block"))
    stop("hank_euler2_residual: this is a discrete-adjustment block ",
         "(hank_het2d_block); its mixed policies do not satisfy the smooth ",
         "FOCs (P-weighted mixtures satisfy neither branch's condition). ",
         "There is no het2d FOC-residual oracle; use the block's own gates ",
         "(budget identities, reductions, ND-vs-fake-news).")
  need <- c("b", "a", "c", "b_grid", "a_grid", "Pi", "beta", "eis",
            "rb", "ra", "chi0", "chi1", "chi2")
  if (!is.list(hh) || !all(need %in% names(hh)))
    stop("hank_euler2_residual: 'hh' must be a solved two-asset household ",
         "(from hank_egm2_solve) or a hank_het2_block.")
  n_e <- dim(hh$c)[1L]; n_b <- dim(hh$c)[2L]; n_a <- dim(hh$c)[3L]

  ## COLLATERAL (D1, brief 19 section 9.12): under theta_coll > 0 the object's
  ## policies live in gap coordinates x = b + theta*a, in which the FOCs keep
  ## their constant-floor FORM with two re-weighted coefficients: the illiquid
  ## arrival envelope carries rho_a = (1+ra) - theta(1+rb), and the illiquid
  ## outlay costs (1-theta) liquid units.  Both reduce bit-for-bit to the
  ## pre-collateral expressions at theta = 0 (pre-D1 objects lack the field;
  ## default it to 0).
  th <- if (is.null(hh$theta_coll)) 0 else hh$theta_coll
  rho_a <- (1 + hh$ra) - th * (1 + hh$rb)

  ## Reconstruct the continuation marginal values from the POLICIES at the
  ## chosen (b', a'), per future income state ep, then take the
  ## Pi-expectation. See the roxygen for why policies are interpolated and
  ## the envelope transforms applied analytically (independence + accuracy).
  Wb_at <- array(0, c(n_e, n_b, n_a)); Wa_at <- Wb_at
  for (e in seq_len(n_e)) {
    bq <- as.numeric(hh$b[e, , ])          # today's chosen (b', a') from (e, x)
    aq <- as.numeric(hh$a[e, , ])
    EVb <- numeric(n_b * n_a); EVa <- numeric(n_b * n_a)
    for (ep in seq_len(n_e)) {
      cp  <- .hank2_bilin(matrix(hh$c[ep, , ], n_b, n_a),
                          hh$b_grid, hh$a_grid, bq, aq)
      app <- .hank2_bilin(matrix(hh$a[ep, , ], n_b, n_a),
                          hh$b_grid, hh$a_grid, bq, aq)
      ucp <- pmax(cp, tiny_floor())^(-1 / hh$eis)
      ## envelope at the QUERY state: tomorrow's illiquid state is today's
      ## choice a', tomorrow's illiquid choice is the interpolated a''
      Psi2p <- .hank_psi(app, aq, hh$ra, hh$chi0, hh$chi1, hh$chi2)$Psi2
      EVb <- EVb + hh$Pi[e, ep] * (1 + hh$rb) * ucp
      EVa <- EVa + hh$Pi[e, ep] * (rho_a - Psi2p) * ucp
    }
    Wb_at[e, , ] <- hh$beta * EVb
    Wa_at[e, , ] <- hh$beta * EVa
  }

  uc <- hh$c^(-1 / hh$eis)
  A_cur <- .hank2_bcast_a(hh$a_grid, n_e, n_b, n_a)
  Psi1 <- .hank_psi(hh$a, A_cur, hh$ra, hh$chi0, hh$chi1, hh$chi2)$Psi1

  liquid   <- uc - Wb_at
  illiquid <- Wa_at - uc * ((1 - th) + Psi1)

  b_con <- hh$b <= (hh$b_grid[1L] + constraint_tol)
  a_con <- hh$a <= (hh$a_grid[1L] + constraint_tol)
  liquid[b_con]   <- NA_real_
  illiquid[a_con] <- NA_real_
  list(max_abs_liquid   = max(abs(liquid),   na.rm = TRUE),
       max_abs_illiquid = max(abs(illiquid), na.rm = TRUE),
       liquid = liquid, illiquid = illiquid)
}


#' Solve the two-asset household problem to a steady-state policy
#'
#' Iterates \code{\link{.hank_egm2_step}} at constant prices until the liquid
#' and illiquid policies stop moving.  The two-asset counterpart of
#' \code{\link{hank_egm_solve}}: a liquid asset \code{b} (return \code{rb},
#' floor \code{b_grid[1]}) and an illiquid asset \code{a} (return \code{ra})
#' that costs \eqn{\Psi(a', a)} to move (\code{\link{.hank_psi}}).  With
#' \code{ra > rb} this produces the wealthy-hand-to-mouth households a
#' one-asset model cannot: rich in illiquid wealth, constrained in liquid.
#'
#' Ported from and validated against SSJ's \code{hh_twoasset.py} (see
#' \code{briefs/19-twoasset-hank-scope.md}).
#'
#' @param b_grid Numeric: increasing LIQUID grid (see
#'   \code{\link{hank_asset_grid}}); \code{b_grid[1]} is the liquid floor.
#' @param a_grid Numeric: increasing ILLIQUID grid.
#' @param y Numeric length-\code{n_e}: labour income by income state.
#' @param rb,ra Numeric: liquid and illiquid returns.
#' @param beta,eis Discount factor and elasticity of intertemporal substitution.
#' @param chi0,chi1,chi2 Adjustment-cost parameters: \code{chi0 > 0},
#'   \code{chi1 >= 0}, \code{chi2 > 1}.  Reference calibration:
#'   \code{chi0 = 0.25}, \code{chi1 = 6.5}, \code{chi2 = 2}.
#' @param Pi Numeric \code{n_e x n_e}: row-stochastic income transition matrix.
#' @param k_grid Optional DECREASING multiplier grid for the liquid-constrained
#'   branch.  Default \code{NULL} builds SSJ's
#'   \code{rev(hank_asset_grid(amax = k_max, n = n_k))}.  It MUST be decreasing:
#'   endogenous liquid assets fall in the multiplier, so a decreasing grid is
#'   what makes them increase along the axis, as step 6's interpolation
#'   requires.
#' @param n_k,k_max Integer/numeric: size and top of the default multiplier grid
#'   (defaults 50 / 1; SSJ's own test config uses \code{n_k = 4}).
#' @param tol,maxit Convergence tolerance on \eqn{\max(|db'|, |da'|)} and
#'   iteration cap.
#' @param Vb_init,Va_init Optional \code{n_e x n_b x n_a} initial marginal
#'   values.  Default \code{NULL} uses SSJ's \code{hh_init} guesses.
#' @param backend Character: \code{"cpp"} (default) or \code{"R"}.  All
#'   validation happens in this wrapper, so both backends reject identical
#'   inputs; \code{"R"} is the reference path.
#' @param theta_coll Scalar LTV in \eqn{[0, 1)} on end-of-period illiquid
#'   collateral (D1): households may borrow up to
#'   \eqn{b' \ge b_{grid}[1] - \theta_{coll}\, a'}.  Implemented by the GAP
#'   COORDINATE \eqn{x = b + \theta a}: when \code{theta_coll > 0},
#'   \code{b_grid} is the grid of \eqn{x} (liquid wealth PLUS collateral
#'   capacity, floor = the unsecured limit), the returned \code{b}/\code{Vb}
#'   policies live in \eqn{x}, and the TRUE liquid position is the additional
#'   returned array \code{b_liq = b - theta_coll * a}.  Collateral runs on the
#'   R reference path only (the compiled kernel is roadmap).  See brief 19
#'   section 9.12 for why the naive \eqn{b}-coordinate implementation is
#'   unsound (the binding boundary must be a grid line).
#' @param threads Worker threads for the COMPILED backend's backward step, or
#'   \code{NULL} (default) to resolve from
#'   \code{getOption("dynhr.hank3_threads")} and then a machine-derived default
#'   (see \code{\link{hank_resolve_threads}}, shared with the three-asset
#'   family); \code{1} forces the serial path.  The kernel spawns its pool ONCE
#'   per solve and parks it on a barrier between backward iterations, so the
#'   per-step overhead is nil even though a two-asset step is only ~1 ms.  Every
#'   parallel write is disjoint and step 2's contraction over income stays on
#'   the main thread, so the result is BIT-IDENTICAL at every thread count --
#'   pinned by \code{identical()} in \code{test-hank-egm2-threads.R}.  Ignored
#'   by \code{backend = "R"} and by the collateral path (which forces it).
#'
#' @return A list with \code{Vb}, \code{Va}, \code{b}, \code{a}, \code{c},
#'   \code{chi} (each \code{n_e x n_b x n_a}), \code{b_liq} (true liquid
#'   policy; identical to \code{b} at \code{theta_coll = 0}),
#'   \code{iterations}, \code{converged}, and the echoed calibration.
#' @seealso \code{\link{hank_het2_block}}, \code{\link{hank_egm_solve}} (the
#'   one-asset solver)
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' bg  <- hank_asset_grid(10, 20, 0)
#' ag  <- hank_asset_grid(40, 20, 0)
#' hh  <- hank_egm2_solve(bg, ag, y = inc$e, rb = 0.005, ra = 0.02,
#'                        beta = 0.95, eis = 0.5, chi0 = 0.25, chi1 = 6.5,
#'                        chi2 = 2, Pi = inc$Pi)
#' range(hh$b)
#' @export
hank_egm2_solve <- function(b_grid, a_grid, y, rb, ra, beta, eis,
                            chi0, chi1, chi2, Pi,
                            k_grid = NULL, n_k = 50L, k_max = 1,
                            tol = 1e-10, maxit = 5000L,
                            Vb_init = NULL, Va_init = NULL,
                            backend = getOption("dynhr.hank_backend", "cpp"),
                            theta_coll = 0, threads = NULL) {
  backend <- match.arg(backend, c("R", "cpp"))
  if (!is.numeric(theta_coll) || length(theta_coll) != 1L ||
      !is.finite(theta_coll) || theta_coll < 0 || theta_coll >= 1)
    stop("hank_egm2_solve: 'theta_coll' must be a scalar in [0, 1) (the LTV ",
         "on end-of-period illiquid collateral; 0 disables it).")
  ## collateral lives on the R reference path only; the validated cpp kernel
  ## is untouched (roadmap)
  if (theta_coll > 0) backend <- "R"
  ## ---- validation (runs for BOTH backends) --------------------------------
  chk_grid <- function(g, nm) {
    if (!is.numeric(g) || length(g) < 2L || !all(is.finite(g)))
      stop("hank_egm2_solve: '", nm, "' must be a finite numeric grid of ",
           "length >= 2.")
    if (any(diff(g) <= 0))
      stop("hank_egm2_solve: '", nm, "' must be strictly increasing.")
  }
  chk_grid(b_grid, "b_grid"); chk_grid(a_grid, "a_grid")
  n_b <- length(b_grid); n_a <- length(a_grid); n_e <- length(y)
  if (!is.numeric(y) || !all(is.finite(y)))
    stop("hank_egm2_solve: 'y' must be a finite numeric vector.")
  .hank_check_markov(Pi, n_e, caller = "hank_egm2_solve")
  chk_scalar <- function(x, nm) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x))
      stop("hank_egm2_solve: '", nm, "' must be a finite numeric scalar.")
  }
  for (nm in c("rb", "ra", "beta", "eis", "chi0", "chi1", "chi2"))
    chk_scalar(get(nm), nm)
  if (chi0 <= 0)
    stop("hank_egm2_solve: 'chi0' must be > 0 (it is the denominator shift ",
         "that keeps the adjustment cost finite at a = 0).")
  if (chi1 < 0) stop("hank_egm2_solve: 'chi1' must be >= 0.")
  if (chi2 <= 1)
    stop("hank_egm2_solve: 'chi2' must be > 1 (the adjustment cost must be ",
         "strictly convex in a'; chi2 = 1 is the kinked/linear cost, which ",
         "this solver's illiquid FOC crossing search does not handle -- see ",
         "briefs/19-twoasset-hank-scope.md).")
  if (beta <= 0 || beta >= 1) stop("hank_egm2_solve: 'beta' must be in (0, 1).")
  if (eis <= 0) stop("hank_egm2_solve: 'eis' must be > 0.")
  if (!is.finite(tol) || tol <= 0)
    stop("hank_egm2_solve: 'tol' must be a finite positive scalar.")
  if (!is.finite(maxit) || maxit < 1L)
    stop("hank_egm2_solve: 'maxit' must be a finite positive integer.")

  ## Multiplier grid: built decreasing, or validated decreasing when supplied.
  if (is.null(k_grid)) {
    if (!is.finite(n_k) || n_k < 2L)
      stop("hank_egm2_solve: 'n_k' must be an integer >= 2.")
    if (!is.finite(k_max) || k_max <= 0)
      stop("hank_egm2_solve: 'k_max' must be a positive scalar.")
    k_grid <- rev(hank_asset_grid(amax = k_max, n = as.integer(n_k), amin = 0))
  } else {
    if (!is.numeric(k_grid) || length(k_grid) < 2L || !all(is.finite(k_grid)))
      stop("hank_egm2_solve: 'k_grid' must be a finite numeric vector of ",
           "length >= 2.")
    if (any(diff(k_grid) >= 0))
      stop("hank_egm2_solve: 'k_grid' must be strictly DECREASING. Endogenous ",
           "liquid assets fall in the multiplier kappa, so only a decreasing ",
           "kappa grid makes them increase along the axis -- which the ",
           "constrained branch's ascending-knot interpolation requires. Use ",
           "rev(hank_asset_grid(amax = k_max, n = n_k)).")
    if (min(k_grid) < 0)
      stop("hank_egm2_solve: 'k_grid' must be non-negative (it is the ",
           "multiplier on the liquid constraint).")
  }

  ## Feasibility at the worst cell: the poorest household, sitting at both grid
  ## floors, must be able to stay there with positive consumption.  Otherwise
  ## c <= 0 forces a NaN marginal value that only surfaces later as a cryptic
  ## interpolation failure (same guard philosophy as hank_egm_solve's min_coh).
  ## In gap coordinates staying put in (x, a) IS staying put in (b, a), and
  ## the interest on the collateralized slice nets to (ra - theta*rb)*a; at
  ## theta_coll = 0 the expression is the pre-collateral one bit-for-bit.
  psi_floor <- .hank_psi(a_grid[1L], a_grid[1L], ra, chi0, chi1, chi2)$Psi
  c_floor <- min(y) + rb * b_grid[1L] +
    (ra - theta_coll * rb) * a_grid[1L] - psi_floor
  if (!is.finite(c_floor) || c_floor <= 0)
    stop("hank_egm2_solve: the poorest household at both grid floors cannot ",
         "consume a positive amount while staying there (min(y) + rb*b_grid[1]",
         " + (ra - theta_coll*rb)*a_grid[1] - Psi = ", format(c_floor),
         " <= 0). Raise the income floor or the grid floors.")

  ## SSJ hh_init guesses (levels raised to -1/eis), broadcast over income states.
  Bb <- .hank2_bcast_mid(b_grid, n_e, n_b, n_a)
  Aa <- .hank2_bcast_a(a_grid, n_e, n_b, n_a)
  ## The SSJ bases are affine in (b, a) and positive only for a NON-BORROWING
  ## liquid grid: 0.6 + 1.1*b + a goes negative once b < -0.545 (at a = 0).  A
  ## negative base raised to a non-integer -1/eis is NaN, so chk_init below
  ## rejects the default init and the block cannot be built at all; and at an
  ## INTEGER -1/eis (eis = 0.5 gives exactly -2) R instead returns a finite
  ## POSITIVE number from a nonsense negative "resource" level, which is worse
  ## -- the bad init is silent.  Floor the base at the worst-case affordable
  ## consumption (c_floor, already proven positive above; capped at 0.5 so the
  ## floor can never bite on a non-borrowing grid, where the bases are >= 0.6
  ## and >= 0.5 respectively).  Byte-identical for every b_min >= 0, a_min >= 0
  ## calibration.  Same guard philosophy as hank_egm_solve's default Va.
  init_floor <- min(0.5, c_floor)
  if (is.null(Va_init))
    Va_init <- pmax(0.6 + 1.1 * Bb + Aa, init_floor)^(-1 / eis)
  if (is.null(Vb_init))
    Vb_init <- pmax(0.5 + Bb + 1.2 * Aa, init_floor)^(-1 / eis)
  chk_init <- function(V, nm) {
    if (!is.numeric(V) || !identical(dim(V), c(n_e, n_b, n_a)) ||
        !all(is.finite(V)))
      stop("hank_egm2_solve: '", nm, "' must be a finite numeric ",
           n_e, " x ", n_b, " x ", n_a, " array.")
  }
  chk_init(Va_init, "Va_init"); chk_init(Vb_init, "Vb_init")

  ## Psi1(a', a) depends only on (a_grid, ra, chi*), so build it once here
  ## rather than per iteration (SSJ's marginal_cost_grid hetinput).
  Psi1_grid <- .hank_psi(matrix(a_grid, n_a, n_a),
                         matrix(a_grid, n_a, n_a, byrow = TRUE),
                         ra, chi0, chi1, chi2)$Psi1

  ## Timed from AFTER validation, matching hank_egm3_solve/hank_egm_solve, so
  ## the three families' elapsed_solve fields mean the same thing.
  .t0 <- proc.time()[["elapsed"]]
  ## The resolved count is recorded, not the `threads` ARGUMENT: NULL is the
  ## common case and reporting NULL would tell a manifest reader nothing about
  ## what actually ran. On the R branch the answer is honestly 1.
  n_thr <- if (backend == "cpp") hank_resolve_threads(threads) else 1L
  ## The R branch measures the policy gap its own convergence test reads; the
  ## compiled loop tests internally and returns only the policies, so its gap
  ## is genuinely unavailable -- NA rather than a re-derived DIFFERENT quantity
  ## reported under the same name.
  gap_final <- NA_real_

  if (backend == "cpp") {
    ## Resolved with the shared HANK resolver deliberately: one resolution
    ## order, one option, one reporter for the whole family. Resolving only on
    ## the cpp branch keeps the R reference path free of any thread machinery.
    out <- hank_egm2_solve_cpp(Vb_init, Va_init, b_grid, a_grid, k_grid, y,
                               rb, ra, beta, eis, chi0, chi1, chi2, Pi,
                               Psi1_grid, tol, as.integer(maxit), n_thr)
  } else {
    Vb <- Vb_init; Va <- Va_init
    b_old <- NULL; a_old <- NULL
    converged <- FALSE; it <- 0L; step <- NULL
    for (it in seq_len(as.integer(maxit))) {
      step <- .hank_egm2_step(Vb, Va, b_grid, a_grid, k_grid, y, rb, ra,
                              beta, eis, chi0, chi1, chi2, Pi,
                              Psi1_grid = Psi1_grid, theta_coll = theta_coll)
      Vb <- step$Vb; Va <- step$Va
      if (!is.null(b_old)) {
        gap <- max(max(abs(step$b - b_old)), max(abs(step$a - a_old)))
        ## NaN-aware: a non-finite iterate must never report convergence (the
        ## cpp loop breaks on the same condition; the shared check below then
        ## raises the identical error for both backends).
        if (!is.finite(gap)) break
        gap_final <- gap
        if (gap < tol) { converged <- TRUE; break }
      }
      b_old <- step$b; a_old <- step$a
    }
    out <- c(step, list(iterations = it, converged = converged))
  }
  .elapsed <- proc.time()[["elapsed"]] - .t0

  ## ---- shared post-checks (identical for both backends) -------------------
  if (!all(is.finite(out$b)) || !all(is.finite(out$a)) ||
      !all(is.finite(out$c)))
    stop("hank_egm2_solve: non-finite policy iterate after ", out$iterations,
         " iterations -- the calibration is infeasible or the grids are too ",
         "coarse.")

  ## The `tiny` consumption floor is a TRANSIENT guard: if it is still binding
  ## at the fixed point, the household cannot afford the portfolio it chose,
  ## the returned policies do not solve the stated problem, and the marginal
  ## values are meaningless (c = 1e-12 with eis = 0.5 gives Vb ~ 1e24). The
  ## one-asset solver rules this out up front with its min_coh check, but here
  ## feasibility depends on the ENDOGENOUS portfolio -- an adjustment cost can
  ## exceed income at cells no pre-check can identify -- so it must be caught
  ## at convergence instead. Failing loudly beats returning converged = TRUE
  ## with a 1e24 marginal value that silently poisons every downstream
  ## Jacobian.
  if (min(out$c) <= tiny_floor())
    stop("hank_egm2_solve: consumption hit the numerical floor (",
         format(tiny_floor()), ") at ", sum(out$c <= tiny_floor()), " of ",
         length(out$c), " cells at the converged policy, so the returned ",
         "policy does not solve the stated problem. ",
         "MOST OFTEN THIS IS RESOLUTION, NOT ECONOMICS: this solver forms c ",
         "as a budget residual from separately interpolated a' and b', so a ",
         "coarse illiquid grid relative to its span accumulates enough ",
         "interpolation error at the richest cells to push c negative. Try ",
         "raising n_a FIRST (a_grid currently has ", length(a_grid),
         " points spanning [", format(a_grid[1L]), ", ",
         format(a_grid[length(a_grid)]), "]; widening the span REQUIRES ",
         "refining with it, since the geomspace top gaps grow). ",
         "Only if refining does not help is the calibration itself ",
         "infeasible -- then the adjustment cost is overwhelming income ",
         "(Psi grows like chi1*|da|^chi2 / ((1+ra)a + chi0)^(chi2-1), so ",
         "raising chi2 above 2 at an unchanged chi1 blows it up): reduce ",
         "chi1 or raise chi0. See briefs/19-twoasset-hank-scope.md F16.")

  ## True liquid policy: b = x - theta*a. At theta_coll = 0 the x and b
  ## coordinates coincide, and b_liq IS out$b (the same object, no copy).
  b_liq <- if (theta_coll > 0) out$b - theta_coll * out$a else out$b

  c(out, list(b_liq = b_liq,
              b_grid = b_grid, a_grid = a_grid, k_grid = k_grid, y = y,
              rb = rb, ra = ra, beta = beta, eis = eis,
              chi0 = chi0, chi1 = chi1, chi2 = chi2,
              theta_coll = theta_coll, Pi = Pi,
              ## Run metadata hank_het2_manifest() reads. `backend` is what
              ## ACTUALLY ran, `threads` the RESOLVED count (1 on the R path),
              ## and last_policy_gap is NA under "cpp" for the reason given at
              ## its initialisation above.
              backend = backend, threads = n_thr, elapsed = .elapsed,
              last_policy_gap = gap_final))
}
