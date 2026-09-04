## R/hank-distribution2.R
## --------------------------------------------------------------------------
## Young's-method distribution machinery for the TWO-ASSET (liquid/illiquid)
## household: the joint forward operator over the (e, b, a) cell space, its
## flattening conventions, and mass-weighted aggregation.
##
## CELL ORDER (the single convention for everything two-asset; see also
## R/hank-egm2.R):
##   index(e, b, a) = (e-1)*n_b*n_a + (b-1)*n_a + a
## i.e. income SLOWEST, liquid b in the middle, ILLIQUID a FASTEST. This
## extends the one-asset rule (income slow / asset fast, see
## R/hank-distribution.R) and matches SSJ's numpy (z, b, a) row-major layout.
## At n_a = 1 it degenerates exactly to the one-asset order over (e, b).
##
## Policies/marginals are n_e x n_b x n_a ARRAYS. Because R is column-major and
## e is the leading axis, .hank2_arr_to_vec() is a single aperm -- and the
## (e, mid)-by-a unfoldings used by the EGM step are free reshapes.
##
## The joint lottery is the PRODUCT of two independent 1-D lotteries (one per
## asset), so each cell scatters to 4 asset destinations x n_e income
## destinations. It reuses the validated one-asset .hank_lottery() per axis
## rather than reimplementing the bracketing/clamping (that helper is already
## shape-agnostic: it reads dim(a_pol)).
##
## The resulting Lambda is an ordinary sparse row-stochastic operator over the
## joint cell space, so hank_stationary_dist() -- including its compiled CSC
## kernel -- works on it UNCHANGED. Only the block constructor's hot path needs
## a fused two-asset kernel.
## --------------------------------------------------------------------------


#' Flatten an (e, b, a) array to the package's two-asset distribution order
#'
#' \code{index(e, b, a) = (e-1)*n_b*n_a + (b-1)*n_a + a} (income slowest,
#' illiquid \code{a} fastest).  The two-asset counterpart of
#' \code{\link{.hank_mat_to_vec}}.
#'
#' @param X Numeric \code{n_e x n_b x n_a} array.
#' @return Numeric length-\code{n_e*n_b*n_a} vector.
#' @keywords internal
.hank2_arr_to_vec <- function(X) as.numeric(aperm(X, c(3L, 2L, 1L)))


#' Inverse of \code{\link{.hank2_arr_to_vec}}
#'
#' @param v Numeric length-\code{n_e*n_b*n_a} vector in distribution order.
#' @param n_e,n_b,n_a Integer dimensions.
#' @return Numeric \code{n_e x n_b x n_a} array.
#' @keywords internal
.hank2_vec_to_arr <- function(v, n_e, n_b, n_a)
  aperm(array(v, c(n_a, n_b, n_e)), c(3L, 2L, 1L))


#' Mass-weighted aggregate over the two-asset cell space
#'
#' \eqn{\sum_x d(x)\, g(x)} over \code{(e, b, a)} cells.  The two-asset sibling
#' of \code{\link{hank_aggregate}}, and it inherits that function's contract:
#' shapes must match cell-for-cell, and mismatches are an ERROR rather than a
#' silent recycle (adversarial review 2026-07-13, P2).
#'
#' @param d Distribution: a length-\code{n_e*n_b*n_a} vector in the package's
#'   two-asset order, or an \code{n_e x n_b x n_a} array.
#' @param x Values to average: same shape rules as \code{d}.
#'
#' @return The scalar aggregate.
#' @seealso \code{\link{hank_aggregate}} (one-asset),
#'   \code{\link{hank_het2_block}}
#' @export
hank_aggregate2 <- function(d, x) {
  flat <- function(z, nm) {
    if (is.array(z) && length(dim(z)) == 3L) return(.hank2_arr_to_vec(z))
    if (is.numeric(z) && is.null(dim(z))) return(as.numeric(z))
    stop("hank_aggregate2(): '", nm, "' must be a numeric vector in the ",
         "two-asset cell order or an n_e x n_b x n_a array.")
  }
  ## When both are arrays the dims must agree exactly; a transposed array has
  ## the same cell count but a different meaning.
  if (is.array(d) && length(dim(d)) == 3L &&
      is.array(x) && length(dim(x)) == 3L && !identical(dim(d), dim(x)))
    stop("hank_aggregate2(): 'd' and 'x' must have identical dimensions ",
         "(got ", paste(dim(d), collapse = " x "), " and ",
         paste(dim(x), collapse = " x "), ").")
  dv <- flat(d, "d"); xv <- flat(x, "x")
  if (length(dv) != length(xv))
    stop("hank_aggregate2(): 'd' (", length(dv), " cells) and 'x' (",
         length(xv), " cells) must cover the same cell space; refusing to ",
         "recycle.")
  if (!all(is.finite(dv)) || !all(is.finite(xv)))
    stop("hank_aggregate2(): 'd' and 'x' must be finite.")
  sum(dv * xv)
}


#' Build the sparse joint forward operator for a two-asset household
#'
#' The two-asset counterpart of \code{\link{hank_forward_operator}}: the
#' row-stochastic \code{(n_e*n_b*n_a) x (n_e*n_b*n_a)} transition matrix over
#' the joint \code{(e, b, a)} cell space, given both savings policies.  Each
#' cell's mass is split by the PRODUCT of two independent Young lotteries (one
#' per asset), giving 4 asset destinations, and then spread across income
#' states by \code{Pi}.  As in the one-asset case \code{Lambda[from, to]} and
#' distributions push forward as \code{t(Lambda) \%*\% d}.
#'
#' The result is an ordinary sparse row-stochastic operator, so
#' \code{\link{hank_stationary_dist}} (and its compiled kernel) applies to it
#' unchanged.
#'
#' @param b_pol,a_pol Numeric \code{n_e x n_b x n_a} arrays: the liquid and
#'   illiquid policies.  Values outside their grids are clamped by the lottery
#'   (a lottery must not create negative mass) -- which is how the deliberate
#'   absence of a policy clamp in \code{\link{.hank_egm2_step}} is absorbed.
#' @param b_grid,a_grid Numeric: the increasing liquid / illiquid grids.
#' @param Pi Numeric \code{n_e x n_e}: row-stochastic income transition matrix.
#'
#' @return A sparse \code{dgCMatrix}, row-stochastic, over the joint cell space.
#' @seealso \code{\link{hank_forward_operator}} (one-asset),
#'   \code{\link{hank_stationary_dist}}, \code{\link{hank_het2_block}}
#' @export
hank_forward_operator2 <- function(b_pol, a_pol, b_grid, a_grid, Pi) {
  chk_pol <- function(P, nm) {
    if (!is.array(P) || length(dim(P)) != 3L || !is.numeric(P) ||
        !all(is.finite(P)))
      stop("hank_forward_operator2(): '", nm, "' must be a finite numeric ",
           "n_e x n_b x n_a array.")
  }
  chk_pol(b_pol, "b_pol"); chk_pol(a_pol, "a_pol")
  if (!identical(dim(b_pol), dim(a_pol)))
    stop("hank_forward_operator2(): 'b_pol' and 'a_pol' must have identical ",
         "dimensions (got ", paste(dim(b_pol), collapse = " x "), " and ",
         paste(dim(a_pol), collapse = " x "), ").")
  n_e <- dim(b_pol)[1L]; n_b <- dim(b_pol)[2L]; n_a <- dim(b_pol)[3L]
  chk_grid <- function(g, n, nm) {
    if (!is.numeric(g) || length(g) != n || !all(is.finite(g)) ||
        (n > 1L && any(diff(g) <= 0)))
      stop("hank_forward_operator2(): '", nm, "' must be a finite, strictly ",
           "increasing numeric vector of length ", n, ".")
  }
  chk_grid(b_grid, n_b, "b_grid"); chk_grid(a_grid, n_a, "a_grid")
  .hank_check_markov(Pi, n_e, caller = "hank_forward_operator2")

  ## Per-axis lotteries. .hank_lottery is shape-agnostic (it uses dim(a_pol)),
  ## so it applies to the (e*b) x a and (e*a) x b unfoldings directly; work on
  ## the (n_e*n_b) x n_a and (n_e*n_a) x n_b views and reassemble in cell order.
  lb <- .hank_lottery(matrix(aperm(b_pol, c(1L, 3L, 2L)), n_e * n_a, n_b),
                      b_grid)
  la <- .hank_lottery(matrix(a_pol, n_e * n_b, n_a), a_grid)
  ## Put both back on the (e, b, a) cell layout, flattened in cell order.
  ib <- .hank2_arr_to_vec(aperm(array(lb$i, c(n_e, n_a, n_b)), c(1L, 3L, 2L)))
  pb <- .hank2_arr_to_vec(aperm(array(lb$p, c(n_e, n_a, n_b)), c(1L, 3L, 2L)))
  ia <- .hank2_arr_to_vec(array(la$i, c(n_e, n_b, n_a)))
  pa <- .hank2_arr_to_vec(array(la$p, c(n_e, n_b, n_a)))

  n_cell <- n_e * n_b * n_a
  cell   <- seq_len(n_cell)
  e_of   <- rep(seq_len(n_e), each = n_b * n_a)   # income state of each cell

  ## 4 joint asset destinations per (cell, e') with product weights.
  dest <- list(
    list(bo = 0L, ao = 0L, w = pb * pa),
    list(bo = 1L, ao = 0L, w = (1 - pb) * pa),
    list(bo = 0L, ao = 1L, w = pb * (1 - pa)),
    list(bo = 1L, ao = 1L, w = (1 - pb) * (1 - pa))
  )

  max_nnz <- n_cell * n_e * 4L
  from <- integer(max_nnz); to <- integer(max_nnz); val <- numeric(max_nnz)
  pos <- 0L
  for (ep in seq_len(n_e)) {
    pr <- Pi[, ep][e_of]                 # Pi[e(cell), ep] for every cell
    if (all(pr == 0)) next
    base_to <- (ep - 1L) * n_b * n_a
    for (d in dest) {
      idx <- (pos + 1L):(pos + n_cell)
      from[idx] <- cell
      to[idx]   <- base_to + (ib + d$bo - 1L) * n_a + (ia + d$ao)
      val[idx]  <- pr * d$w
      pos <- pos + n_cell
    }
  }
  keep <- seq_len(pos)
  Matrix::sparseMatrix(i = from[keep], j = to[keep], x = val[keep],
                       dims = c(n_cell, n_cell))
}


## MATRIX-FREE TWO-ASSET FORWARD PUSH: t(Lambda2(b_pol, a_pol)) %*% D, without
## ever building Lambda2.  The two-asset counterpart of .hank_forward_push()
## (R/hank-distribution.R), which carries the full rationale.
##
## WHY. The two-asset fake-news sweep (.hank_curly_sweep2) builds TWO
## n_cell x n_cell sparse operators per date purely to matvec each one into
## D_ss.  Measured on the installed -O2 build, hank_forward_operator2() is
## 94-97% of the cost of the (build + matvec) pair it feeds across the whole
## size range tested (n_cell = 240 to 2250) -- the assembly, not the product,
## IS the two-asset sweep's forward cost.  The three-asset path has been
## matrix-free since it landed (.hank_forward_direction3 ->
## hank_forward_direction3_cpp); the one-asset path moved in 0.9.0.0025; this
## closes the two-asset gap.
##
## THE IDENTITY is exactly what hank_forward_operator2() writes into its
## triplets, contracted rather than materialized.  With
##   Lambda2[(e,b,a), (e',b',a')] = Pi[e,e'] * lot_b(b->b') * lot_a(a->a'),
## the product factorizes into (i) a PRODUCT-lottery scatter over the four
## joint (b', a') destinations WITHIN each income state, using the same two
## independent 1-D lotteries the sparse builder uses, then (ii) one n_e x n_e
## mixing matmul over income.  So it is an algebraic reassociation: the same
## terms summed in a different ORDER, agreeing to floating-point round-off
## (gated at 1e-14 relative in test-hank-forward-push.R).
##
## Boundary/degenerate cases need no handling here for the same reason they
## need none in the sparse builder: .hank_lottery() clamps i to [1, n-1] and p
## to [0, 1] on each axis, so an off-grid policy puts its whole mass on one
## bracketing node and the other node receives an exact zero.
##
## `D` may be a length-n_cell vector in the package's two-asset cell order or
## an n_e x n_b x n_a array; the return is always a vector in cell order.
##
## BACKEND (0.9.0.0039).  Profiled on the installed -O2 build (n_e = 3,
## n_b = n_a = 50, T_h = 50, 3 inputs x 3 outputs), this function was ~0.600 s
## of a 0.794 s two-asset Jacobian build -- 76% -- called 2*T_h*n_inputs = 300
## times at ~2.0 ms each, with everything else on that path already compiled.
## The cost is not the arithmetic but materializing `key` and `w` at 4*n_cell
## elements each for rowsum().  `backend = "cpp"` moves the LOTTERIES and the
## SCATTER into hank_forward_push2_scatter_cpp(), which accumulates straight
## into `acc`; the income-mixing crossprod(Pi, .) stays here in BLAS either way.
## The two backends are BIT-IDENTICAL (identical(), asserted in
## test-hank-forward-push.R), because the kernel repeats rowsum(reorder=FALSE)'s
## corner-major summation order exactly -- see src/hank_push2.cpp.  The R path
## remains the reference spec and also serves grids with < 2 nodes, which the
## kernel refuses (as .hank_lottery() effectively does).
.hank_forward_push2 <- function(b_pol, a_pol, b_grid, a_grid, Pi, D,
                                backend = getOption("dynhr.hank_backend",
                                                    "cpp")) {
  n_e <- dim(b_pol)[1L]; n_b <- dim(b_pol)[2L]; n_a <- dim(b_pol)[3L]
  n_ba <- n_b * n_a

  Dv <- if (is.array(D) && length(dim(D)) == 3L) .hank2_arr_to_vec(D) else
    as.numeric(D)

  if (identical(backend, "cpp") && n_b >= 2L && n_a >= 2L) {
    acc <- hank_forward_push2_scatter_cpp(as.numeric(b_pol), as.numeric(a_pol),
                                          as.numeric(b_grid),
                                          as.numeric(a_grid), Dv,
                                          n_e, n_b, n_a)
  } else {
    ## Per-axis lotteries, on the same unfoldings hank_forward_operator2() uses.
    lb <- .hank_lottery(matrix(aperm(b_pol, c(1L, 3L, 2L)), n_e * n_a, n_b),
                        b_grid)
    la <- .hank_lottery(matrix(a_pol, n_e * n_b, n_a), a_grid)
    ib <- .hank2_arr_to_vec(aperm(array(lb$i, c(n_e, n_a, n_b)),
                                  c(1L, 3L, 2L)))
    pb <- .hank2_arr_to_vec(aperm(array(lb$p, c(n_e, n_a, n_b)),
                                  c(1L, 3L, 2L)))
    ia <- .hank2_arr_to_vec(array(la$i, c(n_e, n_b, n_a)))
    pa <- .hank2_arr_to_vec(array(la$p, c(n_e, n_b, n_a)))

    ## Destination cell of the (lower b', lower a') corner, in cell order; the
    ## other three corners are +n_a (upper b'), +1 (upper a'), +n_a+1 (both).
    ## The income-block offset is the source state's own block: the scatter
    ## happens BEFORE the income mixing.
    off  <- rep((seq_len(n_e) - 1L) * n_ba, each = n_ba)
    base <- off + (ib - 1L) * n_a + ia
    key  <- c(base, base + n_a, base + 1L, base + n_a + 1L)
    w    <- c(Dv * pb * pa, Dv * (1 - pb) * pa,
              Dv * pb * (1 - pa), Dv * (1 - pb) * (1 - pa))

    ## Scatter-add: many source cells land on the same destination, so this
    ## cannot be an indexed assignment.
    acc <- numeric(n_e * n_ba)
    s <- rowsum(w, key, reorder = FALSE)
    acc[as.integer(rownames(s))] <- s[, 1L]
  }

  ## acc is indexed (e-1)*n_ba + j with the source income state e; mix.
  as.numeric(t(crossprod(Pi, matrix(acc, n_e, n_ba, byrow = TRUE))))
}


## Matrix-free counterpart of hank_forward_operator2d()'s
## t(P Lambda_A + (1-P) Lambda_N) %*% D.  Transposing the branch mixture turns
## the diagonal pre-multipliers into a REWEIGHTING OF THE SOURCE MASS:
##   t(diag(P) L_A + diag(1-P) L_N) D = t(L_A) (P * D) + t(L_N) ((1-P) * D),
## i.e. two ordinary pushes on split mass.  `sol` carries b_A/a_A/b_N/a_N/P as
## in hank_forward_operator2d().  `backend` is passed straight through to both
## branch pushes, so the discrete-choice path inherits the compiled scatter.
.hank_forward_push2d <- function(sol, b_grid, a_grid, Pi, D,
                                 backend = getOption("dynhr.hank_backend",
                                                     "cpp")) {
  Pv <- .hank2_arr_to_vec(sol$P)
  Dv <- if (is.array(D) && length(dim(D)) == 3L) .hank2_arr_to_vec(D) else
    as.numeric(D)
  .hank_forward_push2(sol$b_A, sol$a_A, b_grid, a_grid, Pi, Pv * Dv,
                      backend = backend) +
    .hank_forward_push2(sol$b_N, sol$a_N, b_grid, a_grid, Pi, (1 - Pv) * Dv,
                        backend = backend)
}


#' Sequence-space DISTRIBUTION Jacobian of the two-asset block via the
#' fake-news algorithm
#'
#' Two-asset counterpart of \code{\link{hank_het_dist_jacobian}} (one-asset)
#' and \code{\link{hank_het3_dist_jacobian}} (three-asset): extends
#' \code{\link{hank_het2_jacobian}}'s fake-news algorithm to expose the full
#' distributional response \eqn{J^D[t, s, ] = dD_t/dI_s} (the change in the
#' \code{(n_e*n_b*n_a)}-vector cross-sectional distribution at date \code{t}
#' induced by an anticipated shock to aggregate input \code{i} at date
#' \code{s}), instead of aggregating it into scalar outputs \code{B}/\code{A}/
#' \code{C}/\code{CHI}.
#'
#' Reuses the identical backward sweep (\code{\link{.hank_curly_sweep2}},
#' hence \code{curlyD}) as \code{\link{hank_het2_jacobian}}. Distributions
#' push forward under the transpose of the steady-state joint Young operator;
#' \code{.hank_forward_push2} already IS that transpose contracted
#' matrix-free (its documentation states the convention: distributions push
#' forward as \code{t(Lambda) \%*\% d}), so, unlike the three-asset sibling's
#' \code{.hank_forward_apply3(..., transpose = TRUE)}, no separate transpose
#' flag is needed here -- calling \code{.hank_forward_push2()} directly on the
#' steady-state policies IS the distribution-push direction, by construction
#' of that function. The distribution fake-news matrix cumulates by
#' repeatedly applying that push to \code{curlyD[, s]}, rather than by
#' dotting against an expectation vector (the aggregate-output equivalent of
#' that projection, used by \code{\link{hank_het2_jacobian}} via
#' \code{block$Lambda \%*\% E}, the untransposed direction).
#'
#' TIMING: exactly as for the one- and three-asset blocks, \code{D_t} is the
#' distribution ENTERING period \code{t} (a predetermined state), so
#' \code{D_1 = D_ss} always and row \code{t = 1} of \code{J^D} is identically
#' zero for every shock date \code{s}. \code{curlyD[, s]} is the response of
#' the policy USED in period \code{s} (\code{s = 1} is the direct current-
#' period shock; \code{s >= 2} anticipation terms propagate via the joint
#' \code{(Vb, Va)} derivative -- see \code{\link{.hank_curly_sweep2}}), which
#' the forward operator turns into a distribution change one calendar period
#' later, at \code{t = s + 1}. So the whole cumulation is the aggregate-
#' Jacobian recursion (\code{\link{hank_het2_jacobian}}'s diagonal
#' cumulation) shifted down by one row.
#'
#' Unlike the three-asset block's singleton-\code{f_grid} reduction (whose
#' \code{px} column is an exact zero, see \code{.hank3_px_is_inert}), the
#' two-asset block has no analogous EXACT-ZERO aggregate input. It DOES have
#' one that needs different sweep plumbing: \code{theta_coll} (the collateral
#' coordinate). The AGGREGATE Jacobian's shared backward sweep
#' (\code{\link{.hank_curly_sweep2}}) applies a \code{dD1} coordinate-rebasing
#' correction at \code{s = 2} that is tuned to cancel correctly only once
#' CONTRACTED against a steady output policy -- exactly what
#' \code{\link{hank_het2_jacobian}}'s aggregation does (its own
#' \code{theta_coll} column is validated to machine-ND precision). The raw,
#' uncontracted per-cell distribution response is measurably wrong at the
#' own-shock diagonal if it reuses that same corrected term (2026-07-31
#' finding). Fixed 2026-08-05 by having \code{.hank_curly_sweep2} additionally
#' expose \code{curlyD_raw}, an UNCORRECTED counterpart that omits the
#' contraction-only cancellation; this function uses \code{curlyD_raw} for
#' \code{"theta_coll"} (identical to \code{curlyD} for every other input, so
#' the fix is zero-cost and bit-identical elsewhere). See
#' \code{\link{.hank_curly_sweep2}}'s \code{@return} for the derivation and
#' \code{test-hank-theta-coll-jacobian.R} for the FN-vs-ND oracle.
#'
#' @inheritParams hank_het2_jacobian
#' @param inputs Character subset of \code{c("rb", "ra", "w", "Tr",
#'   "theta_coll")} plus the block's transition-probability inputs, or
#'   \code{NULL} (default) for \code{c("rb", "ra", "w")}.
#'
#' @return Named list \code{JD[[input]]}, each a 3-D array of dimension
#'   \code{T_h x T_h x (n_e*n_b*n_a)} with \code{JD[[i]][t, s, ] = dD_t/dI_s}.
#' @seealso \code{\link{hank_het2_dist_jacobian_nd}} (the numerical oracle
#'   this is validated against), \code{\link{hank_het_dist_jacobian}}
#'   (one-asset), \code{\link{hank_het3_dist_jacobian}} (three-asset),
#'   \code{\link{hank_het2_jacobian}} (the aggregate two-asset Jacobian
#'   sharing this function's backward sweep)
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 2)
#' blk <- hank_het2_block(hank_asset_grid(40, 8, 0), hank_asset_grid(60, 6, 0),
#'                        inc$Pi, inc$e, beta = 0.95, eis = 0.5,
#'                        rb = 0.005, ra = 0.02, w = 1, chi0 = 0.25,
#'                        chi1 = 6.5, chi2 = 2, n_k = 8L)
#' JD <- hank_het2_dist_jacobian(blk, T_h = 3, inputs = c("rb", "w"))
#' dim(JD$rb)
#' @export
hank_het2_dist_jacobian <- function(block, T_h,
                                    inputs = NULL,
                                    delta_in = 1e-5, delta_va = 1e-6,
                                    delta_d = 1e-6,
                                    backend = getOption("dynhr.hank_backend",
                                                        "cpp"),
                                    threads = NULL) {
  if (!inherits(block, "hank_het2_block"))
    stop("hank_het2_dist_jacobian: block must be hank_het2_block")
  if (!is.numeric(T_h) || length(T_h) != 1L || T_h < 1 || !is.finite(T_h))
    stop("hank_het2_dist_jacobian: T_h must be positive")
  backend <- match.arg(backend, c("R", "cpp"))
  threads <- hank_resolve_threads(threads)
  if (is.null(inputs)) inputs <- c("rb", "ra", "w")
  inputs <- .hank_het2_check_inputs(block, inputs)

  n_cell <- length(block$D)
  P <- function(x) .hank_forward_push2(block$b, block$a, block$b_grid,
                                       block$a_grid, block$Pi, x,
                                       backend = backend)

  JD <- setNames(lapply(inputs, function(i) array(0, c(T_h, T_h, n_cell))),
                inputs)
  if (T_h < 2L) {
    ## row t=1 (D_ss, fixed) is the only row -- all-zero EXCEPT theta_coll's
    ## JD[1, 1, ] = dD1 (see the theta_coll-only note below).
    if ("theta_coll" %in% inputs) {
      sweep <- .hank_curly_sweep2(block, T_h, "theta_coll", character(0),
                                  delta_in, delta_va, delta_d,
                                  backend = backend, threads = threads)
      JD[["theta_coll"]][1L, 1L, ] <- sweep$dD1
    }
    return(JD)
  }

  for (i in inputs) {
    ## --- Step 1: backward sweep -> curlyD[, s] (curlyY not needed here) ---
    sweep  <- .hank_curly_sweep2(block, T_h, i, character(0),
                                 delta_in, delta_va, delta_d,
                                 backend = backend, threads = threads)
    ## theta_coll (D1 lift, 2026-08-05): TWO roles need TWO different sweep
    ## terms at column s = 2, because the "- dD1" correction in sweep$curlyD
    ## cancels a double count that only appears once the s=2 column is ADDED
    ## to a prior row -- and row t = 2 (the base row, filled directly with no
    ## addition) is not that; only tt >= 3's cumulation
    ## (JD[tt,s,] = JD[tt-1,s-1,] + FD[tt][,s]) is.
    ##   - curlyD_direct: fills JD[2, s, ] directly (no addition) -> needs the
    ##     UNCORRECTED sweep$curlyD_raw, confirmed against ND to ~1e-6.
    ##   - curlyD_seed: seeds the P-push recursion that builds FD[[3]],
    ##     FD[[4]], ... for tt >= 3 -> needs the CORRECTED sweep$curlyD (same
    ##     one hank_het2_jacobian's aggregate cumulation uses), because
    ##     JD[tt,s,] for tt >= 3 ADDS JD[tt-1,s-1,] (which, at s=2, is
    ##     JD[2,1,] = curlyD[,1], already carrying its own dD1-derived
    ##     contribution) on top of the pushed s=2 term -- using the raw
    ##     (uncorrected) s=2 term there double-counts, exactly the bookkeeping
    ##     .hank_curly_sweep2's Roxygen documents for the aggregate case,
    ##     which turns out to recur here too. For every other input
    ##     curlyD_direct == curlyD_seed == sweep$curlyD bit-for-bit (no cost,
    ##     no behavior change).
    curlyD_direct <- if (i == "theta_coll") sweep$curlyD_raw else sweep$curlyD
    curlyD_seed   <- sweep$curlyD

    ## theta_coll only: row t = 1 is NOT the fixed D_ss zero row every other
    ## input has. An unanticipated shock at s = 1 re-expresses D_1 itself in
    ## the shifted gap coordinate x = b + theta_1*a, so JD[1, 1, ] = dD1 --
    ## the same raw date-1 coordinate-rebasing vector .hank_curly_sweep2
    ## folds (y-weighted) into curlyY[[o]][1] for the aggregate Jacobian.
    ## Every other (t, s) in row 1 stays zero: only the CONTEMPORANEOUS shock
    ## touches the entering distribution (confirmed against
    ## hank_het2_dist_jacobian_nd).
    if (i == "theta_coll") JD[[i]][1L, 1L, ] <- sweep$dD1

    ## --- Step 3: distribution fake-news F^D, indexed by CALENDAR date t ---
    ## FD[[2]][, s] = curlyD_seed[, s]  (SEED for propagation only; row t=2
    ##                                    of JD is filled from curlyD_direct
    ##                                    below, not from FD[[2]])
    ## FD[[t]][, s] = P(FD[[t-1]][, s])   for t >= 3
    ## (FD[[1]] would be t=1, always zero -- omitted; loop starts at t=2.)
    FD <- vector("list", T_h)
    FD[[2L]] <- curlyD_seed
    for (tt in seq_len(T_h - 2L) + 2L) {              # tt = 3 .. T_h, empty if T_h < 3
      prev <- FD[[tt - 1L]]
      cur  <- matrix(0, n_cell, T_h)
      for (s in seq_len(T_h)) cur[, s] <- P(prev[, s])
      FD[[tt]] <- cur
    }

    ## --- Step 4: diagonal cumulation, vector-valued per (t, s) ---
    ## Mirrors hank_het2_jacobian's aggregate assembly exactly, just shifted
    ## down one row (t=1 row is the fixed, unresponsive D_ss and stays
    ## all-zero; the recursion proper starts at t=2).
    JD[[i]][2L, 1L, ] <- curlyD_direct[, 1L]          # JD[2, 1, ] = curlyD_direct[, 1]
    for (s in seq_len(T_h - 1L) + 1L)                 # s = 2 .. T_h
      JD[[i]][2L, s, ] <- curlyD_direct[, s]          # JD[2,s,]=curlyD_direct[,s] (t-1=1 row is 0)
    for (tt in seq_len(T_h - 2L) + 2L) {               # tt = 3 .. T_h, empty if T_h < 3
      JD[[i]][tt, 1L, ] <- FD[[tt]][, 1L]             # JD[t, 1, ] = FD[t][, 1]
      for (s in seq_len(T_h - 1L) + 1L) {              # s = 2 .. T_h
        JD[[i]][tt, s, ] <- JD[[i]][tt - 1L, s - 1L, ] + FD[[tt]][, s]
      }
    }
  }
  JD
}


#' Brute-force numerical-differentiation DISTRIBUTION Jacobian of the
#' two-asset block
#'
#' Reference sequence-space distribution Jacobian \eqn{J^D[t,s,] = dD_t/dI_s}
#' for the two-asset block, computed exactly like
#' \code{\link{hank_het2_jacobian_nd}} but keeping the full \code{Dpath}
#' (already returned unconditionally by \code{\link{hank_td2_nonlinear}})
#' instead of aggregating it into \code{B}/\code{A}/\code{C}/\code{CHI}.
#' Used to validate \code{\link{hank_het2_dist_jacobian}}.
#'
#' @inheritParams hank_het2_jacobian_nd
#' @param inputs Character subset of \code{c("rb", "ra", "w", "Tr",
#'   "theta_coll")} plus the block's transition-probability inputs, or
#'   \code{NULL} (default) for \code{c("rb", "ra", "w")}. This brute-force
#'   oracle computes the true \code{dD} response for ANY admissible input by
#'   perturbing and re-solving, which is what validates
#'   \code{\link{hank_het2_dist_jacobian}}'s fake-news columns
#'   (\code{"theta_coll"} included, since the 2026-08-05 date-0 lift).
#'
#' @return Named list \code{JD_nd[[input]]}, each a 3-D array of dimension
#'   \code{T_h x T_h x (n_e*n_b*n_a)} with \code{JD_nd[[i]][t, s, ] =
#'   dD_t/dI_s} (central difference).
#' @seealso \code{\link{hank_het2_dist_jacobian}} (the fake-news distribution
#'   Jacobian this validates), \code{\link{hank_het_dist_jacobian_nd}}
#'   (one-asset), \code{\link{hank_het3_dist_jacobian_nd}} (three-asset)
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 2)
#' blk <- hank_het2_block(hank_asset_grid(40, 8, 0), hank_asset_grid(60, 6, 0),
#'                        inc$Pi, inc$e, beta = 0.95, eis = 0.5,
#'                        rb = 0.005, ra = 0.02, w = 1, chi0 = 0.25,
#'                        chi1 = 6.5, chi2 = 2, n_k = 8L)
#' JD_nd <- hank_het2_dist_jacobian_nd(blk, T_h = 3, inputs = c("rb", "w"),
#'                                    delta = 3e-6)
#' dim(JD_nd$rb)
#' @keywords internal
#' @export
hank_het2_dist_jacobian_nd <- function(block, T_h, inputs = NULL,
                                       delta = 1e-5,
                                       backend = getOption("dynhr.hank_backend",
                                                           "cpp"),
                                       threads = NULL) {
  if (!inherits(block, "hank_het2_block"))
    stop("hank_het2_dist_jacobian_nd: block must be hank_het2_block")
  if (!is.numeric(T_h) || length(T_h) != 1L || T_h < 1 || !is.finite(T_h))
    stop("hank_het2_dist_jacobian_nd: T_h must be positive")
  backend <- match.arg(backend, c("R", "cpp"))
  threads <- hank_resolve_threads(threads)
  if (is.null(inputs)) inputs <- c("rb", "ra", "w")
  inputs <- .hank_het2_check_inputs(block, inputs)

  n_cell <- length(block$D)
  base <- list(rb = rep(block$rb, T_h), ra = rep(block$ra, T_h),
              w = rep(block$w, T_h), Tr = rep(.hank_block_tr(block), T_h),
              theta_coll = rep(if (is.null(block$theta_coll)) 0
                               else block$theta_coll, T_h))

  run <- function(p, pi_paths) hank_td2_nonlinear(
    block, rb_path = p$rb, ra_path = p$ra, w_path = p$w, T_h = T_h,
    pi_input_paths = pi_paths, Tr_path = p$Tr, theta_path = p$theta_coll,
    backend = backend, threads = threads)

  JD_nd <- setNames(lapply(inputs, function(i) array(0, c(T_h, T_h, n_cell))),
                    inputs)
  for (i in inputs) {
    is_pi <- !(i %in% c("rb", "ra", "w", "Tr", "theta_coll"))
    for (s in seq_len(T_h)) {
      p <- m <- base; pip <- pim <- NULL
      if (is_pi) {
        x0 <- rep(block$Pi_inputs[[i]], T_h)
        xp <- x0; xp[s] <- xp[s] + delta
        xm <- x0; xm[s] <- xm[s] - delta
        pip <- setNames(list(xp), i); pim <- setNames(list(xm), i)
      } else {
        p[[i]][s] <- p[[i]][s] + delta
        m[[i]][s] <- m[[i]][s] - delta
      }
      op <- run(p, pip); om <- run(m, pim)
      dD <- (op$Dpath - om$Dpath) / (2 * delta)   # (n_cell x T_h), col t
      for (tt in seq_len(T_h)) JD_nd[[i]][tt, s, ] <- dD[, tt]
    }
  }
  JD_nd
}
