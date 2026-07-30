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
