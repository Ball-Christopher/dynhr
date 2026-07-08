## R/hank-distribution.R
## --------------------------------------------------------------------------
## Young's (2010) non-stochastic simulation for the cross-sectional
## distribution of a one-asset heterogeneous-agent household block.
##
## The distribution is a histogram (mass vector) on the fixed (e, a) grid.  A
## household with savings policy a'(a, e) that lands between two fixed
## gridpoints a_k <= a' < a_{k+1} has its mass split ("lottery") in proportion
##   p     to a_k       where p = (a_{k+1} - a') / (a_{k+1} - a_k)
##   1 - p to a_{k+1}
## exactly conserving total mass.  Combined with the income transition Pi this
## defines a linear, sparse forward operator Lambda on the (e, a) grid such that
##   d_{t+1} = Lambda' d_t.
##
## STATE ORDERING: the flattened distribution vector d has length n_e * n_a with
## index(e, a) = (e - 1) * n_a + a  (income state OUTER/slow, asset INNER/fast).
## A distribution supplied/returned as a matrix is n_e x n_a with this row-major
## flattening.
##
## The sparse Lambda built here is reused by the fake-news sequence-space
## Jacobian (a later increment): expectation vectors iterate Lambda and the
## distributional shock responses push mass through Lambda'.
## --------------------------------------------------------------------------


#' Lottery weights for Young's non-stochastic simulation
#'
#' For each \code{(e, a)} gridpoint, finds the lower asset gridpoint index the
#' savings policy lands on and the mass fraction assigned to it.
#'
#' @param a_pol Numeric \code{n_e x n_a}: savings policy \eqn{a'(a, e)}.
#' @param a_grid Numeric length-\code{n_a}: fixed asset grid (increasing).
#'
#' @return A list with \code{i} (\code{n_e x n_a} integer: lower gridpoint index
#'   in \code{1..n_a-1}) and \code{p} (\code{n_e x n_a}: mass fraction on the
#'   lower point, in \code{[0, 1]}).
#' @keywords internal
.hank_lottery <- function(a_pol, a_grid) {
  n_a <- length(a_grid)
  ## Lower bracketing index, clamped so i and i+1 are both valid.
  i <- findInterval(a_pol, a_grid)
  i[i < 1L]        <- 1L
  i[i > n_a - 1L]  <- n_a - 1L
  dim(i) <- dim(a_pol)
  aL <- matrix(a_grid[i],      nrow(a_pol), ncol(a_pol))
  aU <- matrix(a_grid[i + 1L], nrow(a_pol), ncol(a_pol))
  p  <- (aU - a_pol) / (aU - aL)
  ## Clamp for policies outside the grid (constraint / top).
  p[p < 0] <- 0; p[p > 1] <- 1
  list(i = i, p = p)
}


#' Build the sparse Young's-method forward transition operator Lambda
#'
#' Constructs the \code{(n_e*n_a) x (n_e*n_a)} row-stochastic transition matrix
#' \code{Lambda} where \code{Lambda[from, to]} is the probability of moving from
#' state \code{from = (e,a)} to the next state \code{(e_next, a_next)} in one period: the
#' savings lottery on assets composed with the income transition \code{Pi}.  The
#' distribution updates as \code{d_next = t(Lambda) \%*\% d}.
#'
#' @param a_pol Numeric \code{n_e x n_a}: savings policy.
#' @param a_grid Numeric length-\code{n_a}: fixed asset grid.
#' @param Pi Numeric \code{n_e x n_e}: income transition matrix.
#'
#' @return A sparse \code{dgCMatrix} (\code{Matrix} package), row-stochastic.
#' @export
hank_forward_operator <- function(a_pol, a_grid, Pi) {
  n_e <- nrow(a_pol); n_a <- ncol(a_pol)
  lot <- .hank_lottery(a_pol, a_grid)

  ## Build (from, to, value) triplets.  For each (e,a) and each e' with
  ## Pi[e,e'] > 0, two asset destinations (lower/upper), unless a boundary
  ## degenerates them.  Preallocate generously.
  max_nnz <- n_e * n_a * n_e * 2L
  from <- integer(max_nnz); to <- integer(max_nnz); val <- numeric(max_nnz)
  pos  <- 0L

  ## from index for all (e,a): (e-1)*n_a + a
  for (e in seq_len(n_e)) {
    base_from <- (e - 1L) * n_a
    ivec <- lot$i[e, ]; pvec <- lot$p[e, ]
    for (ep in seq_len(n_e)) {
      pr <- Pi[e, ep]
      if (pr == 0) next
      base_to <- (ep - 1L) * n_a
      ## lower destinations
      idx <- (pos + 1L):(pos + n_a)
      from[idx] <- base_from + seq_len(n_a)
      to[idx]   <- base_to + ivec
      val[idx]  <- pr * pvec
      pos <- pos + n_a
      ## upper destinations
      idx <- (pos + 1L):(pos + n_a)
      from[idx] <- base_from + seq_len(n_a)
      to[idx]   <- base_to + ivec + 1L
      val[idx]  <- pr * (1 - pvec)
      pos <- pos + n_a
    }
  }
  Matrix::sparseMatrix(i = from[seq_len(pos)], j = to[seq_len(pos)],
                       x = val[seq_len(pos)],
                       dims = c(n_e * n_a, n_e * n_a))
}


#' Stationary cross-sectional distribution
#'
#' Iterates \code{d_{t+1} = t(Lambda) \%*\% d} to the invariant distribution.
#'
#' @param Lambda Sparse forward operator from \code{\link{hank_forward_operator}}.
#' @param d0 Optional initial distribution (length \code{n_e*n_a}); default
#'   uniform.
#' @param tol Convergence tolerance on \code{max|d change|}.
#' @param maxit Maximum iterations.
#'
#' @return A list with \code{d} (stationary distribution vector, sums to 1),
#'   \code{iterations}, and \code{converged}.
#' @export
hank_stationary_dist <- function(Lambda, d0 = NULL, tol = 1e-13,
                                 maxit = 200000L) {
  n <- nrow(Lambda)
  d <- if (is.null(d0)) rep(1 / n, n) else d0 / sum(d0)
  Lt <- Matrix::t(Lambda)
  converged <- FALSE
  it <- 0L

  ## The power iteration calls the mat-vec `Lt %*% d` up to `maxit` times
  ## (typically ~hundreds). `Lt` is an S4 (Matrix package) sparse matrix, so
  ## `Lt %*% d` re-resolves the `%*%` generic (getClass/getClassDef/
  ## .getClassesFromCache) on EVERY iteration -- profiled at ~9% of solve
  ## time on the HANK GE hot path. `.hank_matvec_method()` resolves the
  ## applicable compiled method ONCE (memoized across calls, keyed by
  ## `class(Lt)`, which is always the same concrete class -- "dgCMatrix" --
  ## for every `Lambda` this package builds) and the loop calls that
  ## resolved closure directly, skipping `standardGeneric` dispatch on every
  ## iteration. It is the EXACT SAME compiled routine `Lt %*% d` would use,
  ## so this is a pure dispatch-overhead removal: the iterate sequence is
  ## bit-identical to the un-cached S4 path (verified in
  ## test-hank-hotpath-desr4.R), not merely numerically close. Attempts at
  ## a dense conversion or a hand-rolled sparse mat-vec were measured to be
  ## SLOWER (dense: O(n^2) FLOPs vs. the very sparse n_e*n_a structure;
  ## hand-rolled R accumulation: R-loop/vector overhead exceeds Matrix's
  ## compiled routine) -- see the session notes for benchmarks. If method
  ## resolution ever fails (e.g. a future Matrix version reshapes the
  ## generic), `.hank_matvec_method()` returns `NULL` and we fall back to
  ## the original `%*%` generic, so correctness never depends on the cache.
  mm <- .hank_matvec_method(Lt)
  matvec <- if (is.null(mm)) `%*%` else mm

  for (it in seq_len(maxit)) {
    d_new <- as.numeric(matvec(Lt, d))
    if (max(abs(d_new - d)) < tol) { d <- d_new; converged <- TRUE; break }
    d <- d_new
  }
  d <- d / sum(d)
  list(d = d, iterations = it, converged = converged)
}

## Package-private cache: resolved `%*%` S4 method, keyed by the concrete
## class of the (sparse) left operand. Lazily populated on first use and
## reused for the lifetime of the session -- every `Lambda` this package
## builds (via `hank_forward_operator`) is a `dgCMatrix`, so in practice this
## resolves once per session and every subsequent HANK GE solve reuses it.
.hank_dispatch_cache <- new.env(parent = emptyenv())

#' Resolve (and cache) the compiled \code{\%*\%} method for \code{class(Lt) x numeric}
#'
#' @param Lt A (sparse) \code{Matrix} object; only its class is used as the
#'   cache key.
#' @return The resolved \code{MethodDefinition} closure, or \code{NULL} if
#'   resolution fails for any reason (caller should fall back to the plain
#'   \code{\%*\%} generic in that case).
#' @keywords internal
.hank_matvec_method <- function(Lt) {
  key <- class(Lt)[1L]
  mm <- .hank_dispatch_cache[[key]]
  if (is.null(mm)) {
    mm <- tryCatch(
      methods::selectMethod("%*%", methods::signature(x = key, y = "numeric")),
      error = function(e) NULL)
    if (!is.null(mm)) .hank_dispatch_cache[[key]] <- mm
  }
  mm
}


#' Flatten an \code{n_e x n_a} matrix to distribution (row-major) order
#'
#' The distribution vector uses index \code{(e-1)*n_a + a}; base R's
#' \code{as.numeric()} on a matrix is column-major, so per-agent policy matrices
#' must be flattened row-major to align with the distribution.
#'
#' @param X An \code{n_e x n_a} matrix.
#' @return A length-\code{n_e*n_a} vector in distribution order.
#' @keywords internal
.hank_mat_to_vec <- function(X) as.numeric(t(X))

#' Reshape a distribution-order vector back to an \code{n_e x n_a} matrix
#' @param v A length-\code{n_e*n_a} vector in distribution order.
#' @param n_e,n_a Grid dimensions.
#' @return An \code{n_e x n_a} matrix.
#' @keywords internal
.hank_vec_to_mat <- function(v, n_e, n_a) matrix(v, n_e, n_a, byrow = TRUE)


#' Aggregate a per-agent quantity over the distribution
#'
#' Both arguments are coerced to distribution (row-major) order before the
#' mass-weighted sum, so a matrix per-agent quantity (\code{n_e x n_a}) aligns
#' correctly with the distribution vector regardless of R's column-major default.
#'
#' @param d Distribution: a length-\code{n_e*n_a} vector (distribution order) or
#'   an \code{n_e x n_a} matrix.
#' @param x Per-agent quantity, same shape as \code{d}.
#' @return The scalar mass-weighted mean \code{sum(d * x)}.
#' @export
hank_aggregate <- function(d, x) {
  dv <- if (is.matrix(d)) .hank_mat_to_vec(d) else as.numeric(d)
  xv <- if (is.matrix(x)) .hank_mat_to_vec(x) else as.numeric(x)
  sum(dv * xv)
}
