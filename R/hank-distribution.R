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
  if (!is.matrix(a_pol) || !is.numeric(a_pol) || !all(is.finite(a_pol)))
    stop("hank_forward_operator(): 'a_pol' must be a finite numeric ",
         "n_e x n_a matrix.")
  n_e <- nrow(a_pol); n_a <- ncol(a_pol)
  if (!is.numeric(a_grid) || length(a_grid) != n_a ||
      !all(is.finite(a_grid)) || (n_a > 1L && any(diff(a_grid) <= 0)))
    stop("hank_forward_operator(): 'a_grid' must be a finite, strictly ",
         "increasing numeric vector of length ncol(a_pol) (", n_a, ").")
  .hank_check_markov(Pi, n_e, caller = "hank_forward_operator")
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
#' @param backend Character: \code{"cpp"} (default) or \code{"R"}. \code{"cpp"}
#'   runs the power iteration in a compiled kernel directly off \code{Lambda}'s
#'   CSC slots (\code{d_next = t(Lambda) \%*\% d}, without ever forming the
#'   transpose), for the same numerics as the R path. Defaults to
#'   \code{getOption("dynhr.hank_backend", "cpp")}; the R reference path
#'   remains available via \code{backend = "R"} or
#'   \code{options(dynhr.hank_backend = "R")}.
#'
#' @return A list with \code{d} (stationary distribution vector, sums to 1),
#'   \code{iterations}, and \code{converged}.
#' @export
hank_stationary_dist <- function(Lambda, d0 = NULL, tol = 1e-13,
                                 maxit = 200000L,
                                 backend = getOption("dynhr.hank_backend", "cpp")) {
  backend <- match.arg(backend, c("R", "cpp"))
  ## Input contract (adversarial review 2026-07-13, P1): an unvalidated d0
  ## previously reached the C++ kernel unchecked -- a zero-mass d0 became
  ## NaN/NaN after normalization and the kernel's NaN-swallowing max-diff
  ## test returned converged = TRUE on iteration 1; a short d0 indexed past
  ## the end of a std::vector (undefined behavior). Both backends now reject
  ## invalid inputs identically, before any normalization.
  n <- nrow(Lambda)
  if (is.null(n) || n != ncol(Lambda))
    stop("hank_stationary_dist(): 'Lambda' must be a square matrix.")
  lam_vals <- tryCatch(
    if (methods::is(Lambda, "sparseMatrix")) Lambda@x else
      as.numeric(as.matrix(Lambda)),
    error = function(e) NULL)
  if (!is.null(lam_vals) && !all(is.finite(lam_vals)))
    stop("hank_stationary_dist(): 'Lambda' has non-finite entries.")
  if (!(is.numeric(tol) && length(tol) == 1L && is.finite(tol) && tol > 0))
    stop("hank_stationary_dist(): 'tol' must be a finite positive scalar.")
  if (!(is.numeric(maxit) && length(maxit) == 1L && is.finite(maxit) &&
        maxit >= 1))
    stop("hank_stationary_dist(): 'maxit' must be a finite positive count.")
  if (!is.null(d0)) {
    if (!is.numeric(d0) || length(d0) != n)
      stop("hank_stationary_dist(): 'd0' must be a numeric vector of ",
           "length nrow(Lambda) (", n, "), got length ", length(d0), ".")
    if (!all(is.finite(d0)) || any(d0 < 0))
      stop("hank_stationary_dist(): 'd0' must be finite and non-negative.")
    if (sum(d0) <= 0)
      stop("hank_stationary_dist(): 'd0' must have strictly positive total ",
           "mass (sum = ", format(sum(d0)), ").")
  }
  d <- if (is.null(d0)) rep(1 / n, n) else d0 / sum(d0)

  if (backend == "cpp") {
    Lc  <- methods::as(Lambda, "CsparseMatrix")
    res <- hank_stationary_dist_lambda_cpp(Lc@p, Lc@i, Lc@x, n, d,
                                           tol, as.integer(maxit))
    return(list(d = res$d, iterations = as.integer(res$iterations),
                converged = as.logical(res$converged)))
  }

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
    md <- max(abs(d_new - d))
    d <- d_new
    if (is.na(md)) break                     # NaN iterate: fail loudly, fast
    if (md < tol) { converged <- TRUE; break }
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


## MATRIX-FREE FORWARD PUSH: t(Lambda(a_pol)) %*% D, without ever building
## Lambda.
##
## WHY. The fake-news backward sweep (.hank_curly_sweep) needs the
## distributional response to a perturbed savings policy, which it obtained by
## calling hank_forward_operator() twice per date and multiplying each sparse
## operator into D_ss. Measured (T_h = 200, n_a = 100, n_e = 3, installed -O2):
## hank_forward_operator 0.253 ms vs 0.015 ms for the product it feeds, so
## **86% of the sweep was assembling n_e*n_a x n_e*n_a sparse matrices for a
## matvec that never needs one** -- and the sweep is 94% of a het-block
## structural FD tap. This is the "lazy Lambda" follow-up flagged when the C++
## EGM backend landed (0.9.0.0004) and never taken.
##
## The identity is the same one hank_forward_operator() encodes, contracted
## rather than materialized: Lambda[(e,a), (e',a')] = Pi[e,e'] * lottery(a->a'),
## so
##   (t(Lambda) D)[(e',a')] = sum_e Pi[e,e'] * sum_a D[e,a] * lottery(a->a'),
## i.e. an asset-lottery scatter WITHIN each income state, then one n_e x n_e
## mixing matmul. Both factors are exactly what the triplet builder writes into
## `val` (pr * pvec and pr * (1 - pvec)), so this is an algebraic
## reassociation, not a different discretization: it agrees with the sparse
## path to floating-point round-off (asserted at 1e-15 in
## test-hank-forward-push.R), differing only in summation ORDER.
##
## Boundary/degenerate cases need no special handling here for the same reason
## they need none there: .hank_lottery() clamps i to [1, n_a-1] and p to
## [0, 1], so an off-grid policy puts its whole mass on one bracketing node and
## the other node receives an exact zero.
.hank_forward_push <- function(a_pol, a_grid, Pi, D) {
  n_e <- nrow(a_pol); n_a <- ncol(a_pol)
  Dm <- if (is.matrix(D)) D else .hank_vec_to_mat(D, n_e, n_a)
  lot <- .hank_lottery(a_pol, a_grid)
  off <- (seq_len(n_e) - 1L) * n_a                  # income-state block offset
  key <- c(off + lot$i, off + lot$i + 1L)           # recycles down the columns
  w   <- c(Dm * lot$p, Dm * (1 - lot$p))
  ## Scatter-add: keys repeat (many source cells land on the same node), so
  ## this cannot be an indexed assignment.
  acc <- numeric(n_e * n_a)
  s <- rowsum(w, key, reorder = FALSE)
  acc[as.integer(rownames(s))] <- s[, 1L]
  ## acc is indexed (e-1)*n_a + a', i.e. distribution order; mix over income.
  .hank_mat_to_vec(crossprod(Pi, .hank_vec_to_mat(acc, n_e, n_a)))
}


#' Aggregate a per-agent quantity over the distribution
#'
#' Both arguments are coerced to distribution (row-major) order before the
#' mass-weighted sum, so a matrix per-agent quantity (\code{n_e x n_a}) aligns
#' correctly with the distribution vector regardless of R's column-major default.
#'
#' @param d Distribution: a length-\code{n_e*n_a} vector (distribution order) or
#'   an \code{n_e x n_a} matrix.
#' @param x Per-agent quantity, same shape as \code{d}. Mixing a distribution
#'   VECTOR with a policy MATRIX (or vice versa) is allowed when the total
#'   cell counts agree -- the standard package idiom
#'   \code{hank_aggregate(block$D, block$a)} -- but incompatible sizes are an
#'   error: R's silent recycling previously let an exactly-dividing mismatch
#'   (e.g. \code{length(d) = 4}, \code{length(x) = 2}) return a plausible,
#'   numerically wrong scalar.
#' @return The scalar mass-weighted mean \code{sum(d * x)}.
#' @export
hank_aggregate <- function(d, x) {
  if (is.matrix(d) && is.matrix(x) && !identical(dim(d), dim(x)))
    stop("hank_aggregate(): 'd' (", nrow(d), " x ", ncol(d), ") and 'x' (",
         nrow(x), " x ", ncol(x), ") must have identical dimensions.")
  dv <- if (is.matrix(d)) .hank_mat_to_vec(d) else as.numeric(d)
  xv <- if (is.matrix(x)) .hank_mat_to_vec(x) else as.numeric(x)
  if (length(dv) != length(xv))
    stop("hank_aggregate(): 'd' (", length(dv), " cells) and 'x' (",
         length(xv), " cells) must cover the same number of cells -- ",
         "refusing to recycle.")
  sum(dv * xv)
}


#' Is a Markov transition matrix reducible (more than one communicating class)?
#'
#' At \code{hank_employment_income3}'s degenerate nesting point
#' (\code{p_un = p_nu = 0}, default \code{f_ne = s_en = 0}) the E/U/N
#' employment chain has \code{N} as a CLOSED, unreachable class: the chain is
#' reducible, has more than one invariant distribution, and a uniform-seeded
#' power iteration (\code{\link{hank_stationary_dist}}'s default) strands
#' whatever mass the init put in \code{N} there forever -- see
#' \code{\link{hank_het_block}}'s \code{dist_init} argument and its
#' reducibility guard.
#'
#' Detected via the transitive closure of the "positive one-step transition
#' probability" digraph: \code{Pi} is irreducible iff every state can reach
#' every other state (the digraph, with self-loops added, is strongly
#' connected). \code{n} is always tiny here (the number of income states), so
#' a dense boolean transitive closure by repeated squaring
#' (\code{O(n^3 log n)}, exact -- no eigenvalue tolerance games) is cheap.
#' Entries are compared to exact \code{0}, not a tolerance, so a chain with a
#' tiny-but-nonzero rate (e.g. \code{p_un = p_nu = 1e-9}) is correctly
#' IRREDUCIBLE (no false positive) -- the discontinuity at exactly zero is
#' the point (see the paper bug report this guards against).
#'
#' @param Pi Numeric \code{n x n} transition matrix (rows sum to 1;
#'   non-negative entries assumed -- callers validate that separately via
#'   \code{\link{.hank_check_markov}}).
#' @return \code{TRUE} if \code{Pi} is reducible (more than one communicating
#'   class), \code{FALSE} if it is irreducible (a single communicating class
#'   covering every state).
#' @keywords internal
.hank_pi_reducible <- function(Pi) {
  n <- nrow(Pi)
  if (n <= 1L) return(FALSE)
  A <- Pi > 0
  diag(A) <- TRUE
  ## Reachability-in-<=k-steps via repeated boolean squaring: after k
  ## squarings, R encodes reachability within 2^k steps of the self-looped
  ## graph, so k = ceiling(log2(n)) steps (2^k >= n) certainly covers the
  ## longest possible shortest path (<= n - 1 edges).
  R <- A
  k <- max(1L, ceiling(log2(n)))
  for (i in seq_len(k)) R <- (R %*% R) > 0
  ## Irreducible iff EVERY state reaches every other state, i.e. R is the
  ## all-TRUE matrix (a single communicating class spanning all n states).
  !all(R)
}


#' Resolve a \code{dist_init} argument to a full initial distribution vector
#'
#' Shared by \code{\link{hank_het_block}} (and any future het-block
#' constructor needing the same knob): \code{dist_init} may be
#' \describe{
#'   \item{\code{NULL}}{No seed requested; returns \code{NULL} (callers keep
#'     their own default, e.g. \code{hank_stationary_dist}'s uniform init).}
#'   \item{a length-\code{n_e} vector}{A MARGINAL over income states, spread
#'     UNIFORMLY across the asset grid -- \code{d0[e, a] = dist_init[e] /
#'     n_a} for every \code{a}. This is the shape
#'     \code{hank_employment_income3()$pi_m}-like objects (or, for the
#'     reducibility fix, the FULL combined-state stationary vector
#'     \code{inc$pi}, which is length \code{n_e} in the block's sense once
#'     \code{e} is the combined employment x productivity index) naturally
#'     have.}
#'   \item{an \code{n_e x n_a} matrix, or a length-\code{n_e*n_a} vector}{A
#'     full initial distribution in the block's own \code{(e, a)} shape
#'     (distribution order for the vector form -- see the file header).}
#' }
#' Validated finite, non-negative, and strictly positive total mass (the same
#' contract \code{\link{hank_stationary_dist}} enforces on its own \code{d0},
#' checked here too so a bad \code{dist_init} fails at the block boundary
#' with a block-shaped message rather than deep inside the solver).
#'
#' @param dist_init The raw \code{dist_init} argument (\code{NULL} or
#'   numeric vector/matrix).
#' @param n_e,n_a Grid dimensions.
#' @param caller Function name for error messages.
#' @return \code{NULL}, or a length-\code{n_e*n_a} numeric vector in
#'   distribution order (unnormalized; \code{hank_stationary_dist} normalizes
#'   its \code{d0}).
#' @keywords internal
.hank_dist_init_d0 <- function(dist_init, n_e, n_a, caller) {
  if (is.null(dist_init)) return(NULL)
  if (!is.numeric(dist_init) || !all(is.finite(dist_init)))
    stop(caller, ": 'dist_init' must be a finite numeric vector or matrix.")
  if (any(dist_init < 0))
    stop(caller, ": 'dist_init' must be non-negative (min = ",
         format(min(dist_init)), ").")
  if (sum(dist_init) <= 0)
    stop(caller, ": 'dist_init' must have strictly positive total mass ",
         "(sum = ", format(sum(dist_init)), ").")
  if (is.matrix(dist_init)) {
    if (nrow(dist_init) != n_e || ncol(dist_init) != n_a)
      stop(caller, ": 'dist_init' matrix must be ", n_e, " x ", n_a,
           " (n_e x n_a), got ", nrow(dist_init), " x ", ncol(dist_init), ".")
    return(.hank_mat_to_vec(dist_init))
  }
  if (length(dist_init) == n_e * n_a) return(as.numeric(dist_init))
  if (length(dist_init) == n_e)
    return(.hank_mat_to_vec(matrix(dist_init / n_a, n_e, n_a)))
  stop(caller, ": 'dist_init' must be a length-", n_e,
       " marginal over income states, a length-", n_e * n_a,
       " full distribution vector, or an ", n_e, " x ", n_a,
       " matrix; got length ", length(dist_init), ".")
}


#' Validate a Markov transition matrix at a public HANK boundary
#'
#' Shared input contract (adversarial review 2026-07-13, P2): \code{Pi} must
#' be a square numeric \code{n x n} matrix, finite, entrywise non-negative
#' (within \code{tol}), with every row summing to 1 (within \code{tol}).
#' Inner iteration loops stay validation-free; public constructors/solvers
#' call this once at their entry boundary.
#'
#' @param Pi Candidate transition matrix.
#' @param n Required dimension.
#' @param caller Function name for the error message.
#' @param tol Numerical tolerance for non-negativity and row sums.
#' @keywords internal
.hank_check_markov <- function(Pi, n, caller, tol = 1e-8) {
  if (!is.matrix(Pi) || !is.numeric(Pi) || nrow(Pi) != n || ncol(Pi) != n)
    stop(caller, "(): 'Pi' must be a numeric ", n, " x ", n,
         " matrix (got ",
         if (is.matrix(Pi)) paste0(nrow(Pi), " x ", ncol(Pi)) else
           paste0("a ", class(Pi)[1L]), ").")
  if (!all(is.finite(Pi)))
    stop(caller, "(): 'Pi' has non-finite entries.")
  if (any(Pi < -tol))
    stop(caller, "(): 'Pi' has negative entries (min = ",
         format(min(Pi)), ") -- not a transition matrix.")
  rs <- rowSums(Pi)
  if (any(abs(rs - 1) > tol))
    stop(caller, "(): 'Pi' rows must sum to 1 (max |rowsum - 1| = ",
         format(max(abs(rs - 1))), ") -- not row-stochastic.")
  invisible(TRUE)
}
