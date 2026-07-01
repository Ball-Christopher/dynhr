## R/global-cheb-utils.R
## --------------------------------------------------------------------------
## Chebyshev polynomial utilities for global/projection solution methods.
##
## Functions:
##   cheb_nodes(n)            -- Chebyshev-Gauss nodes on [-1,1]
##   cheb_basis(x, p, d)      -- complete Chebyshev basis matrix
##   cheb_normalize(x, lo, hi) -- rescale [lo,hi] -> [-1,1]
##   cheb_denormalize(u, lo, hi) -- rescale [-1,1] -> [lo,hi]
##   gauss_hermite(n)         -- Gauss-Hermite quadrature for N(0,1)
## --------------------------------------------------------------------------

## Chebyshev-Gauss (NOT Gauss-Lobatto) nodes on [-1,1].
## n nodes: cos(pi*(2*j-1)/(2*n)) for j=1..n, which are the roots of T_n(x).
## These are the standard interior nodes used for Chebyshev collocation.
##
## @param n  number of nodes (integer >= 1)
## @return   numeric vector of length n in descending order (from ~1 to ~-1)
cheb_nodes <- function(n) {
  j <- seq_len(n)
  cos(pi * (2 * j - 1) / (2 * n))
}

## Evaluate Chebyshev polynomials T_0, T_1, ..., T_p at each point in x.
##
## Uses the three-term recurrence: T_0=1, T_1=x, T_k=2x*T_{k-1}-T_{k-2}.
##
## @param x  numeric vector of points in [-1,1]
## @param p  maximum degree (integer >= 0)
## @return   matrix of dim length(x) x (p+1)
cheb_poly_1d <- function(x, p) {
  n <- length(x)
  T_mat <- matrix(0.0, nrow = n, ncol = p + 1L)
  T_mat[, 1L] <- 1.0          # T_0
  if (p >= 1L) T_mat[, 2L] <- x   # T_1
  for (k in seq_len(p - 1L) + 1L) {   # k = 2..p (1-indexed col k+1)
    T_mat[, k + 1L] <- 2.0 * x * T_mat[, k] - T_mat[, k - 1L]
  }
  T_mat
}

## Complete Chebyshev polynomial basis for d-dimensional state vectors.
##
## Evaluates all multi-variate monomials T_{i1}(x1)*T_{i2}(x2)*...*T_{id}(xd)
## with i1+i2+...+id <= p.  The ordering is lexicographic over the index
## tuples in non-decreasing total-degree order (i.e., constant first).
##
## @param x_mat  n x d matrix; each row is a state point in [-1,1]^d
## @param p      maximum total degree (integer >= 0)
## @return       n x n_basis matrix (n_basis = C(p+d, d))
cheb_basis <- function(x_mat, p) {
  if (is.vector(x_mat)) x_mat <- matrix(x_mat, ncol = 1L)
  n <- nrow(x_mat)
  d <- ncol(x_mat)

  ## Pre-evaluate 1-D Chebyshev polys for each dimension
  T_list <- vector("list", d)
  for (j in seq_len(d)) {
    T_list[[j]] <- cheb_poly_1d(x_mat[, j], p)  # n x (p+1) matrix
  }

  ## Enumerate all multi-index tuples (i1,...,id) with sum <= p
  ## Build via recursive enumeration
  idx_list <- .cheb_multiindex(d, p)  # list of integer vectors of length d

  n_basis <- length(idx_list)
  Phi <- matrix(1.0, nrow = n, ncol = n_basis)

  for (b in seq_len(n_basis)) {
    idx <- idx_list[[b]]   # degrees per dimension (0-indexed)
    for (j in seq_len(d)) {
      Phi[, b] <- Phi[, b] * T_list[[j]][, idx[j] + 1L]
    }
  }
  Phi
}

## Generate all d-dimensional multi-index tuples with total degree <= p.
## Returns a list of integer vectors (0-indexed degrees).
.cheb_multiindex <- function(d, p) {
  if (d == 1L) {
    return(lapply(0L:p, function(k) k))
  }
  result <- list()
  ## Enumerate by total degree 0..p
  for (tot in 0L:p) {
    ## partitions of tot into d non-negative parts
    parts <- .integer_partitions_d(tot, d)
    result <- c(result, parts)
  }
  result
}

## Generate all d-tuples of non-negative integers summing to exactly s.
## Returns a list of integer vectors.
.integer_partitions_d <- function(s, d) {
  if (d == 1L) return(list(s))
  result <- list()
  for (k in 0L:s) {
    sub <- .integer_partitions_d(s - k, d - 1L)
    for (part in sub) {
      result <- c(result, list(c(k, part)))
    }
  }
  result
}

## Rescale x from [lo, hi] to [-1, 1].
## @param x   numeric (scalar or vector)
## @param lo  lower bound of original domain
## @param hi  upper bound of original domain
## @return    numeric in [-1, 1]
cheb_normalize <- function(x, lo, hi) {
  2.0 * (x - lo) / (hi - lo) - 1.0
}

## Rescale u from [-1, 1] back to [lo, hi].
## @param u   numeric in [-1, 1]
## @param lo  lower bound of target domain
## @param hi  upper bound of target domain
## @return    numeric in [lo, hi]
cheb_denormalize <- function(u, lo, hi) {
  lo + (u + 1.0) * (hi - lo) / 2.0
}

## Gauss-Hermite quadrature nodes and weights for the standard Normal N(0,1).
##
## The nodes zeta_k and weights w_k satisfy:
##   E[f(X)] = integral f(x) phi(x) dx  ≈  sum_k w_k * f(zeta_k)
## where phi is the N(0,1) density.
##
## This uses the "probabilist" scaling so that sum(w_k) = 1 and the
## approximation is exact for polynomials of degree <= 2*n - 1.
##
## Hand-coded rules for n = 1..7 (sufficient for practical use).
## For n > 7, falls back to an eigenvalue-based computation.
##
## @param n  number of quadrature nodes (integer in 1..20)
## @return   list(nodes = numeric(n), weights = numeric(n))
gauss_hermite <- function(n) {
  ## Pre-computed rules (nodes in increasing order, weights summing to 1)
  ## Source: Abramowitz & Stegun Table 25.10 (physicists) converted to
  ## probabilists' scaling: zeta = sqrt(2)*xi, w_prob = w_phys/sqrt(pi)
  tables <- list(
    `1` = list(
      nodes   = 0.0,
      weights = 1.0
    ),
    `3` = list(
      nodes   = c(-1.7320508075688772, 0.0, 1.7320508075688772),
      weights = c(1/6, 2/3, 1/6)
    ),
    `5` = list(
      nodes   = c(-2.8569700138728056, -1.3556261799742674, 0.0,
                   1.3556261799742674,  2.8569700138728056),
      weights = c(0.011257411327720691, 0.22207592200561264, 0.5333333333333333,
                  0.22207592200561264,  0.011257411327720691)
    ),
    `7` = list(
      nodes   = c(-3.7504397177257425, -2.366759410734541,  -1.1544053548537484,
                   0.0,
                   1.1544053548537484,  2.366759410734541,   3.7504397177257425),
      weights = c(0.0009717812450995192, 0.054515582819664546, 0.4256072526101278,
                  0.8102646175568073,
                  0.4256072526101278,   0.054515582819664546, 0.0009717812450995192)
    )
  )
  ## Re-normalize weights so they sum to 1
  if (as.character(n) %in% names(tables)) {
    tbl <- tables[[as.character(n)]]
    w <- tbl$weights / sum(tbl$weights)
    return(list(nodes = tbl$nodes, weights = w))
  }

  ## Generic eigenvalue-based computation (Golub-Welsch algorithm)
  ## Jacobi matrix for Hermite: tridiagonal with beta_k = sqrt(k) on off-diag
  ## Physicists' polynomials: H_k, beta_k = sqrt(k/2)
  ## Probabilists' scaling: xi = sqrt(2) * x_phys
  betas <- sqrt(seq_len(n - 1L) / 2.0)
  J <- diag(n) * 0.0
  if (n > 1L) {
    J[cbind(seq_len(n - 1L), seq_len(n - 1L) + 1L)] <- betas
    J[cbind(seq_len(n - 1L) + 1L, seq_len(n - 1L))] <- betas
  }
  eig <- eigen(J, symmetric = TRUE)
  nodes_phys <- eig$values
  ## Sort in increasing order
  ord <- order(nodes_phys)
  nodes_phys <- nodes_phys[ord]
  ## Probabilists' nodes: xi_k = sqrt(2) * zeta_k
  nodes <- sqrt(2.0) * nodes_phys
  ## Weights: (first component of eigenvector)^2 * sqrt(pi) normalised to sum=1
  w_raw <- eig$vectors[1L, ord]^2 * sqrt(pi)
  weights <- w_raw / sum(w_raw)

  list(nodes = nodes, weights = weights)
}
