## R/dynamic-perturbation-quadrature.R
## --------------------------------------------------------------------------
## Quadrature rules for stochastic dynamic perturbation (Mennuni et al. 2025).
##
## Provides Gauss-Hermite product-rule quadrature (and a Stroud monomial
## shortcut for the 1-D case) to evaluate E_t[f(u)], u ~ N(0, Sigma_e).
##
## The core function .dp_quadrature_nodes() returns a list:
##   nodes   -- (n_nodes x n_shock) matrix of shock vectors u_k
##   weights -- length-n_nodes vector of positive weights summing to 1
##
## Under the returned nodes/weights, for any function g(u):
##   E[g(u)] ~= sum_k weights[k] * g(nodes[k,])
##
## Reference moments (exact for n_gh >= ceil((max_moment+1)/2) points):
##   E[1]         = 1
##   E[u_i]       = 0
##   E[u_i u_j]   = Sigma_e[i,j]
##   E[u_i^4]     = 3 * Sigma_e[i,i]^2  (for each i)
##   E[u_i^2 u_j^2] = Sigma_e[i,i]*Sigma_e[j,j] + 2*Sigma_e[i,j]^2  (i != j)
## --------------------------------------------------------------------------

#' Gauss-Hermite quadrature nodes and weights for u ~ N(0, Sigma_e)
#'
#' Returns quadrature nodes and weights for computing E[g(u)] where
#' u ~ N(0, Sigma_e).  Uses a product Gauss-Hermite rule in the rotated
#' coordinates of Sigma_e's Cholesky factor; the rule integrates polynomials
#' of degree up to 2*n_gh-1 exactly.
#'
#' @param Sigma_e  Symmetric positive-semidefinite shock covariance matrix
#'   (n_shock x n_shock).  Must be the full Sigma_e (not a Cholesky factor).
#'   Pass \code{matrix(sigma^2, 1, 1)} for a scalar shock with std sigma.
#' @param n_gh     Number of Gauss-Hermite points per dimension (default 5).
#'   The product rule uses \code{n_gh^n_shock} total nodes.  For n_shock > 3,
#'   consider using \code{monomial = TRUE} instead to limit node count.
#' @param monomial Logical; if \code{TRUE}, use the Stroud monomial-2 rule
#'   (2*n_shock nodes, exact for monomials up to degree 3: E[1], E[u_i],
#'   E[u_i u_j], E[u_i^3]).  For n_shock >= 4 this is much cheaper than the
#'   product rule while still capturing the risk-adjustment terms to O(sigma^2).
#'   Default \code{FALSE}.
#' @return A list with:
#'   \describe{
#'     \item{\code{nodes}}{(n_nodes x n_shock) numeric matrix of shock vectors.}
#'     \item{\code{weights}}{Length-n_nodes positive weights summing to 1.}
#'     \item{\code{n_shock}}{Number of shocks.}
#'     \item{\code{n_nodes}}{Number of quadrature nodes.}
#'     \item{\code{Sigma_e}}{The Sigma_e matrix used.}
#'   }
#' @noRd
.dp_quadrature_nodes <- function(Sigma_e, n_gh = 5L, monomial = FALSE) {
  n_shock <- nrow(Sigma_e)
  stopifnot(ncol(Sigma_e) == n_shock, n_shock >= 1L)

  if (monomial) {
    return(.dp_monomial2_nodes(Sigma_e))
  }

  ## Gauss-Hermite nodes/weights on the standard normal (mean 0, var 1/2)
  ## using the Golub-Welsch algorithm (Golub & Welsch 1969).
  ## We want E[g(z)] = int g(z) phi(z) dz where phi is N(0,1).
  ## GH quadrature for the weight w(z) = exp(-z^2) integrates
  ## g~(z) = g(z*sqrt(2)) * (1/sqrt(pi)).
  ## Rescaling: standard GH gives nodes t_k and weights w_k such that
  ##   int f(t) exp(-t^2) dt ~= sum_k w_k f(t_k).
  ## For E_{Z~N(0,1)}[g(Z)] we use x_k = sqrt(2)*t_k, w_k/sqrt(pi).
  gh <- .gh_1d(n_gh)
  x1d <- gh$x   # nodes on N(0,1)
  w1d <- gh$w   # weights summing to 1

  ## Cholesky decomposition of Sigma_e for coordinate transform
  ## u = L %*% z, z ~ N(0, I) => u ~ N(0, Sigma_e)
  L <- tryCatch(
    t(chol(Sigma_e)),
    error = function(e) {
      # Near-singular: add tiny jitter
      t(chol(Sigma_e + diag(1e-14, n_shock)))
    }
  )

  if (n_shock == 1L) {
    ## 1-D case: nodes are L[1,1] * x1d
    sigma <- L[1, 1]
    nodes   <- matrix(sigma * x1d, ncol = 1L)
    weights <- w1d
  } else {
    ## n-D product rule: tensor product of 1-D nodes
    ## Build all combinations of indices in 1..n_gh for each dimension
    idx_list <- rep(list(seq_len(n_gh)), n_shock)
    idx_grid <- as.matrix(expand.grid(idx_list))   # n_gh^n_shock x n_shock
    n_nodes  <- nrow(idx_grid)

    ## Nodes in the standard N(0,I) space
    z_nodes <- matrix(x1d[idx_grid], nrow = n_nodes, ncol = n_shock)

    ## Weights: product of 1-D weights
    w_nodes <- apply(matrix(w1d[idx_grid], nrow = n_nodes, ncol = n_shock),
                     1L, prod)

    ## Transform to N(0, Sigma_e): u_k = L %*% z_k
    nodes   <- t(L %*% t(z_nodes))   # n_nodes x n_shock
    weights <- w_nodes / sum(w_nodes) # normalize (already should sum to 1)
  }

  list(
    nodes   = nodes,
    weights = weights,
    n_shock = n_shock,
    n_nodes = nrow(nodes),
    Sigma_e = Sigma_e
  )
}

#' 1-D Gauss-Hermite nodes and weights for Z ~ N(0,1)
#' Returns nodes x and weights w with sum(w) = 1.
#' Uses the tridiagonal eigenvalue method (Golub-Welsch).
#' @noRd
.gh_1d <- function(n) {
  n <- as.integer(n)
  if (n == 1L) {
    return(list(x = 0, w = 1))
  }
  ## Jacobi matrix for Hermite polynomials: off-diagonal = sqrt(i/2)
  i <- seq_len(n - 1L)
  J <- diag(n) * 0
  J[cbind(i, i + 1L)] <- sqrt(i / 2)
  J[cbind(i + 1L, i)] <- sqrt(i / 2)
  ev  <- eigen(J, symmetric = TRUE)
  x   <- ev$values           # nodes for Hermite w(t) = exp(-t^2)
  w   <- ev$vectors[1L, ]^2 * sqrt(pi)  # weights (unnormalized)

  ## Convert from Hermite weight exp(-t^2) to N(0,1) weight:
  ##   x_k -> sqrt(2) * x_k (for the N(0,1) density from N(0,1/2))
  ## Actually: int f(z) phi(z) dz = int f(sqrt(2) t) exp(-t^2) dt / sqrt(pi)
  ## So x_std = sqrt(2) * t_k, w_std = w_k / sqrt(pi)
  x_std <- sqrt(2) * x
  w_std <- w / sqrt(pi)
  w_std <- w_std / sum(w_std)   # ensure exact normalization

  list(x = x_std, w = w_std)
}

#' Stroud monomial-2 rule for u ~ N(0, Sigma_e): 2*n_shock nodes.
#' Integrates monomials up to total degree 3 exactly.
#' Nodes: +/- sqrt(n_shock) * e_i (in L-rotated coordinates), weight 1/(2*n_shock).
#' @noRd
.dp_monomial2_nodes <- function(Sigma_e) {
  n_shock <- nrow(Sigma_e)
  L <- tryCatch(
    t(chol(Sigma_e)),
    error = function(e) t(chol(Sigma_e + diag(1e-14, n_shock)))
  )

  ## Unit vectors +/- in R^n_shock, scaled by sqrt(n_shock)
  r      <- sqrt(as.numeric(n_shock))
  I_n    <- diag(n_shock)
  z_pos  <- r * I_n   # n_shock x n_shock
  z_neg  <- -r * I_n

  z_all  <- rbind(z_pos, z_neg)   # 2*n_shock x n_shock
  nodes  <- t(L %*% t(z_all))     # 2*n_shock x n_shock
  n_nodes <- 2L * n_shock
  weights <- rep(1 / n_nodes, n_nodes)

  list(
    nodes   = nodes,
    weights = weights,
    n_shock = n_shock,
    n_nodes = n_nodes,
    Sigma_e = Sigma_e
  )
}

#' Validate a quadrature specification and return a quadrature object.
#'
#' Accepts several forms for the \code{quadrature} argument of
#' \code{dynamic_path_perturbation()}:
#' \describe{
#'   \item{\code{NULL}}{No quadrature; returns \code{NULL} (deterministic path).}
#'   \item{integer scalar}{Number of GH points per dimension (product rule).}
#'   \item{character "monomial"}{Stroud monomial-2 rule.}
#'   \item{list with \code{nodes} and \code{weights}}{Custom rule (passed through).}
#' }
#'
#' @param quadrature  Quadrature specification (see above).
#' @param Sigma_e     Shock covariance matrix (required when quadrature is
#'   not already a pre-built rule).
#' @return Either \code{NULL} (deterministic) or a list with fields
#'   \code{nodes}, \code{weights}, \code{n_nodes}, \code{n_shock}.
#' @noRd
.dp_resolve_quadrature <- function(quadrature, Sigma_e) {
  if (is.null(quadrature)) return(NULL)

  if (is.list(quadrature) && !is.null(quadrature$nodes) &&
      !is.null(quadrature$weights)) {
    ## Pre-built rule: validate and pass through
    stopifnot(
      is.matrix(quadrature$nodes),
      is.numeric(quadrature$weights),
      nrow(quadrature$nodes) == length(quadrature$weights)
    )
    q <- quadrature
    if (is.null(q$n_nodes)) q$n_nodes <- nrow(q$nodes)
    if (is.null(q$n_shock)) q$n_shock <- ncol(q$nodes)
    return(q)
  }

  if (identical(quadrature, "monomial") || identical(quadrature, "monomial2")) {
    return(.dp_quadrature_nodes(Sigma_e, monomial = TRUE))
  }

  if (is.numeric(quadrature) && length(quadrature) == 1L) {
    n_gh <- as.integer(quadrature)
    if (n_gh < 1L) stop("quadrature n_gh must be >= 1")
    return(.dp_quadrature_nodes(Sigma_e, n_gh = n_gh))
  }

  stop("'quadrature' must be NULL, an integer (n_gh points), 'monomial', ",
       "or a list(nodes, weights).")
}
