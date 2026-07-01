## R/metric-monge.R
## --------------------------------------------------------------------------
## Stage 3a (manifold-MCMC roadmap): Monge-metric factory for smMALA.
##
## The Monge metric is built from the log-posterior gradient alone -- no
## Hessian, no Fisher information, no third derivatives.  It encodes local
## geometry through gradient magnitude: where the log-density is steep the
## metric stretches, damping step sizes in that direction.
##
## Theory (Hartmann, Girolami & Klami, AISTATS 2022, arXiv:2202.00755):
##   Embed the target density as a Monge patch (graph of log pi) and use the
##   induced Riemannian metric:
##
##     G(theta) = I + alpha^2 * g g'   (rank-1 update of identity)
##
##   where g = grad log pi(theta).  Closed forms via Sherman--Morrison:
##
##     G_inv = I - (alpha^2 / (1 + alpha^2 |g|^2)) * g g'
##     logdet G = log(1 + alpha^2 |g|^2)
##
##   and the Cholesky L = chol(G) can be computed cheaply (rank-1 Cholesky
##   update of the identity).
##
## Usage in smMALA:
##   The factory monge_metric_fn(grad_fn, alpha) returns a closure
##   function(theta) -> list(G, G_inv, L, logdet) compatible with the
##   `metric_fn` slot of dynhr_mala().  No new integrator is needed -- the
##   existing smMALA MH machinery handles the position-dependent correction
##   (the logdet term in .mala_log_q_pd).
##
## Non-finite gradient guard:
##   When grad_fn(theta) returns non-finite values (off-support), the closure
##   returns the identity metric instead of propagating NaN/Inf into the
##   metric matrices.  dynhr_mala() will also auto-reject the proposal if
##   the gradient is non-finite at the proposed position, so this guard is
##   belt-and-suspenders.
## --------------------------------------------------------------------------


#' Monge-metric factory for simplified MMALA (smMALA)
#'
#' Returns a closure \code{function(theta) -> list(G, G_inv, L, logdet)}
#' that evaluates the Monge metric
#' \eqn{G(\theta) = I + \alpha^2\, g g^T} at any position \eqn{\theta},
#' where \eqn{g = \nabla\log\pi(\theta)}.
#'
#' The Monge metric is a rank-1 update of the identity built from the
#' gradient only -- no Hessian, no Fisher information, no third derivatives.
#' Sherman--Morrison gives exact closed forms for \eqn{G^{-1}},
#' \eqn{\log|G|}, and the Cholesky \eqn{L = \text{chol}(G)}.
#'
#' The returned closure is compatible with the \code{metric_fn} argument of
#' \code{\link{dynhr_mala}} (Stage 2 / simplified MMALA interface).
#'
#' @param grad_fn  Gradient function, \code{function(theta) -> numeric vector}
#'   of length \eqn{d}.  Should return \code{NA} or \code{NaN} off-support
#'   (the closure returns the identity metric for non-finite gradient values).
#' @param alpha  Softness parameter controlling how strongly the gradient
#'   magnitude rescales the metric (default 1).  At \code{alpha = 0} the
#'   metric is the identity (plain MALA).  Larger \code{alpha} gives a more
#'   pronounced position-dependent geometry.
#'
#' @return A closure \code{function(theta) -> list(G, G_inv, L, logdet)}:
#'   \describe{
#'     \item{G}{\eqn{d \times d} PD metric matrix (identity + rank-1 term).}
#'     \item{G_inv}{Exact Sherman--Morrison inverse.}
#'     \item{L}{Upper Cholesky factor of \code{G} (rank-1 Cholesky update).}
#'     \item{logdet}{\eqn{\log|G| = \log(1 + \alpha^2 |g|^2)}.}
#'   }
#'   At positions where the gradient is non-finite, returns the identity
#'   metric: \code{G = G_inv = I}, \code{logdet = 0}.
#'
#' @references
#' Hartmann, Girolami & Klami (2022).  "Lagrangian Manifold Monte Carlo on
#' Monge Patches."  \emph{AISTATS 2022}, PMLR v151. arXiv:2202.00755.
#'
#' @noRd
monge_metric_fn <- function(grad_fn, alpha = 1) {
  stopifnot(is.function(grad_fn))
  stopifnot(is.numeric(alpha), length(alpha) == 1L, is.finite(alpha), alpha >= 0)

  function(theta) {
    d <- length(theta)

    ## Evaluate gradient at current position
    g <- grad_fn(theta)

    ## Guard: return identity metric when gradient is non-finite (off-support)
    if (length(g) != d || any(!is.finite(g))) {
      I_d <- diag(d)
      return(list(G = I_d, G_inv = I_d, L = I_d, logdet = 0))
    }

    ## alpha = 0: G = I exactly (plain MALA, no position-dependence)
    if (alpha == 0) {
      I_d <- diag(d)
      return(list(G = I_d, G_inv = I_d, L = I_d, logdet = 0))
    }

    ## Squared gradient norm |g|^2
    g2 <- sum(g^2)
    a2 <- alpha^2

    ## ---- G = I + alpha^2 g g' -------------------------------------------
    ## Use outer product form; for d-dimensional g this is a d x d rank-1
    ## update of the identity.
    gg <- outer(g, g)           ## d x d; gg[i,j] = g[i]*g[j]
    G  <- diag(d) + a2 * gg

    ## ---- G_inv via Sherman--Morrison: I - (a2/(1+a2|g|^2)) g g' ----------
    ## When a2*g2 is very small, the correction is essentially zero; the
    ## formula is numerically stable since 1 + a2*g2 >= 1 always.
    sm_coef <- a2 / (1 + a2 * g2)
    G_inv   <- diag(d) - sm_coef * gg

    ## ---- logdet G = log(1 + alpha^2 |g|^2) --------------------------------
    ## Matrix-determinant lemma: det(I + u v') = 1 + v'u.
    ## Here v = u = alpha*g, so det(G) = 1 + alpha^2 |g|^2.
    logdet <- log1p(a2 * g2)

    ## ---- Cholesky of G (rank-1 update of identity) -------------------------
    ## G = I + alpha^2 g g' = I + (alpha g)(alpha g)'.
    ## Rank-1 Cholesky update of I:  chol(I + u u') where u = alpha * g.
    ## The standard result for a rank-1 update of I gives an upper triangular L
    ## satisfying t(L) %*% L = G.
    ## We use chol() directly -- for small d this is cheap, and G is
    ## guaranteed PD (G = I + rank-1 SPD term with positive coefficient).
    L <- tryCatch(
      chol(G),
      error = function(e) {
        ## Paranoia: force symmetry and add a tiny nugget if chol() fails
        G_sym <- (G + t(G)) / 2
        chol(G_sym + diag(1e-14 * (1 + a2 * g2), d))
      }
    )

    list(G = G, G_inv = G_inv, L = L, logdet = logdet)
  }
}
