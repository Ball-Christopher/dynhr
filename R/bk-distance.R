## R/bk-distance.R
## --------------------------------------------------------------------------
## Differentiating the Blanchard-Kahn criterion (P2 gap #12 / pathological-DSGE
## RESEARCH_AGENDA Thread 4b). At any parameter point theta the reduced
## companion-form pencil (D, E) has generalized eigenvalues lambda solving
##   E x = lambda D x,
## and BK determinacy flips exactly when an eigenvalue crosses the unit circle
## |lambda| = 1. bk_distance() returns, for the eigenvalue nearest that wall:
##   * the SIGNED distance-to-wall d = |lambda_c| - 1 (>0 explosive side),
##   * its NORMAL grad_theta |lambda_c| (the wall's direction in theta-space),
## both from the ANALYTIC generalized-eigenvalue derivative
##   d lambda = (y^H (dE - lambda dD) x) / (y^H D x),
## with dD, dE from a cheap finite difference of the PENCIL ONLY (no QZ solve,
## no eigenvalue re-ordering) via solve_perturbation's pencil_only path.
##
## Uses (per the agenda): (a) an analytic Gaussian proposal-mass-beyond-wall
## integral replacing the Monte-Carlo determinacy-collision rate
## (bk_collision_prob()); (b) the raw ingredients for a determinacy-respecting
## reparameterization. The honest tension (documented there and by the paper's
## own dead-end section): wall-AWARE sampling did not pay -- the value here is
## geometry as a DIAGNOSTIC and theory input, not another sampler.
## --------------------------------------------------------------------------


#' Build the reduced companion-form Blanchard-Kahn pencil at a parameter point
#'
#' Thin wrapper that drives \code{solve_perturbation}'s order-1 machinery only
#' as far as the reduced companion pencil \eqn{(D, E)} (no QZ decomposition),
#' so pencil derivatives can be finite-differenced without paying for -- or
#' depending on the eigenvalue ordering of -- a full solve.
#'
#' @param model,compiled A parsed \code{dynhr_mod} and its
#'   \code{\link{compile_model}} output.
#' @param params Named parameter vector (defaults to \code{model$param_values}).
#' @param ss Optional precomputed steady state (\code{$values} or a bare
#'   numeric vector). When \code{NULL}, \code{\link{solve_steady}} is called.
#' @return A list with \code{D}, \code{E} (the \eqn{p \times p} pencil),
#'   \code{n_minus} (predetermined-block size) and \code{p}; or \code{NULL} if
#'   the model has no dynamic (forward/mixed) block.
#' @keywords internal
.bk_pencil <- function(model, compiled, params = NULL, ss = NULL) {
  if (is.null(params)) params <- model$param_values
  if (is.null(ss)) {
    ssr <- solve_steady(compiled, params, endo_names = model$var_names,
                        exo_names = model$varexo_names, verbose = FALSE)
    if (!isTRUE(ssr$converged)) return(NULL)
    ss <- ssr$values
  } else if (is.list(ss) && !is.null(ss$values)) {
    ss <- ss$values
  }
  sys <- extract_system_matrices(compiled, ss, params)
  out <- tryCatch(
    .solve_from_system(sys, model, compiled, ss, params, verbose = FALSE,
                       pencil_only = TRUE),
    error = function(e) NULL)
  out
}


#' Generalized eigenvalues, right and left eigenvectors of a pencil
#'
#' Solves \eqn{E x = \lambda D x} and \eqn{y^H E = \lambda y^H D}. Returns
#' FINITE eigenvalues only (infinite generalized eigenvalues, \eqn{\beta \approx
#' 0}, are the pencil's static/redundant directions and never the BK wall).
#'
#' @keywords internal
.bk_geigen <- function(D, E, finite_tol = 1e-9) {
  ge <- geigen::geigen(E, D, symmetric = FALSE)   # E v = (alpha/beta) D v
  lam <- ge$alpha / ge$beta
  V   <- ge$vectors                                # right eigenvectors (columns)
  ## Left eigenvectors: right eigenvectors of the transposed pencil
  ## E^H w = conj(lambda) D^H w.
  geL <- geigen::geigen(Conj(t(E)), Conj(t(D)), symmetric = FALSE)
  lamL <- geL$alpha / geL$beta
  W    <- geL$vectors
  finite <- is.finite(lam) & (abs(ge$beta) > finite_tol * max(abs(ge$alpha), 1))
  list(lambda = lam, V = V, lambdaL = lamL, W = W, finite = finite)
}


#' Match each right eigenpair to its left eigenvector by eigenvalue
#'
#' The left problem \eqn{E^H w = \bar\lambda\, D^H w} (solved by
#' \code{.bk_geigen} on the conjugate-transposed pencil) returns eigenvalues
#' \eqn{\bar\lambda}, so the left eigenvector for \eqn{\lambda_i} is matched
#' against \eqn{\bar\lambda_i} -- essential for complex-conjugate pairs, where
#' matching against \eqn{\lambda_i} itself picks the WRONG partner and makes
#' \eqn{y^H D x} collapse to ~0.
#' @keywords internal
.bk_match_left <- function(lambda_i, lamL, W) {
  j <- which.min(abs(lamL - Conj(lambda_i)))
  W[, j]
}


#' Distance to the Blanchard-Kahn determinacy wall and its normal in theta-space
#'
#' At the parameter point \code{params}, finds the generalized eigenvalue of
#' the reduced companion pencil whose modulus is closest to the unit circle
#' (the \dQuote{crossing} eigenvalue \eqn{\lambda_c}) and returns the signed
#' distance-to-wall \eqn{d = |\lambda_c| - 1} together with its gradient
#' \eqn{\nabla_\theta |\lambda_c|} over the requested parameters.
#'
#' The eigenvalue gradient is analytic: for \eqn{E x = \lambda D x} with right
#' eigenvector \eqn{x} and left eigenvector \eqn{y},
#' \deqn{\frac{\partial \lambda}{\partial \theta_k} =
#'   \frac{y^H (\partial_k E - \lambda\, \partial_k D)\, x}{y^H D\, x},}
#' and \eqn{\partial_k |\lambda| = \mathrm{Re}(\bar\lambda\, \partial_k\lambda)/
#' |\lambda|}. Only the PENCIL is finite-differenced (via the
#' \code{pencil_only} path), so \eqn{\partial_k D, \partial_k E} are cheap and
#' free of eigenvalue-ordering discontinuities; the eigen-solve happens once,
#' at the base point.
#'
#' @param model,compiled A parsed \code{dynhr_mod} and its compiled form.
#' @param params Named parameter vector (defaults to \code{model$param_values}).
#' @param param_names Character vector of parameters to differentiate with
#'   respect to (defaults to \code{names(params)} intersected with the model's
#'   parameters). The returned \code{normal}/\code{grad_lambda} are named by
#'   these.
#' @param ss Optional precomputed steady state (see \code{.bk_pencil}). NOTE:
#'   held FIXED across the finite-difference pencil perturbations -- this
#'   measures the pencil's EXPLICIT parameter dependence. For parameters that
#'   also move the steady state, set \code{resolve_ss = TRUE} to re-solve the
#'   steady state at each perturbed parameter (the full total derivative).
#' @param resolve_ss Logical; when \code{TRUE} (default) each finite-difference
#'   pencil is built at a freshly-solved steady state (total derivative). When
#'   \code{FALSE}, the base steady state is reused (partial derivative; faster,
#'   correct for parameters that do not enter the steady state such as most
#'   shock/policy coefficients).
#' @param h Relative finite-difference step for the pencil (default
#'   \code{1e-6}).
#'
#' @return An object of class \code{dynhr_bk_distance}: a list with
#'   \item{distance}{signed \eqn{|\lambda_c| - 1} (positive = explosive side).}
#'   \item{lambda_c}{the crossing generalized eigenvalue (complex).}
#'   \item{modulus}{\eqn{|\lambda_c|}.}
#'   \item{normal}{named numeric \eqn{\nabla_\theta |\lambda_c|} (the wall
#'     normal; also \code{grad_distance}, identical since the wall is at
#'     \eqn{|\lambda| = 1}).}
#'   \item{grad_lambda}{named complex \eqn{\partial_\theta \lambda_c}.}
#'   \item{determinate}{logical Blanchard-Kahn verdict at \code{params}.}
#'   \item{n_unstable,n_forward}{the BK counts behind \code{determinate}.}
#' @seealso \code{\link{bk_collision_prob}}
#' @examples
#' \dontrun{
#' m  <- parse_mod(test_model_path("nk_small"))
#' cm <- compile_model(m)
#' bk <- bk_distance(m, cm, param_names = c("psi1", "psi2"))
#' bk$distance         # how far from the determinacy wall
#' bk$normal           # which way the wall lies in (psi1, psi2)
#' }
#' @export
bk_distance <- function(model, compiled, params = NULL, param_names = NULL,
                        ss = NULL, resolve_ss = TRUE, h = 1e-6) {
  if (is.null(params)) params <- model$param_values
  model_pars <- names(model$param_values)
  if (is.null(param_names)) {
    param_names <- intersect(names(params), model_pars)
  } else {
    bad <- setdiff(param_names, model_pars)
    if (length(bad))
      stop("bk_distance(): unknown parameter(s): ", paste(bad, collapse = ", "))
  }

  base_ss <- ss
  if (is.null(base_ss)) {
    ssr <- solve_steady(compiled, params, endo_names = model$var_names,
                        exo_names = model$varexo_names, verbose = FALSE)
    if (!isTRUE(ssr$converged))
      stop("bk_distance(): steady state did not converge at `params`.")
    base_ss <- ssr$values
  } else if (is.list(base_ss) && !is.null(base_ss$values)) {
    base_ss <- base_ss$values
  }

  pen <- .bk_pencil(model, compiled, params, base_ss)
  if (is.null(pen))
    stop("bk_distance(): model has no dynamic (forward/mixed) block -- the ",
         "Blanchard-Kahn pencil is empty; determinacy is trivial.")
  D <- pen$D; E <- pen$E

  eg <- .bk_geigen(D, E)
  fin <- eg$finite
  if (!any(fin))
    stop("bk_distance(): no finite generalized eigenvalues in the pencil.")
  lam_fin <- eg$lambda[fin]

  ## Crossing eigenvalue: modulus nearest the unit circle among finite roots.
  ic <- which(fin)[which.min(abs(abs(lam_fin) - 1))]
  lambda_c <- eg$lambda[ic]
  x <- eg$V[, ic]
  y <- .bk_match_left(lambda_c, eg$lambdaL, eg$W)
  denom <- as.complex(sum(Conj(y) * (D %*% x)))
  if (Mod(denom) < .Machine$double.eps^0.5)
    stop("bk_distance(): degenerate left/right eigenvector overlap (y^H D x ~ ",
         "0) at the crossing eigenvalue -- eigenvalue is (near) defective; ",
         "the analytic derivative is not defined here.")

  ## BK verdict (finite unstable count vs forward-looking count).
  mods <- Mod(eg$lambda[fin])
  n_unstable <- sum(mods > 1 + 1e-8)
  n_forward  <- pen$p - pen$n_minus
  determinate <- (n_unstable == n_forward)

  ## --- analytic gradient via FD of the pencil only -------------------------
  grad_lambda <- setNames(complex(length(param_names)), param_names)
  for (k in param_names) {
    step <- h * max(abs(params[[k]]), 1)
    pp <- params; pp[[k]] <- pp[[k]] + step
    pm <- params; pm[[k]] <- pm[[k]] - step
    if (resolve_ss) {
      penp <- .bk_pencil(model, compiled, pp, NULL)
      penm <- .bk_pencil(model, compiled, pm, NULL)
    } else {
      penp <- .bk_pencil(model, compiled, pp, base_ss)
      penm <- .bk_pencil(model, compiled, pm, base_ss)
    }
    if (is.null(penp) || is.null(penm) ||
        !identical(dim(penp$D), dim(D)) || !identical(dim(penm$D), dim(D))) {
      ## variable classification changed under perturbation (kink); skip.
      grad_lambda[[k]] <- NA_complex_
      next
    }
    dD <- (penp$D - penm$D) / (2 * step)
    dE <- (penp$E - penm$E) / (2 * step)
    grad_lambda[[k]] <- sum(Conj(y) * ((dE - lambda_c * dD) %*% x)) / denom
  }

  modulus <- Mod(lambda_c)
  normal <- Re(Conj(lambda_c) * grad_lambda) / modulus   # d|lambda_c|/dtheta

  structure(list(
    distance = modulus - 1, lambda_c = lambda_c, modulus = modulus,
    normal = normal, grad_distance = normal, grad_lambda = grad_lambda,
    determinate = determinate, n_unstable = n_unstable, n_forward = n_forward,
    param_names = param_names),
    class = "dynhr_bk_distance")
}


#' @export
print.dynhr_bk_distance <- function(x, ...) {
  cat("Blanchard-Kahn distance-to-wall\n")
  cat(sprintf("  crossing eigenvalue : %.6g%+.6gi  (|lambda| = %.6f)\n",
              Re(x$lambda_c), Im(x$lambda_c), x$modulus))
  cat(sprintf("  signed distance     : %+.6e  (%s side)\n",
              x$distance, if (x$distance > 0) "explosive" else "stable"))
  cat(sprintf("  determinate         : %s  (%d unstable / %d forward)\n",
              x$determinate, x$n_unstable, x$n_forward))
  cat("  wall normal (grad |lambda_c|):\n")
  for (k in x$param_names)
    cat(sprintf("    %-16s % .6e\n", k, x$normal[[k]]))
  invisible(x)
}


#' Analytic Gaussian proposal-mass-beyond-the-determinacy-wall
#'
#' Linearizes the determinacy wall \eqn{|\lambda_c(\theta)| = 1} at
#' \code{params} and returns the probability that a Gaussian proposal
#' \eqn{\mathcal N(\theta, \Sigma)} lands on the OTHER side of the wall (a
#' determinacy \dQuote{collision}) -- the analytic replacement for a
#' Monte-Carlo collision rate. With signed distance \eqn{d = |\lambda_c| - 1}
#' and normal \eqn{g = \nabla_\theta |\lambda_c|}, the wall is at signed
#' distance \eqn{d} in the direction \eqn{g}; the proposal's projection onto
#' \eqn{g} is Gaussian with sd \eqn{\sigma_g = \sqrt{g^\top \Sigma g}}, so the
#' crossing probability is \eqn{\Phi(-|d| / \sigma_g)}.
#'
#' @param bk A \code{\link{bk_distance}} result at the proposal centre.
#' @param Sigma Proposal covariance (\eqn{P \times P}, ordered by
#'   \code{bk$param_names}) or a scalar variance (isotropic).
#' @return A list with \code{prob} (collision probability), \code{sigma_g}
#'   (proposal sd along the wall normal), and \code{z} (\eqn{-|d|/\sigma_g}).
#' @seealso \code{\link{bk_distance}}
#' @export
bk_collision_prob <- function(bk, Sigma) {
  if (!inherits(bk, "dynhr_bk_distance"))
    stop("bk_collision_prob(): `bk` must be a bk_distance() result.")
  g <- bk$normal[bk$param_names]
  if (anyNA(g))
    stop("bk_collision_prob(): the wall normal has NA entries (a variable ",
         "classification kink under perturbation) -- collision probability ",
         "is undefined.")
  P <- length(g)
  if (length(Sigma) == 1L) {
    var_g <- as.numeric(Sigma) * sum(g^2)
  } else {
    Sigma <- as.matrix(Sigma)
    if (!all(dim(Sigma) == P))
      stop("bk_collision_prob(): `Sigma` must be ", P, "x", P,
           " (ordered by bk$param_names) or a scalar.")
    var_g <- as.numeric(t(g) %*% Sigma %*% g)
  }
  sigma_g <- sqrt(var_g)
  z <- -abs(bk$distance) / sigma_g
  list(prob = stats::pnorm(z), sigma_g = sigma_g, z = z)
}


#' Analytic parameter-Jacobian of the solution pencil's generalized eigenvalues
#'
#' The generalized eigenvalues \eqn{\lambda} of the reduced companion pencil
#' \eqn{E x = \lambda D x} ARE the Blanchard-Kahn spectrum (the same values in
#' \code{solve_perturbation}'s \code{dr$eigenvalues}). This returns every
#' finite eigenvalue together with its analytic derivative
#' \eqn{\partial \lambda_i / \partial \theta_k} over the requested parameters,
#' via the generalized-eigenvalue-derivative identity
#' \deqn{\partial_k \lambda_i =
#'   \frac{y_i^H (\partial_k E - \lambda_i \partial_k D) x_i}{y_i^H D\, x_i},}
#' with the pencil derivatives \eqn{\partial_k D, \partial_k E} from a
#' finite difference of the PENCIL ONLY (no QZ decomposition, no eigenvalue
#' re-ordering). This is the analytically-available \dQuote{spectrum}
#' half of differentiating the generalized Schur / QZ solve (Tier-18 A2
#' research remainder): the eigenVALUE Jacobian is exact and cheap; the
#' eigenVECTOR / stable-subspace (Schur-vector) Jacobian -- the genuinely hard
#' and, per the pathological-DSGE RESEARCH_AGENDA, DEMOTED half -- is not
#' needed for the first-order solution adjoint (\code{.solution_adjoint},
#' which reverses the generalized-Sylvester fixed point directly).
#'
#' @inheritParams bk_distance
#' @return An object of class \code{dynhr_pencil_spectrum}:
#'   \item{eigenvalues}{complex vector of the finite generalized eigenvalues.}
#'   \item{modulus}{their moduli.}
#'   \item{jacobian}{a \code{length(eigenvalues) x length(param_names)} complex
#'     matrix \eqn{\partial \lambda_i / \partial \theta_k} (row = eigenvalue,
#'     column = parameter; \code{NA} in a column where a variable-classification
#'     kink changed the pencil dimension under perturbation).}
#'   \item{modulus_jacobian}{the real \eqn{\partial |\lambda_i| / \partial
#'     \theta_k} matrix (the determinacy-relevant part).}
#'   \item{param_names}{the differentiated parameters (column order).}
#' @seealso \code{\link{bk_distance}} (the crossing-eigenvalue specialization).
#' @export
solution_pencil_spectrum <- function(model, compiled, params = NULL,
                                     param_names = NULL, ss = NULL,
                                     resolve_ss = TRUE, h = 1e-6) {
  if (is.null(params)) params <- model$param_values
  model_pars <- names(model$param_values)
  if (is.null(param_names)) {
    param_names <- intersect(names(params), model_pars)
  } else {
    bad <- setdiff(param_names, model_pars)
    if (length(bad))
      stop("solution_pencil_spectrum(): unknown parameter(s): ",
           paste(bad, collapse = ", "))
  }

  base_ss <- ss
  if (is.null(base_ss)) {
    ssr <- solve_steady(compiled, params, endo_names = model$var_names,
                        exo_names = model$varexo_names, verbose = FALSE)
    if (!isTRUE(ssr$converged))
      stop("solution_pencil_spectrum(): steady state did not converge.")
    base_ss <- ssr$values
  } else if (is.list(base_ss) && !is.null(base_ss$values)) {
    base_ss <- base_ss$values
  }

  pen <- .bk_pencil(model, compiled, params, base_ss)
  if (is.null(pen))
    stop("solution_pencil_spectrum(): model has no dynamic block (empty pencil).")
  D <- pen$D; E <- pen$E
  eg <- .bk_geigen(D, E)
  fin <- which(eg$finite)
  if (!length(fin))
    stop("solution_pencil_spectrum(): no finite generalized eigenvalues.")
  lam <- eg$lambda[fin]

  ## precompute left/right eigenvectors and denominators
  denom <- complex(length(fin)); Xr <- vector("list", length(fin))
  Yl <- vector("list", length(fin))
  for (a in seq_along(fin)) {
    i <- fin[a]; x <- eg$V[, i]
    y <- .bk_match_left(eg$lambda[i], eg$lambdaL, eg$W)
    Xr[[a]] <- x; Yl[[a]] <- y
    denom[a] <- sum(Conj(y) * (D %*% x))
  }

  J <- matrix(NA_complex_, length(fin), length(param_names),
              dimnames = list(NULL, param_names))
  for (k in param_names) {
    step <- h * max(abs(params[[k]]), 1)
    pp <- params; pp[[k]] <- pp[[k]] + step
    pm <- params; pm[[k]] <- pm[[k]] - step
    if (resolve_ss) {
      penp <- .bk_pencil(model, compiled, pp, NULL)
      penm <- .bk_pencil(model, compiled, pm, NULL)
    } else {
      penp <- .bk_pencil(model, compiled, pp, base_ss)
      penm <- .bk_pencil(model, compiled, pm, base_ss)
    }
    if (is.null(penp) || is.null(penm) ||
        !identical(dim(penp$D), dim(D)) || !identical(dim(penm$D), dim(D)))
      next                                    # kink; leave column NA
    dD <- (penp$D - penm$D) / (2 * step)
    dE <- (penp$E - penm$E) / (2 * step)
    for (a in seq_along(fin)) {
      if (Mod(denom[a]) < .Machine$double.eps^0.5) next
      J[a, k] <- sum(Conj(Yl[[a]]) *
                     ((dE - lam[a] * dD) %*% Xr[[a]])) / denom[a]
    }
  }

  modulus <- Mod(lam)
  modJ <- Re(sweep(Conj(lam) * J, 1L, modulus, "/"))

  structure(list(eigenvalues = lam, modulus = modulus, jacobian = J,
                 modulus_jacobian = modJ, param_names = param_names),
            class = "dynhr_pencil_spectrum")
}


#' @export
print.dynhr_pencil_spectrum <- function(x, ...) {
  cat(sprintf("Solution pencil spectrum: %d finite generalized eigenvalue(s)\n",
              length(x$eigenvalues)))
  ord <- order(abs(x$modulus - 1))
  for (a in ord) {
    cat(sprintf("  lambda = %9.5g%+9.5gi  |lambda| = %.5f%s\n",
                Re(x$eigenvalues[a]), Im(x$eigenvalues[a]), x$modulus[a],
                if (x$modulus[a] > 1 + 1e-8) "  (unstable)" else ""))
  }
  cat("  d|lambda|/dtheta by parameter:", paste(x$param_names, collapse = ", "),
      "\n")
  invisible(x)
}
