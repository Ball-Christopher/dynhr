## R/hank-survey.R
## --------------------------------------------------------------------------
## Standard survey-statistics calibration: reweight a design-weighted sample
## so that known auxiliary totals are reproduced exactly (or, for the
## entropy/raking variant, in the limit of Newton convergence).
##
## This is a SELF-CONTAINED, generic survey-statistics utility -- it depends
## on no HANK object and knows nothing about DSGE/HANK internals.  It exists
## to eventually bridge raw weighted household survey cross-sections to HANK
## macro aggregates (e.g. reweighting a cross-section so its weighted wealth
## distribution matches known population marginals before feeding it into a
## het-block), but that bridge is future work; here we only implement and
## validate the calibration machinery itself.
##
## Two calibration distances, both classic in the survey literature:
##
##   - "greg"    Chi-square (linear) distance -> closed-form linear
##               calibration a la Deville & Sarndal (1992, JASA) / the GREG
##               estimator of Sarndal, Swensson & Wretman (1992).  Weights
##               can in principle go negative or to zero; when the columns
##               of Z are stratum indicators this collapses exactly to
##               poststratification.
##
##   - "entropy" Kullback-Leibler distance -> multiplicative (raking /
##               entropy-balancing) calibration a la Deville & Sarndal
##               (1992, Sec. 2) / Hainmueller (2012, Pol. Analysis).  Weights
##               are always strictly positive; solved by Newton iteration on
##               the dual (Lagrange-multiplier) parameters.
##
## Both methods solve the same constrained problem -- minimize a distance
## between d_i and w_i subject to Sum_i w_i z_i = X -- they differ only in
## the distance function, which is what determines the functional form of
## w_i(d_i, z_i, lambda).
## --------------------------------------------------------------------------


#' Calibrate design weights to known population control totals (GREG / raking)
#'
#' Reweights a design-weighted sample \code{(d_i, z_i)} so that the
#' calibrated weights \code{w_i*} reproduce known auxiliary/control totals
#' exactly: \code{Sum_i w_i* z_i = X}. This is the standard survey-statistics
#' calibration-estimation problem (Deville & Sarndal 1992; Sarndal, Swensson &
#' Wretman 1992 Ch. 6 for the linear/GREG case; Hainmueller 2012 for the
#' entropy-balancing / raking case).  Both methods minimize a distance between
#' the design weights \code{d} and the calibrated weights \code{w*} subject to
#' the calibration constraint; they differ only in the choice of distance.
#'
#' \strong{method = "greg"} (chi-square distance): the linear-calibration
#' case has a closed-form solution (Deville & Sarndal 1992, eq. 2.1-2.5). With
#' \code{T = Sum_i d_i z_i t(z_i)} (the \code{p x p} weighted cross-product of
#' \code{Z}) and total gap \code{X - Sum_i d_i z_i}, the g-weight for unit i is
#' \deqn{g_i = 1 + (X - \Sigma_i d_i z_i)^{T} T^{-1} z_i}
#' and the calibrated weight is \code{w_i* = d_i g_i}. This is an exact,
#' single linear-algebra step -- no iteration. Calibrated weights can in
#' principle be negative if \code{d} is small or \code{Z} spans extreme
#' values; when the columns of \code{Z} are disjoint 0/1 stratum indicators
#' (poststratification), the GREG solution collapses to the classical
#' poststratified weight (stratum control total / stratum sample count) for
#' every unit in the stratum, and is always positive in that special case.
#'
#' \strong{method = "entropy"} (Kullback-Leibler distance): weights take the
#' exponential-tilt form \code{w_i* = d_i * exp(t(z_i) lambda)}, which is always
#' strictly positive.  There is no closed form for \code{lambda} in general
#' (except in the pure-indicator poststratification case, where it again
#' collapses to poststratified weights); \code{lambda} solves the p equations
#' \code{Sum_i w_i*(lambda) z_i = X} by the Newton method, using the exact
#' gradient \code{Sum_i w_i* z_i} and Hessian \code{Sum_i w_i* z_i t(z_i)} of the
#' dual objective at each step.  Iteration stops when
#' \code{max(abs(Sum_i w_i* z_i - X)) < tol} or after \code{maxit} steps; if
#' it has not converged by then, \code{converged} is set \code{FALSE} and a
#' \code{warning()} is raised (the function does not error, since a
#' near-converged but imperfect calibration can still be a usable diagnostic).
#'
#' @param d Numeric length-N vector of design weights, \code{d_i = 1 / pi_i}
#'   (the number of population units unit i represents under the sampling
#'   design). Must be strictly positive.
#' @param Z Numeric \code{N x p} matrix of auxiliary/calibration variables:
#'   row i is \code{z(x_i)}, the calibration-variable vector for unit i (e.g.
#'   an intercept column, stratum indicator columns, or continuous
#'   covariates).
#' @param X Numeric length-p vector of known population control totals: the
#'   calibration constraint solved for is \code{Sum_i w_i* z_i = X}.
#' @param method Calibration distance: \code{"greg"} (chi-square / linear,
#'   closed-form) or \code{"entropy"} (Kullback-Leibler / raking, Newton
#'   iteration). Partial matching via \code{\link{match.arg}}.
#' @param tol Convergence tolerance on \code{max(abs(Sum_i w_i* z_i - X))} for
#'   \code{method = "entropy"}. Ignored (constraint is exact to machine
#'   precision) for \code{method = "greg"}.
#' @param maxit Maximum Newton iterations for \code{method = "entropy"}.
#'   Ignored for \code{method = "greg"}.
#'
#' @return A list with elements:
#'   \item{w}{Length-N calibrated weights \code{w_i*}.}
#'   \item{g}{Length-N g-weights \code{w_i* / d_i} (the multiplicative
#'     adjustment applied to each design weight).}
#'   \item{lambda}{Length-p vector: for \code{method = "entropy"}, the
#'     converged (or best-effort, if \code{!converged}) dual/Lagrange
#'     multiplier vector \code{lambda} solving
#'     \code{Sum_i d_i exp(t(z_i) lambda) z_i = X}. For \code{method = "greg"},
#'     the corresponding linear multiplier
#'     \code{t(X - Sum_i d_i z_i) T^{-1}} (a length-p vector), i.e. the object
#'     such that \code{g_i = 1 + t(lambda) z_i}, returned for symmetry with the
#'     entropy case and for downstream diagnostics.}
#'   \item{converged}{Logical. Always \code{TRUE} for \code{method = "greg"}
#'     (the linear system is solved exactly, up to the \code{rcond} guard
#'     below). For \code{method = "entropy"}, \code{TRUE} iff the Newton
#'     iteration met \code{tol} within \code{maxit} steps.}
#'   \item{method}{The method used (character scalar), echoed back.}
#'   \item{vhat}{A \code{p x p} plug-in variance-style estimate of the
#'     variance of the calibrated total \code{Sum_i w_i* z_i}, computed as
#'     \code{Sum_i (w_i*)^2 (z_i - zbar*) t(z_i - zbar*)}, where \code{zbar*}
#'     is the calibrated-weighted mean of \code{Z}
#'     (\code{Sum_i w_i* z_i / Sum_i w_i*}). This is a simple design-weighted
#'     plug-in dispersion measure of the calibration variables around their
#'     calibrated mean -- it is NOT a full linearized/replicate-variance
#'     estimator of the calibration estimator itself (which would require
#'     the original sampling design, e.g. stratification/clustering, that
#'     this function does not observe). It is intended only as a
#'     loss-metric building block for later use, not as a publication-grade
#'     variance estimate.}
#'
#' @references
#' Deville, J.-C. and Sarndal, C.-E. (1992). "Calibration Estimators in
#' Survey Sampling." \emph{Journal of the American Statistical Association},
#' 87(418), 376-382.
#'
#' Sarndal, C.-E., Swensson, B. and Wretman, J. (1992). \emph{Model Assisted
#' Survey Sampling}. Springer. (Ch. 6: the generalized regression, GREG,
#' estimator.)
#'
#' Hainmueller, J. (2012). "Entropy Balancing for Causal Effects: A
#' Multivariate Reweighting Method to Produce Balanced Samples in
#' Observational Studies." \emph{Political Analysis}, 20(1), 25-46.
#'
#' @examples
#' ## Two strata, design weights deliberately wrong; control totals are the
#' ## true stratum population counts (60, 40). Both methods recover the
#' ## classical poststratified weight within each stratum.
#' Z <- cbind(strat1 = c(1, 1, 0, 0, 0), strat2 = c(0, 0, 1, 1, 1))
#' d <- c(10, 10, 5, 5, 5)
#' X <- c(60, 40)
#' cal <- hank_calibrate_weights(d, Z, X, method = "greg")
#' cal$w   # 30, 30 (=60/2) in stratum 1; 40/3 each in stratum 2
#'
#' @export
hank_calibrate_weights <- function(d, Z, X,
                                    method = c("greg", "entropy"),
                                    tol = 1e-10, maxit = 100L) {
  method <- match.arg(method)

  d <- as.numeric(d)
  Z <- as.matrix(Z)
  X <- as.numeric(X)
  N <- length(d)
  p <- ncol(Z)

  if (nrow(Z) != N)
    stop(sprintf(
      "hank_calibrate_weights(): nrow(Z) (%d) must equal length(d) (%d).",
      nrow(Z), N))
  if (length(X) != p)
    stop(sprintf(
      "hank_calibrate_weights(): ncol(Z) (%d) must equal length(X) (%d).",
      p, length(X)))
  if (any(!is.finite(d)) || any(d <= 0))
    stop("hank_calibrate_weights(): all entries of d must be finite and > 0.")

  ## Weighted cross-product T = Sum_i d_i z_i z_i' (p x p), and the total gap
  ## X - Sum_i d_i z_i.  Both methods need T (entropy: as the Newton Hessian
  ## at lambda = 0 is not required, but the greg branch needs T^{-1} and the
  ## rcond guard below applies to that same matrix shape/conditioning issue
  ## that would also afflict a degenerate Newton Hessian in the entropy case;
  ## checking it once up front fails loud for both methods on a rank-deficient
  ## Z, mirroring the rcond guard in hank_model_irf(), R/hank-model.R).
  Tmat <- crossprod(Z, d * Z)
  rc <- rcond(Tmat)
  if (!is.finite(rc) || rc < 1e-10)
    stop(sprintf(paste0(
      "hank_calibrate_weights(): the calibration cross-product matrix ",
      "t(Z) %%*%% diag(d) %%*%% Z is numerically singular (rcond = %.2e); ",
      "Z is likely rank-deficient (e.g. a duplicated or collinear column). ",
      "Remove redundant calibration variables and retry."), rc))

  gap <- X - as.numeric(crossprod(Z, d))

  if (method == "greg") {
    lambda <- as.numeric(solve(Tmat, gap))          # length-p multiplier
    g <- 1 + as.numeric(Z %*% lambda)
    w <- d * g
    converged <- TRUE
  } else {
    ## Entropy / raking: Newton iteration on lambda solving
    ##   Sum_i d_i exp(z_i' lambda) z_i = X
    ## Start at lambda = 0 (w = d), gradient = Sum_i d_i z_i - X = -gap,
    ## Hessian = Sum_i d_i z_i z_i' = Tmat at lambda = 0.
    lambda <- rep(0, p)
    converged <- FALSE
    for (it in seq_len(maxit)) {
      eta <- as.numeric(Z %*% lambda)
      w <- d * exp(eta)
      resid <- as.numeric(crossprod(Z, w)) - X       # Sum_i w_i z_i - X
      if (max(abs(resid)) < tol) {
        converged <- TRUE
        break
      }
      H <- crossprod(Z, w * Z)                        # Sum_i w_i z_i z_i'
      rc_h <- rcond(H)
      if (!is.finite(rc_h) || rc_h < 1e-10)
        stop(sprintf(paste0(
          "hank_calibrate_weights(): the entropy Newton Hessian became ",
          "numerically singular (rcond = %.2e) at iteration %d; Z is ",
          "likely rank-deficient (e.g. a duplicated or collinear column)."),
          rc_h, it))
      step <- solve(H, resid)
      lambda <- lambda - step
    }
    ## Final weights at the (possibly non-converged) lambda.
    eta <- as.numeric(Z %*% lambda)
    w <- d * exp(eta)
    if (!converged) {
      resid <- as.numeric(crossprod(Z, w)) - X
      warning(sprintf(paste0(
        "hank_calibrate_weights(): entropy/raking Newton iteration did not ",
        "converge within maxit = %d steps (max|residual| = %.3e, tol = ",
        "%.3e). Returning the best-effort weights from the final ",
        "iteration; treat 'converged = FALSE' as a signal to raise maxit ",
        "or inspect Z/X for near-infeasibility."),
        maxit, max(abs(resid)), tol))
    }
  }

  ## Plug-in variance-style handle: dispersion of Z around its
  ## calibrated-weighted mean, weighted by (w_i*)^2. See @return docs above
  ## for exactly what this is (and is not).
  wsum <- sum(w)
  zbar <- as.numeric(crossprod(w, Z)) / wsum
  Zc <- sweep(Z, 2, zbar, "-")
  vhat <- crossprod(Zc, (w^2) * Zc)

  list(w = w, g = w / d, lambda = lambda, converged = converged,
       method = method, vhat = vhat)
}
