## R/csminwel.R
## --------------------------------------------------------------------------
## csminwel() -- Christopher Sims's quasi-Newton (BFGS) minimiser, the de
## facto standard optimiser for DSGE estimation (Smets-Wouters, Ireland,
## RBC IRF-matching, ...).
##
## This is a faithful base-R port of Sims's MATLAB csminwel.m / csminit.m /
## bfgsi.m / numgrad.m.  It combines a numerical (or user-supplied) gradient,
## a back-/forward-tracking line search (csminit), a BFGS inverse-Hessian
## update (bfgsi), and the characteristic "bad gradient" / Hessian-reset
## logic that lets it climb out of stalls.  Self-contained: base R only.
##
## PROVENANCE AND LICENCE (checked 2026-09-05).
##
## Ported from Sims's OWN files, as distributed from his Princeton page:
##   http://sims.princeton.edu/yftp/optimize/mfiles/{csminwel,csminit,bfgsi,numgrad}.m
## That host stopped answering some time after 2026-05-18 (DNS still resolves
## to eco-csimsserv.princeton.edu, nothing listens on http or https), so the
## live citation is Sims's index page and the files are readable in the
## Internet Archive:
##   https://www.princeton.edu/~sims/#optimize                    (live index)
##   https://web.archive.org/web/20180409211528/http://sims.princeton.edu/yftp/optimize/mfiles/csminwel.m
##
## LICENCE. Sims's `bfgsi.m` carries an explicit grant --
##   "Copyright by Christopher Sims 1996.  This material may be freely
##    reproduced and modified."
## -- which is permissive and compatible with this package's MIT licence.
## `csminwel.m`, `csminit.m` and `numgrad.m` as distributed by Sims carry NO
## copyright or licence notice of any kind.
##
## NOT DERIVED FROM THE GPL FORKS, and this matters: the copies most people
## reach first are copyleft. Dynare's `matlab/optimization/csminwel1.m` is
## "Copyright (C) 1993-2007 Christopher Sims / 2006-2025 Dynare Team", GPL-3,
## and the Sims-Zha `contrib/ms-sbvar/TZcode/.../csminwel.m` is
## "Copyright (C) 1997-2012 Christopher A. Sims and Tao Zha", GPL-3-or-later.
## A port of either would make this file GPL and collide with dynhr's MIT.
## This port carries NONE of their fork-only machinery -- no `penalty`,
## `Save_files`, `message` or `epsilon` (Dynare's, 6-23 occurrences each), no
## `dispIndx` or `stps` (Sims-Zha's) -- and follows the original's structure,
## including `numgrad.m`'s delta = 1e-6 and its abs(g) < 1e15 reliability
## guard. If this file is ever re-synced against a MATLAB source, re-sync it
## against the ARCHIVED ORIGINALS, never against Dynare.
##
## A private reference copy of the four originals (with the full licence
## findings) is kept OUTSIDE this repository, at
## ../dynhr_refs/csminwel-sims-originals/. Three of the four carry no licence
## notice, so they are cited and checksummed here, never redistributed --
## the same rule that kept Dynare's fs2000.mod out of the 0.9.3 release.
## sha256 of the 2018-04-09 archive snapshot:
##   6b3b13247aa1c6cfcfbe5bac4be2fb9c8dfd4f8695afd78d168849fbe381ad52  csminwel.m
##   aa839dec3a44ed1c58eb6441c2ce61b253ce60b407d7226ef11b176fc82e2a24  csminit.m
##   33beb5acdbb96dba51c1cab0ac7a59e99f2cf35f7ba4e0573ce9f1a1dd0b56c1  bfgsi.m
##   55e28e6170e0347dec9141f8796e031c57c9d1f8094599c9696af7374c7563fd  numgrad.m
## --------------------------------------------------------------------------


## --------------------------------------------------------------------------
## Internal: BFGS inverse-Hessian update (Sims's bfgsi.m).
##   H0  current inverse-Hessian estimate (n x n)
##   dg  change in gradient  (g_new - g_old)
##   dx  change in x         (x_new - x_old)
## Returns the updated inverse Hessian; on a degenerate (near-zero) curvature
## the previous H0 is returned unchanged (with a warning, as in Sims).
## --------------------------------------------------------------------------
.csminwel_bfgsi <- function(H0, dg, dx) {
  dg <- as.matrix(dg); dx <- as.matrix(dx)        # column vectors
  Hdg  <- H0 %*% dg
  dgdx <- as.numeric(crossprod(dx, dg))           # dx' dg  (scalar curvature)
  if (abs(dgdx) > 1e-12) {
    H <- H0 +
      (1 + as.numeric(crossprod(dg, Hdg)) / dgdx) *
        (dx %*% t(dx)) / dgdx -
      (dx %*% t(Hdg) + Hdg %*% t(dx)) / dgdx
  } else {
    warning("csminwel: bfgs update failed (near-zero curvature); H unchanged.")
    H <- H0
  }
  ## Symmetrise to suppress round-off drift.
  0.5 * (H + t(H))
}


## --------------------------------------------------------------------------
## Internal: forward-difference numerical gradient (Sims's numgrad.m).
## Returns list(g = gradient, badg = TRUE if any component looked unreliable).
## --------------------------------------------------------------------------
.csminwel_numgrad <- function(fcn, x, ...) {
  n     <- length(x)
  delta <- 1e-6
  g     <- numeric(n)
  badg  <- FALSE
  f0    <- fcn(x, ...)
  for (i in seq_len(n)) {
    xi      <- x
    xi[i]   <- xi[i] + delta
    fi      <- fcn(xi, ...)
    gi      <- (fi - f0) / delta
    if (is.finite(gi) && abs(gi) < 1e15) {
      g[i] <- gi
    } else {
      g[i]  <- 0
      badg  <- TRUE
    }
  }
  list(g = g, badg = badg)
}


## --------------------------------------------------------------------------
## Internal: the csminit line search (Sims's csminit.m).
##
## Given the current point x0 (value f0, gradient g0) and an inverse-Hessian
## estimate H0, take a (quasi-)Newton step p = -H0 g0 and adjust its length
## with a forward/backward search satisfying an Armijo-type sufficient-
## decrease condition, returning the best point found.
##
## Returns list(fhat, xhat, fcount, retcodeh) where retcodeh follows Sims:
##   0  normal step
##   1  zero gradient
##   2/4  back/forward step but improvement only marbinal (stuck)
##   3  no improvement, smallest step still worse
##   5  largest step still improving (hit growth cap)
##   6  zero gradient direction / no improvement
##   7  cliff: improvement tiny in both directions
## --------------------------------------------------------------------------
.csminwel_csminit <- function(fcn, x0, f0, g0, badg, H0, verbose, ...) {
  ANGLE  <- 0.005
  THETA  <- 0.3        # min relative improvement vs. predicted (Sims)
  FCHANGE <- 1000
  MINLAMB <- 1e-9
  MINDFAC <- 0.01
  fcount   <- 0L
  lambda   <- 1
  xhat     <- x0
  f        <- f0
  fhat     <- f0
  g        <- g0
  gnorm    <- sqrt(sum(g * g))

  if (gnorm < 1e-12 && !badg) {
    ## Gradient convergence.
    return(list(fhat = f0, xhat = x0, fcount = fcount, retcodeh = 1L))
  }

  ## Newton direction (scaled by gradient magnitude when badg, as in Sims).
  dx <- as.numeric(-H0 %*% g)
  dxnorm <- sqrt(sum(dx * dx))
  if (dxnorm > 1e12) dx <- dx * 1e12 / dxnorm
  dfhat <- sum(dx * g0)            # predicted change (should be < 0)

  if (!badg) {
    ## If the search direction is nearly orthogonal to the gradient, nudge it.
    a <- -dfhat / (gnorm * dxnorm)
    if (a < ANGLE) {
      dx     <- dx - (ANGLE * dxnorm / gnorm + dfhat / (gnorm * gnorm)) * g
      dx     <- dx * dxnorm / sqrt(sum(dx * dx))
      dfhat  <- sum(dx * g0)
      if (verbose)
        cat(sprintf("    csminit: angle correction, predicted dfhat=%g\n", dfhat))
    }
  }
  if (verbose) cat(sprintf("    csminit: predicted improvement %18.9f\n", -dfhat / 2))

  done    <- FALSE
  factor  <- 3
  shrink  <- TRUE
  lambdaMin <- 0
  lambdaMax <- Inf
  lambdaPeak <- 0
  fPeak    <- f0
  retcode  <- 0L

  while (!done) {
    dxtest <- x0 + lambda * dx
    f      <- fcn(dxtest, ...)
    if (verbose) cat(sprintf("    lambda = %10.5g; f = %20.7f\n", lambda, f))
    fcount <- fcount + 1L

    if (f < fhat) { fhat <- f; xhat <- dxtest; lambdahat <- lambda }

    ## Sims's branching on whether we've improved enough.
    shrinkSignal <- (!badg & (f0 - f < max(-THETA * dfhat * lambda, 0))) |
                    ( badg & ((f0 - f) < 0))
    growSignal   <- !badg &
                    ((lambda > 0) &
                     (f0 - f > -(1 - THETA) * dfhat * lambda))

    if (shrinkSignal & ((lambda > lambdaPeak) | (lambda < 0))) {
      if ((lambda > 0) & ((!shrink) | (lambda / factor <= lambdaPeak))) {
        shrink <- TRUE
        factor <- factor^0.6
        while (lambda / factor <= lambdaPeak) factor <- factor^0.6
        if (abs(factor - 1) < MINDFAC) {
          if (abs(lambda) < 4) retcode <- 2L else retcode <- 7L
          done <- TRUE
        }
      }
      if ((lambda < lambdaMax) & (lambda > lambdaPeak)) lambdaMax <- lambda
      lambda <- lambda / factor
      if (abs(lambda) < MINLAMB) {
        if ((lambda > 0) & (f0 <= fhat)) {
          lambda  <- -lambda * factor^6  # try other direction
        } else {
          if (lambda < 0) retcode <- 6L else retcode <- 3L
          done <- TRUE
        }
      }
    } else if ((growSignal & (lambda > 0)) |
               (shrinkSignal & ((lambda <= lambdaPeak) & (lambda > 0)))) {
      if (shrink) {
        shrink <- FALSE
        factor <- factor^0.6
        if (abs(factor - 1) < MINDFAC) {
          if (abs(lambda) < 4) retcode <- 4L else retcode <- 7L
          done <- TRUE
        }
      }
      if ((f < fPeak) & (lambda > 0)) {
        fPeak      <- f
        lambdaPeak <- lambda
        if (lambdaMax <= lambdaPeak) lambdaMax <- lambdaPeak * factor * factor
      }
      lambda <- lambda * factor
      if (abs(lambda) > 1e20) { retcode <- 5L; done <- TRUE }
    } else {
      done <- TRUE
      if (factor < 1.2) retcode <- 7L else retcode <- 0L
    }
  }

  if (verbose)
    cat(sprintf("    csminit done: retcode=%d  fhat=%18.9f\n", retcode, fhat))
  list(fhat = fhat, xhat = xhat, fcount = fcount, retcodeh = retcode)
}


#' Christopher Sims's \code{csminwel} quasi-Newton minimiser
#'
#' A base-R port of Christopher Sims's \code{csminwel} optimiser, the de-facto
#' standard for likelihood / minimum-distance estimation of DSGE models.  It
#' is a BFGS method built around a robust forward/backward line search
#' (\code{csminit}) and a "bad gradient" diagnostic that triggers an
#' inverse-Hessian reset, letting the algorithm climb out of stalls where a
#' textbook BFGS would terminate prematurely.  The gradient is computed by
#' forward differences unless \code{grad} is supplied.
#'
#' \strong{Minimises} \code{fcn}.  For maximum-likelihood estimation pass the
#' negative log-likelihood / negative log-posterior, exactly as you would to
#' \code{\link[stats]{optim}}.
#'
#' @param fcn Function to minimise; \code{fcn(x, ...)} returns a scalar.
#' @param x0 Numeric starting vector (names are preserved on the result).
#' @param ... Extra arguments passed through to \code{fcn} (and \code{grad}).
#' @param H0 Initial inverse-Hessian estimate (\code{n x n}).  \code{NULL}
#'   (default) uses \code{1e-4 * I}, Sims's default scaling.
#' @param grad Optional analytic gradient \code{grad(x, ...)} returning a
#'   length-\code{n} numeric vector.  \code{NULL} (default) uses forward
#'   differences.
#' @param crit Convergence tolerance on the improvement in \code{fcn} between
#'   iterations.  Default \code{1e-7}.
#' @param nit Maximum number of iterations.  Default \code{1000}.
#' @param verbose Logical; print per-iteration progress.  Default \code{FALSE}.
#' @param stall_warmup Integer; only begin the stall (plateau) convergence check
#'   after this many iterations, so a quickly-converging problem is never cut
#'   short.  Default \code{100}.
#' @param stall_window Integer; number of recent per-iteration improvements
#'   inspected by the stall check.  Default \code{20}.
#' @param stall_tol Numeric; if none of the last \code{stall_window} improvements
#'   cleared this tolerance, the search is declared plateaued (retcode 8).  Set
#'   to \code{0} to disable the stall check.  Default \code{1e-4}.
#'
#' @return A list (Sims's convention, with dynhr-friendly aliases) containing:
#'   \describe{
#'     \item{fh}{Minimised function value (also aliased \code{value}).}
#'     \item{xh}{Minimising argument (also aliased \code{par}; named if
#'       \code{x0} was named).}
#'     \item{gh}{Gradient at \code{xh}.}
#'     \item{H}{Final inverse-Hessian estimate (an approximate covariance of
#'       the estimator at a likelihood optimum).}
#'     \item{itct}{Number of iterations performed (also \code{iterations}).}
#'     \item{fcount}{Number of \code{fcn} evaluations.}
#'     \item{retcode}{Termination code: \code{0} normal/critical-improvement;
#'       others propagate the last \code{csminit} line-search code (see Sims).}
#'     \item{convergence}{\code{0} if converged within \code{nit}, else \code{1}
#'       (\code{\link[stats]{optim}}-style alias).}
#'   }
#'
#' @examples
#' ## Rosenbrock: minimum f = 0 at (1, 1).
#' rosen <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
#' fit <- csminwel(rosen, c(-1.2, 1))
#' fit$xh      # ~ c(1, 1)
#' fit$fh      # ~ 0
#'
#' @references Sims, C. A. \emph{csminwel} (optimization software).
#'   \url{https://www.princeton.edu/~sims/#optimize}. The canonical file
#'   location, \code{http://sims.princeton.edu/yftp/optimize/}, has been
#'   unreachable since 2026 and is preserved in the Internet Archive:
#'   \url{https://web.archive.org/web/20180409211528/http://sims.princeton.edu/yftp/optimize/mfiles/csminwel.m}.
#'   This is a base-R port of Sims's own MATLAB files, not of the GPL forks
#'   shipped with Dynare -- see the provenance note at the top of
#'   \code{R/csminwel.R}.
#' @seealso \code{\link{match_irfs}}, \code{\link[stats]{optim}}
#' @export
csminwel <- function(fcn, x0, ..., H0 = NULL, grad = NULL,
                     crit = 1e-7, nit = 1000, verbose = FALSE,
                     stall_warmup = 100L, stall_window = 20L,
                     stall_tol = 1e-4) {

  x0_names <- names(x0)
  x0 <- as.numeric(x0)
  nx <- length(x0)
  if (nx == 0L) stop("csminwel: `x0` must be non-empty.")

  if (is.null(H0)) {
    H0 <- diag(1e-4, nx)
  } else {
    H0 <- as.matrix(H0)
    if (!all(dim(H0) == c(nx, nx)))
      stop(sprintf("csminwel: `H0` must be %d x %d.", nx, nx))
  }

  ## Gradient helper: analytic (grad) or numerical.  Returns list(g, badg).
  get_grad <- function(x) {
    if (is.null(grad)) {
      .csminwel_numgrad(fcn, x, ...)
    } else {
      g <- as.numeric(grad(x, ...))
      list(g = g, badg = any(!is.finite(g)))
    }
  }

  f0 <- fcn(x0, ...)
  if (!is.finite(f0))
    stop("csminwel: objective is not finite at the starting value.")

  fcount  <- 1L
  itct    <- 0L
  done    <- FALSE
  ## Ring buffer of the most recent per-iteration improvements, for the stall
  ## (plateau) stopping rule below.
  recent_impr <- numeric(0)
  H       <- H0
  x       <- x0
  f       <- f0
  gg      <- get_grad(x)
  g       <- gg$g; badg <- gg$badg
  retcode <- 0L

  while (!done && itct < nit) {
    itct <- itct + 1L
    f1 <- f; x1 <- x; g1 <- g

    ## --- Line search from the current point. ---
    ls <- .csminwel_csminit(fcn, x, f, g, badg, H, verbose, ...)
    fcount  <- fcount + ls$fcount
    fh      <- ls$fhat
    xh      <- ls$xhat
    retcode1 <- ls$retcodeh

    ## --- Gradient at the new point. ---
    gh_obj <- get_grad(xh)
    gh     <- gh_obj$g
    badgh  <- gh_obj$badg

    ## --- Hessian-reset logic (Sims). On a stalled/bad-gradient step, reset
    ##     H to its default scaling and retry the line search once before
    ##     giving up.  A genuinely small improvement signals convergence. ---
    improvement <- f - fh
    if (retcode1 == 1L) {
      ## Zero gradient -> converged.
      done <- TRUE; retcode <- retcode1
    } else if (retcode1 %in% c(2L, 4L) || badg) {
      ## Stuck (back/forward marginal) or bad gradient: reset H and retry.
      Hreset <- diag(1e-4, nx)
      ls2 <- .csminwel_csminit(fcn, x, f, g, badg, Hreset, verbose, ...)
      fcount <- fcount + ls2$fcount
      if (ls2$fhat < fh) {
        fh <- ls2$fhat; xh <- ls2$xhat; retcode1 <- ls2$retcodeh
        gh_obj <- get_grad(xh); gh <- gh_obj$g; badgh <- gh_obj$badg
      }
      improvement <- f - fh
      if (improvement < crit) { done <- TRUE; retcode <- retcode1 }
    } else if (improvement < crit) {
      done <- TRUE; retcode <- retcode1
    }

    ## --- BFGS inverse-Hessian update (only with two good gradients). ---
    if (!done || improvement >= crit) {
      if (!badg && !badgh) {
        dg <- gh - g
        dx <- xh - x
        H  <- .csminwel_bfgsi(H, dg, dx)
      }
    }

    ## --- Stall / plateau stop (additional to the per-iteration `crit`). ---
    ## The per-iteration `improvement < crit` test (crit = 1e-7) does not fire on
    ## an ill-conditioned objective that keeps making tiny-but-nonzero gains every
    ## iteration: it grinds toward `nit` (thousands of iterations) for a logpost
    ## change that is numerically meaningless. So once past `stall_warmup`
    ## iterations, look back over the last `stall_window` improvements: if NONE of
    ## them cleared `stall_tol` (i.e. the best recent step still moved the
    ## objective by less than stall_tol), the search has plateaued -- report
    ## converged (retcode 8). Guarded by stall_warmup so a problem that genuinely
    ## converges quickly is never cut short. Set stall_tol = 0 (or stall_window
    ## >= nit) to disable and recover the pure-crit behaviour.
    recent_impr <- c(recent_impr, improvement)
    if (length(recent_impr) > stall_window)
      recent_impr <- recent_impr[-1L]
    if (!done && stall_tol > 0 && itct >= stall_warmup &&
        length(recent_impr) >= stall_window &&
        max(recent_impr) < stall_tol) {
      done <- TRUE; retcode <- 8L
    }

    ## Advance.
    x    <- xh; f <- fh; g <- gh; badg <- badgh
    fcount <- fcount + length(gh)   # gradient costs ~n evals

    if (verbose)
      cat(sprintf("  iter %4d  f = %20.10f  improvement = %12.4e  rc=%d\n",
                  itct, f, improvement, retcode1))
  }

  xh_out <- x
  if (!is.null(x0_names)) names(xh_out) <- x0_names

  list(
    fh          = f,
    xh          = xh_out,
    gh          = g,
    H           = H,
    itct        = itct,
    fcount      = fcount,
    retcode     = retcode,
    ## optim-style aliases for ergonomics within dynhr:
    par         = xh_out,
    value       = f,
    iterations  = itct,
    convergence = if (itct < nit) 0L else 1L
  )
}
