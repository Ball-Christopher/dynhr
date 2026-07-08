## R/pruned-grad-chain-order3.R
## --------------------------------------------------------------------------
## D1: semi-analytic order-3 pruned gradient (adjoint filter + FD-of-assembly).
##
## Full analytic order-3 gradient is out of scope (no solution_derivatives_
## order3 sibling with d(ghxxx)/dtheta etc., and the fold (A, Cu, stationary
## Sxi) needs a derivative-Lyapunov pass -- see the brief's scope report,
## reproduced in the completion report). This file implements the
## SEMI-ANALYTIC middle path instead:
##
##   d(loglik)/dtheta_j = sum over the 9 SSM inputs M of <bar_M, dM/dtheta_j>_F
##
## where bar_M are the input-gradients from ONE
## .pruned_kf_correlated_adjoint() pass (the order-3 filter is exactly the
## order-2 correlated-noise KF, just fed the fold-adjusted inputs
## Tlin_eff/ZZ/QQ/HH/SS/mu0/Sxi0 -- see .order3_pruned_kf_inputs(), the
## extract-assembly refactor in R/pruned-state-space-order3.R), and
## dM/dtheta_j comes from CENTRAL FD of the ASSEMBLY ONLY (steady state +
## order-3 solve + augmented system + fold + stationary moments -- NO filter
## pass, i.e. NOT a re-run of the O(T) Kalman recursion). Cost per gradient:
## 1 filter+adjoint pass + 2*n_params assembly-only evaluations, versus the
## FD-of-forward fallback's (n_params+1) * (assembly + O(T) filter). This
## wins whenever the filter's T-loop dominates the (fixed-cost) assembly.
##
## Convention: gradients w.r.t. every matrix input are dL/dX_ij treating all
## entries as independent (matches .pruned_kf_correlated_adjoint()'s
## convention exactly), so the two layers compose via a plain Frobenius inner
## product with no extra symmetrization Jacobian -- identical convention to
## R/pruned-grad-chain.R's order-2 chain.
## --------------------------------------------------------------------------


#' Assemble the order-3 pruned-SS KF inputs at a given theta (assembly only,
#' NO filter pass)
#'
#' Solves the model to order 3 at \code{theta}, builds the order-3 AFVRR
#' augmented system, and returns the nine \code{.pruned_kf_correlated}/adjoint
#' inputs via \code{.order3_pruned_kf_inputs()} -- the SAME assembly
#' \code{pruned_ss_loglik3()}/\code{make_log_posterior_pruned3()} use, so the
#' FD taken over this function's outputs is FD of exactly the production
#' assembly path, not a re-derivation.
#'
#' @param theta Named numeric parameter vector (prior-space).
#' @param model,compiled Model + compiled objects.
#' @param par_names Character vector, names of \code{theta}'s entries.
#' @param obs_vars Character vector of observed endogenous variable names.
#' @param me_variance Scalar measurement-error variance (default 0).
#' @return \code{NULL} on any solve/assembly failure, else a list with the
#'   nine KF inputs (\code{Tlin, ZZ, d_y, c_drift, QQ, HH, SS, mu0, Sxi0}).
#' @keywords internal
.pgo3_assemble <- function(theta, model, compiled, par_names, obs_vars,
                           me_variance = 0) {
  names(theta) <- par_names
  params <- .apply_theta_to_params(model, theta)
  ss_result <- tryCatch(
    solve_steady_state(model, compiled, params, verbose = FALSE),
    error = function(e) NULL)
  if (is.null(ss_result) || !isTRUE(ss_result$converged)) return(NULL)
  params <- ss_result$params %||% params

  dr3 <- tryCatch(
    solve_perturbation(model, compiled, ss_result$ss, params,
                       order = 3L, verbose = FALSE),
    error = function(e) NULL)
  if (is.null(dr3) || !isTRUE(dr3$bk_satisfied)) return(NULL)

  pss3 <- tryCatch(pruned_state_space3(dr3, model, params), error = function(e) NULL)
  if (is.null(pss3)) return(NULL)

  inp <- tryCatch(
    .order3_pruned_kf_inputs(pss3$sys, pss3$ys, obs_vars, pss3$endo_names,
                             me_variance = me_variance),
    error = function(e) NULL)
  if (is.null(inp)) return(NULL)

  list(Tlin = inp$Tlin, ZZ = inp$ZZ, d_y = inp$d_y, c_drift = inp$c_drift,
       QQ = inp$QQ, HH = inp$HH, SS = inp$SS, mu0 = inp$mu0, Sxi0 = inp$Sxi0)
}

## The nine KF-input names shared by .pgo3_assemble()'s return list and
## .pruned_kf_correlated_adjoint()'s grad list (single source of truth for
## the loop in .pgo3_grad_chain() below).
.PGO3_INPUT_NAMES <- c("Tlin", "ZZ", "d_y", "c_drift", "QQ", "HH", "SS",
                       "mu0", "Sxi0")

#' Central FD of d(assembly)/dtheta_j, assembly-only (no filter pass)
#'
#' @return \code{NULL} if either perturbed assembly fails (non-finite
#'   gradient contribution for this parameter -- caller falls back to
#'   FD-of-forward), else a list of the nine d(input)/dtheta_j matrices/
#'   vectors, same names/shapes as \code{.pgo3_assemble()}'s output.
#' @keywords internal
.pgo3_assembly_fd <- function(theta, nm, model, compiled, par_names, obs_vars,
                              me_variance = 0, h_rel = 1e-5) {
  h <- h_rel * max(abs(theta[[nm]]), 1e-3)
  tp <- theta; tp[[nm]] <- tp[[nm]] + h
  tm <- theta; tm[[nm]] <- tm[[nm]] - h

  ap <- .pgo3_assemble(tp, model, compiled, par_names, obs_vars, me_variance)
  am <- .pgo3_assemble(tm, model, compiled, par_names, obs_vars, me_variance)
  if (is.null(ap) || is.null(am)) return(NULL)

  d <- setNames(vector("list", length(.PGO3_INPUT_NAMES)), .PGO3_INPUT_NAMES)
  for (nmi in .PGO3_INPUT_NAMES) {
    d[[nmi]] <- (ap[[nmi]] - am[[nmi]]) / (2 * h)
  }
  d
}

#' End-to-end semi-analytic d(order-3 pruned loglik)/d(theta) for ALL params
#'
#' Runs ONE \code{.pruned_kf_correlated_adjoint()} pass at the base theta to
#' get \code{loglik} and the nine input-gradients \code{bar_M}, then for each
#' parameter takes central FD of the ASSEMBLY ONLY (\code{.pgo3_assembly_fd},
#' no filter pass) and contracts: \code{sum(bar_M * dM/dtheta_j)} summed over
#' the nine inputs.
#'
#' @param theta Named numeric parameter vector.
#' @param model,compiled Model + compiled objects.
#' @param par_names Character vector of parameter names to differentiate.
#' @param Y n_obs x T data matrix.
#' @param obs_vars Character vector of observed endogenous variable names.
#' @param me_variance Scalar measurement-error variance (default 0).
#' @return \code{NULL} if the base assembly fails (caller falls back to prior-
#'   only / FD-of-forward for everything), else a list with \code{loglik}
#'   (scalar) and \code{grad} (named numeric vector, one entry per
#'   \code{par_names}; \code{NA} for any parameter whose assembly FD failed
#'   in either direction).
#' @keywords internal
.pgo3_grad_chain <- function(theta, model, compiled, par_names, Y, obs_vars,
                             me_variance = 0) {
  names(theta) <- par_names
  base <- .pgo3_assemble(theta, model, compiled, par_names, obs_vars, me_variance)
  if (is.null(base)) return(NULL)

  ad <- tryCatch(
    .pruned_kf_correlated_adjoint(Y, base$Tlin, base$ZZ, base$d_y,
                                  base$c_drift, base$QQ, base$HH, base$SS,
                                  base$mu0, base$Sxi0),
    error = function(e) NULL)
  if (is.null(ad)) return(NULL)

  grad <- setNames(rep(NA_real_, length(par_names)), par_names)
  for (nm in par_names) {
    dM <- tryCatch(
      .pgo3_assembly_fd(theta, nm, model, compiled, par_names, obs_vars,
                        me_variance = me_variance),
      error = function(e) NULL)
    if (is.null(dM)) next
    gj <- 0
    ok <- TRUE
    for (nmi in .PGO3_INPUT_NAMES) {
      term <- sum(ad$grad[[nmi]] * dM[[nmi]])
      if (!is.finite(term)) { ok <- FALSE; break }
      gj <- gj + term
    }
    if (ok) grad[nm] <- gj
  }

  list(loglik = ad$loglik, grad = grad)
}
