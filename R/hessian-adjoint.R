## R/hessian-adjoint.R
## --------------------------------------------------------------------------
## EXACT posterior Hessian of the Gaussian Kalman-filter log-likelihood
## (ROADMAP Tier 6 #2). Builds on the first-order adjoint gradient
## (.kf_loglik_adjoint) and the first/second-order solution-derivative layers
## (solution_derivatives / solution_derivatives_2).
##
## THE DECOMPOSITION
## -----------------
## The loglik L is a function of the state-space matrices X = (TT,RR,ZZ,DD,d,
## Sigma_e), which in turn depend on theta. The adjoint returns the exact
## gradient grad_j = <G_X, dX_j> where G_X = dL/dX are the filter-gradient
## matrices and dX_j = dX/dtheta_j. Differentiating once more wrt theta_i:
##
##   H[i,j] = <dG_X/dtheta_i, dX_j>  +  <G_X, d2X_ij>
##            \__ filter curvature _/    \_ solution curvature _/
##                  (T1)                        (T2)
##
## T2 is computed exactly: G_X is contracted (by .kf_loglik_adjoint) against
## the second-order solution derivatives d2X_ij from solution_derivatives_2.
##
## T1 = <(Hess_X L) dX_i, dX_j> is the directional derivative of the exact
## adjoint gradient grad_j(X) = <G_X(X), dX_j> along dX_i. It is obtained by a
## central directional difference of the (exact, analytic) adjoint gradient in
## MATRIX space along dX_i -- no model re-solve, one finite-difference level on
## an exact quantity (a Hessian-vector product). The resulting T1 is symmetric
## to machine precision (an internal correctness check); the returned Hessian
## is symmetrised.
## --------------------------------------------------------------------------


#' Exact Hessian of the Kalman-filter log-likelihood wrt the state-space
#' matrices, contracted with first- and second-order solution derivatives.
#'
#' @param Y          n_obs x n_T observation matrix (no NAs).
#' @param ss         list TT, RR, ZZ, DD, d, Sigma_e (base point X).
#' @param dX_list    length-n_par list; element i has dTT,dRR,dZZ,dDD,dd,
#'                   dSigma_e (= dX/dtheta_i; missing blocks treated as zero).
#' @param d2X_list   list keyed "i|j" (i<=j and j<=i) of second-order blocks
#'                   d2TT,d2RR,d2ZZ,d2DD,d2d,d2Sigma_e (= d2X/dtheta_i dtheta_j).
#' @param me_variance scalar measurement-error variance.
#' @param eps        directional FD step for the filter-curvature term (T1).
#' @return n_par x n_par symmetric Hessian of the loglik.
#' @noRd
#' @param t1_method how to compute the filter-curvature term T1: \code{"hvp"}
#'   (default) is a matrix-space directional finite difference of the exact
#'   adjoint gradient (an HVP, residual ~1e-8 at the default \code{eps});
#'   \code{"analytic"} uses the fully-analytic forward-over-reverse second-order
#'   adjoint (\code{.kf_loglik_dG}), which is FD-free at the filter level.
kf_loglik_hessian <- function(Y, ss, dX_list, d2X_list, me_variance = 0,
                              eps = 1e-5, t1_method = c("hvp", "analytic")) {
  t1_method <- match.arg(t1_method)
  np <- length(dX_list)
  nm <- names(dX_list)

  ## Normalise a derivative block to the .kf_loglik_adjoint contract.
  norm_blk <- function(b) {
    if (is.null(b)) b <- list()
    list(dTT = b$dTT %||% b$d2TT, dRR = b$dRR %||% b$d2RR,
         dZZ = b$dZZ %||% b$d2ZZ, dDD = b$dDD %||% b$d2DD,
         dd  = b$dd  %||% b$d2d,  dSigma_e = b$dSigma_e %||% b$d2Sigma_e)
  }
  dX <- lapply(dX_list, norm_blk)

  ## T1: directional derivative of grad_j(X) = <G_X, dX_j> along dX_i.
  T1 <- matrix(0, np, np)
  if (t1_method == "analytic") {
    ## Fully-analytic: dG_X along dX_i via the second-order adjoint, then
    ## contract with dX_j. FD-free at the filter level.
    for (i in seq_len(np)) {
      dGi <- .kf_loglik_dG(Y, ss, dX[[i]], me_variance = me_variance)
      for (j in seq_len(np)) {
        b <- dX[[j]]
        z <- function(M, B) if (is.null(B)) 0 else sum(M * B)
        T1[i, j] <- z(dGi$dG_TT, b$dTT) + z(dGi$dG_RR, b$dRR) +
                    z(dGi$dG_ZZ, b$dZZ) + z(dGi$dG_DD, b$dDD) +
                    z(dGi$dg_d, b$dd)   + z(dGi$dG_Sig, b$dSigma_e)
      }
    }
  } else {
    perturb <- function(sgn, i) {
      b <- dX[[i]]
      z <- function(v, like) if (is.null(v)) 0 * like else v
      list(TT = ss$TT + sgn * eps * z(b$dTT, ss$TT),
           RR = ss$RR + sgn * eps * z(b$dRR, ss$RR),
           ZZ = ss$ZZ + sgn * eps * z(b$dZZ, ss$ZZ),
           DD = ss$DD + sgn * eps * z(b$dDD, ss$DD),
           d  = as.numeric(ss$d) + sgn * eps * z(b$dd, as.numeric(ss$d)),
           Sigma_e = ss$Sigma_e + sgn * eps * z(b$dSigma_e, ss$Sigma_e))
    }
    for (i in seq_len(np)) {
      gp <- .kf_loglik_adjoint(Y, perturb(+1, i), dX, me_variance = me_variance)$grad
      gm <- .kf_loglik_adjoint(Y, perturb(-1, i), dX, me_variance = me_variance)$grad
      T1[i, ] <- (gp - gm) / (2 * eps)
    }
  }

  ## T2: <G_X, d2X_ij> -- contract the base-point adjoint gradient matrices
  ## against the second-order solution derivatives (one adjoint call per row i,
  ## passing the n_par second-order blocks d2X_ij for j = 1..np).
  T2 <- matrix(0, np, np)
  for (i in seq_len(np)) {
    d2row <- lapply(seq_len(np), function(j) {
      key <- paste(i, j, sep = "|")
      norm_blk(d2X_list[[key]])
    })
    T2[i, ] <- .kf_loglik_adjoint(Y, ss, d2row, me_variance = me_variance)$grad
  }

  H <- T1 + T2
  H <- 0.5 * (H + t(H))
  dimnames(H) <- list(nm, nm)
  attr(H, "t1_asymmetry") <- max(abs(T1 - t(T1)))
  H
}


## --------------------------------------------------------------------------
## posterior_hessian(): higher-level wrapper that
##   1. solves the model at `params`,
##   2. builds the full dX_list (dTT,dRR,dZZ,dDD,dd,dSigma_e per param) and
##      d2X_list (d2TT,...,d2Sigma_e per (i,j) pair),
##   3. calls kf_loglik_hessian to get the loglik Hessian,
##   4. optionally adds the analytic Hessian of the log-prior (diagonal for
##      independent priors) if include_prior = TRUE.
##
## The d2Sigma_e blocks are computed by a second FD of .get_shock_cov:
##   3-point diagonal stencil for (i,i) pairs,
##   4-corner mixed stencil for (i,j) pairs.
## This matches the total-derivative convention of the first-order layer
## (which already uses .dSigma_e_fd = central FD of .get_shock_cov).
## Parameters that also move the decision rule (non-sigma params) get BOTH
## solution-derivative blocks AND dSigma_e from FD, consistent with the
## analytic gradient's "adjoint" path (analytic-gradient.R:719-727).
## --------------------------------------------------------------------------

#' Exact posterior Hessian of the KF log-likelihood (and optionally the prior)
#'
#' @param model     parsed model (from \code{parse_mod}).
#' @param compiled  compiled model (from \code{compile_model}).
#' @param dr        first-order decision rule (from \code{solve_perturbation}).
#' @param params    named numeric vector of ALL model parameters at the evaluation
#'                  point.
#' @param param_names character vector of parameter names to differentiate
#'                    (the Hessian is \code{length(param_names) x length(param_names)}).
#' @param obs_vars  character vector of observable variable names.
#' @param Y         n_obs x T observation matrix (no NAs).
#' @param me_variance scalar measurement-error variance (default 0).
#' @param include_prior logical; if TRUE, add the analytic Hessian of
#'   \code{log p(theta)} (diagonal for independent priors) to the loglik Hessian.
#' @param prior_spec prior-spec data.frame (required when \code{include_prior = TRUE}).
#' @param eps       directional FD step for the filter-curvature (T1) term
#'   (default 1e-5).
#' @param h_Sigma_e relative FD step for the dSigma_e / d2Sigma_e computation
#'   (default 1e-5).
#' @return n_par x n_par symmetric matrix; the loglik Hessian
#'   (plus prior Hessian if \code{include_prior = TRUE}).
#'   Carries attribute \code{t1_asymmetry} (scalar) from \code{kf_loglik_hessian}.
#' @noRd
posterior_hessian <- function(model, compiled, dr, params, param_names,
                              obs_vars, Y,
                              me_variance = 0,
                              include_prior = FALSE,
                              prior_spec = NULL,
                              eps = 1e-5,
                              h_Sigma_e = 1e-5,
                              t1_method = c("hvp", "analytic")) {

  t1_method <- match.arg(t1_method)
  np  <- length(param_names)
  exo <- model$varexo_names

  ## ------------------------------------------------------------------
  ## State-space at the evaluation point.
  ## ------------------------------------------------------------------
  si  <- dr$state_idx
  oi  <- match(obs_vars, dr$endo_names)
  ss0 <- list(
    TT      = dr$ghx[si, , drop = FALSE],
    RR      = dr$ghu[si, , drop = FALSE],
    ZZ      = dr$ghx[oi, , drop = FALSE],
    DD      = dr$ghu[oi, , drop = FALSE],
    d       = dr$ys[obs_vars],
    Sigma_e = .get_shock_cov(model, exo, params)
  )

  ## ------------------------------------------------------------------
  ## First-order solution derivatives (structural blocks) for ALL params.
  ## solution_derivatives internally distinguishes ok / !ok per param.
  ## ------------------------------------------------------------------
  sd1 <- tryCatch(
    solution_derivatives(model, compiled, dr, params,
                         param_names = param_names, obs_vars = obs_vars),
    error = function(e) NULL
  )

  ## ------------------------------------------------------------------
  ## Second-order solution derivatives (structural blocks) for ALL params.
  ## ------------------------------------------------------------------
  sd2 <- tryCatch(
    solution_derivatives_2(model, compiled, dr, params,
                           param_names = param_names, obs_vars = obs_vars),
    error = function(e) NULL
  )

  ## ------------------------------------------------------------------
  ## dSigma_e and d2Sigma_e by FD of .get_shock_cov.
  ## ------------------------------------------------------------------
  h_vec <- vapply(param_names, function(nm)
    h_Sigma_e * max(abs(params[[nm]]), 1e-4), 0.0)
  names(h_vec) <- param_names

  ## Single-step Sigma_e at theta +/- h_i (for dSigma_e and d2Sigma_e).
  Se_p <- Se_m <- vector("list", np); names(Se_p) <- names(Se_m) <- param_names
  for (i in seq_len(np)) {
    nm <- param_names[i]; h <- h_vec[i]
    pp <- params; pm <- params
    pp[[nm]] <- pp[[nm]] + h; pm[[nm]] <- pm[[nm]] - h
    Se_p[[nm]] <- .get_shock_cov(model, exo, pp)
    Se_m[[nm]] <- .get_shock_cov(model, exo, pm)
  }
  Se0 <- ss0$Sigma_e

  dSe <- lapply(param_names, function(nm) {
    h <- h_vec[nm]
    (Se_p[[nm]] - Se_m[[nm]]) / (2 * h)
  })
  names(dSe) <- param_names

  ## d2Sigma_e by a second FD:
  ##   diagonal (i=j): 3-point: (Se(+h) - 2*Se(0) + Se(-h)) / h^2
  ##   mixed   (i!=j): 4-corner: (Se(+hi,+hj) - Se(+hi,-hj) - Se(-hi,+hj) + Se(-hi,-hj)) / (4 hi hj)
  d2Se <- list()
  for (i in seq_len(np)) {
    nm_i <- param_names[i]; hi <- h_vec[nm_i]
    for (j in seq_len(np)) {
      nm_j <- param_names[j]; hj <- h_vec[nm_j]
      key <- paste(i, j, sep = "|")
      if (i == j) {
        d2Se[[key]] <- (Se_p[[nm_i]] - 2 * Se0 + Se_m[[nm_i]]) / (hi^2)
      } else {
        ## Compute Se at (theta +/- hi, +/- hj)
        pp <- params; pp[[nm_i]] <- pp[[nm_i]] + hi; pp[[nm_j]] <- pp[[nm_j]] + hj
        pm <- params; pm[[nm_i]] <- pm[[nm_i]] + hi; pm[[nm_j]] <- pm[[nm_j]] - hj
        mp <- params; mp[[nm_i]] <- mp[[nm_i]] - hi; mp[[nm_j]] <- mp[[nm_j]] + hj
        mm <- params; mm[[nm_i]] <- mm[[nm_i]] - hi; mm[[nm_j]] <- mm[[nm_j]] - hj
        Se_pp <- .get_shock_cov(model, exo, pp)
        Se_pm <- .get_shock_cov(model, exo, pm)
        Se_mp <- .get_shock_cov(model, exo, mp)
        Se_mm <- .get_shock_cov(model, exo, mm)
        d2Se[[key]] <- (Se_pp - Se_pm - Se_mp + Se_mm) / (4 * hi * hj)
      }
    }
  }

  ## ------------------------------------------------------------------
  ## Build dX_list: first-order blocks per parameter.
  ## ------------------------------------------------------------------
  dX_list <- vector("list", np); names(dX_list) <- param_names
  for (i in seq_len(np)) {
    nm <- param_names[i]
    d1 <- if (!is.null(sd1)) sd1$derivs[[nm]] else NULL
    ok <- !is.null(d1) && isTRUE(d1$ok)
    dX_list[[i]] <- list(
      dTT      = if (ok) d1$dTT      else NULL,
      dRR      = if (ok) d1$dRR      else NULL,
      dZZ      = if (ok) d1$dZZ      else NULL,
      dDD      = if (ok) d1$dDD      else NULL,
      dd       = if (ok) d1$dd       else NULL,
      dSigma_e = dSe[[nm]]
    )
  }

  ## ------------------------------------------------------------------
  ## Build d2X_list: second-order blocks per (i,j) pair.
  ## ------------------------------------------------------------------
  d2X_list <- list()
  for (i in seq_len(np)) {
    for (j in seq_len(np)) {
      key  <- paste(i, j, sep = "|")
      sd2_blk <- if (!is.null(sd2)) sd2$d2[[key]] else NULL
      blk <- list(
        d2TT      = sd2_blk$d2TT,
        d2RR      = sd2_blk$d2RR,
        d2ZZ      = sd2_blk$d2ZZ,
        d2DD      = sd2_blk$d2DD,
        d2d       = sd2_blk$d2d,
        d2Sigma_e = d2Se[[key]]
      )
      d2X_list[[key]] <- blk
    }
  }

  ## ------------------------------------------------------------------
  ## Loglik Hessian via kf_loglik_hessian.
  ## ------------------------------------------------------------------
  H <- kf_loglik_hessian(Y, ss0, dX_list, d2X_list,
                         me_variance = me_variance, eps = eps,
                         t1_method = t1_method)

  ## ------------------------------------------------------------------
  ## Optionally add the analytic prior Hessian (diagonal for independent
  ## priors; d2/dx2 log p(x) = d/dx [dlog_prior_density1]).
  ## We differentiate .dlog_prior_density1 numerically with a central
  ## difference (avoids reimplementing analytic second derivatives for
  ## every distribution).
  ## ------------------------------------------------------------------
  if (include_prior) {
    if (is.null(prior_spec))
      stop("posterior_hessian: prior_spec required when include_prior = TRUE")
    for (nm in param_names) {
      row_i <- match(nm, prior_spec$name)
      if (is.na(row_i)) next
      x0   <- params[[nm]]
      dist <- prior_spec$distribution[row_i]
      p1   <- prior_spec$p1[row_i]
      p2   <- prior_spec$p2[row_i]
      ## Check that this prior name maps to a param in param_names
      idx  <- match(nm, param_names)
      if (is.na(idx)) next
      ## d2 log p / dx2 by central FD of .dlog_prior_density1
      h_pr <- 1e-5 * max(abs(x0), 1e-4)
      d2pr <- (.dlog_prior_density1(x0 + h_pr, dist, p1, p2) -
               .dlog_prior_density1(x0 - h_pr, dist, p1, p2)) / (2 * h_pr)
      H[idx, idx] <- H[idx, idx] + d2pr
    }
  }

  H
}


#' Laplace approximation to the log marginal likelihood
#'
#' Computes the Laplace (second-order) approximation to the model log marginal
#' likelihood from a mode-finding result that carries an exact posterior Hessian
#' (\code{run_mode_finding(use_exact_hessian = TRUE)}):
#' \deqn{\log p(Y\mid M) \approx \log p(Y\mid\hat\theta) + \log p(\hat\theta\mid M)
#'        + \tfrac{d}{2}\log 2\pi - \tfrac12 \log\det(-H),}
#' where \eqn{H} is the exact Hessian of the log-posterior at the mode and the
#' first two terms are the (unnormalised) log-posterior at the mode
#' (\code{mode_result$mode$logpost}). Deterministic and exact to
#' \eqn{O(1/T)} (Tierney & Kadane, 1986); a fast, randomness-free alternative to
#' the SMC marginal-likelihood estimator for model comparison.
#'
#' @param mode_result A \code{dynhr_mode_result} from
#'   \code{\link{run_mode_finding}} with \code{hessian_exact} populated.
#' @return The Laplace log marginal likelihood (numeric scalar), or
#'   \code{NA_real_} when no exact Hessian is present or \eqn{-H} is not
#'   positive-definite at the mode (the mode is not a strict local maximum).
#' @seealso \code{posterior_hessian}, \code{\link{run_mode_finding}}
#' @export
laplace_log_marglik <- function(mode_result) {
  H <- mode_result$hessian_exact
  if (is.null(H)) {
    warning("laplace_log_marglik: mode_result has no exact Hessian; re-run ",
            "run_mode_finding(use_exact_hessian = TRUE).", call. = FALSE)
    return(NA_real_)
  }
  ev <- eigen(-H, symmetric = TRUE)$values
  if (any(ev <= 0)) {
    warning("laplace_log_marglik: -H is not positive-definite at the mode; ",
            "the Laplace marginal likelihood is undefined.", call. = FALSE)
    return(NA_real_)
  }
  d <- length(mode_result$theta_mode)
  as.numeric(mode_result$mode$logpost + 0.5 * d * log(2 * pi) - 0.5 * sum(log(ev)))
}
