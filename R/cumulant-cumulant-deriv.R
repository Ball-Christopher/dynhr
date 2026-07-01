## R/cumulant-cumulant-deriv.R
## --------------------------------------------------------------------------
## Sensitivity of the order-3 and order-4 cumulant moment vectors to
## structural parameters, for use in the implicit-differentiation gradient
## of the cumulant log-likelihood (R/cumulant-gradient.R).
##
## Entry point: cumulant_moment_derivs_3_4()
##
## Given the base decision-rule quantities, this file computes:
##   d_c3_obs[k] — derivative of the n_obs x n_obs^2 cumulant matrix wrt θ_k
##   d_c4_obs[k] — derivative of the n_obs x n_endo^3 kurtosis matrix wrt θ_k
##
## STATUS (2026-07-01, updated)
## --------------------
## Both order-3 and order-4 derivatives are central differences of the forward
## functions (compute_third_cumulant() / compute_fourth_cumulant()), but the
## perturbed decision rule is now obtained by an ANALYTIC-DIRECTION step
## (.dr_along_direction): the order-1/2 solution derivatives from
## solution_derivatives_order2() (dG/dH/dys/d_ghxx/d_ghxu/d_ghuu/d_ghss) are
## used to shift ghx/ghu/ghxx/.../ys by eps*(d/dtheta) WITHOUT a perturbation
## re-solve.  Because the closed-form cumulants are smooth in those inputs, the
## directional derivative along the analytic direction equals d(cumulant)/dtheta
## exactly (validated to ~1e-6 vs the legacy re-solving FD; see
## test-cumulant-grad-directional.R).  This is more accurate than perturbing
## theta and re-solving, and avoids the per-parameter steady-state + order-2
## perturbation solve (a real saving when re-solving dominates, e.g. larger
## estimated models; a wash for tiny models where the forward cumulant
## dominates).  If the solution derivatives are unavailable the helpers fall
## back to the legacy re-solving FD via .resolve_order2().
##
## A fully closed-form analytic sensitivity of the CHAIN term (differentiating
## the C211/C2211 Lyapunov tensors -- d_C211/d_C2211 are derived+validated in
## .claude/orchestration/cumulant-oracle/{dC211_dev,dC2211_dev}.R -- and the
## chain assembly) was scoped and DELIBERATELY NOT IMPLEMENTED: profiling
## (2026-07-01) shows the order-4 forward cost is 99.8% the trace q-form and
## only 0.2% the closed-form chain, so a symbolic chain derivative removes ~0%
## of the gradient runtime while adding a large, error-prone tensor
## differentiation (the recurring source of self-consistent-but-wrong gradient
## bugs in this module). The trace term dominates and its truncated-window
## tr((MS)^4) quadrature has no cheap symbolic derivative; the only real speedup
## lever is a closed-form trace (needs the full Andreasen-FV-RR 4th-moment
## hierarchy / C2222 -- a standalone research build, deferred: the trace is
## ~1.3% of kappa4 and its truncation error is inside the MC4 gate). Order-3 is
## already pure closed-form (~2.5 ms) with no trace. Net: the directional FD
## path below is the right design -- correct, and re-solve-avoiding.
##
## LEGACY NOTE: previously both orders re-solved at theta +- h via
## .resolve_order2() (defined in R/cumulant-gradient.R).
##
## History / why FD and not a closed-form tensor-Lyapunov sensitivity:
##  - An earlier version of this file differentiated compute_third_cumulant()'s
##    OLD (pre-fix) Lyapunov RHS `2*Sigma_x*H*Sigma_x + tr(H*Sigma_x)*Sigma_x`
##    (hxx family only) and an OLD single-permutation observation projection.
##    When compute_third_cumulant() was corrected (commit 809cd8d: hx*Sigma_x
##    feedthrough, all 3 permutations, plus the hxu/huu driving families; see
##    R/cumulant-likelihood.R's .third_cumulant_rhs()), that analytic
##    derivative became stale and test-cumulant-grad-order34.R's direct
##    d_c3_obs-vs-FD check failed with rel_err up to ~2.8.
##  - compute_fourth_cumulant() was separately rewritten from a crude
##    O(sigma^2)~0 ratio heuristic (kurtosis_i = 3*(Var_x2/Var_total)^2) to the
##    EXACT marginal excess kurtosis via a linear-plus-quadratic-form /
##    generalized-chi-square construction in the (truncated) Gaussian
##    innovation history (.fourth_cumulant_qform_marginal(); see
##    R/cumulant-likelihood.R). That representation involves an O(n_lag) sum
##    of Kronecker-structured quadratic forms; differentiating it in closed
##    form is a substantial undertaking on its own.
##
## Rather than maintain a second, easily-divergent analytic differentiation
## of each forward formula, both orders now use central finite differences of
## the (corrected, exact) forward functions directly — this keeps the
## gradient automatically consistent with whatever compute_third_cumulant() /
## compute_fourth_cumulant() compute, at the cost of two extra steady-state +
## order-2-perturbation solves per parameter per order (4 extra solves total
## for orders 3+4 combined). Orders 1-2 (mean, variance) remain analytic via
## solution_derivatives_order2() / implicit differentiation in
## R/cumulant-gradient.R.
##
## Reference: Mutschler (2015), Sections 3.1-3.2; agent-F-tensor-lyapunov.md
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## Helper: central-FD derivative of compute_third_cumulant()$c3_obs wrt a
## single parameter pnm, evaluated at theta +- h (one-sided if only one
## perturbed solve succeeds, zero if both fail).
##
## @param model,compiled  dynhr_mod / dynhr_compiled, needed to re-solve at
##   theta +- h.
## @param dr2         DecisionRules2 at the base point (used only for dims
##   and as the .resolve_order2() template).
## @param params      full named parameter vector at the base point
## @param pnm         name of the parameter being differentiated
## @param h           absolute FD step for this parameter
## @param c3_obs_base base-point c3_obs (n_endo x n_endo^2), used as the
##   one-sided fallback if a perturbed solve fails.
## @return n_endo x n_endo^2 matrix: d(c3_obs)/d(theta_pnm)
## ---------------------------------------------------------------------------
## Perturb a DecisionRules2 along the ANALYTIC solution-derivative direction.
## Returns dr2 with ghx/ghu/ghxx/ghxu/ghuu/ghss/ys shifted by eps*(d/dtheta),
## using the order-1 (o1: dG,dH,dys) and order-2 (o2: d_ghxx,d_ghxu,d_ghuu,
## d_ghss) solution derivatives from solution_derivatives_order2().  No
## perturbation re-solve: the directional derivative of the (smooth) closed-form
## cumulant along this analytic direction equals d(cumulant)/dtheta exactly
## (validated to ~1e-6 vs the re-solving FD; faster and more accurate, since it
## reuses the analytic order-1/2 solution derivatives instead of an extra solve).
## @noRd
.dr_along_direction <- function(dr2, o1, o2, eps) {
  d <- dr2
  if (!is.null(o1$dG))      d$ghx  <- dr2$ghx  + eps * o1$dG
  if (!is.null(o1$dH))      d$ghu  <- dr2$ghu  + eps * o1$dH
  if (!is.null(o2$d_ghxx))  d$ghxx <- dr2$ghxx + eps * o2$d_ghxx
  if (!is.null(o2$d_ghxu))  d$ghxu <- dr2$ghxu + eps * o2$d_ghxu
  if (!is.null(o2$d_ghuu))  d$ghuu <- dr2$ghuu + eps * o2$d_ghuu
  if (!is.null(o2$d_ghss))  d$ghss <- dr2$ghss + eps * o2$d_ghss
  if (!is.null(o1$dys) && length(o1$dys) == length(dr2$ys))
    d$ys <- dr2$ys + eps * o1$dys
  d
}

.d_cumulant_order3 <- function(model, compiled, dr2, params, pnm, h,
                                c3_obs_base, o1 = NULL, o2 = NULL) {
  n_endo <- length(dr2$endo_names)

  params_p <- params; params_p[[pnm]] <- params[[pnm]] + h
  params_m <- params; params_m[[pnm]] <- params[[pnm]] - h

  ## Analytic-direction directional derivative (no re-solve) when the solution
  ## derivatives are available; otherwise fall back to the re-solving FD.
  if (!is.null(o1) && !is.null(o2)) {
    dr_p <- .dr_along_direction(dr2, o1, o2,  h)
    dr_m <- .dr_along_direction(dr2, o1, o2, -h)
  } else {
    dr_p <- .resolve_order2(model, compiled, params_p, dr2)
    dr_m <- .resolve_order2(model, compiled, params_m, dr2)
  }

  c3_p <- if (!is.null(dr_p)) {
    tryCatch(compute_third_cumulant(dr_p, model, params_p)$c3_obs,
             error = function(e) NULL)
  } else NULL
  c3_m <- if (!is.null(dr_m)) {
    tryCatch(compute_third_cumulant(dr_m, model, params_m)$c3_obs,
             error = function(e) NULL)
  } else NULL

  if (!is.null(c3_p) && !is.null(c3_m)) return((c3_p - c3_m) / (2 * h))
  if (!is.null(c3_p)) return((c3_p - c3_obs_base) / h)
  if (!is.null(c3_m)) return((c3_obs_base - c3_m) / h)
  matrix(0, n_endo, n_endo * n_endo)
}


## ---------------------------------------------------------------------------
## Helper: central-FD derivative of compute_fourth_cumulant()$c4_obs wrt a
## single parameter pnm (same pattern as .d_cumulant_order3 above).
##
## @param c4_obs_base  base-point c4_obs (n_endo x n_endo^3), used as the
##   one-sided fallback if a perturbed solve fails.
## @return n_endo x n_endo^3 matrix: d(c4_obs)/d(theta_pnm)
## ---------------------------------------------------------------------------
.d_cumulant_order4 <- function(model, compiled, dr2, params, pnm, h,
                                c4_obs_base, o1 = NULL, o2 = NULL) {
  n_endo <- length(dr2$endo_names)

  params_p <- params; params_p[[pnm]] <- params[[pnm]] + h
  params_m <- params; params_m[[pnm]] <- params[[pnm]] - h

  if (!is.null(o1) && !is.null(o2)) {
    dr_p <- .dr_along_direction(dr2, o1, o2,  h)
    dr_m <- .dr_along_direction(dr2, o1, o2, -h)
  } else {
    dr_p <- .resolve_order2(model, compiled, params_p, dr2)
    dr_m <- .resolve_order2(model, compiled, params_m, dr2)
  }

  c4_p <- if (!is.null(dr_p)) {
    tryCatch(compute_fourth_cumulant(dr_p, model, params_p)$c4_obs,
             error = function(e) NULL)
  } else NULL
  c4_m <- if (!is.null(dr_m)) {
    tryCatch(compute_fourth_cumulant(dr_m, model, params_m)$c4_obs,
             error = function(e) NULL)
  } else NULL

  if (!is.null(c4_p) && !is.null(c4_m)) {
    return((c4_p - c4_m) / (2 * h))
  }
  if (!is.null(c4_p)) {
    return((c4_p - c4_obs_base) / h)
  }
  if (!is.null(c4_m)) {
    return((c4_obs_base - c4_m) / h)
  }
  matrix(0, n_endo, n_endo * n_endo * n_endo)
}


#' Finite-difference derivatives of the order-3 and order-4 cumulant moment
#' blocks
#'
#' Computes d_c3_obs and d_c4_obs per parameter via central finite differences
#' of \code{compute_third_cumulant()} / \code{compute_fourth_cumulant()}. This
#' is called from \code{.cumulant_loglik_grad_implicit()} so that orders 3-4
#' track whatever those forward functions compute (see the module-level
#' "STATUS" note above for why this isn't a closed-form analytic derivative).
#'
#' @param dr2        DecisionRules2 at the base point
#' @param model      dynhr_mod
#' @param params     Named numeric parameter vector (base point)
#' @param o2d        Output of \code{solution_derivatives_order2()} (only
#'   \code{o2d$param_names} is used, to keep parameter ordering consistent
#'   with the orders 1-2 path)
#' @param obs_idx    Integer: indices of obs_vars in dr2$endo_names
#' @param obs_vars   Character: observed variable names
#' @param Sigma_x    n_s × n_s base state covariance (unused directly here,
#'   kept for interface compatibility / potential future analytic terms)
#' @param Sigma_y_full  n_endo × n_endo base total variance (unused directly)
#' @param c3_base    Output of \code{compute_third_cumulant()} at base
#' @param c4_base    Output of \code{compute_fourth_cumulant()} at base
#' @param compiled   dynhr_compiled, needed to re-solve at theta +- h for the
#'   finite-difference orders 3-4. If \code{NULL}, both derivatives are
#'   returned as all-zero (with a warning) rather than failing.
#' @param h_rel      Relative FD step size (default 1e-4).
#'
#' @return Named list (by param_names) of lists:
#'   \item{d_c3_obs_sub}{n_obs × n_obs^2 derivative of c3 moment block}
#'   \item{d_c4_obs_sub}{n_obs × n_endo^3 derivative of c4 moment block}
#'   \item{ok}{logical}
#' @noRd
cumulant_moment_derivs_3_4 <- function(dr2, model, params, o2d, obs_idx,
                                        obs_vars, Sigma_x, Sigma_y_full,
                                        c3_base, c4_base,
                                        compiled = NULL, h_rel = 1e-4) {
  param_names <- o2d$param_names
  np          <- length(param_names)

  ## Orders 3-4 are central finite differences of the (corrected) forward
  ## compute_third_cumulant()/compute_fourth_cumulant(), so they need a compiled
  ## model to re-solve at theta +- h.  If the caller did not supply one,
  ## recompile from `model` (cheap; max_order = 2L) so the FD path still works
  ## instead of degrading to a zero derivative.
  if (is.null(compiled)) {
    compiled <- tryCatch(
      compile_model(model, verbose = FALSE, max_order = 2L),
      error = function(e) NULL)
  }

  endo_names  <- dr2$endo_names
  n_endo      <- length(endo_names)
  n_obs       <- length(obs_vars)

  c3_obs_base <- c3_base$c3_obs   # n_endo × n_endo^2
  c4_obs_base <- c4_base$c4_obs   # n_endo × n_endo^3

  result <- vector("list", np)
  names(result) <- param_names

  for (k in seq_len(np)) {
    pnm  <- param_names[k]
    pval <- params[[pnm]]
    h_k  <- h_rel * (abs(pval) + 1e-8)

    if (is.null(compiled)) {
      warning("cumulant_moment_derivs_3_4: 'compiled' not supplied; orders ",
              "3-4 finite-difference derivatives set to zero for ", pnm,
              call. = FALSE)
      result[[pnm]] <- list(
        d_c3_obs_sub = matrix(0, n_obs, n_obs * n_obs),
        d_c4_obs_sub = matrix(0, n_obs, n_endo * n_endo * n_endo),
        ok           = TRUE
      )
      next
    }

    ## Per-parameter analytic solution derivatives (order-1 dG/dH/dys + order-2
    ## d_ghxx/...) enabling the no-re-solve directional derivative.  NULL -> the
    ## helpers fall back to the re-solving FD.
    o1_k <- tryCatch(o2d$first$derivs[[pnm]], error = function(e) NULL)
    o2_k <- tryCatch(o2d$derivs[[pnm]],       error = function(e) NULL)
    if (!is.null(o1_k) && !isTRUE(o1_k$ok)) o1_k <- NULL
    if (!is.null(o2_k) && !isTRUE(o2_k$ok)) o2_k <- NULL

    ## ---- Order 3 ----
    d_c3_obs_full <- tryCatch(
      .d_cumulant_order3(model, compiled, dr2, params, pnm, h_k, c3_obs_base,
                         o1 = o1_k, o2 = o2_k),
      error = function(e) {
        warning(sprintf("cumulant_moment_derivs_3_4: order-3 deriv failed for %s: %s",
                        pnm, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(d_c3_obs_full)) {
      result[[pnm]] <- list(
        d_c3_obs_sub = matrix(NA_real_, n_obs, n_obs * n_obs),
        d_c4_obs_sub = matrix(NA_real_, n_obs, n_endo * n_endo * n_endo),
        ok           = FALSE
      )
      next
    }

    ## Subset d_c3_obs to observable rows/cols (mirrors .cumulant_loglik)
    ## c3_obs_only[a, dst_col] <- c3_model_raw[a, src_col]
    ##   where src_col = (obs_idx[a]-1)*n_endo + obs_idx[b], dst_col = (a-1)*n_obs + b
    d_c3_obs_raw <- d_c3_obs_full[obs_idx, , drop = FALSE]   # n_obs × n_endo^2
    d_c3_obs_sub <- matrix(0, n_obs, n_obs * n_obs)
    for (a in seq_len(n_obs)) {
      for (b in seq_len(n_obs)) {
        src_col <- (obs_idx[a] - 1L) * n_endo + obs_idx[b]
        dst_col <- (a - 1L) * n_obs + b
        d_c3_obs_sub[a, dst_col] <- d_c3_obs_raw[a, src_col]
      }
    }

    ## ---- Order 4 ----
    d_c4_obs_full <- tryCatch(
      .d_cumulant_order4(model, compiled, dr2, params, pnm, h_k, c4_obs_base,
                         o1 = o1_k, o2 = o2_k),
      error = function(e) {
        warning(sprintf("cumulant_moment_derivs_3_4: order-4 deriv failed for %s: %s",
                        pnm, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(d_c4_obs_full)) {
      d_c4_obs_sub <- matrix(0, n_obs, n_endo * n_endo * n_endo)
    } else {
      ## Row-subset only (mirrors .cumulant_loglik / .build_moment_vector)
      d_c4_obs_sub <- d_c4_obs_full[obs_idx, , drop = FALSE]   # n_obs × n_endo^3
    }

    result[[pnm]] <- list(
      d_c3_obs_sub = d_c3_obs_sub,
      d_c4_obs_sub = d_c4_obs_sub,
      ok           = TRUE
    )
  }

  result
}
