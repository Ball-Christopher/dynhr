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


## ===========================================================================
## Order-3-semi: reverse (adjoint) of the third-cumulant tensor-Lyapunov solve
## ---------------------------------------------------------------------------
## Forward (.solve_third_cross_cumulant): given hx (n_s x n_s) and rhs
## (n_s x n_s^2), solve the tensor-Lyapunov
##
##     (I - hx^{⊗3}) vec(C3) = 1/2 vec(rhs)                              (*)
##
## with the (j,k)-symmetrization dropped here (the adjoint targets the core
## linear solve; the caller symmetrizes rhs upstream as in the forward).
##
## REVERSE: given a cotangent bar_C3 (n_s x n_s^2) = dL/dC3, produce
##   (a) bar_rhs   (n_s x n_s^2) = dL/d(rhs)
##   (b) bar_hx    (n_s x n_s)   = dL/d(hx)
## via ONE transposed tensor-Lyapunov solve
##
##     (I - hx^{⊗3})' vec(M) = vec(bar_C3),                              (**)
##
## then bar_rhs = 1/2 M, and bar_hx is the sum of the three tensor-mode
## contractions of  d[(hx^{⊗3}) vec(C3)]/d hx  applied with cotangent M.
##
## For the mode-wise hx derivative: with C3 as the 3-tensor T[i,j,k]
## (i = mode-1 row, (j,k) flattened col), the operator applies hx to each
## mode: (hx^{⊗3} T)[a,b,c] = sum_{i,j,k} hx[a,i] hx[b,j] hx[c,k] T[i,j,k].
## Reverse-accumulating the cotangent M[a,(b,c)] into hx (three modes):
##   mode1: bar_hx[a,i] += sum_{j,k,b,c} M[a,(b,c)] hx[b,j] hx[c,k] T[i,(j,k)]
##   mode2: bar_hx[b,j] += sum_{i,k,a,c} M[a,(b,c)] hx[a,i] hx[c,k] T[i,(j,k)]
##   mode3: bar_hx[c,k] += sum_{i,j,a,b} M[a,(b,c)] hx[a,i] hx[b,j] T[i,(j,k)]
## Because the forward RHS is (**) with the SAME operator (I - hx^{⊗3}), the
## bar_hx contribution from (*) is the NEGATIVE of the mode contractions of
## the operator hx^{⊗3} applied to the SOLUTION C3 (the "I" part has no hx).
## ---------------------------------------------------------------------------

#' Apply (I - hx^{⊗3})' to a tensor given as an n_s x n_s^2 matrix.
#' Uses (hx^{⊗3})' = (hx')^{⊗3}: transpose each mode, i.e. apply t(hx) to
#' mode 1 (left-multiply) and t(hx) ⊗ t(hx) to modes 2,3 (via .apply_kron2).
#' @noRd
.tensor_lyap3_apply_T <- function(hx, X) {
  ## (hx^{⊗3}) X : mode1 hx %*% X ; modes 2,3 via .apply_kron2(hx, .)
  HX <- hx %*% X
  HX <- .apply_kron2(hx, HX)
  X - HX          # (I - hx^{⊗3}) X  ; for transpose pass hx = t(hx0)
}

#' Reverse-mode adjoint of .solve_third_cross_cumulant (core linear solve).
#'
#' @param hx      n_s x n_s state transition (base point).
#' @param c3      n_s x n_s^2 forward solution C3 (from
#'   .solve_third_cross_cumulant BEFORE the final (j,k)-symmetrization, or the
#'   symmetric one — the transposed solve is agnostic; bar_hx uses c3 as the
#'   operand tensor).
#' @param bar_c3  n_s x n_s^2 cotangent dL/dC3.
#' @return list(bar_rhs = n_s x n_s^2, bar_hx = n_s x n_s).
#' @noRd
.solve_third_cross_cumulant_adjoint <- function(hx, c3, bar_c3) {
  n_s <- nrow(hx)
  if (n_s == 0L)
    return(list(bar_rhs = matrix(0, 0, 0), bar_hx = matrix(0, 0, 0)))

  ## The forward symmetrizes C3 in (j,k) AFTER the solve: C3 = Sym(solve).
  ## Sym is self-adjoint, so the cotangent on the raw solve is Sym(bar_c3).
  for (i in seq_len(n_s)) {
    B <- matrix(bar_c3[i, ], n_s, n_s)
    bar_c3[i, ] <- as.numeric((B + t(B)) * 0.5)
  }

  ## (a) transposed tensor-Lyapunov solve (**): (I - hx^{⊗3})' M = bar_c3
  ## Eigen-solve in the eigenbasis of hx' (== conj-eigenbasis of hx), mirroring
  ## the forward denominator 1 - lam_i lam_j lam_k.
  eigT <- eigen(t(hx))
  V   <- eigT$vectors
  lam <- eigT$values
  Vi  <- solve(V)
  b_tfm <- Vi %*% bar_c3
  b_tfm <- .apply_kron2(Vi, b_tfm)
  M_tfm <- matrix(0, n_s, n_s * n_s)
  for (i in seq_len(n_s)) {
    row_b <- matrix(b_tfm[i, ], n_s, n_s)
    for (j in seq_len(n_s)) for (k in seq_len(n_s)) {
      denom <- 1 - lam[i] * lam[j] * lam[k]
      row_b[j, k] <- if (abs(denom) > 1e-14) row_b[j, k] / denom else 0
    }
    M_tfm[i, ] <- as.numeric(row_b)
  }
  M <- V %*% M_tfm
  M <- .apply_kron2(V, M)
  M <- Re(M)

  ## (b) bar_rhs = 1/2 M  (rhs enters (*) as 1/2 vec(rhs))
  bar_rhs <- 0.5 * M

  ## (c) bar_hx: the operator side (*) is (I - hx^{⊗3}) C3 = 1/2 rhs. With M the
  ## adjoint state, bar on the residual r(hx,C3) = (I - hx^{⊗3})C3 is (-M) for
  ## the implicit-function reverse (bar_hx = - d r/d hx contracted with M via
  ## the standard adjoint sign for A(hx) C3 = b => bar_hx = -M (dA/dhx) C3).
  ## dA/dhx acts only through -hx^{⊗3}; the three modes each contribute.
  bar_hx <- matrix(0, n_s, n_s)
  ## Reshape helpers: T3[i,j,k] flattened as c3[i, (k-1)*n_s + j]? The n_s x n_s^2
  ## layout stores row i as as.numeric(matrix(n_s,n_s)) = column-major over (j,k)
  ## with j fastest. So col index = (k-1)*n_s + j, T3[i,j,k] = c3[i,(k-1)*n_s+j].
  ## Likewise M[a,(b,c)] = M[a,(c-1)*n_s+b].
  T3 <- array(0, dim = c(n_s, n_s, n_s))
  Ma <- array(0, dim = c(n_s, n_s, n_s))
  for (i in seq_len(n_s)) {
    T3[i, , ] <- matrix(c3[i, ], n_s, n_s)   # [j,k]
    Ma[i, , ] <- matrix(M[i, ], n_s, n_s)    # [b,c]
  }
  ## sign: A = I - hx^{⊗3}; residual (A C3 - b); bar_hx = -M · dA/dhx · C3
  ## = +M · d(hx^{⊗3})/dhx · C3.  Mode contractions (positive):
  for (a in seq_len(n_s)) for (i in seq_len(n_s)) {
    ## mode1: bar_hx[a,i] += sum_{b,c,j,k} M[a,b,c] hx[b,j] hx[c,k] T3[i,j,k]
    s <- 0
    for (b in seq_len(n_s)) for (cc in seq_len(n_s)) {
      hb <- hx[b, ]; hc <- hx[cc, ]
      s <- s + Ma[a, b, cc] * as.numeric(t(hb) %*% T3[i, , ] %*% hc)
    }
    bar_hx[a, i] <- bar_hx[a, i] + s
  }
  ## modes 2 and 3 by symmetry of the operator structure (b<->j, c<->k):
  for (b in seq_len(n_s)) for (j in seq_len(n_s)) {
    s2 <- 0; s3 <- 0
    for (a in seq_len(n_s)) for (cc in seq_len(n_s)) {
      ha <- hx[a, ]; hc <- hx[cc, ]
      ## mode2: bar_hx[b,j] += sum_{a,c,i,k} M[a,b,c] hx[a,i] hx[c,k] T3[i,j,k]
      s2 <- s2 + sum(vapply(seq_len(n_s), function(i)
        Ma[a, b, cc] * ha[i] * as.numeric(t(T3[i, j, ]) %*% hc), 0))
      ## mode3: bar_hx[c,k] handled in its own loop below
    }
    bar_hx[b, j] <- bar_hx[b, j] + s2
  }
  for (cc in seq_len(n_s)) for (k in seq_len(n_s)) {
    s3 <- 0
    for (a in seq_len(n_s)) for (b in seq_len(n_s)) {
      ha <- hx[a, ]; hb <- hx[b, ]
      s3 <- s3 + sum(vapply(seq_len(n_s), function(i)
        Ma[a, b, cc] * ha[i] * as.numeric(t(T3[i, , k]) %*% hb), 0))
    }
    bar_hx[cc, k] <- bar_hx[cc, k] + s3
  }

  list(bar_rhs = bar_rhs, bar_hx = Re(bar_hx))
}


## ===========================================================================
## Consumer 1: reverse-mode (adjoint_solution) cumulant gradient, O(1) in P.
## ---------------------------------------------------------------------------
## The forward chain differentiated here (see R/cumulant-likelihood.R):
##   solution blocks (ghx, ghu, ghxx, ghxu, ghuu, ghss, Sigma_e)
##     -> stationary moments Sigma_x (Lyapunov), mean, Sigma_y (order 1-2)
##     -> third-cumulant state tensor C211 (tensor-Lyapunov)
##     -> observable third cumulant c3_obs (B/C/D/E projection)
##   -> moment vector m(theta) -> quadratic-form loglik.
##
## Given the moment-vector cotangent bar_m (= dL/dm), this reverses the chain
## to block cotangents on the SEVEN blocks above, then ONE .solution_adjoint
## (first-order ghx/ghu/ys channel), ONE .solution_adjoint_order2 (the four
## order-2 blocks), and a per-parameter dSigma_e contraction turn those into
## the structural-parameter gradient -- all O(1) in the number of solves wrt P
## (the two adjoint kernels share their factorizations across parameters; only
## a cheap Frobenius contraction and one dSigma_e FD run per parameter).
##
## Order 4 (kurtosis) is NOT reversed in closed form here: its forward is a
## truncated-window generalized-chi-square trace + a closed-form chain whose
## symbolic reverse is a large, error-prone build (see the module-level STATUS
## note). When orders include 4, that block's cotangent contribution is taken
## by the exact FD-of-forward fallback in the caller. Orders 1-3 go fully
## through the reverse path.
## ---------------------------------------------------------------------------

#' Reverse of the discrete Lyapunov solve  P = A P A' + Q.
#'
#' Given cotangent bar_P (n x n) on the solution P (== Sigma_x), returns
#'   bar_Q  (n x n): dL/dQ,  solving the transposed Lyapunov bar_Q = A' bar_Q A + bar_P
#'   bar_A  (n x n): dL/dA = (bar_Q + bar_Q') A P     (P symmetric)
#' @noRd
.lyap_solve_adjoint <- function(A, P, bar_P) {
  n <- nrow(A)
  if (n == 0L) return(list(bar_Q = matrix(0, 0, 0), bar_A = matrix(0, 0, 0)))
  ## bar_Q solves M' vec(bar_Q) = vec(bar_P) with M = I - kron(A, A);
  ## M' = I - kron(A', A'), i.e. the Lyapunov with A -> A'.
  bar_Q <- .solve_lyapunov(t(A), bar_P)
  ## dP = A dP A' + (dA P A' + A P dA' + dQ); adjoint wrt A:
  ##   bar_A = (bar_Q + t(bar_Q)) A P   (using P = P')
  bar_A <- (bar_Q + t(bar_Q)) %*% A %*% P
  list(bar_Q = bar_Q, bar_A = bar_A)
}


#' Reverse-mode of compute_third_cumulant()'s observable projection + the
#' third-cumulant tensor-Lyapunov, into block cotangents.
#'
#' Mirrors compute_third_cumulant() EXACTLY (B/C/D/E terms, all permutations),
#' reversing each explicit projection and the .solve_third_cross_cumulant
#' tensor-Lyapunov (via .solve_third_cross_cumulant_adjoint).
#'
#' @param dr        DecisionRules2 at the base point.
#' @param model,params  as compute_third_cumulant.
#' @param bar_c3_obs n_endo x n_endo^2 cotangent dL/d(c3_obs).
#' @param fwd       optional cached compute_third_cumulant(dr, model, params).
#' @return list of cotangents on the blocks (all in the "full endo-row" layout
#'   used by the DR object): bar_ghx, bar_ghu, bar_ghxx, bar_ghxu, bar_ghuu,
#'   bar_Sigma_x (n_s x n_s), bar_Sigma_e (n_exo x n_exo), bar_hx (n_s x n_s),
#'   bar_hu (n_s x n_exo). Sigma_x/Sigma_e/hx/hu cotangents are the DIRECT ones;
#'   the caller folds bar_Sigma_x through the Lyapunov reverse.
#' @noRd
.compute_third_cumulant_adjoint <- function(dr, model, params, bar_c3_obs,
                                            fwd = NULL) {
  endo      <- dr$endo_names
  exo       <- dr$exo_names
  state_idx <- dr$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)
  n_s       <- length(state_idx)

  ghx  <- dr$ghx;   ghu  <- dr$ghu
  ghxx <- dr$ghxx;  ghuu <- dr$ghuu;  ghxu <- dr$ghxu

  shock_stderr <- .get_shock_stderr(model, exo, params)
  Sigma_e <- diag(shock_stderr^2, n_exo)

  hx  <- ghx[state_idx, , drop = FALSE]
  hu  <- ghu[state_idx, , drop = FALSE]
  hxx <- ghxx[state_idx, , drop = FALSE]
  hxu <- if (!is.null(ghxu)) ghxu[state_idx, , drop = FALSE] else NULL
  huu <- if (!is.null(ghuu)) ghuu[state_idx, , drop = FALSE] else NULL

  Sigma_x <- .state_covariance(hx, hu, Sigma_e)

  Z  <- ghx[, seq_len(n_s), drop = FALSE]   # n_endo x n_s
  Zu <- ghu                                  # n_endo x n_exo

  ## Forward intermediates (recompute cheaply; identical to the forward fn).
  rhs    <- .third_cumulant_rhs(hxx, Sigma_x, hxu = hxu, huu = huu,
                                hu = hu, hx = hx, Sigma_e = Sigma_e)
  c3_211 <- .solve_third_cross_cumulant(hx, rhs)   # n_s x n_s^2
  ZC3    <- Z %*% c3_211

  ## ---- cotangent accumulators on the blocks -------------------------------
  bar_Z    <- matrix(0, n_endo, n_s)      # -> ghx state cols
  bar_Zu   <- matrix(0, n_endo, n_exo)    # -> ghu
  bar_ghxx <- matrix(0, n_endo, n_s * n_s)
  bar_ghxu <- if (!is.null(ghxu)) matrix(0, n_endo, n_exo * n_s) else NULL
  bar_ghuu <- if (!is.null(ghuu)) matrix(0, n_endo, n_exo * n_exo) else NULL
  bar_Sigma_x <- matrix(0, n_s, n_s)
  bar_Sigma_e <- matrix(0, n_exo, n_exo)
  bar_ZC3  <- matrix(0, n_endo, n_s * n_s)   # cotangent on ZC3 = Z %*% c3_211
  bar_c3_211 <- matrix(0, n_s, n_s * n_s)    # cotangent on the state tensor

  ## Pre-materialise the forward per-observable matrices to reverse C/D/E terms.
  ZSig    <- Z %*% Sigma_x                      # n_endo x n_s  (C-term)
  ZuSe    <- Zu %*% Sigma_e                      # n_endo x n_exo (E-term)

  ## Reverse each observable row i, mirroring the forward accumulation order.
  for (i in seq_len(n_endo)) {
    bi <- matrix(bar_c3_obs[i, ], n_endo, n_endo)   # cotangent block (j,k)

    ## ============ B-term (x2 state, via C211) ============
    ## Perm 1: c3_obs[i,] += vec(Z %*% mat(ZC3[i,]) %*% t(Z))
    ## Perms 2/3: M2[jj,] = zi %*% mat(ZC3[jj,]) %*% t(Z); M3 = t(M2)
    ## reverse perm1:  A1 = Z' bi Z  (n_s x n_s) -> cotangent on mat(ZC3[i,])
    A1 <- t(Z) %*% bi %*% Z
    bar_ZC3[i, ] <- bar_ZC3[i, ] + as.numeric(A1)
    ## d(Z %*% M %*% t(Z))/dZ contributes to bar_Z:
    Mi <- matrix(ZC3[i, ], n_s, n_s)
    bar_Z <- bar_Z + bi %*% Z %*% t(Mi) + t(bi) %*% Z %*% Mi
    ## reverse perms 2 & 3:  forward added as.numeric(M2)+as.numeric(M3), M3=t(M2)
    ## so cotangent on M2 is bi + t(bi).
    barM2 <- bi + t(bi)                              # n_endo x n_endo
    zi <- Z[i, , drop = FALSE]                       # 1 x n_s
    for (jj in seq_len(n_endo)) {
      ## M2[jj,] = zi %*% mat(ZC3[jj,]) %*% t(Z) ; cotangent barM2[jj,] (len n_endo)
      g <- barM2[jj, ]                               # length n_endo
      Mjj <- matrix(ZC3[jj, ], n_s, n_s)
      ## d/d(mat(ZC3[jj,])) of zi %*% Mjj %*% t(Z) with cotangent g:
      ##   = t(zi) %*% g' %*% Z   (n_s x n_s), flatten
      bar_ZC3[jj, ] <- bar_ZC3[jj, ] + as.numeric(t(zi) %*% t(g) %*% Z)
      ## d/dZ  (two channels: through zi = Z[i,] and through the trailing t(Z))
      bar_Z[i, ] <- bar_Z[i, ] + as.numeric(Mjj %*% t(Z) %*% g)
      bar_Z <- bar_Z + outer(g, as.numeric(zi %*% Mjj))
    }

    ## ============ C-term (ghxx quadratic) ============
    ## Forward: A_m = ZSig %*% H_m %*% t(ZSig)  (n_endo x n_endo), ZSig = Z Sig_x.
    ## Perm1: c3_obs[i,] += vec(A_i)         -> cotangent on A_i is bi.
    ## Perms2/3: P2[jj,]=A_jj[i,]; P3[,jj]=A_jj[i,]; forward adds P2+P3.
    ##   => cotangent on row i of A_jj is bi[jj,] (P2) + bi[,jj] (P3).
    ## Accumulate ZSig cotangents here (local), reverse ZSig once at the end.
    bar_ZSig <- matrix(0, n_endo, n_s)
    for (jj in seq_len(n_endo)) {
      ## cotangent on the (full) matrix A_jj:
      ##   perm1 (only jj==i): += bi
      ##   perms2/3: only row i receives (bi[jj,] + bi[,jj])
      barA <- matrix(0, n_endo, n_endo)
      if (jj == i) barA <- barA + bi
      barA[i, ] <- barA[i, ] + bi[jj, ] + bi[, jj]
      if (all(barA == 0)) next
      Hjj <- matrix(ghxx[jj, ], n_s, n_s)
      ## A_jj = ZSig Hjj ZSig'
      bar_ghxx[jj, ] <- bar_ghxx[jj, ] + as.numeric(t(ZSig) %*% barA %*% ZSig)
      bar_ZSig <- bar_ZSig + barA %*% ZSig %*% t(Hjj) + t(barA) %*% ZSig %*% Hjj
    }
    ## reverse ZSig = Z %*% Sigma_x -> bar_Z, bar_Sigma_x
    bar_Z <- bar_Z + bar_ZSig %*% t(Sigma_x)
    bar_Sigma_x <- bar_Sigma_x + t(Z) %*% bar_ZSig

    ## ============ D-term (ghxu cross) ============
    if (!is.null(ghxu)) {
      ## B_m = Zu %*% Sigma_e %*% Xu_m %*% Sigma_x %*% t(Z)   (n_endo x n_endo)
      ## Perm1: c3_obs[i,] += vec(B_i + t(B_i))
      ## Perms2/3: D2[jj,]=B_jj[,i]+B_jj[i,]; D3[,jj]=B_jj[,i]+B_jj[i,]
      ## reverse perm1: cotangent on B_i is (bi + t(bi))
      barB_i <- bi + t(bi)
      bar_B <- vector("list", n_endo)
      bar_B[[i]] <- barB_i
      for (jj in seq_len(n_endo)) {
        if (jj == i) next
        bar_B[[jj]] <- matrix(0, n_endo, n_endo)
      }
      ## perms 2/3: cotangent on (B_jj[,i]+B_jj[i,]) is bi[jj,] (from D2) +
      ## bi[,jj] (from D3). Distribute onto B_jj columns/rows.
      for (jj in seq_len(n_endo)) {
        w <- bi[jj, ] + bi[, jj]                     # length n_endo
        if (is.null(bar_B[[jj]])) bar_B[[jj]] <- matrix(0, n_endo, n_endo)
        bar_B[[jj]][, i] <- bar_B[[jj]][, i] + w      # from B_jj[,i]
        bar_B[[jj]][i, ] <- bar_B[[jj]][i, ] + w      # from B_jj[i,]
      }
      ## reverse each B_m = Zu Se Xu_m Sig Z' into blocks.
      for (m in seq_len(n_endo)) {
        bm <- bar_B[[m]]
        if (all(bm == 0)) next
        Xu_m <- matrix(ghxu[m, ], n_exo, n_s, byrow = TRUE)
        ## Let P = Zu %*% Sigma_e (n_endo x n_exo), Q = Sigma_x %*% t(Z) (n_s x n_endo)
        P <- Zu %*% Sigma_e
        Qm <- Sigma_x %*% t(Z)
        ## B_m = P %*% Xu_m %*% Qm
        ## bar_Xu_m = P' bm Qm'
        bar_Xu <- t(P) %*% bm %*% t(Qm)               # n_exo x n_s
        ## store row-major into bar_ghxu[m,]
        bar_ghxu[m, ] <- bar_ghxu[m, ] + as.numeric(t(bar_Xu))
        ## bar_P += bm Qm' Xu_m'  ; P = Zu Sigma_e
        bar_P <- bm %*% t(Qm) %*% t(Xu_m)             # n_endo x n_exo
        bar_Zu <- bar_Zu + bar_P %*% Sigma_e
        bar_Sigma_e <- bar_Sigma_e + t(Zu) %*% bar_P
        ## bar_Qm += Xu_m' P' bm ; Qm = Sigma_x Z'
        bar_Qm <- t(Xu_m) %*% t(P) %*% bm             # n_s x n_endo
        bar_Sigma_x <- bar_Sigma_x + bar_Qm %*% Z
        bar_Z <- bar_Z + t(bar_Qm) %*% Sigma_x        # Q = Sig Z' -> bar_Z via Sig' bar_Qm? handle t
      }
    }

    ## ============ E-term (ghuu quadratic in e) ============
    if (!is.null(ghuu)) {
      ## E_m = ZuSe %*% Uu_m %*% t(ZuSe)  (n_endo x n_endo), ZuSe = Zu Sigma_e
      ## Perm1: c3_obs[i,] += vec(E_i)
      ## Perms2/3: E2[jj,]=E_jj[i,]; E3[,jj]=E_jj[i,]
      barE_i <- bi
      Uu_i <- matrix(ghuu[i, ], n_exo, n_exo, byrow = TRUE)
      bar_Uu_i <- t(ZuSe) %*% barE_i %*% ZuSe
      bar_ghuu[i, ] <- bar_ghuu[i, ] + as.numeric(t(bar_Uu_i))
      bar_ZuSe <- matrix(0, n_endo, n_exo)
      bar_ZuSe <- bar_ZuSe + barE_i %*% ZuSe %*% t(Uu_i) +
                             t(barE_i) %*% ZuSe %*% Uu_i
      for (jj in seq_len(n_endo)) {
        grow <- bi[jj, ] + bi[, jj]                  # cotangent on E_jj[i,]
        Uu_jj <- matrix(ghuu[jj, ], n_exo, n_exo, byrow = TRUE)
        zusei <- ZuSe[i, , drop = FALSE]             # 1 x n_exo
        bar_Uu_jj <- t(zusei) %*% t(grow) %*% ZuSe   # n_exo x n_exo
        bar_ghuu[jj, ] <- bar_ghuu[jj, ] + as.numeric(t(bar_Uu_jj))
        bar_ZuSe[i, ] <- bar_ZuSe[i, ] + as.numeric(Uu_jj %*% t(ZuSe) %*% grow)
        bar_ZuSe <- bar_ZuSe + outer(grow, as.numeric(zusei %*% Uu_jj))
      }
      ## reverse ZuSe = Zu %*% Sigma_e
      bar_Zu <- bar_Zu + bar_ZuSe %*% Sigma_e
      bar_Sigma_e <- bar_Sigma_e + t(Zu) %*% bar_ZuSe
    }
  }

  ## ---- reverse ZC3 = Z %*% c3_211 -> bar_Z, bar_c3_211 --------------------
  bar_Z <- bar_Z + bar_ZC3 %*% t(c3_211)
  bar_c3_211 <- bar_c3_211 + t(Z) %*% bar_ZC3

  ## ---- reverse the tensor-Lyapunov solve (C211) ---------------------------
  ## .solve_third_cross_cumulant(hx, rhs) with the (j,k) symmetrization handled
  ## inside the adjoint. Returns bar_rhs and bar_hx (state block).
  tl <- .solve_third_cross_cumulant_adjoint(hx, c3_211, bar_c3_211)
  bar_rhs <- tl$bar_rhs
  bar_hx  <- tl$bar_hx                       # n_s x n_s

  ## ---- reverse .third_cumulant_rhs -> bar_hxx, bar_hxu, bar_huu, hx/hu/Se --
  rev_rhs <- .third_cumulant_rhs_adjoint(hxx, Sigma_x, hxu, huu, hu, hx,
                                         Sigma_e, bar_rhs)
  bar_hxx <- rev_rhs$bar_hxx                  # n_s x n_s^2 (state rows)
  bar_hxu <- rev_rhs$bar_hxu                  # n_s x (n_exo*n_s) or NULL
  bar_huu <- rev_rhs$bar_huu                  # n_s x n_exo^2 or NULL
  bar_hx  <- bar_hx + rev_rhs$bar_hx
  bar_hu  <- rev_rhs$bar_hu                   # n_s x n_exo
  bar_Sigma_x <- bar_Sigma_x + rev_rhs$bar_Sigma_x
  bar_Sigma_e <- bar_Sigma_e + rev_rhs$bar_Sigma_e

  ## ---- assemble full endo-row block cotangents ----------------------------
  ## bar_Z is on ghx[, 1:n_s]; bar_Zu is on ghu; bar_hx/bar_hu/bar_hxx/... are
  ## on the state ROWS of the full blocks.
  bar_ghx <- matrix(0, n_endo, ncol(ghx))
  bar_ghx[, seq_len(n_s)] <- bar_Z
  ## hx cotangent adds into ghx[state_idx, 1:n_s]
  bar_ghx[state_idx, seq_len(n_s)] <- bar_ghx[state_idx, seq_len(n_s)] + bar_hx

  bar_ghu <- bar_Zu
  bar_ghu[state_idx, ] <- bar_ghu[state_idx, , drop = FALSE] + bar_hu

  ## ghxx/ghxu/ghuu: the direct (C/D/E) cotangents are already full endo-row;
  ## the state-row (B-term via C211 driving) cotangents add onto state rows.
  bar_ghxx[state_idx, ] <- bar_ghxx[state_idx, , drop = FALSE] + bar_hxx
  if (!is.null(bar_ghxu) && !is.null(bar_hxu))
    bar_ghxu[state_idx, ] <- bar_ghxu[state_idx, , drop = FALSE] + bar_hxu
  if (!is.null(bar_ghuu) && !is.null(bar_huu))
    bar_ghuu[state_idx, ] <- bar_ghuu[state_idx, , drop = FALSE] + bar_huu

  list(bar_ghx = bar_ghx, bar_ghu = bar_ghu,
       bar_ghxx = bar_ghxx, bar_ghxu = bar_ghxu, bar_ghuu = bar_ghuu,
       bar_Sigma_x = bar_Sigma_x, bar_Sigma_e = bar_Sigma_e)
}


#' Reverse of .third_cumulant_rhs into cotangents on the state-block inputs.
#' Forward: rhs[i,] = 2*(hxS Hxx_i hxS') + [2*(W_i + W_i') for hxu]
#'          + [2*(huSe Uu_i huSe') for huu], hxS = hx Sigma_x, huSe = hu Sigma_e,
#'          W_i = huSe Xu_i hxS'. All matrices n_s-square blocks per state row i.
#' @noRd
.third_cumulant_rhs_adjoint <- function(hxx, Sigma_x, hxu, huu, hu, hx,
                                        Sigma_e, bar_rhs) {
  n_s   <- nrow(Sigma_x)
  n_exo <- nrow(Sigma_e)
  hxS   <- hx %*% Sigma_x                    # n_s x n_s
  huSe  <- hu %*% Sigma_e                    # n_s x n_exo

  bar_hxx <- matrix(0, n_s, n_s * n_s)
  bar_hxu <- if (!is.null(hxu)) matrix(0, n_s, n_exo * n_s) else NULL
  bar_huu <- if (!is.null(huu)) matrix(0, n_s, n_exo * n_exo) else NULL
  bar_hxS  <- matrix(0, n_s, n_s)
  bar_huSe <- matrix(0, n_s, n_exo)
  bar_hx   <- matrix(0, n_s, n_s)
  bar_hu   <- matrix(0, n_s, n_exo)
  bar_Sigma_x <- matrix(0, n_s, n_s)
  bar_Sigma_e <- matrix(0, n_exo, n_exo)

  for (i in seq_len(n_s)) {
    Bi <- matrix(bar_rhs[i, ], n_s, n_s)     # cotangent on rhs row i (n_s x n_s)

    ## hxx term: R = 2 hxS Hi hxS'
    Hi <- matrix(hxx[i, ], n_s, n_s)
    bar_Hi <- 2 * (t(hxS) %*% Bi %*% hxS)
    bar_hxx[i, ] <- bar_hxx[i, ] + as.numeric(bar_Hi)
    ## bar_hxS from R = 2 hxS Hi hxS' : 2*(Bi hxS Hi' + Bi' hxS Hi)
    bar_hxS <- bar_hxS + 2 * (Bi %*% hxS %*% t(Hi) + t(Bi) %*% hxS %*% Hi)

    ## hxu term: R = 2 (W_i + W_i'), W_i = huSe Xu_i hxS'
    if (!is.null(hxu)) {
      Xu_i <- matrix(hxu[i, ], n_exo, n_s, byrow = TRUE)
      ## cotangent on W_i is 2*(Bi + Bi')
      barW <- 2 * (Bi + t(Bi))               # n_s x n_s
      ## W_i = huSe Xu_i hxS'
      bar_Xu <- t(huSe) %*% barW %*% hxS      # n_exo x n_s
      bar_hxu[i, ] <- bar_hxu[i, ] + as.numeric(t(bar_Xu))
      bar_huSe <- bar_huSe + barW %*% hxS %*% t(Xu_i)
      bar_hxS  <- bar_hxS  + t(barW) %*% huSe %*% Xu_i
    }

    ## huu term: R = 2 huSe Uu_i huSe'
    if (!is.null(huu)) {
      Uu_i <- matrix(huu[i, ], n_exo, n_exo, byrow = TRUE)
      bar_Uu <- 2 * (t(huSe) %*% Bi %*% huSe)  # n_exo x n_exo
      bar_huu[i, ] <- bar_huu[i, ] + as.numeric(t(bar_Uu))
      bar_huSe <- bar_huSe + 2 * (Bi %*% huSe %*% t(Uu_i) + t(Bi) %*% huSe %*% Uu_i)
    }
  }

  ## reverse hxS = hx Sigma_x
  bar_hx      <- bar_hx      + bar_hxS %*% t(Sigma_x)
  bar_Sigma_x <- bar_Sigma_x + t(hx) %*% bar_hxS
  ## reverse huSe = hu Sigma_e
  bar_hu      <- bar_hu      + bar_huSe %*% t(Sigma_e)
  bar_Sigma_e <- bar_Sigma_e + t(hu) %*% bar_huSe

  list(bar_hxx = bar_hxx, bar_hxu = bar_hxu, bar_huu = bar_huu,
       bar_hx = bar_hx, bar_hu = bar_hu,
       bar_Sigma_x = bar_Sigma_x, bar_Sigma_e = bar_Sigma_e)
}
