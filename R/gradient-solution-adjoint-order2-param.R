## R/gradient-solution-adjoint-order2-param.R
## --------------------------------------------------------------------------
## ANALYTIC, d2X-FREE second-parameter-derivative (Hessian) of the frozen-bars
## solution-adjoint gradient (Tier 18 A2, "adjoint_solution" T2 route for
## posterior_hessian).
##
## posterior_hessian's T2 term is T2[i,j] = <G_X, d2X_ij>, the Hessian of
##   phi(theta) = <G_X, X(theta)>
## with the Kalman-filter adjoint bar matrices G_X FROZEN at the base point
## (they depend only on Y, ss, me_variance -- NOT on theta). Since
## grad phi = <G_X, dX/dtheta> is EXACTLY one \code{.solution_adjoint} call
## (gradient-solution-adjoint.R) with \code{bars = G_X}, T2 is exactly the
## JACOBIAN, in theta_i, of that call's analytic gradient formula:
##
##   grad_j(theta) = <bar_df_plus,    df_plus_j>  + <bar_df_zero, df_zero_j>
##                 + <bar_df_minus_S, df_minus_j[,state_idx]>
##                 + <bar_df_u,       df_exo_j>    + <bar_ys, dys_j>          (*)
##
## where bar_df_plus, bar_df_zero, bar_df_minus_S, bar_df_u are FUNCTIONS OF
## THETA (through A, GS, SG, G, H -- the base-point solution objects -- via
## two transposed linear solves), and bar_ys is a CONSTANT (it is built solely
## from the frozen bars$g_d and the fixed obs_vars index set, with no
## dependence on A/G/H/theta at all).
##
## Differentiating (*) in theta_i by the product rule and using
## d(bar_ys)/dtheta_i = 0:
##
##   T2[i,j] = <d(bar_df_plus)/dtheta_i,    df_plus_j> + <bar_df_plus,    d2f_plus_ij>
##           + <d(bar_df_zero)/dtheta_i,    df_zero_j> + <bar_df_zero,    d2f_zero_ij>
##           + <d(bar_df_minus_S)/dtheta_i, df_minus_j[,state_idx]>
##                                          + <bar_df_minus_S, d2f_minus_ij[,state_idx]>
##           + <d(bar_df_u)/dtheta_i,       df_exo_j>   + <bar_df_u,       d2f_exo_ij>
##           + <bar_ys, d2ys_ij>                                              (**)
##
## The SECOND total primitive derivatives (d2f_plus_ij, d2f_zero_ij,
## d2f_minus_ij, d2f_exo_ij, d2ys_ij) are obtained EXACTLY the way
## solution_derivatives_2 (gradient-solution-deriv-2.R) does: the analytic
## 5-term second-total-derivative contraction (\code{.sd2_build_ctx} /
## \code{.sd2_d2J_pair}) when the compiled model carries the order-2 parameter
## codegen, else the 3-point diagonal / 4-corner mixed FD stencil of
## \code{.sd2_primitives_at}. We never form d2G_ij / d2H_ij / d2TT_ij etc (the
## O(n_endo x n_state) second-order SOLUTION blocks solution_derivatives_2
## computes) -- only the much smaller second-order PRIMITIVE blocks
## (n_endo x n_endo / n_endo x n_exo), contracted directly against the
## base-point bar matrices. That is the d2X-free saving.
##
## THE d(bar)/dtheta_i CHAIN (standard adjoint-of-linear-solve, applied twice)
## ----------------------------------------------------------------------
## Recall (gradient-solution-adjoint.R) the two transposed solves:
##   W = A^{-T} V_H              (V_H CONSTANT; bar_df_u = -W, bar_dA = -W H')
##   A' L + f_plus' L SG' = V_G  (V_G = V_G_const + f_plus' bar_dA S'; L solved
##                                by the SAME dedicated real-Schur Sylvester
##                                factorization; bar_RHS = -L)
## Differentiating (with dA_i = df_plus_i GS + f_plus dGS_i + df_zero_i, the
## SAME A_i formula solution_derivatives_2 uses; dGS_i = dG_i S, dSG_i = S dG_i;
## dG_i, dH_i from the validated FIRST-order layer, solution_derivatives):
##
##   A' dW_i = -dA_i' W                       (dW_i: ONE extra transposed
##                                              A-solve, reusing At_qr)
##   dbar_df_u_i = -dW_i,   dbar_dA_i = -dW_i H' - W dH_i'
##
##   dV_G_i = df_plus_i' bar_dA S' + f_plus' dbar_dA_i S'
##   A' dL_i + f_plus' dL_i SG' =
##       dV_G_i - dA_i' L - df_plus_i' L SG' - f_plus' L dSG_i'   (dL_i: ONE
##                                              extra transposed Sylvester
##                                              solve, reusing sylv_fac_T)
##   dbar_RHS_i = -dL_i
##
##   dbar_df_plus_i = dbar_dA_i (GS)' + bar_dA (dGS_i)'
##                  + dbar_RHS_i SG' G' + bar_RHS (dSG_i)' G' + bar_RHS SG' dG_i'
##   dbar_df_zero_i = dbar_dA_i + dbar_RHS_i G' + bar_RHS dG_i'
##   dbar_df_minus_S_i = dbar_RHS_i
##
## Cost: the SAME O(1) factorizations (At_qr, sylv_fac_T) built once at the
## base point are reused for every i (one extra cheap solve per i, no new
## factorization) -- exactly the "standard adjoint-of-linear-solve" pattern.
## The O(np^2) cost lives entirely in the second-primitive stencil (shared
## with solution_derivatives_2's own cost profile), not in any linear solve.
##
## Certified against t2_method = "contract_once" (the exact d2X-based ground
## truth) by test-hessian-t2-adjoint-solution.R, including the flattest
## (smallest-magnitude) eigenvalues.
## --------------------------------------------------------------------------


#' Analytic Hessian of the frozen-bars solution-adjoint gradient:
#' T2[i,j] = d^2/dtheta_i dtheta_j <G_X, X(theta)>, G_X FROZEN at the base
#' point (\code{bars}). FD-free at the solution-adjoint layer (two transposed
#' linear solves per parameter, reusing the base-point factorizations); the
#' second-order PRIMITIVE derivatives are analytic or FD exactly as
#' \code{solution_derivatives_2} computes them, never materialising any
#' second-order SOLUTION block (d2G/d2H/d2TT/...).
#'
#' @param model,compiled,dr,params,obs_vars  as in \code{solution_derivatives}.
#' @param param_names  structural parameters to differentiate (the returned
#'   Hessian is \code{length(param_names) x length(param_names)}).
#' @param bars  list(G_TT, G_RR, G_ZZ, G_DD, g_d) -- the frozen Kalman-filter
#'   adjoint bar matrices (same contract as \code{.solution_adjoint}'s
#'   \code{bars} argument; any element may be NULL, treated as zero).
#' @param h_rel  FD step for the FIRST-order primitive fallback (only used
#'   when the analytic first-order codegen is unavailable).
#' @param h_rel2 relative FD step for the SECOND-order primitive stencils
#'   (default 1e-4, matching \code{solution_derivatives_2}).
#' @param use_analytic NULL (option-driven, default TRUE) / TRUE / FALSE --
#'   gates the FIRST-order analytic primitive path (same contract as
#'   \code{solution_derivatives}/\code{.solution_adjoint}). The SECOND-order
#'   primitive path has its own independent analytic gate (mirroring
#'   \code{solution_derivatives_2}: requires \code{param_deriv2_ok} AND
#'   \code{static_param2_built}).
#' @return list(T2 = matrix, ok = named logical (per-parameter first-order
#'   primitive success), used_analytic_first = logical, used_analytic_second
#'   = logical)
#' @noRd
.solution_adjoint_hessian <- function(model, compiled, dr, params, param_names,
                                      obs_vars, bars, h_rel = 1e-6,
                                      h_rel2 = 1e-4, use_analytic = NULL) {

  if (!all(param_names %in% names(params))) {
    missing_p <- param_names[!param_names %in% names(params)]
    stop(sprintf(".solution_adjoint_hessian: parameter(s) not found in `params`: %s",
                 paste(missing_p, collapse = ", ")))
  }

  endo <- dr$endo_names; n_endo <- length(endo)
  exo  <- dr$exo_names;  n_exo  <- length(exo)
  state_idx <- dr$state_idx; n_state <- length(state_idx)
  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop(".solution_adjoint_hessian: observed variables not found: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))
  n_obs <- length(obs_idx)
  np <- length(param_names)

  G <- dr$ghx; H <- dr$ghu; ys <- dr$ys

  ## ---- base system (mirrors .solution_adjoint / solution_derivatives_2 exactly)
  S <- matrix(0, nrow = n_state, ncol = n_endo)
  if (n_state > 0) S[cbind(seq_len(n_state), state_idx)] <- 1
  sys0 <- extract_system_matrices(compiled, ys, params)
  f_plus <- sys0$f_plus; f_zero <- sys0$f_zero
  f_minus <- sys0$f_minus; f_u <- sys0$f_exo
  SG <- if (n_state > 0) S %*% G else matrix(0, 0, 0)
  GS <- if (n_state > 0) G %*% S else matrix(0, n_endo, n_endo)
  A  <- f_plus %*% GS + f_zero

  zmat <- function(x, nr, nc) {
    if (is.null(x)) return(matrix(0, nr, nc))
    x <- as.matrix(x)
    stopifnot(nrow(x) == nr, ncol(x) == nc)
    x
  }
  G_TT <- zmat(bars$G_TT, n_state, n_state)
  G_RR <- zmat(bars$G_RR, n_state, n_exo)
  G_ZZ <- zmat(bars$G_ZZ, n_obs,  n_state)
  G_DD <- zmat(bars$G_DD, n_obs,  n_exo)
  g_d  <- if (is.null(bars$g_d)) numeric(n_obs) else as.numeric(bars$g_d)
  stopifnot(length(g_d) == n_obs)

  ## ---- reverse of the block extraction (CONSTANT in theta) -----------------
  V_G_const <- matrix(0, n_endo, n_state)
  V_H       <- matrix(0, n_endo, n_exo)
  if (n_state > 0) {
    V_G_const[state_idx, ] <- V_G_const[state_idx, , drop = FALSE] + G_TT
    V_H[state_idx, ]       <- V_H[state_idx, , drop = FALSE] + G_RR
  }
  V_G_const[obs_idx, ] <- V_G_const[obs_idx, , drop = FALSE] + G_ZZ
  V_H[obs_idx, ]       <- V_H[obs_idx, , drop = FALSE] + G_DD
  bar_ys <- setNames(numeric(n_endo), endo)
  bar_ys[obs_vars] <- bar_ys[obs_vars] + g_d

  ## ---- base-point transposed A-solve: W = A^{-T} V_H -----------------------
  At_qr <- qr(t(A))
  W <- qr.solve(At_qr, V_H)
  bar_df_u <- -W
  bar_dA   <- -W %*% t(H)

  bar_df_plus <- bar_dA %*% t(GS)
  bar_df_zero <- bar_dA

  V_G <- V_G_const
  if (n_state > 0) V_G <- V_G_const + t(f_plus) %*% bar_dA %*% t(S)

  ## ---- base-point transposed Sylvester solve: A'L + f_plus'L SG' = V_G -----
  tA <- t(A); tfp <- t(f_plus)
  tSG <- if (n_state > 0) t(SG) else matrix(0, 0, 0)
  sylv_fac_T <- NULL; L <- NULL
  bar_RHS <- matrix(0, n_endo, n_state)
  if (n_state > 0) {
    sylv_fac_T <- .gen_sylvester_k1_factor(tA, tfp, tSG)
    L <- if (!is.null(sylv_fac_T))
      .gen_sylvester_k1_solve(sylv_fac_T, V_G)
    else
      .solve_kron_compact(tA, tfp, tSG, k = 1L, RHS = V_G)
    resid <- max(abs(tA %*% L + tfp %*% L %*% tSG - V_G))
    if (!is.finite(resid) || resid > 1e-6 * max(1, max(abs(V_G))))
      stop(sprintf(".solution_adjoint_hessian: transposed Sylvester residual %.3e", resid))
    bar_RHS <- -L
    bar_df_plus <- bar_df_plus + bar_RHS %*% tSG %*% t(G)
    bar_df_zero <- bar_df_zero + bar_RHS %*% t(G)
  }
  bar_df_minus_S <- bar_RHS

  ## ---- per-parameter FIRST-order solution + total-primitive derivatives ----
  if (is.null(use_analytic))
    use_analytic <- isTRUE(getOption("dynhr.use_analytic_primitives", TRUE))
  analytic <- NULL
  if (use_analytic && .can_use_analytic_primitive_deriv(compiled)) {
    dys_all <- .analytic_dys(compiled, ys, params)
    dprim   <- if (!is.null(dys_all))
      .analytic_dprimitives(compiled, ys, params, dys_all) else NULL
    if (!is.null(dys_all) && !is.null(dprim))
      analytic <- list(dys = dys_all, dprim = dprim)
  }

  ## dG_k, dH_k -- the SAME validated first-order layer solution_derivatives_2
  ## uses for its G1/H1 (needed here for the dA_i / dGS_i / dSG_i chain).
  sd1 <- solution_derivatives(model, compiled, dr, params, param_names,
                              obs_vars = endo, h_rel = h_rel)

  ok <- setNames(rep(FALSE, np), param_names)
  prim1 <- vector("list", np); names(prim1) <- param_names
  for (k in seq_len(np)) {
    pnm <- param_names[k]
    d1 <- sd1$derivs[[pnm]]
    if (is.null(d1) || !isTRUE(d1$ok)) next
    pr <- tryCatch({
      if (!is.null(analytic)) {
        c(list(dys = setNames(as.numeric(analytic$dys[, pnm]),
                              rownames(analytic$dys))),
          analytic$dprim[[pnm]])
      } else {
        .solution_adjoint_fd_prims(pnm, model, compiled, ys, params, h_rel)
      }
    }, error = function(e) NULL)
    if (is.null(pr) || is.null(pr$df_plus)) next
    prim1[[pnm]] <- list(dG = d1$dG, dH = d1$dH,
                         df_plus = pr$df_plus, df_zero = pr$df_zero,
                         df_minus = pr$df_minus, df_exo = pr$df_exo,
                         dys = pr$dys)
    ok[pnm] <- TRUE
  }

  ## ---- SECOND-order primitive machinery (mirrors solution_derivatives_2) ---
  analytic_ok2 <- isTRUE(getOption("dynhr.use_analytic_primitives", TRUE)) &&
    .can_use_analytic_primitive_deriv(compiled) &&
    isTRUE(compiled$dynamic$param_deriv2_ok) &&
    isTRUE(compiled$static$static_param2_built)

  ad2 <- NULL; pidx <- NULL; d2ys_full <- NULL
  if (analytic_ok2) {
    dys_full   <- if (!is.null(analytic)) analytic$dys else .analytic_dys(compiled, ys, params)
    dprim_full <- if (!is.null(analytic)) analytic$dprim else
      if (!is.null(dys_full)) .analytic_dprimitives(compiled, ys, params, dys_full) else NULL
    d2ys_full  <- if (!is.null(dys_full)) .analytic_d2ys(compiled, ys, params, dys_full) else NULL
    ad2 <- if (!is.null(dys_full) && !is.null(d2ys_full))
      .sd2_build_ctx(compiled, ys, params, dys_full, d2ys_full) else NULL
    if (is.null(dys_full) || is.null(dprim_full) || is.null(d2ys_full) || is.null(ad2)) {
      analytic_ok2 <- FALSE
    } else {
      pidx <- match(param_names, compiled$model$param_names)
    }
  }

  hvec <- vapply(param_names, function(p) max(h_rel2 * abs(params[[p]]), 1e-5), 0)
  names(hvec) <- param_names

  ## Single-sided perturbed primitives for the FD fallback (3-pt diagonal +
  ## 4-corner mixed stencils), computed ONCE and cached -- same construction
  ## solution_derivatives_2 uses (.sd2_primitives_at).
  prim_p <- prim_m <- NULL
  if (!analytic_ok2) {
    prim_p <- prim_m <- vector("list", np); names(prim_p) <- names(prim_m) <- param_names
    for (k in seq_len(np)) {
      if (!ok[param_names[k]]) next
      pnm <- param_names[k]; h <- hvec[k]
      tp <- params; tp[[pnm]] <- tp[[pnm]] + h
      tm <- params; tm[[pnm]] <- tm[[pnm]] - h
      prim_p[[k]] <- tryCatch(.sd2_primitives_at(tp, model, compiled, ys), error = function(e) NULL)
      prim_m[[k]] <- tryCatch(.sd2_primitives_at(tm, model, compiled, ys), error = function(e) NULL)
    }
  }

  ## Second-total-derivative block for pair (i,j): list(d2fp, d2f0, d2fm, d2fu, d2ys).
  d2primf <- function(i, j) {
    if (analytic_ok2) {
      a <- pidx[i]; b <- pidx[j]
      dpr2 <- .sd2_d2J_pair(ad2, a, b)
      return(list(d2fp = dpr2$df_plus, d2f0 = dpr2$df_zero, d2fm = dpr2$df_minus,
                  d2fu = dpr2$df_exo, d2ys = d2ys_full[, a, b]))
    }
    pi_name <- param_names[i]; pj_name <- param_names[j]
    hi <- hvec[i]; hj <- hvec[j]
    if (is.null(prim_p[[i]]) || is.null(prim_m[[i]]) ||
        is.null(prim_p[[j]]) || is.null(prim_m[[j]])) return(NULL)
    if (i == j) {
      list(d2fp = (prim_p[[i]]$f_plus  - 2 * f_plus  + prim_m[[i]]$f_plus)  / hi^2,
           d2f0 = (prim_p[[i]]$f_zero  - 2 * f_zero  + prim_m[[i]]$f_zero)  / hi^2,
           d2fm = (prim_p[[i]]$f_minus - 2 * f_minus + prim_m[[i]]$f_minus) / hi^2,
           d2fu = (prim_p[[i]]$f_exo   - 2 * f_u      + prim_m[[i]]$f_exo)  / hi^2,
           d2ys = (prim_p[[i]]$ys      - 2 * ys       + prim_m[[i]]$ys)     / hi^2)
    } else {
      tpp <- params; tpp[[pi_name]] <- tpp[[pi_name]] + hi; tpp[[pj_name]] <- tpp[[pj_name]] + hj
      tpm <- params; tpm[[pi_name]] <- tpm[[pi_name]] + hi; tpm[[pj_name]] <- tpm[[pj_name]] - hj
      tmp <- params; tmp[[pi_name]] <- tmp[[pi_name]] - hi; tmp[[pj_name]] <- tmp[[pj_name]] + hj
      tmm <- params; tmm[[pi_name]] <- tmm[[pi_name]] - hi; tmm[[pj_name]] <- tmm[[pj_name]] - hj
      Fpp <- tryCatch(.sd2_primitives_at(tpp, model, compiled, ys), error = function(e) NULL)
      Fpm <- tryCatch(.sd2_primitives_at(tpm, model, compiled, ys), error = function(e) NULL)
      Fmp <- tryCatch(.sd2_primitives_at(tmp, model, compiled, ys), error = function(e) NULL)
      Fmm <- tryCatch(.sd2_primitives_at(tmm, model, compiled, ys), error = function(e) NULL)
      if (is.null(Fpp) || is.null(Fpm) || is.null(Fmp) || is.null(Fmm)) return(NULL)
      den <- 4 * hi * hj
      mix <- function(fld) (Fpp[[fld]] - Fpm[[fld]] - Fmp[[fld]] + Fmm[[fld]]) / den
      list(d2fp = mix("f_plus"), d2f0 = mix("f_zero"), d2fm = mix("f_minus"),
           d2fu = mix("f_exo"), d2ys = mix("ys"))
    }
  }

  ## ---- per-i derivative of the base-point bar matrices (forward-over-reverse) --
  T2 <- matrix(0, np, np, dimnames = list(param_names, param_names))

  for (i in seq_len(np)) {
    if (!ok[param_names[i]]) next
    pi <- prim1[[param_names[i]]]
    Gi <- pi$dG; Hi <- pi$dH
    dfp_i <- pi$df_plus; df0_i <- pi$df_zero

    dGS_i <- if (n_state > 0) Gi %*% S else matrix(0, n_endo, n_endo)
    dSG_i <- if (n_state > 0) S %*% Gi else matrix(0, 0, 0)
    dA_i  <- dfp_i %*% GS + f_plus %*% dGS_i + df0_i

    dW_i <- qr.solve(At_qr, -t(dA_i) %*% W)
    dbar_df_u_i <- -dW_i
    dbar_dA_i   <- -dW_i %*% t(H) - W %*% t(Hi)

    dbar_df_plus_i <- dbar_dA_i %*% t(GS) + bar_dA %*% t(dGS_i)
    dbar_df_zero_i <- dbar_dA_i

    if (n_state > 0) {
      dV_G_i <- t(dfp_i) %*% bar_dA %*% t(S) + t(f_plus) %*% dbar_dA_i %*% t(S)
      RHS_L_i <- dV_G_i - t(dA_i) %*% L - t(dfp_i) %*% L %*% tSG -
        t(f_plus) %*% L %*% t(dSG_i)
      dL_i <- if (!is.null(sylv_fac_T))
        .gen_sylvester_k1_solve(sylv_fac_T, RHS_L_i)
      else
        .solve_kron_compact(tA, tfp, tSG, k = 1L, RHS = RHS_L_i)
      dbar_RHS_i <- -dL_i
      dbar_df_plus_i <- dbar_df_plus_i +
        dbar_RHS_i %*% tSG %*% t(G) + bar_RHS %*% t(dSG_i) %*% t(G) +
        bar_RHS %*% tSG %*% t(Gi)
      dbar_df_zero_i <- dbar_df_zero_i + dbar_RHS_i %*% t(G) + bar_RHS %*% t(Gi)
      dbar_df_minus_S_i <- dbar_RHS_i
    } else {
      dbar_df_minus_S_i <- matrix(0, n_endo, 0)
    }

    for (j in seq_len(np)) {
      if (!ok[param_names[j]]) next
      pj <- prim1[[param_names[j]]]
      d2p <- d2primf(i, j)
      if (is.null(d2p)) next

      gij <- sum(dbar_df_plus_i * pj$df_plus) + sum(bar_df_plus * d2p$d2fp) +
             sum(dbar_df_zero_i * pj$df_zero) + sum(bar_df_zero * d2p$d2f0) +
             sum(dbar_df_u_i    * pj$df_exo)  + sum(bar_df_u    * d2p$d2fu) +
             sum(bar_ys * d2p$d2ys[endo])
      if (n_state > 0)
        gij <- gij +
          sum(dbar_df_minus_S_i * pj$df_minus[, state_idx, drop = FALSE]) +
          sum(bar_df_minus_S * d2p$d2fm[, state_idx, drop = FALSE])

      T2[i, j] <- gij
    }
  }

  list(T2 = T2, ok = ok, used_analytic_first = !is.null(analytic),
       used_analytic_second = analytic_ok2)
}
