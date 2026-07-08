## R/gradient-solution-adjoint.R
## --------------------------------------------------------------------------
## ADJOINT (reverse-mode) of the first-order perturbation solve + steady-state
## fixed point (Tier 18 A2, phase 1).
##
## The forward implicit layer (gradient-solution-deriv.R) solves, for EVERY
## structural parameter j, one generalized Sylvester equation
##
##   A dG_j + f_plus dG_j SG = -RHS_j,
##   RHS_j = df_plus_j G SG + df_zero_j G + df_minus_S_j            (1)
##
## plus one A-solve for dH_j, then contracts the state-space blocks
## (dTT_j, dRR_j, dZZ_j, dDD_j, dd_j) with the Kalman-filter adjoint's bar
## matrices (G_TT, G_RR, G_ZZ, G_DD, g_d).  Cost: O(p) Sylvester solves.
##
## This file reverses that chain: given the bar matrices ONCE, TWO transposed
## solves produce the adjoints of the smooth primitives, and every parameter's
## gradient is then a Frobenius inner product against its (analytic) primitive
## derivatives -- no per-parameter solve at all.
##
## THE MATH (all inner products Frobenius, <X,Y> = sum(X*Y) = tr(X'Y))
## -------------------------------------------------------------------
## Reverse of the block extraction  dTT = dG[state_idx,], dZZ = dG[obs_idx,],
## dRR = dH[state_idx,], dDD = dH[obs_idx,], dd = dys[obs_vars]:
##
##   V_G[state_idx,] += G_TT,  V_G[obs_idx,] += G_ZZ     (n_endo x n_state)
##   V_H[state_idx,] += G_RR,  V_H[obs_idx,] += G_DD     (n_endo x n_exo)
##   bar_ys[obs_vars] += g_d
##
## Reverse of dH = -A^{-1}(dA H + df_u):  with W = A^{-T} V_H (ONE transposed
## A-solve),
##   bar_df_u = -W,        bar_dA = -W H'.
##
## Reverse of dA = df_plus (G S) + f_plus dG S + df_zero:
##   bar_df_plus += bar_dA (G S)',   bar_df_zero += bar_dA,
##   V_G         += f_plus' bar_dA S'.
##
## Reverse of the Sylvester solve (1):  the forward operator is
## M(X) = A X + f_plus X SG, whose adjoint under <.,.> is
## M*(Y) = A' Y + f_plus' Y SG'.  Solve ONE transposed Sylvester
##
##   A' L + f_plus' L SG' = V_G                                      (2)
##
## (same dedicated real-Schur solver, applied to (A', f_plus', SG')); then
## bar_RHS = -L and, reversing (1),
##
##   bar_df_plus  += bar_RHS SG' G',   bar_df_zero += bar_RHS G',
##   bar_df_minus_S = bar_RHS.
##
## Finally, for every parameter j,
##
##   grad_j = <bar_df_plus,  df_plus_j> + <bar_df_zero, df_zero_j>
##          + <bar_df_minus_S, df_minus_j[, state_idx]>
##          + <bar_df_u, df_u_j> + <bar_ys, dys_j>                   (3)
##
## where the TOTAL primitive derivatives (df_*_j, dys_j) come from the
## ANALYTIC layer (gradient-primitive-deriv.R: compile-time parameter
## Jacobian + model Hessian contracted with the implicit-function-theorem
## dys -- this is the steady-state fixed-point chain of A2) or, as fallback,
## from the same smooth-primitive central FD the forward layer uses.
##
## Consistency: (3) reproduces the forward implicit gradient EXACTLY (same
## primitives, same linear systems, transposed) -- pinned to ~1e-10 by
## test-gradient-solution-adjoint.R.  The Sigma_e / stderr-parameter channel
## is deliberately NOT handled here (it never touches the perturbation
## solve); make_posterior_grad routes it separately.
##
## Phase-1 scope note: second-order solution derivatives (posterior_hessian's
## T2; P2 gap #3) still use the forward layer + FD second primitives --
## the adjoint-of-adjoint and analytic d2f primitives are the remaining A2
## phase-2 items, along with the Rcpp port and the make_posterior_grad
## default flip.
## --------------------------------------------------------------------------


#' Reverse-mode structural-parameter gradient through the perturbation solve
#'
#' Given the Kalman-filter adjoint's bar matrices wrt the state-space system
#' (\code{G_TT}, \code{G_RR}, \code{G_ZZ}, \code{G_DD}, \code{g_d}), returns
#' the gradient of the log-likelihood wrt structural parameters by reversing
#' the first-order perturbation fixed point -- two transposed solves total,
#' then one Frobenius contraction per parameter, instead of the forward
#' layer's one generalized-Sylvester solve per parameter.
#'
#' @param model      dynhr_mod
#' @param compiled   dynhr_compiled
#' @param dr         DecisionRules from \code{solve_perturbation} (order 1)
#' @param params     Named numeric parameter vector (base point)
#' @param param_names Character: structural parameters to differentiate wrt
#' @param obs_vars   Character: observed variable names
#' @param bars       list(G_TT [n_state x n_state], G_RR [n_state x n_exo],
#'   G_ZZ [n_obs x n_state], G_DD [n_obs x n_exo], g_d [n_obs]) -- the
#'   adjoint-KF gradients wrt (TT, RR, ZZ, DD, d).  Any element may be NULL
#'   (treated as zero).
#' @param h_rel      FD step for the primitive fallback (matches the forward
#'   layer's default)
#' @param use_analytic NULL (option-driven, default TRUE) / TRUE / FALSE --
#'   same contract as \code{solution_derivatives}.
#' @return list(grad = named numeric (NA where a parameter's primitives
#'   failed), ok = named logical, used_analytic = logical)
#' @noRd
.solution_adjoint <- function(model, compiled, dr, params, param_names,
                              obs_vars, bars, h_rel = 1e-6,
                              use_analytic = NULL) {

  if (!all(param_names %in% names(params))) {
    missing_p <- param_names[!param_names %in% names(params)]
    stop(sprintf(".solution_adjoint: parameter(s) not found in `params`: %s",
                 paste(missing_p, collapse = ", ")))
  }

  endo      <- dr$endo_names
  n_endo    <- length(endo)
  n_exo     <- length(dr$exo_names)
  state_idx <- dr$state_idx
  n_state   <- length(state_idx)

  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop(".solution_adjoint: observed variables not found: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))
  n_obs <- length(obs_idx)

  G  <- dr$ghx
  H  <- dr$ghu
  ys <- dr$ys

  ## Same base-system preamble as solution_derivatives() (kept in lockstep;
  ## the parity test pins the two layers against each other).
  S <- matrix(0, nrow = n_state, ncol = n_endo)
  if (n_state > 0) S[cbind(seq_len(n_state), state_idx)] <- 1
  sys0 <- extract_system_matrices(compiled, ys, params)
  f_plus <- sys0$f_plus; f_zero <- sys0$f_zero
  f_minus <- sys0$f_minus; f_u <- sys0$f_exo
  SG <- if (n_state > 0) S %*% G else matrix(0, nrow = 0, ncol = 0)
  GS <- if (n_state > 0) G %*% S else matrix(0, nrow = n_endo, ncol = n_endo)
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

  ## ---- reverse of the block extraction ------------------------------------
  V_G <- matrix(0, n_endo, n_state)
  V_H <- matrix(0, n_endo, n_exo)
  if (n_state > 0) {
    V_G[state_idx, ] <- V_G[state_idx, , drop = FALSE] + G_TT
    V_H[state_idx, ] <- V_H[state_idx, , drop = FALSE] + G_RR
  }
  V_G[obs_idx, ] <- V_G[obs_idx, , drop = FALSE] + G_ZZ
  V_H[obs_idx, ] <- V_H[obs_idx, , drop = FALSE] + G_DD
  bar_ys <- setNames(numeric(n_endo), endo)
  bar_ys[obs_vars] <- bar_ys[obs_vars] + g_d

  ## ---- reverse of dH = -A^{-1}(dA H + df_u) --------------------------------
  At_qr <- qr(t(A))
  W <- qr.solve(At_qr, V_H)                      # A^{-T} V_H  [n_endo x n_exo]
  bar_df_u <- -W
  bar_dA   <- -W %*% t(H)                        # [n_endo x n_endo]

  ## ---- reverse of dA = df_plus GS + f_plus dG S + df_zero ------------------
  bar_df_plus <- bar_dA %*% t(GS)
  bar_df_zero <- bar_dA
  if (n_state > 0)
    V_G <- V_G + t(f_plus) %*% bar_dA %*% t(S)

  ## ---- reverse of the generalized Sylvester solve (one transposed solve) ---
  if (n_state > 0) {
    tA <- t(A); tfp <- t(f_plus); tSG <- t(SG)
    sylv_fac_T <- .gen_sylvester_k1_factor(tA, tfp, tSG)
    L <- if (!is.null(sylv_fac_T))
      .gen_sylvester_k1_solve(sylv_fac_T, V_G)
    else
      .solve_kron_compact(tA, tfp, tSG, k = 1L, RHS = V_G)
    ## residual guard: A' L + f_plus' L SG' must equal V_G
    resid <- max(abs(tA %*% L + tfp %*% L %*% tSG - V_G))
    if (!is.finite(resid) || resid > 1e-6 * max(1, max(abs(V_G))))
      stop(sprintf(".solution_adjoint: transposed Sylvester residual %.3e", resid))
    bar_RHS <- -L
    bar_df_plus <- bar_df_plus + bar_RHS %*% tSG %*% t(G)
    bar_df_zero <- bar_df_zero + bar_RHS %*% t(G)
    bar_df_minus_S <- bar_RHS
  } else {
    bar_df_minus_S <- matrix(0, n_endo, 0L)
  }

  ## ---- per-parameter primitive derivatives + contraction -------------------
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

  grad <- setNames(rep(NA_real_, length(param_names)), param_names)
  ok   <- setNames(rep(FALSE, length(param_names)), param_names)

  for (pname in param_names) {
    prim <- tryCatch({
      if (!is.null(analytic)) {
        c(list(dys = setNames(as.numeric(analytic$dys[, pname]),
                              rownames(analytic$dys))),
          analytic$dprim[[pname]])
      } else {
        .solution_adjoint_fd_prims(pname, model, compiled, ys, params, h_rel)
      }
    }, error = function(e) NULL)
    if (is.null(prim) || is.null(prim$df_plus)) next

    gj <- sum(bar_df_plus * prim$df_plus) +
          sum(bar_df_zero * prim$df_zero) +
          sum(bar_df_u    * prim$df_exo)  +
          sum(bar_ys      * prim$dys[endo])
    if (n_state > 0)
      gj <- gj + sum(bar_df_minus_S * prim$df_minus[, state_idx, drop = FALSE])

    grad[pname] <- gj
    ok[pname]   <- TRUE
  }

  list(grad = grad, ok = ok, used_analytic = !is.null(analytic))
}


#' FD fallback for the smooth-primitive TOTAL derivatives of one parameter.
#' Identical construction to .solution_deriv_one's FD branch (central FD of
#' the dynamic Jacobian at the perturbed steady state, NOT of the QZ solve).
#' @noRd
.solution_adjoint_fd_prims <- function(pname, model, compiled, ys, params,
                                       h_rel) {
  theta_j <- params[[pname]]
  h <- max(h_rel * abs(theta_j), 1e-7)
  params_p <- params; params_p[[pname]] <- theta_j + h
  params_m <- params; params_m[[pname]] <- theta_j - h
  ss_p <- solve_steady(compiled, params_p, y0 = ys,
                       endo_names = model$var_names,
                       exo_names = model$varexo_names, verbose = FALSE)
  if (!isTRUE(ss_p$converged))
    stop(sprintf("steady state did not converge at %s + h", pname))
  ss_m <- solve_steady(compiled, params_m, y0 = ys,
                       endo_names = model$var_names,
                       exo_names = model$varexo_names, verbose = FALSE)
  if (!isTRUE(ss_m$converged))
    stop(sprintf("steady state did not converge at %s - h", pname))
  ys_p <- ss_p$values; ys_m <- ss_m$values
  sys_p <- extract_system_matrices(compiled, ys_p,
                                   .ssm_consistent_params(model, params_p))
  sys_m <- extract_system_matrices(compiled, ys_m,
                                   .ssm_consistent_params(model, params_m))
  list(dys      = (ys_p - ys_m) / (2 * h),
       df_plus  = (sys_p$f_plus  - sys_m$f_plus)  / (2 * h),
       df_zero  = (sys_p$f_zero  - sys_m$f_zero)  / (2 * h),
       df_minus = (sys_p$f_minus - sys_m$f_minus) / (2 * h),
       df_exo   = (sys_p$f_exo   - sys_m$f_exo)   / (2 * h))
}
