## R/gradient-solution-deriv-order2.R
## --------------------------------------------------------------------------
## Implicit differentiation of the ORDER-2 perturbation solution.
##
## Computes d(ghxx)/dθ_j, d(ghxu)/dθ_j, d(ghuu)/dθ_j, d(ghss)/dθ_j,
## and d(Sigma_x)/dθ_j by differentiating the second-order Kronecker/linear
## systems analytically.  This mirrors the validated first-order Foundation A
## (gradient-solution-deriv.R) one order up.
##
## THE MATH
## --------
## The order-2 system (solve-perturbation-order2.R, lines ~500-580):
##
##   K_xx vec(ghxx) = -vec(Phi_xx),
##   K_xx = kron(I_{n_s^2}, A_L) + kron(hx' ⊗ hx', fp)
##   A_L  = f0 + fp * ghx * S'   (same as first-order A)
##   Phi_xx = H(T_x, T_x)        (model Hessian contracted with transfer matrices)
##
## Differentiating wrt θ_j:
##   K_xx vec(d ghxx_j) = -vec(dPhi_xx_j) - (dK_xx_j) vec(ghxx)
##
##   dK_xx_j = kron(I, dA_L_j) + kron(d(hx' ⊗ hx')_j, fp) + kron(hx' ⊗ hx', dfp_j)
##   d(hx' ⊗ hx')_j = kron(dhx_j', hx') + kron(hx', dhx_j')
##   dA_L_j = dfp_j * G * S' + fp * dG_j * S' + df0_j   (from first-order layer)
##
## dPhi_xx_j needs d(Phi_xx)/dθ_j. With Phi_xx[e] = vec(T_x' H[e] T_x):
##   dPhi_xx_j[e] = vec(dT_x_j' H[e] T_x + T_x' H[e] dT_x_j + T_x' dH[e]_j T_x)
##   dT_x_j = dT_x/dθ_j via chain rule through ghx/ghu (first-order derivs)
##   dH[e]_j = d(model Hessian)/dθ_j: central FD of hessian2_fn at θ±h
##             (the ONE residual primitive FD, same status as df_plus_j in Fdn A)
##
## For ghxu, ghuu: direct A_L solves given ghxx:
##   A_L * d(ghxu)_j = -(dPhi_xu_j + dfp_j * ghxx * (hu⊗hx) + fp * d(ghxx)_j * (hu⊗hx)
##                        + fp * ghxx * (dhu_j⊗hx + hu⊗dhx_j))   - dA_L_j * ghxu
##   Similarly for ghuu.
##
## For ghss: (A_L + fp) ghss = RHS_ss  =>
##   (A_L + fp) d(ghss)_j = dRHS_ss_j - d(A_L + fp)_j ghss
##
## For d(Sigma_x)/dθ_j: Sigma_x solves hx Sigma_x hx' + B - Sigma_x = 0
##   where B = hu Sigma_e hu'. Differentiating:
##   hx d(Sigma_x) hx' - d(Sigma_x) = -(dhx * Sigma_x * hx' + hx * Sigma_x * dhx'
##                                       + dhu * Sigma_e * hu' + hu * Sigma_e * dhu')
##   = another Lyapunov eq. (same hx, perturbed RHS).
##
## Factorization reuse: K_xx and A_L are built ONCE and all P parameterdirections
## share the same solves. The model-Hessian FD is the only per-parameter FD.
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## Helper: evaluate the model Hessian (as array n_eq x total_cols x total_cols)
## at a given parameter vector, using the symbolic hessian2_fn if available.
## Returns NULL on failure.
## ---------------------------------------------------------------------------
.o2sd_hessian_at <- function(compiled, params, ss) {
  dyn  <- compiled$dynamic
  dy   <- .build_dy_ss_o2(compiled, ss)

  ## Try symbolic first
  H <- .compute_model_hessian_symbolic(compiled, dy, params, ss)
  if (!is.null(H)) return(H)

  ## Fallback: numerical FD
  H <- tryCatch(.compute_model_hessian(dyn, dy, params, ss),
                error = function(e) NULL)
  H
}


## ---------------------------------------------------------------------------
## Helper: ANALYTIC total θ-derivative of the model Hessian, all parameters.
##
## Returns dH_compiled[e, c1, c2, k] = d H_e(c1,c2) / dθ_k evaluated at
## (dy(ȳ), θ, ȳ), in COMPILED column space and EQUATION-row order (i.e. the
## same layout as H_base BEFORE the equation->declaration row permutation):
##
##   dH[e,c1,c2,k] = ∂³F_e/(∂w_c1 ∂w_c2 ∂θ_k)                 (explicit = param_hessian2_fn)
##                 + Σ_{c3} ∂³F_e/(∂w_c1 ∂w_c2 ∂w_c3) · V[c3,k]   (ss chain, hess3 · V)
##
## with V[c3,k] = dys[var(c3),k] (0 for shock columns).  The chain term is the
## single-index contraction of the canonical, fully-symmetric hess3 tensor
## against the column-broadcast V; the orbit expansion matches `.contract_h3`
## exactly (contract the third slot, scatter into the (c1,c2) pair).
##
## Columns are indexed by the compile-time parameter order (compiled$model$
## param_names).  Returns NULL on failure so the caller falls back to FD.
## ---------------------------------------------------------------------------
.o2sd_analytic_dH <- function(compiled, ys, params, dys) {
  if (is.null(dys)) return(NULL)
  dyn <- compiled$dynamic
  if (is.null(dyn$param_hessian2_fn)) return(NULL)

  n_eq       <- dyn$n_eq
  total_cols <- dyn$total_cols
  np         <- ncol(dys)
  dy         <- .build_dy_ss_o2(compiled, ys)

  ## V[c,k] = dys[var(c),k] (0 for shock cols), aligned to compiled columns.
  dcm <- dyn$dyn_col_map
  col_var_idx <- match(dcm$name, compiled$model$var_names)  # NA for exo cols
  V <- matrix(0, total_cols, np)
  for (r in seq_len(nrow(dcm))) {
    vi <- col_var_idx[r]
    if (!is.na(vi)) V[dcm$col[r], ] <- dys[vi, ]
  }

  dH <- array(0, dim = c(n_eq, total_cols, total_cols, np))

  ## --- Explicit term: param_hessian2_fn sparse values (eq,c1,c2,param) ---
  ph2_vals <- tryCatch(dyn$param_hessian2_fn(dy, params, ys),
                       error = function(e) NULL)
  if (is.null(ph2_vals)) return(NULL)
  if (length(ph2_vals) && any(!is.finite(ph2_vals))) return(NULL)
  trip2 <- dyn$param_hess2_triplets
  for (j in seq_along(trip2)) {
    t <- trip2[[j]]; v <- ph2_vals[j]
    dH[t$eq, t$col1, t$col2, t$param] <- dH[t$eq, t$col1, t$col2, t$param] + v
    if (t$col1 != t$col2)
      dH[t$eq, t$col2, t$col1, t$param] <- dH[t$eq, t$col2, t$col1, t$param] + v
  }

  ## --- SSM-computed-parameter chain (Tier 12 #2): the total θ-derivative of the
  ## model Hessian also moves p_c. Capture the PURE ∂³F/∂w²∂p_c columns now (just
  ## the explicit param_hessian2 term, before the hess3 ȳ-chain), to fold into the
  ## free columns after, scaled by the total dp_c/dθ_k. ---
  ssm_dpc <- NULL; dH_pc <- NULL
  if (.ssm_assigns_param(compiled$model)) {
    ch1 <- .ssm_param_chain_derivs(compiled, ys, params)
    if (is.null(ch1)) return(NULL)
    endoN <- compiled$model$var_names; parsN <- compiled$model$param_names
    comp <- ch1$computed; freeN <- setdiff(parsN, comp)
    dH_pc <- dH[, , , match(comp, parsN), drop = FALSE]    # pure ph2 p_c cols (copy)
    ssm_dpc <- matrix(0, length(comp), np, dimnames = list(comp, parsN))
    for (pc in comp) for (k in freeN)
      ssm_dpc[pc, k] <- ch1$dg_dtheta[pc, k] + sum(ch1$dg_dy[pc, endoN] * dys[endoN, k])
  }

  ## --- Steady-state chain term: Σ_{c3} hess3[e,c1,c2,c3] · V[c3,k] ---
  ## hess3 stored canonical c1<=c2<=c3 but fully symmetric: expand the orbit
  ## with `.orbit_3`, contract the THIRD slot of each placement against V and
  ## scatter into the (first,second) pair -- matches `.contract_h3(h3,I,I,V)`.
  h3vals <- tryCatch(dyn$hessian3_fn(dy, params, ys),
                     error = function(e) NULL)
  if (is.null(h3vals)) return(NULL)
  if (length(h3vals) && any(!is.finite(h3vals))) return(NULL)
  trip3 <- dyn$hess3_triplets
  for (j in seq_along(trip3)) {
    val <- h3vals[j]
    if (val == 0) next
    t <- trip3[[j]]; e <- t$eq
    orbit <- .orbit_3(t$col1, t$col2, t$col3)
    for (perm in orbit) {
      a <- perm[1]; b <- perm[2]; c3 <- perm[3]
      dH[e, a, b, ] <- dH[e, a, b, ] + val * V[c3, ]
    }
  }

  ## --- SSM p_c channel fold-in: dH[,,,k] += (dp_c/dθ_k) · (∂³F/∂w²∂p_c) ---
  if (!is.null(ssm_dpc)) {
    parsN <- compiled$model$param_names
    comp  <- rownames(ssm_dpc); freeN <- setdiff(parsN, comp)
    for (jc in seq_along(comp)) for (k in freeN) {
      w <- ssm_dpc[comp[jc], k]
      if (w != 0) {
        ki <- match(k, parsN)
        dH[, , , ki] <- dH[, , , ki] + w * dH_pc[, , , jc]
      }
    }
  }

  if (any(!is.finite(dH))) return(NULL)
  dimnames(dH) <- list(NULL, NULL, NULL, compiled$model$param_names)
  dH
}


## ---------------------------------------------------------------------------
## Helper: compute Phi_xx/Phi_xu/Phi_uu given H, T_x, T_u, n_eq
## (thin wrapper around the existing .compute_phi_matrices)
## ---------------------------------------------------------------------------
.o2sd_phi <- function(H, T_x, T_u, n_eq) {
  .compute_phi_matrices(H, T_x, T_u, n_eq)
}


## ---------------------------------------------------------------------------
## Helper: d(Sigma_x)/dθ_j via a Lyapunov solve.
##   hx dX hx' - dX = -(dhx * Sx * hx' + hx * Sx * dhx' + dhu * Se * hu' + hu * Se * dhu')
## ---------------------------------------------------------------------------
.o2sd_dSigma_x <- function(hx, hu, dhx, dhu, Sigma_x, Sigma_e) {
  dB <- dhx %*% Sigma_x %*% t(hx) +
        hx  %*% Sigma_x %*% t(dhx) +
        dhu %*% Sigma_e %*% t(hu)  +
        hu  %*% Sigma_e %*% t(dhu)
  ## solve hx X hx' - X = -dB  =>  X = solve_lyapunov(hx, dB)
  solve_lyapunov(hx, dB)
}


#' Implicit-differentiation derivatives of the order-2 decision rule
#'
#' Computes \code{d(ghxx)/dθ_j}, \code{d(ghxu)/dθ_j}, \code{d(ghuu)/dθ_j},
#' \code{d(ghss)/dθ_j}, and \code{d(Sigma_x)/dθ_j} for each requested
#' parameter, via implicit differentiation of the second-order perturbation
#' systems.  The coefficient matrices \code{K_xx} and \code{A_L} are
#' factorized once and reused for all P parameters.
#'
#' The only per-parameter FD is of the \emph{model Hessian} at \code{θ ± h}
#' (needed for \code{dPhi_xx_j}).  All other quantities come from the
#' first-order layer via \code{solution_derivatives()}.
#'
#' @param model      dynhr_mod from \code{parse_mod()}.
#' @param compiled   dynhr_compiled with \code{max_order >= 2}.
#' @param dr2        DecisionRules2 (order 2) at the base parameter point.
#' @param params     Named numeric parameter vector (the base point).
#' @param param_names Character: parameters to differentiate wrt.
#' @param h_rel      Relative FD step for the first-order layer (default 1e-6).
#' @param h_hess     Relative FD step for the model-Hessian perturb (default 1e-4).
#' @return A list with:
#'   \item{first}{output of \code{solution_derivatives()} (per-parameter dG/dH/...)}
#'   \item{derivs}{named list (by param_names) of order-2 derivatives, each with
#'     \code{d_ghxx}, \code{d_ghxu}, \code{d_ghuu}, \code{d_ghss}, \code{d_Sigma_x}, \code{ok}}
#'   \item{param_names}{echoed}
#' @noRd
solution_derivatives_order2 <- function(model, compiled, dr2, params, param_names,
                                         h_rel   = 1e-6,
                                         h_hess  = 1e-4) {

  if (!all(param_names %in% names(params)))
    stop("solution_derivatives_order2: unknown parameter(s): ",
         paste(setdiff(param_names, names(params)), collapse = ", "))
  if (!inherits(dr2, "DecisionRules2"))
    stop("solution_derivatives_order2: dr2 must be a DecisionRules2 object")

  ## -----------------------------------------------------------------------
  ## Dimensions and base quantities
  ## -----------------------------------------------------------------------
  endo_names <- dr2$endo_names
  exo_names  <- dr2$exo_names
  state_idx  <- dr2$state_idx
  n          <- length(endo_names)
  n_s        <- length(state_idx)
  n_u        <- length(exo_names)
  np         <- length(param_names)

  ghx    <- dr2$ghx        # n x n_s
  ghu    <- dr2$ghu        # n x n_u
  ghxx   <- dr2$ghxx       # n x n_s^2
  ghxu   <- dr2$ghxu       # n x n_s*n_u
  ghuu   <- dr2$ghuu       # n x n_u^2
  ghss   <- dr2$ghss       # n
  Sigma_e <- dr2$Sigma_e   # n_u x n_u
  ys     <- dr2$ys

  hx  <- ghx[state_idx, , drop = FALSE]   # n_s x n_s
  hu  <- ghu[state_idx, , drop = FALSE]   # n_s x n_u

  ## Selection S: n_s x n  (rows are state-var unit vectors)
  S <- matrix(0, n_s, n)
  if (n_s > 0) S[cbind(seq_len(n_s), state_idx)] <- 1

  ## -----------------------------------------------------------------------
  ## Base system matrices and A_L
  ## -----------------------------------------------------------------------
  sys0 <- extract_system_matrices(compiled, ys, params)
  f0   <- sys0$f_zero
  fp   <- sys0$f_plus
  fm   <- sys0$f_minus
  fu   <- sys0$f_exo

  ## A_L = f0 + fp * ghx * S'  (n x n), using S: n_s x n so S' is n x n_s
  ## and ghx * S' is n x n
  GS  <- ghx %*% S           # n x n  (ghx * S')
  A_L <- f0 + fp %*% GS      # n x n

  ## -----------------------------------------------------------------------
  ## K_xx factorization (built ONCE)
  ## -----------------------------------------------------------------------
  ns2   <- n_s * n_s
  hxt   <- t(hx)                                              # n_s x n_s
  K_xx  <- kronecker(diag(ns2), A_L) + kronecker(hxt %x% hxt, fp)
  K_xx_qr <- qr(K_xx)

  ## A_L factorization (reused for ghxu, ghuu)
  AL_qr <- qr(A_L)

  ## (A_L + fp) factorization (reused for ghss)
  AB_qr <- qr(A_L + fp)

  ## -----------------------------------------------------------------------
  ## Transfer matrices T_x, T_u at base (needed for Phi and dPhi)
  ## -----------------------------------------------------------------------
  tm  <- .build_transfer_matrices(compiled$dynamic, ghx, ghu, state_idx,
                                   hx, hu, endo_names, exo_names)
  T_x <- tm$T_x   # total_cols x n_s
  T_u <- tm$T_u   # total_cols x n_u

  ## Base model Hessian H (n x total_cols x total_cols) — built ONCE.
  ## We also need it for dRHS_ghss later.
  H_base <- .o2sd_hessian_at(compiled, params, ys)
  if (is.null(H_base))
    stop("solution_derivatives_order2: could not compute model Hessian at base")

  ## Reorder Hessian rows from compiled to declaration order (mirrors order2 solver)
  if (!is.null(compiled$model$equations)) {
    eq_to_decl <- .build_eq_to_decl(compiled$model)
    if (all(eq_to_decl > 0L) && !identical(eq_to_decl, seq_len(n))) {
      perm <- order(eq_to_decl)
      H_perm <- array(0, dim = dim(H_base))
      for (k in seq_len(n)) H_perm[k, , ] <- H_base[perm[k], , ]
      H_base <- H_perm
    }
  }

  ## State covariance at base
  Sigma_x <- .state_covariance(hx, hu, Sigma_e)   # n_s x n_s

  ## T_up for ghss derivative (forward-looking block)
  has_lead <- sys0$is_fwd | sys0$is_mixed
  T_up <- .build_T_up(compiled$dynamic, ghu, endo_names, exo_names, has_lead)

  ## -----------------------------------------------------------------------
  ## First-order solution derivatives for ALL requested parameters.
  ## We need dG_j, dH_j, dfp_j, df0_j, dfm_j, dfu_j for each j, plus the total
  ## θ-derivative of the model Hessian dH_mat (for dPhi_xx / dRHS_ss).
  ##
  ## ANALYTIC path (Tier 11 #3): when the model is parameter-differentiable and
  ## the order-2 parameter tensors are compiled, the primitive derivatives and
  ## dH_mat come from the symbolic param-Jacobian / param-Hessian2 + hess3·dys
  ## chain -- no per-parameter steady-state re-solve, no central FD of the model
  ## Hessian.  Falls back to the FD path (per-param solve_steady + central FD)
  ## otherwise, or if any analytic builder returns NULL.
  ## -----------------------------------------------------------------------
  use_analytic <- .can_use_analytic_primitive_deriv(compiled) &&
    isTRUE(compiled$dynamic$param_hess2_built) &&
    isTRUE(compiled$dynamic$hessian3_built)
  ## SSM-computed-parameter models now use the order-2 augmented chain
  ## (.analytic_dprimitives channel 3 + .o2sd_analytic_dH p_c fold-in).

  dprim_an <- NULL    # analytic primitive derivs, list by param_name
  dH_an    <- NULL    # analytic dH_mat, list by param_name (row-permuted)
  if (use_analytic) {
    dys_all <- .analytic_dys(compiled, ys, params)
    dprim_all <- if (!is.null(dys_all))
      .analytic_dprimitives(compiled, ys, params, dys_all) else NULL
    dH_all <- if (!is.null(dys_all))
      .o2sd_analytic_dH(compiled, ys, params, dys_all) else NULL
    if (is.null(dys_all) || is.null(dprim_all) || is.null(dH_all)) {
      use_analytic <- FALSE   # graceful fall back to FD
    } else {
      ## Row permutation applied to H_base (lines below); apply the SAME to the
      ## analytic dH so its rows are in declaration order to match H_base.
      perm_dh <- NULL
      if (!is.null(compiled$model$equations)) {
        e2d <- .build_eq_to_decl(compiled$model)
        if (all(e2d > 0L) && !identical(e2d, seq_len(n))) perm_dh <- order(e2d)
      }
      dprim_an <- vector("list", np); names(dprim_an) <- param_names
      dH_an    <- vector("list", np); names(dH_an)    <- param_names
      for (k in seq_len(np)) {
        pnm <- param_names[k]
        dprim_an[[pnm]] <- dprim_all[[pnm]]
        dHk <- dH_all[, , , pnm]                 # n_eq x total_cols x total_cols
        if (!is.null(perm_dh)) {
          dHk_perm <- array(0, dim = dim(dHk))
          for (e in seq_len(n)) dHk_perm[e, , ] <- dHk[perm_dh[e], , ]
          dHk <- dHk_perm
        }
        dH_an[[pnm]] <- dHk
      }
    }
  }

  ## Use a private FD step (same h_rel) to get the primitive derivatives.
  hvec <- vapply(param_names, function(p) max(h_rel * abs(params[[p]]), 1e-7), 0)
  names(hvec) <- param_names

  ## Compute perturbed system matrices for each param (central FD) -- only the
  ## FD path needs the per-parameter steady-state re-solves.
  prim_p <- prim_m <- vector("list", np)
  names(prim_p) <- names(prim_m) <- param_names
  if (!use_analytic) for (k in seq_len(np)) {
    pnm <- param_names[k]; h <- hvec[k]
    tp  <- params; tp[[pnm]] <- tp[[pnm]] + h
    tm_ <- params; tm_[[pnm]] <- tm_[[pnm]] - h
    ## Perturbed steady states (warm-started)
    ss_p <- tryCatch(
      solve_steady(compiled, tp, y0 = ys, endo_names = model$var_names,
                   exo_names = model$varexo_names, verbose = FALSE),
      error = function(e) list(converged = FALSE))
    ss_m <- tryCatch(
      solve_steady(compiled, tm_, y0 = ys, endo_names = model$var_names,
                   exo_names = model$varexo_names, verbose = FALSE),
      error = function(e) list(converged = FALSE))
    if (!isTRUE(ss_p$converged) || !isTRUE(ss_m$converged)) {
      prim_p[[k]] <- NULL; prim_m[[k]] <- NULL; next
    }
    ## Re-derive SSM-computed parameters so the system matrices use the
    ## consistent (not stale) p_c (no-op for non-SSM-parameter models).
    sys_p <- extract_system_matrices(compiled, ss_p$values, .ssm_consistent_params(model, tp))
    sys_m <- extract_system_matrices(compiled, ss_m$values, .ssm_consistent_params(model, tm_))
    prim_p[[k]] <- c(sys_p, list(ys = ss_p$values))
    prim_m[[k]] <- c(sys_m, list(ys = ss_m$values))
  }

  ## Also call solution_derivatives() for dG, dH per parameter.
  ## We use obs_vars = endo_names (all) so it doesn't error; we only use dG, dH.
  first <- solution_derivatives(model, compiled, dr2, params, param_names,
                                 obs_vars = endo_names, h_rel = h_rel)

  ## -----------------------------------------------------------------------
  ## Per-parameter implicit differentiation
  ## -----------------------------------------------------------------------
  derivs <- vector("list", np)
  names(derivs) <- param_names

  for (k in seq_len(np)) {
    pnm <- param_names[k]
    h   <- hvec[k]

    if (!use_analytic && (is.null(prim_p[[k]]) || is.null(prim_m[[k]]))) {
      warning(sprintf("solution_derivatives_order2: steady state failed at %s; NA", pnm))
      derivs[[pnm]] <- list(
        d_ghxx = matrix(NA_real_, n, ns2),
        d_ghxu = matrix(NA_real_, n, n_s * n_u),
        d_ghuu = matrix(NA_real_, n, n_u^2),
        d_ghss = rep(NA_real_, n),
        d_Sigma_x = matrix(NA_real_, n_s, n_s),
        ok = FALSE
      )
      next
    }

    ## ---- Primitive derivatives: analytic (Tier 11 #3) or central FD ----
    if (use_analytic) {
      dp_an <- dprim_an[[pnm]]
      dfp <- dp_an$df_plus
      df0 <- dp_an$df_zero
      dfm <- dp_an$df_minus
      dfu <- dp_an$df_exo
    } else {
      dfp <- (prim_p[[k]]$f_plus  - prim_m[[k]]$f_plus)  / (2 * h)
      df0 <- (prim_p[[k]]$f_zero  - prim_m[[k]]$f_zero)  / (2 * h)
      dfm <- (prim_p[[k]]$f_minus - prim_m[[k]]$f_minus) / (2 * h)
      dfu <- (prim_p[[k]]$f_exo   - prim_m[[k]]$f_exo)   / (2 * h)
    }

    ## ---- dG_j, dH_j from first-order layer ----
    d1 <- first$derivs[[pnm]]
    if (!isTRUE(d1$ok)) {
      warning(sprintf("solution_derivatives_order2: first-order deriv failed at %s; NA", pnm))
      derivs[[pnm]] <- list(
        d_ghxx = matrix(NA_real_, n, ns2),
        d_ghxu = matrix(NA_real_, n, n_s * n_u),
        d_ghuu = matrix(NA_real_, n, n_u^2),
        d_ghss = rep(NA_real_, n),
        d_Sigma_x = matrix(NA_real_, n_s, n_s),
        ok = FALSE
      )
      next
    }
    dG  <- d1$dG    # n x n_s
    dH  <- d1$dH    # n x n_u
    dhx <- dG[state_idx, , drop = FALSE]   # n_s x n_s
    dhu <- dH[state_idx, , drop = FALSE]   # n_s x n_u

    ## ---- dA_L_j = dfp * G * S' + fp * dG * S' + df0 ----
    dGS  <- dG %*% S   # n x n  (dG * S')
    dA_L <- dfp %*% GS + fp %*% dGS + df0    # n x n

    ## ---- d(hx' ⊗ hx')_j  (note: hxt = t(hx), dhxt = t(dhx)) ----
    dhxt <- t(dhx)   # n_s x n_s
    ## d(hxt ⊗ hxt)/dθ = dhxt ⊗ hxt + hxt ⊗ dhxt
    d_hxtkron <- kronecker(dhxt, hxt) + kronecker(hxt, dhxt)   # n_s^2 x n_s^2

    ## ---- dK_xx_j ----
    dK_xx <- kronecker(diag(ns2), dA_L) +
             kronecker(d_hxtkron, fp)   +
             kronecker(hxt %x% hxt, dfp)   # n*n_s^2 x n*n_s^2

    ## ---- d(model Hessian)/dθ_j: analytic (Tier 11 #3) or central FD ----
    ## dH_mat[e,c1,c2] = ∂³F_e/(∂w_c1 ∂w_c2 ∂θ) + Σ_c3 hess3[e,c1,c2,c3] dys[c3].
    ## The analytic tensor is already row-permuted (declaration order) to match
    ## H_base.  FD fallback: central difference of the model Hessian at θ±h.
    if (use_analytic) {
      dH_mat <- dH_an[[pnm]]
    } else {
      ys_p   <- prim_p[[k]]$ys
      ys_m   <- prim_m[[k]]$ys
      ## Re-derive SSM-computed parameters so the model Hessian uses the
      ## consistent (not stale) p_c (no-op for non-SSM-parameter models).
      params_p <- params; params_p[[pnm]] <- params[[pnm]] + h
      params_m <- params; params_m[[pnm]] <- params[[pnm]] - h
      params_p <- .ssm_consistent_params(model, params_p)
      params_m <- .ssm_consistent_params(model, params_m)

      H_p <- .o2sd_hessian_at(compiled, params_p, ys_p)
      H_m <- .o2sd_hessian_at(compiled, params_m, ys_m)
      if (is.null(H_p) || is.null(H_m)) {
        warning(sprintf("solution_derivatives_order2: Hessian FD failed at %s", pnm))
        derivs[[pnm]] <- list(
          d_ghxx = matrix(NA_real_, n, ns2),
          d_ghxu = matrix(NA_real_, n, n_s * n_u),
          d_ghuu = matrix(NA_real_, n, n_u^2),
          d_ghss = rep(NA_real_, n),
          d_Sigma_x = matrix(NA_real_, n_s, n_s),
          ok = FALSE
        )
        next
      }

      ## Reorder H_p and H_m (same permutation as H_base)
      if (!is.null(compiled$model$equations)) {
        eq_to_decl <- .build_eq_to_decl(compiled$model)
        if (all(eq_to_decl > 0L) && !identical(eq_to_decl, seq_len(n))) {
          perm <- order(eq_to_decl)
          Hp_perm <- array(0, dim = dim(H_p)); Hm_perm <- array(0, dim = dim(H_m))
          for (e in seq_len(n)) {
            Hp_perm[e, , ] <- H_p[perm[e], , ]
            Hm_perm[e, , ] <- H_m[perm[e], , ]
          }
          H_p <- Hp_perm; H_m <- Hm_perm
        }
      }

      dH_mat <- (H_p - H_m) / (2 * h)   # n x total_cols x total_cols
    }

    ## ---- dT_x_j: perturb transfer matrix through ghx/ghu ----
    ## The DCM exo names we need to skip in the dT_x helper:
    exo_set <- exo_names
    ## dT_x[c, s]:
    ##   ll == 0:  dG_j[var, ]
    ##   ll == 1:  d(G * hx)/dθ_j [var, ] = (dG %*% hx + G %*% dhx)[var, ]
    dcm        <- compiled$dynamic$dyn_col_map
    total_cols <- compiled$dynamic$total_cols

    dT_x <- matrix(0, total_cols, n_s)
    dT_u <- matrix(0, total_cols, n_u)   # dT_u needed for dPhi_xu, dPhi_uu

    G_hx  <- ghx %*% hx   # n x n_s
    G_hu  <- ghx %*% hu   # n x n_u  (lead response to shock through state)
    dG_hx <- dG %*% hx + ghx %*% dhx   # n x n_s
    dG_hu <- dG %*% hu + ghx %*% dhu   # n x n_u

    for (kk in seq_len(nrow(dcm))) {
      c_col <- dcm$col[kk]
      nm    <- dcm$name[kk]
      ll    <- dcm$lead_lag[kk]
      if (nm %in% exo_set) {
        ## Exo columns: T_u[c_col, k_exo] = 1; dT_u stays 0 (no param dependence)
        next
      }
      j <- which(endo_names == nm)
      if (length(j) != 1L) next
      if (ll == 0L) {
        dT_x[c_col, ] <- dG[j, ]
        dT_u[c_col, ] <- dH[j, ]
      } else if (ll == 1L) {
        dT_x[c_col, ] <- dG_hx[j, ]
        dT_u[c_col, ] <- dG_hu[j, ]
      }
      ## ll == -1: constant selector, derivative 0
    }

    ## ---- dPhi_xx_j ----
    ##   dPhi_xx[e] = vec(dT_x' H[e] T_x + T_x' H[e] dT_x + T_x' dH[e] T_x)
    dPhi_xx <- matrix(0, n, ns2)
    dPhi_xu <- matrix(0, n, n_s * n_u)
    dPhi_uu <- matrix(0, n, n_u * n_u)
    Tx_t    <- t(T_x)
    Tu_t    <- t(T_u)
    dTx_t   <- t(dT_x)
    dTu_t   <- t(dT_u)
    for (e in seq_len(n)) {
      He  <- H_base[e, , ]
      dHe <- dH_mat[e, , ]
      ## Phi_xx
      if (n_s > 0L) {
        dPhi_xx[e, ] <- as.vector(
          dTx_t %*% He  %*% T_x  +
          Tx_t  %*% He  %*% dT_x +
          Tx_t  %*% dHe %*% T_x
        )
        ## Phi_xu
        if (n_u > 0L) {
          dPhi_xu[e, ] <- as.vector(
            dTx_t %*% He  %*% T_u  +
            Tx_t  %*% He  %*% dT_u +
            Tx_t  %*% dHe %*% T_u
          )
        }
      }
      ## Phi_uu
      if (n_u > 0L) {
        dPhi_uu[e, ] <- as.vector(
          dTu_t %*% He  %*% T_u  +
          Tu_t  %*% He  %*% dT_u +
          Tu_t  %*% dHe %*% T_u
        )
      }
    }

    ## ---- d(ghxx)_j: solve K_xx vec(d_ghxx) = -vec(dPhi_xx) - dK_xx vec(ghxx) ----
    rhs_ghxx <- -as.vector(dPhi_xx) - dK_xx %*% as.vector(ghxx)
    d_ghxx_vec <- qr.solve(K_xx_qr, rhs_ghxx)
    d_ghxx <- matrix(d_ghxx_vec, n, ns2)

    ## ---- d(ghxu)_j ----
    ## A_L d(ghxu) = -(dPhi_xu + dfp * ghxx * (hu⊗hx) + fp * d_ghxx * (hu⊗hx)
    ##                + fp * ghxx * (dhu⊗hx + hu⊗dhx)) - dA_L * ghxu
    rhs_xu <- -(dPhi_xu +
                dfp %*% ghxx %*% (hu %x% hx) +
                fp  %*% d_ghxx %*% (hu %x% hx) +
                fp  %*% ghxx %*% (dhu %x% hx + hu %x% dhx)) -
               dA_L %*% ghxu
    d_ghxu <- qr.solve(AL_qr, rhs_xu)

    ## ---- d(ghuu)_j ----
    ## A_L d(ghuu) = -(dPhi_uu + dfp * ghxx * (hu⊗hu) + fp * d_ghxx * (hu⊗hu)
    ##                + fp * ghxx * (dhu⊗hu + hu⊗dhu)) - dA_L * ghuu
    rhs_uu <- -(dPhi_uu +
                dfp %*% ghxx %*% (hu %x% hu) +
                fp  %*% d_ghxx %*% (hu %x% hu) +
                fp  %*% ghxx %*% (dhu %x% hu + hu %x% dhu)) -
               dA_L %*% ghuu
    d_ghuu <- qr.solve(AL_qr, rhs_uu)

    ## ---- d(ghss)_j ----
    ## ghss solves: (A_L + fp) ghss = RHS_ss
    ##   RHS_ss = -(fp * ghuu * vec(Sigma_e) + H2(T_up, T_up) * vec(Sigma_e))
    ## Differentiating: (A_L + fp) d(ghss) = dRHS_ss - d(A_L + fp) * ghss
    ##   d(A_L + fp) = dA_L + dfp
    ##   dRHS_ss = -(dfp * ghuu * vSe + fp * d_ghuu * vSe
    ##               + d(.bilinear_h2(H_base, T_up, T_up)) * vSe)
    ## The third term requires d(H2(T_up, T_up))/dθ which involves dH_mat and dT_up.
    ## T_up is the "jumper" block — its derivative through ghu is analogous to dT_u.
    ## We build dT_up similarly to dT_u for the current-period ghu response.
    dT_up <- .o2sd_dT_up(compiled$dynamic, dH, ghu, endo_names,
                          exo_names, has_lead)

    vSe <- as.numeric(Sigma_e)
    ## d(fp * ghuu * vSe)/dθ = dfp * ghuu * vSe + fp * d_ghuu * vSe
    term1 <- dfp %*% ghuu %*% vSe + fp %*% d_ghuu %*% vSe

    ## d(.bilinear_h2(H_base, T_up, T_up) * vSe)/dθ
    ## = .bilinear_h2(dH_mat, T_up, T_up) * vSe + 2 * .bilinear_h2(H_base, dT_up, T_up) * vSe
    term2 <- (.bilinear_h2(dH_mat, T_up, T_up)  +
              2 * .bilinear_h2(H_base, dT_up, T_up)) %*% vSe

    dRHS_ss <- -(term1 + term2)

    d_AB_ghss <- (dA_L + dfp) %*% ghss
    d_ghss_vec <- qr.solve(AB_qr, dRHS_ss - d_AB_ghss)
    d_ghss <- as.numeric(d_ghss_vec)

    ## ---- d(Sigma_x)_j ----
    d_Sigma_x <- .o2sd_dSigma_x(hx, hu, dhx, dhu, Sigma_x, Sigma_e)

    ## Name the outputs consistently with dr2 conventions
    names(d_ghss) <- endo_names
    rownames(d_ghxx) <- endo_names
    rownames(d_ghxu) <- endo_names
    rownames(d_ghuu) <- endo_names

    derivs[[pnm]] <- list(
      d_ghxx    = d_ghxx,
      d_ghxu    = d_ghxu,
      d_ghuu    = d_ghuu,
      d_ghss    = d_ghss,
      d_Sigma_x = d_Sigma_x,
      ok        = TRUE
    )
  }

  list(first = first, derivs = derivs, param_names = param_names)
}


## ---------------------------------------------------------------------------
## Helper: build dT_up/dθ_j (derivative of the "jumper" transfer matrix).
##
## T_up[c, k] = ghu[j, k] for jumper compound-columns c (forward-looking vars).
## dT_up[c, k] = dH_j[j, k]  (from the first-order layer dH = d(ghu)/dθ_j).
## ---------------------------------------------------------------------------
.o2sd_dT_up <- function(dyn, dH, ghu, endo_names, exo_names, has_lead) {
  total_cols <- dyn$total_cols
  n_u        <- length(exo_names)
  dT_up      <- matrix(0, total_cols, n_u)

  jm <- .identify_jumpers(dyn, endo_names, has_lead)
  for (i in seq_along(jm$jumper_idx)) {
    j <- jm$jumper_idx[i]
    c <- jm$jumper_compound_c[i]
    dT_up[c, ] <- dH[j, ]
  }
  dT_up
}
