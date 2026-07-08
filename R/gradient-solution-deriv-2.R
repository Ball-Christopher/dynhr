## R/gradient-solution-deriv-2.R
## --------------------------------------------------------------------------
## SECOND-ORDER solution-derivative layer (Foundation A for the exact
## posterior Hessian, ROADMAP Tier 6 #2).
##
## Computes the second total derivatives of the first-order decision rule
##   d2G/dtheta_i dtheta_j,  d2H/dtheta_i dtheta_j,  d2ys/dtheta_i dtheta_j
## (and the state-space blocks d2TT, d2RR, d2ZZ, d2DD, d2d) by differentiating
## the implicit first-order system a SECOND time. The first-order layer
## (gradient-solution-deriv.R) solves, for each parameter j,
##
##   M vec(dG_j) = -vec(RHS_j),
##     M      = (S G)' (x) f_plus + I_{n_state} (x) A,   A = f_plus G S + f_zero
##     RHS_j  = df_plus_j G (S G) + df_zero_j G + df_minus_S_j
##
## Differentiating the fully-expanded first total derivative once more wrt
## theta_i, the SAME coefficient matrix M reappears on the second-derivative
## unknown:
##
##   M vec(d2G_ij) = -vec(K_ij)
##
## where K_ij collects every term that does NOT multiply d2G_ij:
##
##   K_ij =  d2f_plus_ij G (SG) + df_plus_j G_i (SG) + df_plus_j G (S G_i)
##         + df_plus_i G_j (SG) + f_plus G_j (S G_i)
##         + df_plus_i G  (S G_j) + f_plus G_i (S G_j)
##         + d2f_zero_ij G + df_zero_j G_i
##         + df_zero_i G_j
##         + d2f_minus_S_ij
##
## (G_i = dG/dtheta_i etc., from the validated first-order layer; SG = S G,
## S G_i = S dG_i.) Then, from A H = -f_u differentiated twice,
##
##   A H_ij = -( d2f_u_ij + A_ij H + A_i H_j + A_j H_i )
##   A_j    = df_plus_j G S + f_plus G_j S + df_zero_j
##   A_ij   = d2f_plus_ij G S + df_plus_i G_j S + df_plus_j G_i S
##            + f_plus G_ij S + d2f_zero_ij
##
## reusing the SAME QR factorizations of M and A built by the first-order
## layer. The only new numerical inputs are the TOTAL second derivatives of
## the smooth model primitives (d2f_*_ij, d2ys_ij), obtained by finite-
## difference stencils of extract_system_matrices / solve_steady -- exactly
## the smooth-primitive FD the first-order layer already uses for df_*_j,
## just one order higher. The implicit-function structure is otherwise exact.
## --------------------------------------------------------------------------


## Extract the TOTAL model primitives at a parameter point: re-solve the
## steady state (warm-started from ys0) and evaluate the dynamic Jacobian
## blocks there. Returns NULL if the steady state fails to converge.
.sd2_primitives_at <- function(theta_vec, model, compiled, ys0) {
  ss <- solve_steady(compiled, theta_vec, y0 = ys0,
                     endo_names = model$var_names,
                     exo_names = model$varexo_names, verbose = FALSE)
  if (!isTRUE(ss$converged)) return(NULL)
  ## Re-derive steady_state_model-computed parameters so the system matrices use
  ## the consistent (not stale) p_c (no-op for non-SSM-parameter models).
  sys <- extract_system_matrices(compiled, ss$values, .ssm_consistent_params(model, theta_vec))
  list(ys = ss$values, f_plus = sys$f_plus, f_zero = sys$f_zero,
       f_minus = sys$f_minus, f_exo = sys$f_exo)
}


## ---------------------------------------------------------------------------
## Analytic second-total primitive context + per-pair builder (Tier 11 #3, 3b).
##
## .sd2_build_ctx evaluates, ONCE at the base point, the four sparse symbolic
## tensors the analytic second total derivative of the dynamic Jacobian needs:
##   param2_jacobian_fn  d3F/(dw_c dtheta_a dtheta_b)   [explicit, T1]
##   param_hessian2_fn   d3F/(dw_c1 dw_c2 dtheta_k)     [T2/T3]
##   hessian2_fn         d2F/(dw_c1 dw_c2)              [T4, contract with d2ys]
##   hessian3_fn         d3F/(dw_c1 dw_c2 dw_c3)        [T5, contract with dys (x)dys]
## plus the column-broadcast steady-state sensitivities V[c,k]=dys[var(c),k].
##
## .sd2_d2J_pair(ctx,a,b) assembles the n_eq x total_cols second total
## derivative d2(dF/dw)/(dtheta_a dtheta_b) (compiled equation order) via the
## 5-term formula and partitions it (shared .partition_dJ) into the
## df_plus/df_zero/df_minus/df_exo blocks -- the analytic replacement for the
## FD diagonal/4-corner stencils. The hess3 orbit handling matches the verified
## .o2sd_analytic_dH / .contract_h3 convention.
## ---------------------------------------------------------------------------
.sd2_build_ctx <- function(compiled, ys, params, dys, d2ys) {
  dyn <- compiled$dynamic
  if (is.null(dyn$param_hessian2_fn) || is.null(dyn$param2_jacobian_fn) ||
      is.null(dyn$hessian3_fn) || is.null(dyn$hessian2_fn)) return(NULL)

  layout <- .dsys_layout(compiled)
  tc     <- layout$total_cols
  np     <- ncol(dys)
  col_var_idx <- layout$col_var_idx
  endo_cols   <- which(!is.na(col_var_idx))
  dy <- .build_dy_ss_o2(compiled, ys)

  V <- matrix(0, tc, np)
  V[endo_cols, ] <- dys[col_var_idx[endo_cols], , drop = FALSE]

  ph2 <- tryCatch(dyn$param_hessian2_fn(dy, params, ys), error = function(e) NULL)
  h2  <- tryCatch(dyn$hessian2_fn(dy, params, ys),        error = function(e) NULL)
  h3  <- tryCatch(dyn$hessian3_fn(dy, params, ys),        error = function(e) NULL)
  p2j <- tryCatch(dyn$param2_jacobian_fn(dy, params, ys), error = function(e) NULL)
  if (is.null(ph2) || is.null(h2) || is.null(h3) || is.null(p2j)) return(NULL)
  if (any(!is.finite(ph2)) || any(!is.finite(h2)) || any(!is.finite(h3)) ||
      any(!is.finite(p2j))) return(NULL)

  ## SSM-computed-parameter chain data (Tier 12 #2, second order). The second
  ## total derivative of the dynamic Jacobian for such models uses TOTAL
  ## parameter directions W[,a] = e_a + Σ_pc (dp_c/dθ_a) e_pc in the param-tensor
  ## contractions (T1 param2_jac, T2/T3 param_hess2) and adds a new term
  ## (∂²F/∂w∂p_c)·(d² p_c/dθ_a dθ_b) (T6). Validated vs full re-solve (~1.8e-9).
  ssm <- NULL
  if (.ssm_assigns_param(compiled$model)) {
    ch1 <- .ssm_param_chain_derivs(compiled, ys, params)
    ch2 <- .ssm_param_chain_derivs2(compiled, ys, params)
    if (is.null(ch1) || is.null(ch2)) return(NULL)
    pars <- compiled$model$param_names
    endo <- compiled$model$var_names
    comp <- ch1$computed; free <- setdiff(pars, comp)
    dpc <- matrix(0, length(comp), np, dimnames = list(comp, pars))
    for (pc in comp) for (a in free)
      dpc[pc, a] <- ch1$dg_dtheta[pc, a] + sum(ch1$dg_dy[pc, endo] * dys[endo, a])
    W <- matrix(0, np, np, dimnames = list(pars, pars))
    for (a in free) { W[a, a] <- 1; for (pc in comp) W[pc, a] <- W[pc, a] + dpc[pc, a] }
    PJ <- matrix(dyn$param_jacobian_fn(dy, params, ys), layout$n_eq * tc, np)
    ssm <- list(comp = comp, pci = match(comp, pars), W = W, dpc = dpc,
                PJ_pc = PJ[, match(comp, pars), drop = FALSE],
                d2g = ch2$d2, dg_dy = ch1$dg_dy, dys_e = dys[endo, , drop = FALSE],
                endo = endo, free = free, pars = pars)
  }

  list(layout = layout, compiled = compiled, n_eq = layout$n_eq, total_cols = tc,
       V = V, col_var_idx = col_var_idx, endo_cols = endo_cols, d2ys = d2ys, ssm = ssm,
       ph2_trip = dyn$param_hess2_triplets, ph2_vals = ph2,
       h2_trip  = dyn$hess2_triplets,       h2_vals  = h2,
       h3_trip  = dyn$hess3_triplets,       h3_vals  = h3,
       p2j_trip = dyn$param2_jac_triplets,  p2j_vals = p2j)
}

## Total second derivative d²p_c/dθ_a dθ_b for an SSM-computed parameter pc
## (a, b are full param indices). Standard order-2 chain through the SSM block.
.sd2_total_d2pc <- function(ctx, pc, a, b) {
  ssm <- ctx$ssm; endo <- ssm$endo
  pa <- ssm$pars[a]; pb <- ssm$pars[b]
  if (!(pa %in% ssm$free && pb %in% ssm$free)) return(0)
  M <- ssm$d2g[[pc]]                                  # coords x coords (c(endo, free))
  da <- ssm$dys_e[endo, a]; db <- ssm$dys_e[endo, b]
  M[pa, pb] + sum(M[endo, pa] * db) + sum(M[endo, pb] * da) +
    sum(outer(da, db) * M[endo, endo]) +
    sum(ssm$dg_dy[pc, endo] * ctx$d2ys[endo, a, b])
}

.sd2_d2J_pair <- function(ctx, a, b) {
  n_eq <- ctx$n_eq; tc <- ctx$total_cols; V <- ctx$V
  d2J  <- matrix(0, n_eq, tc)
  ssm  <- ctx$ssm

  if (is.null(ssm)) {
    ## --- non-SSM: exact param-index match (unchanged hot path) ---
    amin <- min(a, b); amax <- max(a, b)
    ## T1: param2_jac (row,col,pa,pb), stored canonical pa<=pb.
    for (k in seq_along(ctx$p2j_trip)) {
      t <- ctx$p2j_trip[[k]]
      if (t$pa == amin && t$pb == amax)
        d2J[t$row, t$col] <- d2J[t$row, t$col] + ctx$p2j_vals[k]
    }
    ## T2 (param==a, contract free w-index with V[,b]) + T3 (param==b, V[,a]).
    for (k in seq_along(ctx$ph2_trip)) {
      t <- ctx$ph2_trip[[k]]; v <- ctx$ph2_vals[k]
      if (t$param == a) {
        d2J[t$eq, t$col1] <- d2J[t$eq, t$col1] + v * V[t$col2, b]
        if (t$col1 != t$col2) d2J[t$eq, t$col2] <- d2J[t$eq, t$col2] + v * V[t$col1, b]
      }
      if (t$param == b) {
        d2J[t$eq, t$col1] <- d2J[t$eq, t$col1] + v * V[t$col2, a]
        if (t$col1 != t$col2) d2J[t$eq, t$col2] <- d2J[t$eq, t$col2] + v * V[t$col1, a]
      }
    }
  } else {
    ## --- SSM: TOTAL parameter directions W[,a], W[,b] (Tier 12 #2) ---
    W <- ssm$W
    ## T1: Σ_{p,q} param2_jac[:,:,p,q] W[p,a] W[q,b] (canonical pa<=pb, symmetric).
    for (k in seq_along(ctx$p2j_trip)) {
      t <- ctx$p2j_trip[[k]]; v <- ctx$p2j_vals[k]; if (v == 0) next
      w <- W[t$pa, a] * W[t$pb, b]
      if (t$pa != t$pb) w <- w + W[t$pb, a] * W[t$pa, b]
      if (w != 0) d2J[t$row, t$col] <- d2J[t$row, t$col] + v * w
    }
    ## T2/T3: param_hess2 with TOTAL param direction (weight W[param, a]/[param, b]).
    for (k in seq_along(ctx$ph2_trip)) {
      t <- ctx$ph2_trip[[k]]; v <- ctx$ph2_vals[k]; if (v == 0) next
      wa <- W[t$param, a]; wb <- W[t$param, b]
      if (wa != 0) {
        d2J[t$eq, t$col1] <- d2J[t$eq, t$col1] + v * wa * V[t$col2, b]
        if (t$col1 != t$col2) d2J[t$eq, t$col2] <- d2J[t$eq, t$col2] + v * wa * V[t$col1, b]
      }
      if (wb != 0) {
        d2J[t$eq, t$col1] <- d2J[t$eq, t$col1] + v * wb * V[t$col2, a]
        if (t$col1 != t$col2) d2J[t$eq, t$col2] <- d2J[t$eq, t$col2] + v * wb * V[t$col1, a]
      }
    }
  }

  ## T4: hess2 contracted with the column-broadcast d2ys (V2[c]=d2ys[var(c),a,b]).
  V2 <- numeric(tc)
  ec <- ctx$endo_cols
  V2[ec] <- ctx$d2ys[ctx$col_var_idx[ec], a, b]
  for (k in seq_along(ctx$h2_trip)) {
    t <- ctx$h2_trip[[k]]; v <- ctx$h2_vals[k]
    d2J[t$eq, t$col1] <- d2J[t$eq, t$col1] + v * V2[t$col2]
    if (t$col1 != t$col2) d2J[t$eq, t$col2] <- d2J[t$eq, t$col2] + v * V2[t$col1]
  }

  ## T5: hess3 contracted on TWO indices with V[,a] and V[,b], keeping the
  ## first slot. Same .orbit_3 expansion as .contract_h3(h3, I, V_a, V_b).
  for (k in seq_along(ctx$h3_trip)) {
    val <- ctx$h3_vals[k]; if (val == 0) next
    t <- ctx$h3_trip[[k]]; e <- t$eq
    orbit <- .orbit_3(t$col1, t$col2, t$col3)
    for (perm in orbit) {
      keep <- perm[1]; x <- perm[2]; y <- perm[3]
      d2J[e, keep] <- d2J[e, keep] + val * V[x, a] * V[y, b]
    }
  }

  ## T6 (SSM only): (∂²F/∂w∂p_c)·(d²p_c/dθ_a dθ_b), the computed-parameter analog
  ## of T4 in param space (PJ_pc is the pure dynamic param-Jacobian p_c column).
  if (!is.null(ssm)) {
    for (jc in seq_along(ssm$comp)) {
      d2pc <- .sd2_total_d2pc(ctx, ssm$comp[jc], a, b)
      if (d2pc != 0) d2J <- d2J + matrix(ssm$PJ_pc[, jc] * d2pc, n_eq, tc)
    }
  }

  .partition_dJ(d2J, ctx$layout, ctx$compiled)
}


#' Second-order implicit-differentiation derivatives of the first-order rule
#'
#' @param model,compiled,dr,params,obs_vars  as in \code{solution_derivatives}
#' @param param_names parameters to differentiate wrt (Hessian is over their
#'   pairwise combinations)
#' @param h_rel  relative FD step for the first-order layer (default 1e-6)
#' @param h_rel2 relative FD step for the second-derivative primitive stencils
#'   (default 1e-4; larger because second differences are noisier)
#' @return list with:
#'   \item{first}{output of \code{solution_derivatives} (per-parameter dG/dH/...) }
#'   \item{d2}{named-by-"i|j" list of second-derivative blocks
#'     (d2G, d2H, d2ys, d2TT, d2RR, d2ZZ, d2DD, d2d), symmetric in (i,j)}
#'   \item{param_names}{echoed}
#' @noRd
solution_derivatives_2 <- function(model, compiled, dr, params, param_names,
                                    obs_vars, h_rel = 1e-6, h_rel2 = 1e-4) {

  if (!all(param_names %in% names(params)))
    stop("solution_derivatives_2: unknown parameter(s): ",
         paste(setdiff(param_names, names(params)), collapse = ", "))

  endo <- dr$endo_names; exo <- dr$exo_names
  n_endo <- length(endo); n_exo <- length(exo)
  state_idx <- dr$state_idx; n_state <- length(state_idx)
  obs_idx <- match(obs_vars, endo)
  np <- length(param_names)

  G  <- dr$ghx; H <- dr$ghu; ys <- dr$ys

  S <- matrix(0, n_state, n_endo)
  if (n_state > 0) S[cbind(seq_len(n_state), state_idx)] <- 1
  SG <- if (n_state > 0) S %*% G else matrix(0, 0, n_state)

  ## Base primitives + the SAME M, A factorizations used at first order.
  sys0 <- extract_system_matrices(compiled, ys, params)
  f_plus <- sys0$f_plus; f_zero <- sys0$f_zero
  f_minus <- sys0$f_minus; f_u <- sys0$f_exo
  A <- f_plus %*% (if (n_state > 0) G %*% S else matrix(0, n_endo, n_endo)) + f_zero
  M_qr <- if (n_state > 0)
    qr(kronecker(t(SG), f_plus) + kronecker(diag(n_state), A)) else NULL
  A_qr <- qr(A)

  ## ---- First-order derivatives for every requested parameter ----
  first <- solution_derivatives(model, compiled, dr, params, param_names,
                                obs_vars, h_rel = h_rel)
  G1 <- lapply(param_names, function(p) first$derivs[[p]]$dG)
  H1 <- lapply(param_names, function(p) first$derivs[[p]]$dH)
  names(G1) <- names(H1) <- param_names

  ## ---- Analytic vs finite-difference primitive derivatives ----
  ## When the second-order parameter codegen is present
  ## (compile_model(param_deriv = "second")) the first AND second total
  ## derivatives of the dynamic-Jacobian primitives are computed analytically
  ## (no per-parameter / per-pair steady re-solves); otherwise fall back to the
  ## finite-difference stencils.
  analytic_ok <- isTRUE(getOption("dynhr.use_analytic_primitives", TRUE)) &&
    .can_use_analytic_primitive_deriv(compiled) &&
    isTRUE(compiled$dynamic$param_deriv2_ok) &&
    isTRUE(compiled$static$static_param2_built)
  ## SSM-computed-parameter models now use the order-2 augmented chain
  ## (.analytic_d2ys SSM branch + .sd2_d2J_pair total directions + T6).

  dys_full <- dprim_full <- d2ys_full <- ad2 <- pidx <- NULL
  if (analytic_ok) {
    dys_full   <- .analytic_dys(compiled, ys, params)
    dprim_full <- if (!is.null(dys_full)) .analytic_dprimitives(compiled, ys, params, dys_full) else NULL
    d2ys_full  <- if (!is.null(dys_full)) .analytic_d2ys(compiled, ys, params, dys_full) else NULL
    ad2        <- if (!is.null(dys_full) && !is.null(d2ys_full))
                    .sd2_build_ctx(compiled, ys, params, dys_full, d2ys_full) else NULL
    if (is.null(dys_full) || is.null(dprim_full) || is.null(d2ys_full) || is.null(ad2))
      analytic_ok <- FALSE
  }

  hvec <- vapply(param_names, function(p) max(h_rel2 * abs(params[[p]]), 1e-5), 0)
  names(hvec) <- param_names

  if (analytic_ok) {
    ## First TOTAL primitive derivatives from the analytic first-order layer.
    pidx <- match(param_names, compiled$model$param_names)
    dfp <- lapply(param_names, function(p) dprim_full[[p]]$df_plus)
    df0 <- lapply(param_names, function(p) dprim_full[[p]]$df_zero)
    dfm <- lapply(param_names, function(p) dprim_full[[p]]$df_minus)
    dfu <- lapply(param_names, function(p) dprim_full[[p]]$df_exo)
  } else {
    ## ---- Single-sided perturbed primitives (FD fallback; for diagonal
    ##      stencils + first derivatives df_*_i used inside K_ij) ----
    prim_p <- prim_m <- vector("list", np); names(prim_p) <- names(prim_m) <- param_names
    for (i in seq_len(np)) {
      pi_name <- param_names[i]; h <- hvec[i]
      tp <- params; tp[[pi_name]] <- tp[[pi_name]] + h
      tm <- params; tm[[pi_name]] <- tm[[pi_name]] - h
      prim_p[[i]] <- .sd2_primitives_at(tp, model, compiled, ys)
      prim_m[[i]] <- .sd2_primitives_at(tm, model, compiled, ys)
      if (is.null(prim_p[[i]]) || is.null(prim_m[[i]]))
        stop(sprintf("solution_derivatives_2: steady state failed at %s +/- h", pi_name))
    }
    df <- function(field, i) (prim_p[[i]][[field]] - prim_m[[i]][[field]]) / (2 * hvec[i])
    dfp <- lapply(seq_len(np), function(i) df("f_plus",  i))
    df0 <- lapply(seq_len(np), function(i) df("f_zero",  i))
    dfm <- lapply(seq_len(np), function(i) df("f_minus", i))
    dfu <- lapply(seq_len(np), function(i) df("f_exo",   i))
  }

  Smask <- if (n_state > 0) state_idx else integer(0)

  d2 <- list()
  for (i in seq_len(np)) {
    for (j in i:np) {
      pi_name <- param_names[i]; pj_name <- param_names[j]
      hi <- hvec[i]; hj <- hvec[j]

      ## ---- TOTAL second derivatives of the smooth primitives ----
      if (analytic_ok) {
        ## Analytic 5-term second total derivative + analytic d2ys.
        a <- pidx[i]; b <- pidx[j]
        dpr2 <- .sd2_d2J_pair(ad2, a, b)
        d2fp <- dpr2$df_plus;  d2f0 <- dpr2$df_zero
        d2fm <- dpr2$df_minus; d2fu <- dpr2$df_exo
        d2ys <- d2ys_full[, a, b]
      } else if (i == j) {
        ## 3-point: F_ii = (F(+) - 2 F0 + F(-)) / h^2
        d2fp <- (prim_p[[i]]$f_plus  - 2 * f_plus  + prim_m[[i]]$f_plus)  / hi^2
        d2f0 <- (prim_p[[i]]$f_zero  - 2 * f_zero  + prim_m[[i]]$f_zero)  / hi^2
        d2fm <- (prim_p[[i]]$f_minus - 2 * f_minus + prim_m[[i]]$f_minus) / hi^2
        d2fu <- (prim_p[[i]]$f_exo   - 2 * f_u      + prim_m[[i]]$f_exo)  / hi^2
        d2ys <- (prim_p[[i]]$ys      - 2 * ys       + prim_m[[i]]$ys)     / hi^2
      } else {
        ## 4-corner mixed stencil with steady re-solve at each corner.
        tpp <- params; tpp[[pi_name]] <- tpp[[pi_name]] + hi; tpp[[pj_name]] <- tpp[[pj_name]] + hj
        tpm <- params; tpm[[pi_name]] <- tpm[[pi_name]] + hi; tpm[[pj_name]] <- tpm[[pj_name]] - hj
        tmp <- params; tmp[[pi_name]] <- tmp[[pi_name]] - hi; tmp[[pj_name]] <- tmp[[pj_name]] + hj
        tmm <- params; tmm[[pi_name]] <- tmm[[pi_name]] - hi; tmm[[pj_name]] <- tmm[[pj_name]] - hj
        Fpp <- .sd2_primitives_at(tpp, model, compiled, ys)
        Fpm <- .sd2_primitives_at(tpm, model, compiled, ys)
        Fmp <- .sd2_primitives_at(tmp, model, compiled, ys)
        Fmm <- .sd2_primitives_at(tmm, model, compiled, ys)
        if (is.null(Fpp) || is.null(Fpm) || is.null(Fmp) || is.null(Fmm))
          stop(sprintf("solution_derivatives_2: steady state failed at corner (%s,%s)",
                       pi_name, pj_name))
        den <- 4 * hi * hj
        mix <- function(fld) (Fpp[[fld]] - Fpm[[fld]] - Fmp[[fld]] + Fmm[[fld]]) / den
        d2fp <- mix("f_plus"); d2f0 <- mix("f_zero")
        d2fm <- mix("f_minus"); d2fu <- mix("f_exo"); d2ys <- mix("ys")
      }
      d2fmS <- if (n_state > 0) d2fm[, Smask, drop = FALSE] else matrix(0, n_endo, 0)

      Gi <- G1[[i]]; Gj <- G1[[j]]; Hi <- H1[[i]]; Hj <- H1[[j]]
      dfp_i <- dfp[[i]]; dfp_j <- dfp[[j]]; df0_i <- df0[[i]]; df0_j <- df0[[j]]

      if (n_state > 0) {
        SGi <- S %*% Gi; SGj <- S %*% Gj
        ## K_ij: all terms not multiplying d2G_ij (see header derivation).
        K <- d2fp %*% G %*% SG +
             dfp_j %*% Gi %*% SG + dfp_j %*% G %*% SGi +
             dfp_i %*% Gj %*% SG + f_plus %*% Gj %*% SGi +
             dfp_i %*% G %*% SGj + f_plus %*% Gi %*% SGj +
             d2f0 %*% G + df0_j %*% Gi +
             df0_i %*% Gj +
             d2fmS
        d2G <- matrix(qr.solve(M_qr, -as.vector(K)), n_endo, n_state)
      } else {
        d2G <- matrix(0, n_endo, 0)
      }

      ## A_i, A_j, A_ij and H_ij
      GiS <- if (n_state > 0) Gi %*% S else matrix(0, n_endo, n_endo)
      GjS <- if (n_state > 0) Gj %*% S else matrix(0, n_endo, n_endo)
      d2GS <- if (n_state > 0) d2G %*% S else matrix(0, n_endo, n_endo)
      GS  <- if (n_state > 0) G %*% S else matrix(0, n_endo, n_endo)
      A_i  <- dfp_i %*% GS + f_plus %*% GiS + df0_i
      A_j  <- dfp_j %*% GS + f_plus %*% GjS + df0_j
      A_ij <- d2fp %*% GS + dfp_i %*% GjS + dfp_j %*% GiS + f_plus %*% d2GS + d2f0
      d2H  <- -qr.solve(A_qr, d2fu + A_ij %*% H + A_i %*% Hj + A_j %*% Hi)

      key <- paste(i, j, sep = "|")
      blk <- list(
        d2G = d2G, d2H = d2H, d2ys = setNames(as.numeric(d2ys), endo),
        d2TT = if (n_state > 0) d2G[state_idx, , drop = FALSE] else matrix(0, 0, 0),
        d2RR = if (n_state > 0) d2H[state_idx, , drop = FALSE] else matrix(0, 0, n_exo),
        d2ZZ = d2G[obs_idx, , drop = FALSE],
        d2DD = d2H[obs_idx, , drop = FALSE],
        d2d  = setNames(as.numeric(d2ys)[obs_idx], obs_vars)
      )
      d2[[key]] <- blk
      if (i != j) d2[[paste(j, i, sep = "|")]] <- blk  # symmetric
    }
  }

  list(first = first, d2 = d2, param_names = param_names,
       ## Which second-primitive path produced d2f: "analytic" (param_deriv =
       ## "second" codegen, Tier 11 #3) or "fd" (the stencil fallback). Both
       ## are certified to agree on regular models (nk_small: identical eigen
       ## spectra); callers surface this so a Hessian consumer can tell which
       ## path it got without re-deriving the compile-time gating.
       second_primitives = if (analytic_ok) "analytic" else "fd")
}
