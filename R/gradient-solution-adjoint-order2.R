## R/gradient-solution-adjoint-order2.R
## --------------------------------------------------------------------------
## ADJOINT (reverse-mode) of the ORDER-2 perturbation solve (Tier 18 A2,
## "order-up" analogue of the phase-1 first-order .solution_adjoint).
##
## The forward layer (solution_derivatives_order2, gradient-solution-deriv-
## order2.R) solves, for EVERY structural parameter j:
##
##   K_xx vec(d_ghxx_j) = -vec(dPhi_xx_j) - dK_xx_j vec(ghxx)          (1)
##   A_L  d_ghxu_j      = rhs_xu_j(d_ghxx_j, dA_L_j, ...)              (2)
##   A_L  d_ghuu_j      = rhs_uu_j(d_ghxx_j, dA_L_j, ...)             (3)
##  (A_L+fp) d_ghss_j   = dRHS_ss_j(d_ghuu_j, ...) - (dA_L_j+dfp_j) ghss (4)
##
## i.e. one K_xx-solve + three A_L / (A_L+fp)-solves PER PARAMETER.
##
## This file reverses that DAG.  Given cotangents (bar_ghxx, bar_ghxu,
## bar_ghuu, bar_ghss) on the four solution blocks (= dL/d each block for a
## scalar loss L), it back-propagates in REVERSE topological order
## (ghss -> ghuu,ghxu -> ghxx), accumulating cotangents on the *per-parameter
## primitive inputs* dfp, df0, dfu, dhx, dhu, dG, dH, dvSe and the model-
## Hessian derivative dH_mat.  ONE transposed K_xx solve and a handful of
## transposed A_L / (A_L+fp) solves are shared across all parameters; then
## every parameter's gradient is a Frobenius contraction of the accumulated
## primitive cotangents against that parameter's (analytic or FD) primitive
## derivatives -- the SAME primitives the forward layer consumes.  Result:
## O(1) factorizations+solves, P cheap contractions, EXACT agreement with the
## forward layer (pinned by test-solution-adjoint-order2.R).
##
## Reverse of the linear solves uses the adjoint operator identity: if
## y = A^{-1} b then bar_b = A^{-T} bar_y and bar_A = -bar_b y'.  For the
## generalized Sylvester (1), M(X) = K_xx vec(X) has adjoint M*(Y) solving
## the transposed Kronecker system (same .solve_kron_compact on transposed
## factors).
## --------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## Adjoint of .bilinear_h2: out[e,] = vec(T_a' H2[e] T_b).
## Given bar_out (n_eq x (nca*ncb)), returns the cotangents
##   bar_H2[e,,] += T_a bar_M_e T_b'      (n_eq x ta_rows x tb_rows)
##   bar_Ta      += H2[e] T_b bar_M_e'    summed over e
##   bar_Tb      += H2[e]' T_a bar_M_e    summed over e
## where bar_M_e = matrix(bar_out[e,], nca, ncb).
## ---------------------------------------------------------------------------
.bilinear_h2_adjoint <- function(H2, T_a, T_b, bar_out,
                                 want_H2 = TRUE, want_Ta = TRUE,
                                 want_Tb = TRUE) {
  n_eq <- dim(H2)[1]
  nca  <- ncol(T_a); ncb <- ncol(T_b)
  bar_H2 <- if (want_H2) array(0, dim = dim(H2)) else NULL
  bar_Ta <- if (want_Ta) matrix(0, nrow(T_a), nca) else NULL
  bar_Tb <- if (want_Tb) matrix(0, nrow(T_b), ncb) else NULL
  for (e in seq_len(n_eq)) {
    bM <- matrix(bar_out[e, ], nca, ncb)   # cotangent on T_a' H2_e T_b
    He <- H2[e, , ]
    if (want_H2) bar_H2[e, , ] <- T_a %*% bM %*% t(T_b)
    if (want_Ta) bar_Ta <- bar_Ta + He %*% T_b %*% t(bM)
    if (want_Tb) bar_Tb <- bar_Tb + t(He) %*% T_a %*% bM
  }
  list(bar_H2 = bar_H2, bar_Ta = bar_Ta, bar_Tb = bar_Tb)
}


## ---------------------------------------------------------------------------
## Transposed generalized-Sylvester / Kronecker solve for K_xx' L = rhs.
## K_xx = kron(I_{ns2}, A_L) + kron(hxt %x% hxt, fp), so
## K_xx' = kron(I_{ns2}, A_L') + kron((hxt %x% hxt)', fp').
## We solve it with the same dense QR the forward uses (K_xx is already
## explicitly assembled and factorized ONCE upstream); this helper just
## solves against the transpose.  rhs is the vec of an (n x ns2) matrix.
## ---------------------------------------------------------------------------
.o2adj_kxx_transpose_solve <- function(K_xx_t_qr, rhs_vec) {
  qr.solve(K_xx_t_qr, rhs_vec)
}


#' Reverse-mode order-2 structural-parameter gradient through the perturbation
#' solve.
#'
#' Given cotangents on the four second-order solution blocks
#' (\code{bar_ghxx}, \code{bar_ghxu}, \code{bar_ghuu}, \code{bar_ghss}),
#' returns \code{dL/dθ_j} for a scalar loss whose sensitivity to each block is
#' the corresponding bar, by reversing the order-2 perturbation DAG.  One
#' transposed \code{K_xx} solve and a handful of transposed \code{A_L} /
#' \code{A_L+fp} solves are shared across parameters; then a Frobenius
#' contraction per parameter.
#'
#' @param model,compiled,dr2,params,param_names  as in
#'   \code{solution_derivatives_order2}.
#' @param bars list with any of \code{bar_ghxx} (n x ns2), \code{bar_ghxu}
#'   (n x ns*nu), \code{bar_ghuu} (n x nu^2), \code{bar_ghss} (length n).  NULL
#'   entries treated as zero.
#' @param h_rel,h_hess FD steps (used only on the FD primitive fallback).
#' @return list(grad = named numeric (NA where a param's primitives failed),
#'   ok = named logical, used_analytic = logical)
#' @noRd
.solution_adjoint_order2 <- function(model, compiled, dr2, params, param_names,
                                     bars, h_rel = 1e-6, h_hess = 1e-4) {

  if (!all(param_names %in% names(params)))
    stop(".solution_adjoint_order2: unknown parameter(s): ",
         paste(setdiff(param_names, names(params)), collapse = ", "))
  if (!inherits(dr2, "DecisionRules2"))
    stop(".solution_adjoint_order2: dr2 must be a DecisionRules2 object")

  endo_names <- dr2$endo_names
  exo_names  <- dr2$exo_names
  state_idx  <- dr2$state_idx
  n   <- length(endo_names)
  n_s <- length(state_idx)
  n_u <- length(exo_names)
  np  <- length(param_names)

  ghx  <- dr2$ghx; ghu <- dr2$ghu
  ghxx <- dr2$ghxx; ghxu <- dr2$ghxu; ghuu <- dr2$ghuu; ghss <- dr2$ghss
  Sigma_e <- dr2$Sigma_e
  ys  <- dr2$ys

  hx <- ghx[state_idx, , drop = FALSE]
  hu <- ghu[state_idx, , drop = FALSE]
  ns2 <- n_s * n_s

  S <- matrix(0, n_s, n)
  if (n_s > 0) S[cbind(seq_len(n_s), state_idx)] <- 1

  ## ---- base system + operators (mirror the forward layer exactly) ----------
  sys0 <- extract_system_matrices(compiled, ys, params)
  f0 <- sys0$f_zero; fp <- sys0$f_plus; fu <- sys0$f_exo

  GS  <- ghx %*% S
  A_L <- f0 + fp %*% GS
  hxt <- t(hx)
  K_xx <- kronecker(diag(ns2), A_L) + kronecker(hxt %x% hxt, fp)
  K_xx_t_qr <- qr(t(K_xx))
  AL_t_qr   <- qr(t(A_L))
  AB_t_qr   <- qr(t(A_L + fp))

  tm  <- .build_transfer_matrices(compiled$dynamic, ghx, ghu, state_idx,
                                  hx, hu, endo_names, exo_names)
  T_x <- tm$T_x; T_u <- tm$T_u

  H_base <- .o2sd_hessian_at(compiled, params, ys)
  if (is.null(H_base))
    stop(".solution_adjoint_order2: could not compute model Hessian at base")
  if (!is.null(compiled$model$equations)) {
    eq_to_decl <- .build_eq_to_decl(compiled$model)
    if (all(eq_to_decl > 0L) && !identical(eq_to_decl, seq_len(n))) {
      perm <- order(eq_to_decl)
      H_perm <- array(0, dim = dim(H_base))
      for (k in seq_len(n)) H_perm[k, , ] <- H_base[perm[k], , ]
      H_base <- H_perm
    }
  }

  has_lead <- sys0$is_fwd | sys0$is_mixed
  T_up <- .build_T_up(compiled$dynamic, ghu, endo_names, exo_names, has_lead)
  vSe  <- as.numeric(Sigma_e)

  ## ---- cotangents on solution blocks ---------------------------------------
  zmat <- function(x, nr, nc) {
    if (is.null(x)) return(matrix(0, nr, nc))
    x <- as.matrix(x); stopifnot(nrow(x) == nr, ncol(x) == nc); x
  }
  bar_ghxx <- zmat(bars$bar_ghxx, n, ns2)
  bar_ghxu <- zmat(bars$bar_ghxu, n, n_s * n_u)
  bar_ghuu <- zmat(bars$bar_ghuu, n, n_u * n_u)
  bar_ghss <- if (is.null(bars$bar_ghss)) numeric(n) else as.numeric(bars$bar_ghss)
  stopifnot(length(bar_ghss) == n)

  ## Base bilinear-H2(T_up,T_up) matrix (constant across params).
  bilH_up  <- .bilinear_h2(H_base, T_up, T_up)     # n x nu^2

  ## =========================================================================
  ## REVERSE PASS.  Accumulate cotangents on the per-parameter primitives.
  ## Primitive cotangents accumulated below (base-point matrices):
  ##   bar_dfp   (n x n),  bar_df0 (n x n),  bar_dfu (n x n_u)
  ##   bar_dhx   (n_s x n_s), bar_dhu (n_s x n_u)   -- state-block of dG/dH
  ##   bar_dG    (n x n_s),   bar_dH  (n x n_u)      -- full dG/dH cotangents
  ##   bar_dvSe  (length n_u^2)
  ##   bar_dHmat (n x total_cols x total_cols)       -- model-Hessian derivative
  ## Together with bar_dT_x / bar_dT_u / bar_dT_up which we then fold into
  ## bar_dG / bar_dH via the transfer-matrix chain (dT depends on dG/dH).
  ## =========================================================================
  total_cols <- compiled$dynamic$total_cols
  bar_dfp   <- matrix(0, n, n)
  bar_df0   <- matrix(0, n, n)
  bar_dfu   <- matrix(0, n, n_u)
  bar_dG    <- matrix(0, n, n_s)
  bar_dH    <- matrix(0, n, n_u)
  bar_dvSe  <- numeric(n_u * n_u)
  bar_dHmat <- array(0, dim = c(n, total_cols, total_cols))
  bar_dT_x  <- matrix(0, total_cols, n_s)
  bar_dT_u  <- matrix(0, total_cols, n_u)
  bar_dT_up <- matrix(0, total_cols, n_u)
  ## dA_L cotangent (dA_L = dfp GS + fp dG S + df0); folded at the end.
  bar_dA_L  <- matrix(0, n, n)
  ## d_ghxx cotangent (accumulated from ghxu/ghuu branches + the direct bar).
  bar_dghxx <- bar_ghxx

  ## ---- (4) reverse d_ghss ---------------------------------------------------
  ## d_ghss = (A_L+fp)^{-1} ( dRHS_ss - (dA_L+dfp) ghss ),  dRHS_ss = -(term1+term2)
  if (any(bar_ghss != 0)) {
    br_ss <- qr.solve(AB_t_qr, bar_ghss)            # cotangent on the RHS vector
    ## -(dA_L+dfp) ghss :
    bar_dA_L <- bar_dA_L + (-br_ss) %*% t(ghss)
    bar_dfp  <- bar_dfp  + (-br_ss) %*% t(ghss)
    ## dRHS_ss = -(term1 + term2); so cotangent on (term1+term2) is +br_ss.
    ## term1 = dfp ghuu vSe + fp d_ghuu vSe + fp ghuu dvSe
    guv <- ghuu %*% vSe                              # n
    bar_dfp    <- bar_dfp + (-br_ss) %*% t(guv)
    bar_d_ghuu <- (-t(fp)) %*% br_ss %*% t(vSe)      # n x nu^2 (from fp d_ghuu vSe)
    bar_dvSe   <- bar_dvSe + as.numeric(-t(fp %*% ghuu) %*% br_ss)
    ## term2 = (bil(dHmat,Tup,Tup) + 2 bil(H,dTup,Tup)) vSe + bil(H,Tup,Tup) dvSe
    ##   bil(dHmat,Tup,Tup) vSe : cotangent on the (n x nu^2) matrix is (-br_ss) vSe'
    bar_bilA <- (-br_ss) %*% t(vSe)                  # n x nu^2  (feeds dHmat & dTup)
    ##   split: dHmat channel + 2*dTup channel
    a1 <- .bilinear_h2_adjoint(H_base, T_up, T_up, bar_bilA,
                               want_H2 = FALSE, want_Ta = TRUE, want_Tb = TRUE)
    ## dHmat channel: bar_out = bar_bilA feeds bil(dHmat,...) => cotangent on dHmat
    a1H <- .bilinear_h2_adjoint(H_base, T_up, T_up, bar_bilA,
                                want_H2 = TRUE, want_Ta = FALSE, want_Tb = FALSE)
    bar_dHmat <- bar_dHmat + a1H$bar_H2
    ## 2 * bil(H, dTup, Tup) vSe : symmetric in the two Tup slots at fixed H,
    ## cotangent on dTup is 2 * (adjoint wrt the FIRST T slot of bil(H,dTup,Tup)).
    a2 <- .bilinear_h2_adjoint(H_base, T_up, T_up, 2 * bar_bilA,
                               want_H2 = FALSE, want_Ta = TRUE, want_Tb = FALSE)
    bar_dT_up <- bar_dT_up + a2$bar_Ta
    ##   bil(H,Tup,Tup) dvSe : cotangent on dvSe
    bar_dvSe  <- bar_dvSe + as.numeric(-t(bilH_up) %*% br_ss)
    ## fold d_ghuu cotangent into the ghuu-branch accumulator
    bar_dghuu_from_ss <- bar_d_ghuu
  } else {
    bar_dghuu_from_ss <- matrix(0, n, n_u * n_u)
  }

  ## ---- (3) reverse d_ghuu ---------------------------------------------------
  ## A_L d_ghuu = rhs_uu, rhs_uu = -(dPhi_uu + dfp ghxx (hu⊗hu) + fp d_ghxx (hu⊗hu)
  ##              + fp ghxx (dhu⊗hu + hu⊗dhu)) - dA_L ghuu
  bar_ghuu_tot <- bar_ghuu + bar_dghuu_from_ss
  bar_dhu <- matrix(0, n_s, n_u)   # accumulate dhu here (state block) + into bar_dH later? no: dhu = dH[state,]
  if (any(bar_ghuu_tot != 0) && n_u > 0) {
    br_uu <- qr.solve(AL_t_qr, bar_ghuu_tot)         # cotangent on rhs_uu (n x nu^2)
    ## -dA_L ghuu
    bar_dA_L <- bar_dA_L + (-br_uu) %*% t(ghuu)
    ## rhs_uu inner = -( ... ); cotangent on the inner (...) is +br_uu.
    huhu <- hu %x% hu                                # nu^2 x nu^2  (=(hu⊗hu))
    ## dPhi_uu term: cotangent on dPhi_uu = -br_uu
    bar_dPhi_uu <- -br_uu
    ## dfp ghxx (hu⊗hu):   bar_dfp += (-br_uu) (ghxx (hu⊗hu))'
    bar_dfp <- bar_dfp + (-br_uu) %*% t(ghxx %*% huhu)
    ## fp d_ghxx (hu⊗hu):  bar_d_ghxx += fp' (-br_uu) (hu⊗hu)'
    bar_dghxx <- bar_dghxx + t(fp) %*% (-br_uu) %*% t(huhu)
    ## fp ghxx (dhu⊗hu + hu⊗dhu): cotangent on dhu.
    ## Let C = fp ghxx (n x nu^2). inner = C (dhu⊗hu + hu⊗dhu).
    ## bar wrt (dhu⊗hu + hu⊗dhu) is C' (-br_uu) (nu^2 x nu^2).
    C <- fp %*% ghxx
    bar_kron_uu <- t(C) %*% (-br_uu)                 # ns^2 x nu^2
    ## hu %x% hu: both factors hu (n_s x n_u); left & right adjoints add.
    bar_dhu <- bar_dhu + .adj_kron_left(bar_kron_uu, hu, n_s, n_u) +
                         .adj_kron_right(bar_kron_uu, hu, n_s, n_u)
  } else {
    bar_dPhi_uu <- matrix(0, n, n_u * n_u)
  }

  ## ---- (2) reverse d_ghxu ---------------------------------------------------
  ## A_L d_ghxu = rhs_xu, rhs_xu = -(dPhi_xu + dfp ghxx (hu⊗hx) + fp d_ghxx (hu⊗hx)
  ##              + fp ghxx (dhu⊗hx + hu⊗dhx)) - dA_L ghxu
  bar_dhx <- matrix(0, n_s, n_s)
  if (any(bar_ghxu != 0) && n_u > 0 && n_s > 0) {
    br_xu <- qr.solve(AL_t_qr, bar_ghxu)             # cotangent on rhs_xu (n x ns*nu)
    bar_dA_L <- bar_dA_L + (-br_xu) %*% t(ghxu)
    huhx <- hu %x% hx                                # (ns*nu) x (ns*nu)? (hu⊗hx): nu*ns rows
    bar_dPhi_xu <- -br_xu
    bar_dfp <- bar_dfp + (-br_xu) %*% t(ghxx %*% huhx)
    bar_dghxx <- bar_dghxx + t(fp) %*% (-br_xu) %*% t(huhx)
    C <- fp %*% ghxx
    bar_kron_xu <- t(C) %*% (-br_xu)                 # ns^2 x (nu*ns)
    ## hu %x% hx: left factor hu (n_s x n_u), right factor hx (n_s x n_s).
    ## dhu⊗hx -> adjoint wrt left factor dhu (dims n_s x n_u), Q = hx
    bar_dhu <- bar_dhu + .adj_kron_left(bar_kron_xu, hx, n_s, n_u)
    ## hu⊗dhx -> adjoint wrt right factor dhx (dims n_s x n_s), P = hu
    bar_dhx <- bar_dhx + .adj_kron_right(bar_kron_xu, hu, n_s, n_s)
  } else {
    bar_dPhi_xu <- matrix(0, n, n_s * n_u)
  }

  ## ---- (1) reverse d_ghxx (ONE transposed K_xx solve) ----------------------
  ## K_xx vec(d_ghxx) = -vec(dPhi_xx) - dK_xx vec(ghxx)
  ## bar on vec(d_ghxx) is vec(bar_dghxx).  L = K_xx^{-T} bar.
  bar_dPhi_xx <- matrix(0, n, ns2)
  if (n_s > 0 && any(bar_dghxx != 0)) {
    Lvec <- .o2adj_kxx_transpose_solve(K_xx_t_qr, as.numeric(bar_dghxx))
    ## cotangent on the RHS ( -vec(dPhi_xx) - dK_xx vec(ghxx) ) is Lvec.
    ## => cotangent on dPhi_xx is -matrix(Lvec).
    bar_dPhi_xx <- -matrix(Lvec, n, ns2)
    ## dK_xx vec(ghxx) term: cotangent on dK_xx is (-Lvec) (vec(ghxx))'
    ## dK_xx = kron(I_ns2, dA_L) + kron(d_hxtkron, fp) + kron(hxt⊗hxt, dfp)
    ## We contract via matrix identities instead of forming dK_xx.
    vg  <- as.numeric(ghxx)                          # ns2*n
    negL <- -Lvec
    ## (a) kron(I_ns2, dA_L): dK_xx vec(ghxx)[block b] = dA_L %*% ghxx[,b].
    ##     bar_dA_L += Σ_b (negL_block_b) (ghxx[,b])'
    NL <- matrix(negL, n, ns2)                       # negL reshaped (n x ns2)
    bar_dA_L <- bar_dA_L + NL %*% t(ghxx)            # sum_b outer(NL[,b], ghxx[,b])
    ## (c) kron(hxt⊗hxt, dfp): block structure -> bar_dfp += fp-analogue.
    ##     (hxt⊗hxt) is ns2 x ns2; term contributes
    ##     NL %*% t( ghxx %*% t(hxt⊗hxt) )   [since (M⊗dfp) vec = dfp G M' pattern]
    HH <- hxt %x% hxt                                # ns2 x ns2
    bar_dfp <- bar_dfp + NL %*% t(ghxx %*% t(HH))
    ## (b) kron(d_hxtkron, fp): d_hxtkron = dhxt⊗hxt + hxt⊗dhxt.
    ##     term = fp %*% ghxx %*% t(d_hxtkron).  cotangent flows to d_hxtkron:
    ##     bar_dhxtkron = (fp %*% ghxx)' (NL)  contracted -> then to dhx.
    Cg <- fp %*% ghxx                                # n x ns2
    bar_dhxtkron <- t(Cg) %*% NL                     # ns2 x ns2 (cotangent on d_hxtkron^T pattern)
    ## d_hxtkron = dhxt⊗hxt + hxt⊗dhxt ; term uses t(d_hxtkron), so
    ## bar on t(d_hxtkron) is bar_dhxtkron => bar on d_hxtkron is t(bar_dhxtkron).
    bd <- t(bar_dhxtkron)
    ## dhxt⊗hxt : left factor dhxt (ns x ns), right hxt
    bar_dhxt <- .adj_kron_left(bd, hxt, n_s, n_s) +
                .adj_kron_right(bd, hxt, n_s, n_s)
    bar_dhx <- bar_dhx + t(bar_dhxt)                 # dhxt = t(dhx)
  }

  ## ---- reverse dPhi_{xx,xu,uu} -> dT_x, dT_u, dHmat -------------------------
  ## dPhi_xx[e] = vec(dTx' He Tx + Tx' He dTx + Tx' dHe Tx)  (He = H_base[e])
  ## dPhi_xu[e] = vec(dTx' He Tu + Tx' He dTu + Tx' dHe Tu)
  ## dPhi_uu[e] = vec(dTu' He Tu + Tu' He dTu + Tu' dHe Tu)
  Tx_t <- t(T_x); Tu_t <- t(T_u)
  for (e in seq_len(n)) {
    He <- H_base[e, , ]
    ## --- xx ---
    if (n_s > 0 && any(bar_dPhi_xx[e, ] != 0)) {
      B <- matrix(bar_dPhi_xx[e, ], n_s, n_s)        # cotangent on the (ns x ns) block
      ## d/dTx of (dTx' He Tx): appears as dTx' He Tx -> bar_dTx += He Tx B'
      ## and Tx' He dTx -> bar_dTx += (He' Tx) B
      bar_dT_x <- bar_dT_x + He %*% T_x %*% t(B) + t(He) %*% T_x %*% B
      ## Tx' dHe Tx -> bar_dHe += Tx B Tx'
      bar_dHmat[e, , ] <- bar_dHmat[e, , ] + T_x %*% B %*% Tx_t
    }
    ## --- xu ---
    if (n_s > 0 && n_u > 0 && any(bar_dPhi_xu[e, ] != 0)) {
      B <- matrix(bar_dPhi_xu[e, ], n_s, n_u)
      ## dTx' He Tu -> bar_dTx += He Tu B'
      bar_dT_x <- bar_dT_x + He %*% T_u %*% t(B)
      ## Tx' He dTu -> bar_dTu += He' Tx B
      bar_dT_u <- bar_dT_u + t(He) %*% T_x %*% B
      ## Tx' dHe Tu -> bar_dHe += Tx B Tu'
      bar_dHmat[e, , ] <- bar_dHmat[e, , ] + T_x %*% B %*% Tu_t
    }
    ## --- uu ---
    if (n_u > 0 && any(bar_dPhi_uu[e, ] != 0)) {
      B <- matrix(bar_dPhi_uu[e, ], n_u, n_u)
      bar_dT_u <- bar_dT_u + He %*% T_u %*% t(B) + t(He) %*% T_u %*% B
      bar_dHmat[e, , ] <- bar_dHmat[e, , ] + T_u %*% B %*% Tu_t
    }
  }

  ## ---- reverse dT_x, dT_u, dT_up -> dG, dH ---------------------------------
  ## dT_x[c,]/dT_u[c,] built from dG/dH per the DCM (mirror the forward loop):
  ##   ll==0 endo j: dT_x[c,]=dG[j,], dT_u[c,]=dH[j,]
  ##   ll==1 endo j: dT_x[c,]=(dG hx + G dhx)[j,], dT_u[c,]=(dG hu + G dhu)[j,]
  ## dT_up[c,]=dH[j,] for jumper compounds.
  dcm <- compiled$dynamic$dyn_col_map
  ## dT_up first (feeds bar_dH directly)
  if (any(bar_dT_up != 0)) {
    jm <- .identify_jumpers(compiled$dynamic, endo_names, has_lead)
    for (i in seq_along(jm$jumper_idx)) {
      j <- jm$jumper_idx[i]; c <- jm$jumper_compound_c[i]
      bar_dH[j, ] <- bar_dH[j, ] + bar_dT_up[c, ]
    }
  }
  ## dT_x / dT_u
  for (kk in seq_len(nrow(dcm))) {
    c_col <- dcm$col[kk]; nm <- dcm$name[kk]; ll <- dcm$lead_lag[kk]
    if (nm %in% exo_names) next
    j <- which(endo_names == nm)
    if (length(j) != 1L) next
    if (ll == 0L) {
      bar_dG[j, ] <- bar_dG[j, ] + bar_dT_x[c_col, ]
      bar_dH[j, ] <- bar_dH[j, ] + bar_dT_u[c_col, ]
    } else if (ll == 1L) {
      ## dT_x[c,] = (dG hx + G dhx)[j,]
      ## bar wrt dG[j,] += bar_dT_x[c,] hx'  ; bar wrt dhx += G[j,]' bar_dT_x[c,]
      bar_dG[j, ] <- bar_dG[j, ] + as.numeric(hx %*% bar_dT_x[c_col, ])
      bar_dhx <- bar_dhx + outer(ghx[j, ], bar_dT_x[c_col, ])
      if (n_u > 0) {
        bar_dG[j, ] <- bar_dG[j, ] + as.numeric(hu %*% bar_dT_u[c_col, ])
        bar_dhu <- bar_dhu + outer(ghx[j, ], bar_dT_u[c_col, ])
      }
    }
  }

  ## ---- fold dhx/dhu (state block) into bar_dG/bar_dH -----------------------
  ## dhx = dG[state_idx,], dhu = dH[state_idx,].
  if (n_s > 0) {
    bar_dG[state_idx, ] <- bar_dG[state_idx, , drop = FALSE] + bar_dhx
    if (n_u > 0)
      bar_dH[state_idx, ] <- bar_dH[state_idx, , drop = FALSE] + bar_dhu
  }

  ## ---- fold dA_L = dfp GS + fp dG S + df0 -----------------------------------
  ## bar_dfp += bar_dA_L (GS)'; bar_df0 += bar_dA_L; bar_dG += fp' bar_dA_L S'.
  bar_dfp <- bar_dfp + bar_dA_L %*% t(GS)
  bar_df0 <- bar_df0 + bar_dA_L
  if (n_s > 0)
    bar_dG <- bar_dG + t(fp) %*% bar_dA_L %*% t(S)

  ## =========================================================================
  ## Per-parameter primitive derivatives + contraction.
  ## We reuse the forward layer's primitive assembly (analytic or FD) so the
  ## contraction reproduces the forward result EXACTLY.  dG_j/dH_j come from
  ## the first-order layer; dHmat_j from the analytic tensor or central FD;
  ## dfp/df0/dfu from analytic/FD; dvSe from .o2sd_dSigma_e.
  ## =========================================================================
  prims <- .o2adj_param_primitives(model, compiled, dr2, params, param_names,
                                   H_base, h_rel, h_hess, n, n_s, n_u,
                                   endo_names, exo_names, state_idx)

  grad <- setNames(rep(NA_real_, np), param_names)
  ok   <- setNames(rep(FALSE, np), param_names)
  for (pnm in param_names) {
    pj <- prims$derivs[[pnm]]
    if (is.null(pj) || !isTRUE(pj$ok)) next
    gj <- sum(bar_dfp * pj$dfp) + sum(bar_df0 * pj$df0)
    if (n_u > 0)  gj <- gj + sum(bar_dfu * pj$dfu)
    if (n_s > 0)  gj <- gj + sum(bar_dG * pj$dG)
    if (n_u > 0)  gj <- gj + sum(bar_dH * pj$dH)
    if (length(bar_dvSe)) gj <- gj + sum(bar_dvSe * as.numeric(pj$dvSe))
    gj <- gj + sum(bar_dHmat * pj$dHmat)
    grad[pnm] <- gj
    ok[pnm]   <- TRUE
  }

  list(grad = grad, ok = ok, used_analytic = prims$used_analytic)
}


## ---------------------------------------------------------------------------
## Adjoint of a Kronecker product wrt its LEFT factor.
## For M = kron(P, Q) with P (rp x cp), Q (rq x cq), and cotangent bar_M
## (rp*rq x cp*cq), the cotangent on P is
##   bar_P[i,j] = <bar_M block (i,j) sub, Q> = sum(bar_M[(i-1)rq+(1:rq),
##                (j-1)cq+(1:cq)] * Q).
## Here we pass Q and the DIMENSIONS of P (rp=cp given via np1/np2? we take
## P square-ish through explicit dims).  Signature: bar_M, Q, rp, cp with
## Q's dims inferred; P dims = (rp x cp).
## ---------------------------------------------------------------------------
.adj_kron_left <- function(bar_M, Q, rp, cp) {
  rq <- nrow(Q); cq <- ncol(Q)
  bar_P <- matrix(0, rp, cp)
  for (i in seq_len(rp)) for (j in seq_len(cp)) {
    blk <- bar_M[((i - 1) * rq + 1):(i * rq),
                 ((j - 1) * cq + 1):(j * cq), drop = FALSE]
    bar_P[i, j] <- sum(blk * Q)
  }
  bar_P
}

## Adjoint of kron(P, Q) wrt its RIGHT factor Q, given P and Q's dims (rq x cq).
##   bar_Q = sum_{i,j} P[i,j] * bar_M block (i,j).
.adj_kron_right <- function(bar_M, P, rq, cq) {
  rp <- nrow(P); cp <- ncol(P)
  bar_Q <- matrix(0, rq, cq)
  for (i in seq_len(rp)) for (j in seq_len(cp)) {
    blk <- bar_M[((i - 1) * rq + 1):(i * rq),
                 ((j - 1) * cq + 1):(j * cq), drop = FALSE]
    bar_Q <- bar_Q + P[i, j] * blk
  }
  bar_Q
}


## ---------------------------------------------------------------------------
## Per-parameter primitive derivatives for the order-2 adjoint contraction.
## Mirrors solution_derivatives_order2's primitive assembly: analytic when
## available (param-Jacobian + param-Hessian2 + hess3.dys), else central FD.
## Returns list(derivs = by-param list(dfp,df0,dfu,dG,dH,dvSe,dHmat,ok),
##              used_analytic).
## ---------------------------------------------------------------------------
.o2adj_param_primitives <- function(model, compiled, dr2, params, param_names,
                                    H_base, h_rel, h_hess, n, n_s, n_u,
                                    endo_names, exo_names, state_idx) {
  np <- length(param_names)
  ## dG/dH from the first-order layer (shared with the forward path).
  first <- solution_derivatives(model, compiled, dr2, params, param_names,
                                obs_vars = endo_names, h_rel = h_rel)

  use_analytic <- .can_use_analytic_primitive_deriv(compiled) &&
    isTRUE(compiled$dynamic$param_hess2_built) &&
    isTRUE(compiled$dynamic$hessian3_built)

  dprim_an <- NULL; dH_an <- NULL
  if (use_analytic) {
    dys_all   <- .analytic_dys(compiled, ys = dr2$ys, params = params)
    dprim_all <- if (!is.null(dys_all))
      .analytic_dprimitives(compiled, dr2$ys, params, dys_all) else NULL
    dH_all    <- if (!is.null(dys_all))
      .o2sd_analytic_dH(compiled, dr2$ys, params, dys_all) else NULL
    if (is.null(dys_all) || is.null(dprim_all) || is.null(dH_all)) {
      use_analytic <- FALSE
    } else {
      perm_dh <- NULL
      if (!is.null(compiled$model$equations)) {
        e2d <- .build_eq_to_decl(compiled$model)
        if (all(e2d > 0L) && !identical(e2d, seq_len(n))) perm_dh <- order(e2d)
      }
      dprim_an <- dprim_all
      dH_an <- vector("list", np); names(dH_an) <- param_names
      for (pnm in param_names) {
        dHk <- dH_all[, , , pnm]
        if (!is.null(perm_dh)) {
          dHk_perm <- array(0, dim = dim(dHk))
          for (e in seq_len(n)) dHk_perm[e, , ] <- dHk[perm_dh[e], , ]
          dHk <- dHk_perm
        }
        dH_an[[pnm]] <- dHk
      }
    }
  }

  hvec <- vapply(param_names, function(p) max(h_rel * abs(params[[p]]), 1e-7), 0)
  names(hvec) <- param_names

  derivs <- vector("list", np); names(derivs) <- param_names
  for (k in seq_len(np)) {
    pnm <- param_names[k]; h <- hvec[k]
    d1 <- first$derivs[[pnm]]
    if (!isTRUE(d1$ok)) { derivs[[pnm]] <- list(ok = FALSE); next }

    if (use_analytic) {
      dp <- dprim_an[[pnm]]
      dfp <- dp$df_plus; df0 <- dp$df_zero; dfu <- dp$df_exo
      dHmat <- dH_an[[pnm]]
    } else {
      tp <- params; tp[[pnm]] <- tp[[pnm]] + h
      tm_ <- params; tm_[[pnm]] <- tm_[[pnm]] - h
      ss_p <- tryCatch(solve_steady(compiled, tp, y0 = dr2$ys,
                                    endo_names = model$var_names,
                                    exo_names = model$varexo_names, verbose = FALSE),
                       error = function(e) list(converged = FALSE))
      ss_m <- tryCatch(solve_steady(compiled, tm_, y0 = dr2$ys,
                                    endo_names = model$var_names,
                                    exo_names = model$varexo_names, verbose = FALSE),
                       error = function(e) list(converged = FALSE))
      if (!isTRUE(ss_p$converged) || !isTRUE(ss_m$converged)) {
        derivs[[pnm]] <- list(ok = FALSE); next
      }
      pp <- .ssm_consistent_params(model, tp); pm <- .ssm_consistent_params(model, tm_)
      sp <- extract_system_matrices(compiled, ss_p$values, pp)
      sm <- extract_system_matrices(compiled, ss_m$values, pm)
      dfp <- (sp$f_plus - sm$f_plus) / (2 * h)
      df0 <- (sp$f_zero - sm$f_zero) / (2 * h)
      dfu <- (sp$f_exo  - sm$f_exo)  / (2 * h)
      H_p <- .o2sd_hessian_at(compiled, pp, ss_p$values)
      H_m <- .o2sd_hessian_at(compiled, pm, ss_m$values)
      if (is.null(H_p) || is.null(H_m)) { derivs[[pnm]] <- list(ok = FALSE); next }
      if (!is.null(compiled$model$equations)) {
        e2d <- .build_eq_to_decl(compiled$model)
        if (all(e2d > 0L) && !identical(e2d, seq_len(n))) {
          perm <- order(e2d)
          Hp <- array(0, dim = dim(H_p)); Hm <- array(0, dim = dim(H_m))
          for (e in seq_len(n)) { Hp[e, , ] <- H_p[perm[e], , ]; Hm[e, , ] <- H_m[perm[e], , ] }
          H_p <- Hp; H_m <- Hm
        }
      }
      dHmat <- (H_p - H_m) / (2 * h)
    }

    dvSe <- .o2sd_dSigma_e(model, params, pnm, h, n_u)

    derivs[[pnm]] <- list(
      dfp = dfp, df0 = df0, dfu = dfu,
      dG = d1$dG, dH = d1$dH, dvSe = dvSe, dHmat = dHmat, ok = TRUE)
  }

  list(derivs = derivs, used_analytic = use_analytic)
}
