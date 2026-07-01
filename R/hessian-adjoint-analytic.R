## R/hessian-adjoint-analytic.R
## --------------------------------------------------------------------------
## FULLY-ANALYTIC second-order adjoint of the Kalman filter:
## "forward-over-reverse" = propagate a tangent through the primal + backward
## sweep to compute dG/dX along a given direction dX.
##
## .kf_loglik_dG(Y, ss, dX, me_variance = 0) returns the DIRECTIONAL
## DERIVATIVE of all six adjoint gradient matrices along dX:
##   list(dG_TT, dG_RR, dG_ZZ, dG_DD, dg_d, dG_Sig)
## These are d/dε G_X(X + ε dX)|_{ε=0}, computed analytically.
##
## Convention: primal quantities use the SAME naming as gradient-adjoint-kf.R
## (bar_s, bar_P, bar_F, bar_K, bar_v, bar_Mnum); their tangents carry a "d"
## prefix (dbar_s, dbar_P, dbar_F, dbar_K, dbar_v, dbar_Mnum).
##
## The tangent of the FORWARD stored quantities (ds, dP, dv, dFi, dK, dA, dB)
## is computed by differentiating the same recursion in gradient-tangent-kf.R.
## --------------------------------------------------------------------------

#' TRUE when the compiled forward-over-reverse dG kernel is available.
#' @noRd
.HAS_RCPP_KF_DG <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("kf_loglik_dG_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

#' Directional derivative of the adjoint gradient matrices
#'
#' Runs a combined forward (primal + tangent) and backward (primal + tangent)
#' pass of the Kalman filter, returning the directional derivative of the six
#' gradient matrices along the state-space direction \code{dX}.
#'
#' @param Y          n_obs x n_T observation matrix (no NAs).
#' @param ss         list: TT, RR, ZZ, DD, d, Sigma_e (base point X).
#' @param dX         list with (any of) dTT, dRR, dZZ, dDD, dd, dSigma_e
#'                   (the direction; missing blocks treated as zero).
#' @param me_variance scalar measurement-error variance (parameter-independent).
#'
#' @return list(dG_TT, dG_RR, dG_ZZ, dG_DD, dg_d, dG_Sig) -- the directional
#'   derivative of each gradient matrix along dX.
#'   Returns NULL invisibly on filter failure (non-PD F).
#' @noRd
.kf_loglik_dG <- function(Y, ss, dX, me_variance = 0) {

  TT <- ss$TT; RR <- ss$RR; ZZ <- ss$ZZ; DD <- ss$DD
  d  <- as.numeric(ss$d); Sigma_e <- ss$Sigma_e

  n_state <- nrow(TT)
  n_obs   <- nrow(ZZ)
  n_exo   <- ncol(RR)

  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  if (anyNA(Y)) stop(".kf_loglik_dG: Y must not contain missing values.")
  n_T <- ncol(Y)

  ## ---- Extract direction components (default to zero if absent) -----------
  z_state <- matrix(0, n_state, n_state)
  z_RR    <- matrix(0, n_state, n_exo)
  z_ZZ    <- matrix(0, n_obs, n_state)
  z_DD    <- matrix(0, n_obs, n_exo)
  z_d     <- numeric(n_obs)
  z_Sig   <- matrix(0, n_exo, n_exo)

  dTT  <- if (!is.null(dX$dTT))     dX$dTT     else z_state
  dRR  <- if (!is.null(dX$dRR))     dX$dRR     else z_RR
  dZZ  <- if (!is.null(dX$dZZ))     dX$dZZ     else z_ZZ
  dDD  <- if (!is.null(dX$dDD))     dX$dDD     else z_DD
  dd   <- if (!is.null(dX$dd))      as.numeric(dX$dd)  else z_d
  dSig <- if (!is.null(dX$dSigma_e)) dX$dSigma_e else z_Sig

  ## ---- Fast path: compiled forward-over-reverse kernel --------------------
  if (.HAS_RCPP_KF_DG()) {
    out <- kf_loglik_dG_cpp(Y, TT, RR, ZZ, DD, d, Sigma_e,
                            dTT, dRR, dZZ, dDD, dd, dSig,
                            me_variance, .KF_LL_MIN)
    if (!isTRUE(out$ok)) return(invisible(NULL))
    return(list(dG_TT  = out$dG_TT,  dG_RR  = out$dG_RR,
                dG_ZZ  = out$dG_ZZ,  dG_DD  = out$dG_DD,
                dg_d   = as.numeric(out$dg_d),
                dG_Sig = out$dG_Sig))
  }

  tZZ  <- t(ZZ)
  tdZZ <- t(dZZ)
  QQ   <- tcrossprod(RR %*% Sigma_e, RR)
  HH   <- tcrossprod(DD %*% Sigma_e, DD)
  SS   <- RR %*% Sigma_e %*% t(DD)
  me_diag  <- me_variance * diag(n_obs)
  ll_const <- -0.5 * n_obs * log(2 * pi)

  ## Tangent of noise matrices w.r.t. dX:
  ## dQQ = dRR Sig RR' + RR dSig RR' + RR Sig dRR'
  dQQ <- .sym(tcrossprod(dRR %*% Sigma_e, RR) +
              tcrossprod(RR  %*% dSig,    RR) +
              tcrossprod(RR  %*% Sigma_e, dRR))
  ## dHH = dDD Sig DD' + DD dSig DD' + DD Sig dDD'
  dHH <- .sym(tcrossprod(dDD %*% Sigma_e, DD) +
              tcrossprod(DD  %*% dSig,    DD) +
              tcrossprod(DD  %*% Sigma_e, dDD))
  ## dSS = dRR Sig DD' + RR dSig DD' + RR Sig dDD'
  dSS <- dRR %*% Sigma_e %*% t(DD) +
         RR  %*% dSig    %*% t(DD) +
         RR  %*% Sigma_e %*% t(dDD)

  ## ---- Lyapunov initial conditions and their tangents ---------------------
  P0 <- .solve_lyapunov(TT, QQ)
  if (!all(is.finite(P0))) return(invisible(NULL))

  ## dP0 = solve_lyapunov(TT, RHS), RHS = dTT P0 TT' + TT P0 dTT' + dQQ
  rhs_dP0 <- .sym(dTT %*% P0 %*% t(TT) + TT %*% P0 %*% t(dTT) + dQQ)
  dP0 <- .solve_lyapunov(TT, rhs_dP0)
  if (!all(is.finite(dP0))) return(invisible(NULL))
  dP0 <- .sym(dP0)

  ## ---- Allocate storage for forward pass ----------------------------------
  s_store  <- vector("list", n_T)
  P_store  <- vector("list", n_T)
  v_store  <- vector("list", n_T)
  Fi_store <- vector("list", n_T)
  K_store  <- vector("list", n_T)
  A_store  <- vector("list", n_T)
  B_store  <- vector("list", n_T)

  ds_store  <- vector("list", n_T)
  dP_store  <- vector("list", n_T)
  dv_store  <- vector("list", n_T)
  dFi_store <- vector("list", n_T)
  dK_store  <- vector("list", n_T)
  dA_store  <- vector("list", n_T)
  dB_store  <- vector("list", n_T)

  s  <- numeric(n_state)
  P  <- P0
  ds <- numeric(n_state)
  dP <- dP0

  loglik <- 0.0

  ## ---- Forward pass: primal + tangent at each step ------------------------
  for (t in seq_len(n_T)) {
    s_store[[t]]  <- s
    P_store[[t]]  <- P
    ds_store[[t]] <- ds
    dP_store[[t]] <- dP

    ## Primal innovation covariance
    PZ  <- P %*% tZZ
    Ft  <- .sym(ZZ %*% PZ + HH + me_diag)

    Fc <- tryCatch(chol(Ft), error = function(e) NULL)
    if (is.null(Fc)) return(invisible(NULL))
    Fi  <- chol2inv(Fc)
    ldf <- 2 * sum(log(diag(Fc)))

    v   <- Y[, t] - d - as.numeric(ZZ %*% s)
    Fiv <- as.numeric(Fi %*% v)

    ll_t <- ll_const - 0.5 * (ldf + sum(v * Fiv))
    if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) return(invisible(NULL))
    loglik <- loglik + ll_t

    ## Primal gain, A, B
    K <- (TT %*% PZ + SS) %*% Fi
    A <- TT - K %*% ZZ
    B <- RR - K %*% DD

    ## ---- Tangent of forward quantities ------------------------------------
    ## dv = -dd - dZZ s - ZZ ds
    dv_t <- -dd - as.numeric(dZZ %*% s) - as.numeric(ZZ %*% ds)

    ## dF = dZZ P ZZ' + ZZ dP ZZ' + ZZ P dZZ' + dHH   (symmetrise)
    dF_t <- .sym(tcrossprod(dZZ %*% P,  ZZ) +
                 ZZ %*% tcrossprod(dP,  ZZ) +
                 ZZ %*% tcrossprod(P,  dZZ) +
                 dHH)

    ## dFi = -Fi dF Fi
    dFi_t <- -Fi %*% dF_t %*% Fi

    ## dK = (dTT P ZZ' + TT dP ZZ' + TT P dZZ' + dSS) Fi - K dF Fi
    ## Product rule on K = Mnum * Fi (Mnum = TT P ZZ' + SS):
    ##   dK = dMnum * Fi + Mnum * dFi = dMnum * Fi + K * F * dFi
    ## Since F * dFi = F * (-Fi dF Fi) = -dF Fi (because F*Fi=I),
    ##   dK = dMnum * Fi - K * dF * Fi
    dPZ   <- dP %*% tZZ
    dK_num <- tcrossprod(dTT %*% P, ZZ) +
              TT %*% dPZ +
              TT %*% P %*% tdZZ +
              dSS
    dK_t   <- dK_num %*% Fi - K %*% (dF_t %*% Fi)

    ## dA = dTT - dK ZZ - K dZZ
    dA_t <- dTT - dK_t %*% ZZ - K %*% dZZ

    ## dB = dRR - dK DD - K dDD
    dB_t <- dRR - dK_t %*% DD - K %*% dDD

    ## Store tangents
    v_store[[t]]   <- v
    Fi_store[[t]]  <- Fi
    K_store[[t]]   <- K
    A_store[[t]]   <- A
    B_store[[t]]   <- B
    dv_store[[t]]  <- dv_t
    dFi_store[[t]] <- dFi_t
    dK_store[[t]]  <- dK_t
    dA_store[[t]]  <- dA_t
    dB_store[[t]]  <- dB_t

    ## Advance primal state
    s_new <- as.numeric(TT %*% s) + as.numeric(K %*% v)
    P_new <- .sym(tcrossprod(A %*% P, A) + tcrossprod(B %*% Sigma_e, B))

    ## Advance tangent state:
    ## ds' = dTT s + TT ds + dK v + K dv
    ds_new <- as.numeric(dTT %*% s) + as.numeric(TT %*% ds) +
              as.numeric(dK_t %*% v) + as.numeric(K %*% dv_t)

    ## dP' = dA P A' + A dP A' + A P dA' + dB Sig B' + B dSig B' + B Sig dB'
    AP   <- A %*% P
    BSig <- B %*% Sigma_e
    dP_new <- .sym(tcrossprod(dA_t %*% P,  A) +
                   A %*% tcrossprod(dP,     A) +
                   AP %*% t(dA_t) +
                   tcrossprod(dB_t %*% Sigma_e, B) +
                   B %*% tcrossprod(dSig, B) +
                   BSig %*% t(dB_t))

    s  <- s_new;  P  <- P_new
    ds <- ds_new; dP <- dP_new
  }

  ## ---- Backward pass: primal adjoint + its tangent ------------------------
  ##
  ## The backward sweep computes G_TT,G_RR,G_ZZ,G_DD,g_d,G_Sig.
  ## Simultaneously we differentiate each line by the product rule, carrying:
  ##   - tangents of stored forward quantities: ds_prev, dP_prev, dv, dFi, dK, dA, dB
  ##   - tangents of bar-variables: dbar_s, dbar_P, dbar_F, dbar_K, dbar_v, dbar_Mnum
  ##
  ## For each assignment in the primal:
  ##   result = f(forward_stuff, bar_stuff)
  ## the tangent is:
  ##   dresult = f_forward(d_forward_stuff, bar_stuff) + f_bar(forward_stuff, dbar_stuff)
  ##
  ## Notation: primal gradients G_* are accumulated; their tangents dG_* are
  ## accumulated in parallel.

  bar_s  <- numeric(n_state)
  bar_P  <- matrix(0, n_state, n_state)
  dbar_s <- numeric(n_state)
  dbar_P <- matrix(0, n_state, n_state)

  G_TT  <- matrix(0, n_state, n_state)
  G_RR  <- matrix(0, n_state, n_exo)
  G_ZZ  <- matrix(0, n_obs,   n_state)
  G_DD  <- matrix(0, n_obs,   n_exo)
  g_d   <- numeric(n_obs)
  G_Sig <- matrix(0, n_exo,   n_exo)

  dG_TT  <- matrix(0, n_state, n_state)
  dG_RR  <- matrix(0, n_state, n_exo)
  dG_ZZ  <- matrix(0, n_obs,   n_state)
  dG_DD  <- matrix(0, n_obs,   n_exo)
  dg_d   <- numeric(n_obs)
  dG_Sig <- matrix(0, n_exo,   n_exo)

  for (t in n_T:1) {
    s_prev  <- s_store[[t]]
    P_prev  <- P_store[[t]]
    v       <- v_store[[t]]
    Fi      <- Fi_store[[t]]
    K       <- K_store[[t]]
    A       <- A_store[[t]]
    B       <- B_store[[t]]

    ds_prev  <- ds_store[[t]]
    dP_prev  <- dP_store[[t]]
    dv       <- dv_store[[t]]
    dFi      <- dFi_store[[t]]
    dK       <- dK_store[[t]]
    dA       <- dA_store[[t]]
    dB       <- dB_store[[t]]

    Fiv  <- as.numeric(Fi %*% v)
    dFiv <- as.numeric(dFi %*% v) + as.numeric(Fi %*% dv)

    ## ==================================================================
    ## Step 1: Adjoint of P_t = sym(A P_prev A' + B Sigma_e B')
    ## ==================================================================

    ## Primal: bar_P <- sym(bar_P)
    bar_P  <- .sym(bar_P)
    dbar_P <- .sym(dbar_P)

    ## Primal: bar_A = 2 bar_P A P_prev
    bar_A  <- 2 * bar_P %*% A %*% P_prev
    ## Tangent:  dbar_A = 2 (dbar_P A P_prev + bar_P dA P_prev + bar_P A dP_prev)
    dbar_A <- 2 * (dbar_P %*% A %*% P_prev +
                   bar_P  %*% dA %*% P_prev +
                   bar_P  %*% A  %*% dP_prev)

    ## Primal: bar_B = 2 bar_P B Sigma_e
    bar_B  <- 2 * bar_P %*% B %*% Sigma_e
    ## Tangent: dbar_B = 2 (dbar_P B Sigma_e + bar_P dB Sigma_e + bar_P B dSig)
    dbar_B <- 2 * (dbar_P %*% B  %*% Sigma_e +
                   bar_P  %*% dB %*% Sigma_e +
                   bar_P  %*% B  %*% dSig)

    ## Primal: bar_P_prev_from_AP = A' bar_P A
    bar_P_prev_from_AP  <- t(A) %*% bar_P  %*% A
    ## Tangent: dA' bar_P A + A' dbar_P A + A' bar_P dA
    dbar_P_prev_from_AP <- t(dA) %*% bar_P  %*% A  +
                           t(A)  %*% dbar_P %*% A  +
                           t(A)  %*% bar_P  %*% dA

    ## Primal: G_Sig += B' bar_P B
    G_Sig  <- G_Sig  + t(B)  %*% bar_P  %*% B
    ## Tangent: dG_Sig += dB' bar_P B + B' dbar_P B + B' bar_P dB
    dG_Sig <- dG_Sig + t(dB) %*% bar_P  %*% B  +
                       t(B)  %*% dbar_P %*% B  +
                       t(B)  %*% bar_P  %*% dB

    ## ==================================================================
    ## Step 2: Adjoint of s_t = TT s_prev + K v
    ## ==================================================================

    ## Primal: G_TT += outer(bar_s, s_prev)
    G_TT  <- G_TT  + outer(bar_s,  s_prev)
    ## Tangent: dG_TT += outer(dbar_s, s_prev) + outer(bar_s, ds_prev)
    dG_TT <- dG_TT + outer(dbar_s, s_prev) + outer(bar_s, ds_prev)

    ## Primal: bar_K = outer(bar_s, v)
    bar_K  <- outer(bar_s,  v)
    ## Tangent: dbar_K = outer(dbar_s, v) + outer(bar_s, dv)
    dbar_K <- outer(dbar_s, v) + outer(bar_s, dv)

    ## Primal: bar_v = K' bar_s
    bar_v  <- as.numeric(t(K)  %*% bar_s)
    ## Tangent: dbar_v = dK' bar_s + K' dbar_s
    dbar_v <- as.numeric(t(dK) %*% bar_s) + as.numeric(t(K) %*% dbar_s)

    ## Primal: bar_s_prev = TT' bar_s
    bar_s_prev  <- as.numeric(t(TT)  %*% bar_s)
    ## Tangent: dbar_s_prev = dTT' bar_s + TT' dbar_s
    dbar_s_prev <- as.numeric(t(dTT) %*% bar_s) + as.numeric(t(TT) %*% dbar_s)

    ## ==================================================================
    ## Step 3: Adjoint of ll_t = -0.5*(log|F| + v' Fi v)
    ## ==================================================================

    ## Primal: bar_F = -0.5*(Fi - outer(Fiv, Fiv))
    bar_F  <- -0.5 * (Fi  - outer(Fiv,  Fiv))
    ## Tangent: dbar_F = -0.5*(dFi - outer(dFiv, Fiv) - outer(Fiv, dFiv))
    dbar_F <- -0.5 * (dFi - outer(dFiv, Fiv) - outer(Fiv, dFiv))

    ## Primal: bar_v += -Fiv
    bar_v  <- bar_v  - Fiv
    dbar_v <- dbar_v - dFiv

    ## ==================================================================
    ## Step 4: Adjoint of A = TT - K ZZ
    ## ==================================================================

    ## Primal: G_TT += bar_A
    G_TT  <- G_TT  + bar_A
    dG_TT <- dG_TT + dbar_A

    ## Primal: G_ZZ -= K' bar_A
    G_ZZ  <- G_ZZ  - t(K)  %*% bar_A
    dG_ZZ <- dG_ZZ - t(dK) %*% bar_A - t(K) %*% dbar_A

    ## Primal: bar_K -= bar_A ZZ'
    bar_K  <- bar_K  - bar_A  %*% tZZ
    dbar_K <- dbar_K - dbar_A %*% tZZ - bar_A %*% tdZZ

    ## ==================================================================
    ## Step 5: Adjoint of B = RR - K DD
    ## ==================================================================

    ## Primal: G_RR += bar_B
    G_RR  <- G_RR  + bar_B
    dG_RR <- dG_RR + dbar_B

    ## Primal: G_DD -= K' bar_B
    G_DD  <- G_DD  - t(K)  %*% bar_B
    dG_DD <- dG_DD - t(dK) %*% bar_B - t(K) %*% dbar_B

    ## Primal: bar_K -= bar_B DD'
    bar_K  <- bar_K  - bar_B  %*% t(DD)
    dbar_K <- dbar_K - dbar_B %*% t(DD) - bar_B %*% t(dDD)

    ## ==================================================================
    ## Step 6: Adjoint of K = Mnum Fi (Mnum = TT P_prev ZZ' + SS)
    ## ==================================================================

    ## Primal: bar_F += -K' bar_K Fi   (then symmetrize)
    bar_F  <- bar_F  - t(K)  %*% bar_K  %*% Fi
    dbar_F <- dbar_F - t(dK) %*% bar_K  %*% Fi  -
                       t(K)  %*% dbar_K %*% Fi  -
                       t(K)  %*% bar_K  %*% dFi

    ## Primal: bar_F <- sym(bar_F)
    bar_F  <- .sym(bar_F)
    dbar_F <- .sym(dbar_F)

    ## Primal: bar_Mnum = bar_K Fi
    bar_Mnum  <- bar_K  %*% Fi
    dbar_Mnum <- dbar_K %*% Fi + bar_K %*% dFi

    ## From Mnum = TT P_prev ZZ' + SS:
    ## Primal: G_TT += bar_Mnum ZZ P_prev
    G_TT  <- G_TT  + bar_Mnum  %*% ZZ  %*% P_prev
    dG_TT <- dG_TT + dbar_Mnum %*% ZZ  %*% P_prev  +
                     bar_Mnum  %*% dZZ %*% P_prev  +
                     bar_Mnum  %*% ZZ  %*% dP_prev

    ## Primal: G_ZZ += bar_Mnum' TT P_prev
    G_ZZ  <- G_ZZ  + t(bar_Mnum)  %*% TT  %*% P_prev
    dG_ZZ <- dG_ZZ + t(dbar_Mnum) %*% TT  %*% P_prev  +
                     t(bar_Mnum)  %*% dTT %*% P_prev  +
                     t(bar_Mnum)  %*% TT  %*% dP_prev

    ## Primal: bar_P_prev_from_Mnum = TT' bar_Mnum ZZ
    bar_P_prev_from_Mnum  <- t(TT)  %*% bar_Mnum  %*% ZZ
    dbar_P_prev_from_Mnum <- t(dTT) %*% bar_Mnum  %*% ZZ  +
                             t(TT)  %*% dbar_Mnum %*% ZZ  +
                             t(TT)  %*% bar_Mnum  %*% dZZ

    ## From SS = RR Sig DD' in Mnum:
    ## Primal: G_RR += bar_Mnum DD Sigma_e
    G_RR  <- G_RR  + bar_Mnum  %*% DD  %*% Sigma_e
    dG_RR <- dG_RR + dbar_Mnum %*% DD  %*% Sigma_e  +
                     bar_Mnum  %*% dDD %*% Sigma_e  +
                     bar_Mnum  %*% DD  %*% dSig

    ## Primal: G_DD += bar_Mnum' RR Sigma_e
    G_DD  <- G_DD  + t(bar_Mnum)  %*% RR  %*% Sigma_e
    dG_DD <- dG_DD + t(dbar_Mnum) %*% RR  %*% Sigma_e  +
                     t(bar_Mnum)  %*% dRR %*% Sigma_e  +
                     t(bar_Mnum)  %*% RR  %*% dSig

    ## Primal: G_Sig += RR' bar_Mnum DD
    G_Sig  <- G_Sig  + t(RR)  %*% bar_Mnum  %*% DD
    dG_Sig <- dG_Sig + t(dRR) %*% bar_Mnum  %*% DD  +
                       t(RR)  %*% dbar_Mnum %*% DD  +
                       t(RR)  %*% bar_Mnum  %*% dDD

    ## ==================================================================
    ## Step 7: Adjoint of F = sym(ZZ P_prev ZZ' + HH + me_diag)
    ## ==================================================================

    ## bar_F already symmetrized above.
    ## Primal: bar_P_prev_from_F = ZZ' bar_F ZZ
    bar_P_prev_from_F  <- t(ZZ)  %*% bar_F  %*% ZZ
    dbar_P_prev_from_F <- t(dZZ) %*% bar_F  %*% ZZ  +
                          t(ZZ)  %*% dbar_F %*% ZZ  +
                          t(ZZ)  %*% bar_F  %*% dZZ

    ## Primal: G_ZZ += 2 bar_F ZZ P_prev
    G_ZZ  <- G_ZZ  + 2 * bar_F  %*% ZZ  %*% P_prev
    dG_ZZ <- dG_ZZ + 2 * (dbar_F %*% ZZ  %*% P_prev  +
                           bar_F  %*% dZZ %*% P_prev  +
                           bar_F  %*% ZZ  %*% dP_prev)

    ## From HH = DD Sig DD':
    ## Primal: G_DD += 2 bar_F DD Sigma_e
    G_DD  <- G_DD  + 2 * bar_F  %*% DD  %*% Sigma_e
    dG_DD <- dG_DD + 2 * (dbar_F %*% DD  %*% Sigma_e  +
                           bar_F  %*% dDD %*% Sigma_e  +
                           bar_F  %*% DD  %*% dSig)

    ## Primal: G_Sig += DD' bar_F DD
    G_Sig  <- G_Sig  + t(DD)  %*% bar_F  %*% DD
    dG_Sig <- dG_Sig + t(dDD) %*% bar_F  %*% DD  +
                       t(DD)  %*% dbar_F %*% DD  +
                       t(DD)  %*% bar_F  %*% dDD

    ## ==================================================================
    ## Step 8: Adjoint of v_t = y_t - d - ZZ s_prev
    ## ==================================================================

    ## Primal: g_d -= bar_v
    g_d  <- g_d  - bar_v
    dg_d <- dg_d - dbar_v

    ## Primal: G_ZZ -= outer(bar_v, s_prev)
    G_ZZ  <- G_ZZ  - outer(bar_v,  s_prev)
    dG_ZZ <- dG_ZZ - outer(dbar_v, s_prev) - outer(bar_v, ds_prev)

    ## Primal: bar_s_prev -= ZZ' bar_v
    bar_s_prev  <- bar_s_prev  - as.numeric(t(ZZ)  %*% bar_v)
    dbar_s_prev <- dbar_s_prev - as.numeric(t(dZZ) %*% bar_v) -
                                 as.numeric(t(ZZ)  %*% dbar_v)

    ## ==================================================================
    ## Step 9: Collect bar_P_{t-1}
    ## ==================================================================
    bar_P_prev_new  <- .sym(bar_P_prev_from_AP  + bar_P_prev_from_Mnum  + bar_P_prev_from_F)
    dbar_P_prev_new <- .sym(dbar_P_prev_from_AP + dbar_P_prev_from_Mnum + dbar_P_prev_from_F)

    ## Update carry variables
    bar_s  <- bar_s_prev
    bar_P  <- bar_P_prev_new
    dbar_s <- dbar_s_prev
    dbar_P <- dbar_P_prev_new
  }

  ## ==================================================================
  ## Lyapunov adjoint: bar_QQ = solve_lyapunov(TT', bar_P0)
  ## ==================================================================
  bar_P0  <- .sym(bar_P)
  dbar_P0 <- .sym(dbar_P)

  if (n_state == 1L) {
    ## Scalar: P_0 = QQ / (1 - TT^2)
    denom   <- 1 - TT[1, 1]^2
    bar_QQ  <- matrix(bar_P0[1, 1]  / denom, 1, 1)
    dbar_QQ <- matrix(dbar_P0[1, 1] / denom, 1, 1)
    ## Primal scalar Lyapunov TT contribution:
    ## G_TT[1,1] += bar_P0 * 2 * TT * P0 / denom
    G_TT[1, 1]  <- G_TT[1, 1]  + bar_P0[1, 1]  * 2 * TT[1, 1] * P0[1, 1] / denom
    ## Tangent of scalar TT G_TT term:
    ## d/deps (bar_P0 * 2 * TT * P0 / (1-TT^2))
    ## = dbar_P0 * 2 TT P0/denom + bar_P0 * 2 dTT P0/denom + bar_P0 * 2 TT dP0/denom
    ##   + bar_P0 * 2 TT P0 * d(1/denom)/ddeps
    ## d(1/denom)/d eps = 2 TT dTT / (1-TT^2)^2  -- but encoded via dP0 tangent
    ## Simplest: use the product rule directly on the scalar formula
    ## P0 = QQ/(1-TT^2), dP0 = (dQQ - 2 TT P0 dTT)/(1-TT^2) -- already in dP0
    dG_TT[1, 1] <- dG_TT[1, 1] +
      dbar_P0[1, 1] * 2 * TT[1, 1] * P0[1, 1] / denom +
      bar_P0[1, 1]  * 2 * dTT[1, 1] * P0[1, 1] / denom +
      bar_P0[1, 1]  * 2 * TT[1, 1]  * dP0[1, 1] / denom +
      bar_P0[1, 1]  * 2 * TT[1, 1]  * P0[1, 1] * 2 * TT[1, 1] * dTT[1, 1] / denom^2
  } else {
    ## General: bar_QQ = solve_lyapunov(TT', bar_P0)
    bar_QQ <- tryCatch(.solve_lyapunov(t(TT), bar_P0), error = function(e) NULL)
    if (is.null(bar_QQ) || !all(is.finite(bar_QQ))) return(invisible(NULL))
    bar_QQ <- .sym(bar_QQ)

    ## Tangent of bar_QQ = solve_lyapunov(TT', bar_P0):
    ## The Lyapunov equation is (I - TT'^T ⊗ TT'^T) vec(bar_QQ) = vec(bar_P0).
    ## Differentiating: (I - TT'^T ⊗ TT'^T) vec(dbar_QQ) = vec(dbar_P0) + d[(TT'^T ⊗ TT'^T)] vec(bar_QQ)
    ## The latter is: d[(TT'^T ⊗ TT'^T)] vec(bar_QQ)
    ##   = (dTT'^T ⊗ TT'^T + TT'^T ⊗ dTT'^T) vec(bar_QQ)
    ##   = vec(TT' bar_QQ dTT + dTT' bar_QQ TT)  [Kronecker product formula]
    ## So the RHS for dbar_QQ is:  dbar_P0 + TT' bar_QQ dTT + dTT' bar_QQ TT
    ## = d/deps of [TT' bar_QQ TT = bar_P0], i.e. standard Lyapunov tangent with LHS TT'.
    ## rhs_dbar_QQ = dbar_P0 + d(TT' bar_QQ TT)/deps
    ##             = dbar_P0 + dTT' bar_QQ TT + TT' dbar_QQ TT + TT' bar_QQ dTT  ... circular!
    ## Actually: differentiate the Lyapunov equation  bar_QQ - TT' bar_QQ TT = bar_P0
    ## wrt eps:  dbar_QQ - dTT' bar_QQ TT - TT' dbar_QQ TT - TT' bar_QQ dTT = dbar_P0
    ## => dbar_QQ - TT' dbar_QQ TT = dbar_P0 + dTT' bar_QQ TT + TT' bar_QQ dTT
    ## So: dbar_QQ = solve_lyapunov(TT', rhs)
    ## where rhs = dbar_P0 + dTT' bar_QQ TT + TT' bar_QQ dTT  (symmetrized)
    rhs_dbar_QQ <- .sym(dbar_P0 + t(dTT) %*% bar_QQ %*% TT + t(TT) %*% bar_QQ %*% dTT)
    dbar_QQ <- tryCatch(.solve_lyapunov(t(TT), rhs_dbar_QQ), error = function(e) NULL)
    if (is.null(dbar_QQ) || !all(is.finite(dbar_QQ))) return(invisible(NULL))
    dbar_QQ <- .sym(dbar_QQ)

    ## Primal: G_TT += 2 bar_QQ TT P0
    G_TT  <- G_TT  + 2 * bar_QQ  %*% TT  %*% P0
    ## Tangent: dG_TT += 2 (dbar_QQ TT P0 + bar_QQ dTT P0 + bar_QQ TT dP0)
    dG_TT <- dG_TT + 2 * (dbar_QQ %*% TT  %*% P0  +
                           bar_QQ  %*% dTT %*% P0  +
                           bar_QQ  %*% TT  %*% dP0)
  }

  ## Primal: G_Sig += RR' bar_QQ RR
  G_Sig  <- G_Sig  + t(RR)  %*% bar_QQ  %*% RR
  dG_Sig <- dG_Sig + t(dRR) %*% bar_QQ  %*% RR  +
                     t(RR)  %*% dbar_QQ %*% RR  +
                     t(RR)  %*% bar_QQ  %*% dRR

  ## Primal: G_RR += 2 bar_QQ RR Sigma_e
  G_RR  <- G_RR  + 2 * bar_QQ  %*% RR  %*% Sigma_e
  dG_RR <- dG_RR + 2 * (dbar_QQ %*% RR  %*% Sigma_e  +
                         bar_QQ  %*% dRR %*% Sigma_e  +
                         bar_QQ  %*% RR  %*% dSig)

  list(dG_TT = dG_TT, dG_RR = dG_RR, dG_ZZ = dG_ZZ,
       dG_DD = dG_DD, dg_d  = dg_d,  dG_Sig = dG_Sig)
}
