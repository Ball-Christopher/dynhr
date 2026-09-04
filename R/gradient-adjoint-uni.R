## R/gradient-adjoint-uni.R
## --------------------------------------------------------------------------
## Adjoint (reverse-mode) Kalman-filter likelihood gradient with MISSING-DATA
## support and an optional caller-supplied initial covariance P0.

#' TRUE when the compiled missing-data univariate adjoint is available.
#' @noRd
.HAS_RCPP_KF_ADJOINT_UNI <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("kf_adjoint_uni_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}
##
## This is the dense multivariate adjoint (.kf_loglik_adjoint, R reference)
## extended with PER-PERIOD OBSERVED-ROW SUBSETTING -- the textbook
## multivariate treatment of missing observations: at period t only the
## observed rows O_t of (Y, ZZ, DD, d) enter the innovation, gain, and
## log-likelihood; the Joseph covariance update uses A_t = TT - K_t ZZ[O_t,],
## B_t = RR - K_t DD[O_t,]. A fully-missing period is a pure prediction step
## (K_t = 0). The backward sweep mirrors the dense one row-for-row, scattering
## G_ZZ / G_DD / g_d contributions back to the observed rows O_t only.
##
## When all observations are present this reduces EXACTLY to .kf_loglik_adjoint
## (validated to agree to ~1e-10).
##
## P0: when NULL, the stationary P0 = solve_lyapunov(TT, RR Sigma_e RR') is used
## and the Lyapunov adjoint is included. When supplied (e.g. a big-kappa
## approximate-diffuse P0 for a near/exact unit-root model), P0 is treated as a
## FIXED parameter-independent input -- the Lyapunov adjoint term is omitted and
## P0's gradient contribution is zero by construction.
##
## shock_scale (n_exo x T multiplicative shock-std factors, heteroskedastic
## shocks / SV conditional path): supported by the R reference since
## 0.9.0.0008. Per period the effective shock covariance is
## Se_t = diag(sc_t) Sigma_e diag(sc_t) with sc_t FIXED DATA, so every
## backward-sweep accumulation into the Sigma_e gradient is sandwiched
## elementwise by outer(sc_t, sc_t) (chain rule through Se_t), and every
## date-t use of Sigma_e (HHo/SSo/QQ_t/covariance update, and the Se factors
## feeding G_RR/G_DD) substitutes Se_t. The stationary Lyapunov P0 stays on
## the BASELINE Sigma_e -- matching kalman_filter's lik_init = "stationary"
## convention under shock_scale -- so the Lyapunov adjoint contributions are
## NOT sandwiched. Mirrors the dense adjoint's shock_scale treatment
## (R/gradient-adjoint-kf.R). The compiled fast path (kf_adjoint_uni_cpp)
## does not take shock_scale; the R reference runs instead.
##
## Scope: scalar me_variance; no me_extra (use the dense adjoint for that).
## --------------------------------------------------------------------------

#' Adjoint KF gradient with missing-data + optional supplied P0
#'
#' @param Y          n_obs x n_T matrix; NA entries are skipped per period.
#' @param ss         list TT, RR, ZZ, DD, d, Sigma_e.
#' @param d_ss_list  list length n_par; element j has any of dTT,dRR,dZZ,dDD,
#'                   dd,dSigma_e.
#' @param me_variance scalar measurement-error variance.
#' @param P0         optional n_state x n_state initial state covariance. NULL
#'                   => stationary Lyapunov P0 (gradient included). Supplied =>
#'                   fixed input (no Lyapunov adjoint).
#' @param shock_scale optional n_exo x n_T matrix of multiplicative shock-std
#'                   factors (fixed data; see the file header). NULL = baseline.
#' @return list(loglik, grad).
#' @noRd
.kf_loglik_adjoint_uni <- function(Y, ss, d_ss_list, me_variance = 0,
                                   P0 = NULL, shock_scale = NULL) {
  TT <- ss$TT; RR <- ss$RR; ZZ <- ss$ZZ; DD <- ss$DD
  d  <- as.numeric(ss$d); Sigma_e <- ss$Sigma_e

  n_state <- nrow(TT); n_obs <- nrow(ZZ); n_exo <- ncol(RR)
  n_par <- length(d_ss_list)

  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)

  has_sc <- !is.null(shock_scale)
  if (has_sc && (!is.matrix(shock_scale) || nrow(shock_scale) != n_exo ||
                 ncol(shock_scale) != n_T))
    stop(".kf_loglik_adjoint_uni: shock_scale must be n_exo x n_T.",
         call. = FALSE)
  ## Per-period effective shock covariance (baseline when no scaling).
  Se_at <- if (has_sc) {
    function(t) Sigma_e * outer(shock_scale[, t], shock_scale[, t])
  } else {
    function(t) Sigma_e
  }

  QQ <- tcrossprod(RR %*% Sigma_e, RR)
  ll_2pi <- log(2 * pi)
  fail <- list(loglik = -Inf, grad = rep(NA_real_, n_par))
  p0_supplied <- !is.null(P0)

  if (p0_supplied) {
    if (!all(is.finite(P0))) return(fail)
    P0_use <- P0
  } else {
    P0_use <- .solve_lyapunov(TT, QQ)
    if (!all(is.finite(P0_use))) return(fail)
  }

  ## Sandwich a Sigma_e-gradient accumulation through Se_t = S_t Sigma_e S_t
  ## (chain rule with S_t = diag(sc_t) fixed data); identity when unscaled.
  sand <- function(M, t) {
    if (has_sc) M * outer(shock_scale[, t], shock_scale[, t]) else M
  }

  ## -- Fast path: compiled univariate adjoint (kf_adjoint_uni_cpp) ----------
  ## (not extended for shock_scale -- the R reference below handles it)
  if (.HAS_RCPP_KF_ADJOINT_UNI() && !has_sc) {
    zTT0  <- matrix(0, n_state, n_state)
    zRR0  <- matrix(0, n_state, n_exo)
    zZZ0  <- matrix(0, n_obs, n_state)
    zDD0  <- matrix(0, n_obs, n_exo)
    zSig0 <- matrix(0, n_exo, n_exo)

    dTT_cube    <- array(0, c(n_state, n_state, n_par))
    dRR_cube    <- array(0, c(n_state, n_exo,   n_par))
    dZZ_cube    <- array(0, c(n_obs,   n_state, n_par))
    dDD_cube    <- array(0, c(n_obs,   n_exo,   n_par))
    dd_mat      <- matrix(0, n_obs, n_par)
    dSigma_cube <- array(0, c(n_exo, n_exo, n_par))

    for (j in seq_len(n_par)) {
      dpar <- d_ss_list[[j]]
      if (is.null(dpar)) dpar <- list()
      dTT_cube[, , j]    <- .dpiece2(dpar, "dTT",      zTT0)
      dRR_cube[, , j]    <- .dpiece2(dpar, "dRR",      zRR0)
      dZZ_cube[, , j]    <- .dpiece2(dpar, "dZZ",      zZZ0)
      dDD_cube[, , j]    <- .dpiece2(dpar, "dDD",      zDD0)
      dd_mat[, j]        <- .dpiece2(dpar, "dd",       numeric(n_obs))
      dSigma_cube[, , j] <- .dpiece2(dpar, "dSigma_e", zSig0)
    }

    ## Pass P0_use (already computed above: either supplied or Lyapunov).
    ## p0_supplied tells the C++ kernel whether to run the Lyapunov adjoint.
    out <- kf_adjoint_uni_cpp(Y, TT, RR, ZZ, DD, d, Sigma_e,
                              dTT_cube, dRR_cube, dZZ_cube, dDD_cube, dd_mat,
                              dSigma_cube, me_variance, .KF_LL_MIN,
                              P0_use, p0_supplied)

    if (!isTRUE(out$ok)) {
      return(list(loglik = -Inf, grad = rep(NA_real_, n_par)))
    }
    return(list(loglik = out$loglik, grad = as.numeric(out$grad)))
  }

  ## -- Forward pass with per-period observed-row subsetting -----------------
  s_store  <- vector("list", n_T); P_store <- vector("list", n_T)
  v_store  <- vector("list", n_T); Fi_store <- vector("list", n_T)
  K_store  <- vector("list", n_T); A_store <- vector("list", n_T)
  B_store  <- vector("list", n_T); O_store <- vector("list", n_T)

  s <- numeric(n_state); P <- P0_use; loglik <- 0.0

  for (t in seq_len(n_T)) {
    s_store[[t]] <- s; P_store[[t]] <- P
    O <- which(is.finite(Y[, t]))   # observed rows this period
    O_store[[t]] <- O
    q <- length(O)
    Se_t <- Se_at(t)

    if (q == 0L) {
      ## Pure prediction (K = 0): A = TT, B = RR.
      A_store[[t]] <- NULL; B_store[[t]] <- NULL
      s <- as.numeric(TT %*% s)
      QQ_t <- if (has_sc) tcrossprod(RR %*% Se_t, RR) else QQ
      P <- .sym(tcrossprod(TT %*% P, TT) + QQ_t)
      next
    }

    ZZo <- ZZ[O, , drop = FALSE]; DDo <- DD[O, , drop = FALSE]
    HHo <- tcrossprod(DDo %*% Se_t, DDo)
    SSo <- RR %*% Se_t %*% t(DDo)
    me_o <- me_variance * diag(q)

    PZ <- P %*% t(ZZo)
    Ft <- .sym(ZZo %*% PZ + HHo + me_o)
    Fc <- tryCatch(chol(Ft), error = function(e) NULL)
    if (is.null(Fc)) return(fail)
    Fi <- chol2inv(Fc); ldf <- 2 * sum(log(diag(Fc)))

    v   <- Y[O, t] - d[O] - as.numeric(ZZo %*% s)
    Fiv <- as.numeric(Fi %*% v)
    ll_t <- -0.5 * (q * ll_2pi + ldf + sum(v * Fiv))
    if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) return(fail)
    loglik <- loglik + ll_t

    K <- (TT %*% PZ + SSo) %*% Fi
    A <- TT - K %*% ZZo
    B <- RR - K %*% DDo

    Fi_store[[t]] <- Fi; v_store[[t]] <- v; K_store[[t]] <- K
    A_store[[t]] <- A;   B_store[[t]] <- B

    s <- as.numeric(TT %*% s) + as.numeric(K %*% v)
    P_raw <- tcrossprod(A %*% P, A) + tcrossprod(B %*% Se_t, B)
    ## TRUE measurement-noise law (F3-D): P' += K me_o K'.
    if (me_variance != 0) P_raw <- P_raw + me_variance * tcrossprod(K)
    P <- .sym(P_raw)
  }

  ## -- Backward sweep -------------------------------------------------------
  bar_s <- numeric(n_state); bar_P <- matrix(0, n_state, n_state)
  G_TT <- matrix(0, n_state, n_state); G_RR <- matrix(0, n_state, n_exo)
  G_ZZ <- matrix(0, n_obs, n_state);   G_DD <- matrix(0, n_obs, n_exo)
  g_d  <- numeric(n_obs);              G_Sig <- matrix(0, n_exo, n_exo)

  for (t in n_T:1) {
    s_prev <- s_store[[t]]; P_prev <- P_store[[t]]
    O <- O_store[[t]]; q <- length(O)
    bar_P <- .sym(bar_P)
    Se_t <- Se_at(t)

    if (q == 0L) {
      ## Pure-prediction adjoint: A = TT, B = RR, no measurement terms.
      bar_A <- 2 * bar_P %*% TT %*% P_prev
      bar_B <- 2 * bar_P %*% RR %*% Se_t
      bar_P_prev <- t(TT) %*% bar_P %*% TT
      G_Sig <- G_Sig + sand(t(RR) %*% bar_P %*% RR, t)
      G_TT  <- G_TT + outer(bar_s, s_prev)
      bar_s_prev <- as.numeric(t(TT) %*% bar_s)
      G_TT  <- G_TT + bar_A
      G_RR  <- G_RR + bar_B
      bar_s <- bar_s_prev
      bar_P <- .sym(bar_P_prev)
      next
    }

    ZZo <- ZZ[O, , drop = FALSE]; DDo <- DD[O, , drop = FALSE]
    v  <- v_store[[t]]; Fi <- Fi_store[[t]]; K <- K_store[[t]]
    A  <- A_store[[t]]; B <- B_store[[t]]
    Fiv <- as.numeric(Fi %*% v)

    ## Step 1: P_t = sym(A P_prev A' + B Se_t B')
    bar_A <- 2 * bar_P %*% A %*% P_prev
    bar_B <- 2 * bar_P %*% B %*% Se_t
    bar_P_prev_AP <- t(A) %*% bar_P %*% A
    G_Sig <- G_Sig + sand(t(B) %*% bar_P %*% B, t)

    ## Step 2: s_t = TT s_prev + K v
    G_TT  <- G_TT + outer(bar_s, s_prev)
    bar_K <- outer(bar_s, v)
    ## Adjoint of the ME Joseph term P_t += K me_o K' (me_o is DATA):
    ## d tr(bar_P K me K') / dK = 2 me bar_P K.
    if (me_variance != 0)
      bar_K <- bar_K + 2 * me_variance * (bar_P %*% K)
    bar_v <- as.numeric(t(K) %*% bar_s)
    bar_s_prev <- as.numeric(t(TT) %*% bar_s)

    ## Step 3: ll_t
    bar_F <- -0.5 * (Fi - outer(Fiv, Fiv))
    bar_v <- bar_v - Fiv

    ## Step 4: A = TT - K ZZo
    G_TT <- G_TT + bar_A
    G_ZZ[O, ] <- G_ZZ[O, ] - t(K) %*% bar_A
    bar_K <- bar_K - bar_A %*% t(ZZo)

    ## Step 5: B = RR - K DDo
    G_RR <- G_RR + bar_B
    G_DD[O, ] <- G_DD[O, ] - t(K) %*% bar_B
    bar_K <- bar_K - bar_B %*% t(DDo)

    ## Step 6: K = (TT P_prev ZZo' + SSo) Fi, SSo = RR Se_t DDo'
    bar_F <- bar_F - t(K) %*% bar_K %*% Fi
    bar_F <- .sym(bar_F)
    bar_Mnum <- bar_K %*% Fi
    G_TT <- G_TT + bar_Mnum %*% ZZo %*% P_prev
    G_ZZ[O, ] <- G_ZZ[O, ] + t(bar_Mnum) %*% TT %*% P_prev
    bar_P_prev_Mnum <- t(TT) %*% bar_Mnum %*% ZZo
    G_RR <- G_RR + bar_Mnum %*% DDo %*% Se_t
    G_DD[O, ] <- G_DD[O, ] + t(bar_Mnum) %*% RR %*% Se_t
    G_Sig <- G_Sig + sand(t(RR) %*% bar_Mnum %*% DDo, t)

    ## Step 7: F = sym(ZZo P_prev ZZo' + HHo + me), HHo = DDo Se_t DDo'
    bar_P_prev_F <- t(ZZo) %*% bar_F %*% ZZo
    G_ZZ[O, ] <- G_ZZ[O, ] + 2 * bar_F %*% ZZo %*% P_prev
    G_DD[O, ] <- G_DD[O, ] + 2 * bar_F %*% DDo %*% Se_t
    G_Sig <- G_Sig + sand(t(DDo) %*% bar_F %*% DDo, t)

    ## Step 8: v = y[O] - d[O] - ZZo s_prev
    g_d[O] <- g_d[O] - bar_v
    G_ZZ[O, ] <- G_ZZ[O, ] - outer(bar_v, s_prev)
    bar_s_prev <- bar_s_prev - as.numeric(t(ZZo) %*% bar_v)

    ## Step 9: collect
    bar_s <- bar_s_prev
    bar_P <- .sym(bar_P_prev_AP + bar_P_prev_Mnum + bar_P_prev_F)
  }

  ## -- Lyapunov adjoint (only when P0 computed internally) ------------------
  bar_P0 <- .sym(bar_P)
  if (!p0_supplied) {
    if (n_state == 1L) {
      denom <- 1 - TT[1, 1]^2
      bar_QQ <- matrix(bar_P0[1, 1] / denom, 1, 1)
      G_TT[1, 1] <- G_TT[1, 1] + bar_P0[1, 1] * 2 * TT[1, 1] * P0_use[1, 1] / denom
    } else {
      bar_QQ <- tryCatch(.solve_lyapunov(t(TT), bar_P0), error = function(e) NULL)
      if (is.null(bar_QQ) || !all(is.finite(bar_QQ))) return(fail)
      bar_QQ <- .sym(bar_QQ)
      G_TT <- G_TT + 2 * bar_QQ %*% TT %*% P0_use
    }
    G_Sig <- G_Sig + t(RR) %*% bar_QQ %*% RR
    G_RR  <- G_RR + 2 * bar_QQ %*% RR %*% Sigma_e
  }

  ## -- Contraction ----------------------------------------------------------
  grad <- numeric(n_par)
  for (j in seq_len(n_par)) {
    dpar <- d_ss_list[[j]]
    if (is.null(dpar)) { grad[j] <- 0; next }
    gj <- 0
    if (!is.null(dpar$dTT))      gj <- gj + sum(G_TT  * dpar$dTT)
    if (!is.null(dpar$dRR))      gj <- gj + sum(G_RR  * dpar$dRR)
    if (!is.null(dpar$dZZ))      gj <- gj + sum(G_ZZ  * dpar$dZZ)
    if (!is.null(dpar$dDD))      gj <- gj + sum(G_DD  * dpar$dDD)
    if (!is.null(dpar$dd))       gj <- gj + sum(g_d   * dpar$dd)
    if (!is.null(dpar$dSigma_e)) gj <- gj + sum(G_Sig * dpar$dSigma_e)
    grad[j] <- gj
  }
  list(loglik = loglik, grad = grad)
}
