## R/gradient-adjoint-diffuse.R
## --------------------------------------------------------------------------
## Adjoint (reverse-mode) gradient of the EXACT-DIFFUSE Kalman-filter
## log-likelihood (Durbin-Koopman 2012, ch.5).
##
## Differentiates through:
##   Stage 1: the diffuse-phase recursion (Case A / Case B) and the
##            stationary tail.
##   Stage 2: the Schur-based initialisation (P_inf_0 = U J U' and
##            P_star_0 = U_s Pa_ss U_s'). Both are invariant to within-block
##            Schur rotations, so the only hidden variable is the inter-block
##            coupling Omega_su (Sylvester T_ss Omega_su - Omega_su T_uu = -G_su).
##            The exact init-gradient is assembled by the analytic FORWARD
##            directional derivative of (P_inf_0, P_star_0) contracted with
##            (bar_P_inf_0, bar_P_star_0) -- a one-time O(n^2) Jacobian build.
##
## The full gradient is exact for ALL parameters, including TT-movers and
## general non-triangular TT (U != I), validated vs numDeriv to ~1e-9. NOTE:
## the exact-diffuse loglik is differentiable only in UNIT-ROOT-PRESERVING TT
## directions; a perturbation that moves an eigenvalue off the unit circle
## changes the diffuse classification (nunit) and is genuinely non-smooth. The
## adjoint returns the correct within-regime derivative. In DSGE estimation the
## unit roots are structural (held by the model), so estimated parameters move
## TT within the smooth regime.
##
## State-space convention (identical to gradient-adjoint-kf.R):
##   QQ = RR Sigma_e RR',   HH = DD Sigma_e DD',   SS = RR Sigma_e DD'
##
## Diffuse init:
##   P_inf_0  = A_inf A_inf'   (A_inf = U[:,1:nunit], U from ordered Schur)
##   P_star_0 = U Pa U'  (Pa_{s,s} = lyap(Ts_s, QQ_s); zero elsewhere)
##
## Diffuse Case B (F_inf nonsingular):
##   v      = y - d - ZZ s
##   F_inf  = ZZ P_inf ZZ'        (symmetrised)
##   F_star = ZZ P_star ZZ' + HH + me_diag
##   M_inf  = TT P_inf  ZZ'
##   M_star = TT P_star ZZ' + SS
##   K0     = M_inf F_inf^{-1}
##   ll_t   = -0.5 log|F_inf|
##   s'       = TT s + K0 v
##   P_inf'   = sym(TT P_inf  TT' - K0 F_inf  K0')
##   P_star'  = sym(TT P_star TT' + QQ - K0 M_star' - M_star K0' + K0 F_star K0')
##
## Diffuse Case A (F_inf ~ 0):
##   standard .kf_step on (s, P_star); propagate P_inf <- sym(TT P_inf TT').
##
## Adjoint notation: bar_X = d loglik / d X  (shape same as X).
## --------------------------------------------------------------------------

#' TRUE when the compiled exact-diffuse adjoint is available.
#' @noRd
.HAS_RCPP_KF_DIFFUSE <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("kf_adjoint_diffuse_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

#' Adjoint gradient of the exact-diffuse Kalman-filter log-likelihood
#'
#' Forward pass runs the exact-diffuse Kalman filter.  Backward sweep
#' accumulates gradient matrices (G_TT, G_RR, G_ZZ, G_DD, g_d, G_Sig)
#' and contracts them against d_ss_list.
#'
#' @param Y          n_obs x n_T observation matrix (no NAs during diffuse phase).
#' @param ss         list: TT, RR, ZZ, DD, d, Sigma_e.
#' @param d_ss_list  list of length n_par; each element has (any subset of)
#'                   dTT, dRR, dZZ, dDD, dd, dSigma_e.
#' @param me_variance scalar ME variance (parameter-independent).
#' @param ur_tol     unit-root tolerance passed to .kf_diffuse_P0.
#'
#' @return list(loglik, grad, stage1_ok, stage2_ok, bar_P_inf_0,
#'             bar_P_star_0, nunit, d_diffuse, min_eig_margin,
#'             near_regime_boundary). \code{min_eig_margin} is the distance of
#'             the nearest stable (non-unit) eigenvalue of \code{TT} from the
#'             unit circle; \code{near_regime_boundary} is \code{TRUE} (and a
#'             warning is emitted) when that margin is below \code{100 * ur_tol},
#'             signalling that the fixed-\code{nunit} gradient is near a
#'             non-smooth unit-root regime change.
#' @noRd
.kf_loglik_adjoint_diffuse <- function(Y, ss, d_ss_list, me_variance = 0,
                                       ur_tol = 1e-6) {
  TT      <- ss$TT;  RR <- ss$RR;  ZZ <- ss$ZZ;  DD <- ss$DD
  d_obs   <- as.numeric(ss$d);     Sigma_e <- ss$Sigma_e

  n_state <- nrow(TT)
  n_obs   <- nrow(ZZ)
  n_exo   <- ncol(RR)
  n_par   <- length(d_ss_list)

  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)

  ## ---- Regime-boundary detection (computed R-side regardless of fast path) --
  ## Must be computed before the fast-path return so the warning is always
  ## emitted and min_eig_margin / near_regime_boundary are always returned.
  ev_mag_fp      <- abs(eigen(TT, only.values = TRUE)$values)
  dist_uc_fp     <- abs(ev_mag_fp - 1)
  ## We need nunit to filter stable eigenvalues. Defer until after we have it
  ## (from C++ or from .kf_diffuse_P0 below). Captured in a helper closure:
  .check_regime_boundary <- function(nunit_val) {
    ur_mags  <- ev_mag_fp[dist_uc_fp < ur_tol]
    all_dist <- dist_uc_fp[dist_uc_fp >= ur_tol]
    mem      <- if (length(all_dist)) min(all_dist) else Inf
    caution  <- 100 * ur_tol
    nrb      <- is.finite(mem) && mem < caution
    if (nrb)
      warning(sprintf(paste0(
        ".kf_loglik_adjoint_diffuse: a stable eigenvalue is within %.1e of the ",
        "unit circle (margin %.2e, nunit = %d). The exact-diffuse gradient is ",
        "non-smooth across a unit-root regime change; a small parameter ",
        "perturbation may move this eigenvalue across the |lambda| = 1 boundary ",
        "and change nunit."), caution, mem, nunit_val), call. = FALSE)
    list(min_eig_margin = mem, near_regime_boundary = nrb)
  }

  ## ---- Fast path: compiled exact-diffuse adjoint recursion -----------------
  if (.HAS_RCPP_KF_DIFFUSE()) {
    cpp_out <- kf_adjoint_diffuse_cpp(Y, TT, RR, ZZ, DD, d_obs, Sigma_e,
                                      me_variance, ur_tol)

    if (!isTRUE(cpp_out$ok)) {
      ## Return fail-list with regime-boundary fields set from R-side check.
      ## nunit unknown on failure; use NA and skip warning.
      return(list(loglik = -Inf, grad = rep(NA_real_, n_par),
                  stage1_ok = FALSE, stage2_ok = FALSE,
                  bar_P_inf_0 = NULL, bar_P_star_0 = NULL,
                  nunit = NA_integer_, d_diffuse = NA_integer_,
                  min_eig_margin = NA_real_, near_regime_boundary = NA))
    }

    nunit_cpp <- cpp_out$nunit
    rb_info   <- .check_regime_boundary(nunit_cpp)   # may emit warning

    ## Final contraction (R-side, same as reference implementation).
    G_TT  <- cpp_out$G_TT;  G_RR <- cpp_out$G_RR;  G_ZZ <- cpp_out$G_ZZ
    G_DD  <- cpp_out$G_DD;  g_d  <- as.numeric(cpp_out$g_d)
    G_Sig <- cpp_out$G_Sig

    grad <- numeric(n_par)
    for (j in seq_len(n_par)) {
      dpar <- d_ss_list[[j]]
      if (is.null(dpar)) { grad[j] <- 0; next }
      gj <- 0
      dTT  <- dpar[["dTT"]];      if (!is.null(dTT))  gj <- gj + sum(G_TT  * dTT)
      dRR  <- dpar[["dRR"]];      if (!is.null(dRR))  gj <- gj + sum(G_RR  * dRR)
      dZZ  <- dpar[["dZZ"]];      if (!is.null(dZZ))  gj <- gj + sum(G_ZZ  * dZZ)
      dDD  <- dpar[["dDD"]];      if (!is.null(dDD))  gj <- gj + sum(G_DD  * dDD)
      dd_j <- dpar[["dd"]];       if (!is.null(dd_j)) gj <- gj + sum(g_d   * dd_j)
      dSig <- dpar[["dSigma_e"]]; if (!is.null(dSig)) gj <- gj + sum(G_Sig * dSig)
      grad[j] <- gj
    }

    return(list(
      loglik       = cpp_out$loglik,
      grad         = grad,
      stage1_ok    = isTRUE(cpp_out$stage1_ok),
      stage2_ok    = isTRUE(cpp_out$stage2_ok),
      bar_P_inf_0  = cpp_out$bar_P_inf_0,
      bar_P_star_0 = cpp_out$bar_P_star_0,
      nunit        = nunit_cpp,
      d_diffuse    = cpp_out$d_diffuse,
      min_eig_margin       = rb_info$min_eig_margin,
      near_regime_boundary = rb_info$near_regime_boundary
    ))
  }

  QQ      <- tcrossprod(RR %*% Sigma_e, RR)
  HH      <- tcrossprod(DD %*% Sigma_e, DD)
  SS      <- RR %*% Sigma_e %*% t(DD)
  tZZ     <- t(ZZ)
  me_diag <- me_variance * diag(n_obs)
  HHme    <- HH + me_diag
  ll_const <- -0.5 * n_obs * log(2 * pi)

  fail <- list(loglik = -Inf, grad = rep(NA_real_, n_par),
               stage1_ok = FALSE, stage2_ok = FALSE,
               bar_P_inf_0 = NULL, bar_P_star_0 = NULL,
               nunit = NA_integer_, d_diffuse = NA_integer_,
               min_eig_margin = NA_real_, near_regime_boundary = NA)

  ## ---- Diffuse initialization (Schur decomposition) -----------------------
  P0_list <- .kf_diffuse_P0(TT, QQ, ur_tol)
  nunit    <- P0_list$nunit

  ## Regime-boundary detection (Tier 11 #4) -- reuse pre-computed eigenvalues.
  rb_info_r       <- .check_regime_boundary(nunit)
  min_eig_margin  <- rb_info_r$min_eig_margin
  near_regime_boundary <- rb_info_r$near_regime_boundary
  P_inf_0  <- P0_list$P_inf
  P_star_0 <- P0_list$P_star

  ## ---- Forward pass -------------------------------------------------------
  ## Tolerances mirroring .kf_diffuse_phase
  diffuse_tol <- 1e-10
  conv_tol    <- 1e-8
  cap         <- min(n_T, 100L)

  ## Per-step store (list of per-step records)
  step_store <- vector("list", n_T)

  ## Inline step function for Case A / stationary tail
  .do_stat_step <- function(s_in, P_in, v_in) {
    PZ  <- P_in %*% tZZ
    Ft  <- .sym(ZZ %*% PZ + HHme)
    Fc  <- tryCatch(chol(Ft), error = function(e) NULL)
    if (is.null(Fc)) return(NULL)
    Fi  <- chol2inv(Fc)
    ldf <- 2 * sum(log(diag(Fc)))
    Fiv <- as.numeric(Fi %*% v_in)
    ll  <- ll_const - 0.5 * (ldf + sum(v_in * Fiv))
    if (!is.finite(ll) || ll < .KF_LL_MIN) return(NULL)
    K   <- (TT %*% PZ + SS) %*% Fi
    A   <- TT - K %*% ZZ
    B   <- RR - K %*% DD
    s_n <- as.numeric(TT %*% s_in) + as.numeric(K %*% v_in)
    P_n <- .sym(tcrossprod(A %*% P_in, A) + tcrossprod(B %*% Sigma_e, B))
    list(ll = ll, s = s_n, P = P_n, K = K, Fi = Fi, A = A, B = B, Fiv = Fiv)
  }

  s      <- numeric(n_state)
  P_inf  <- P_inf_0
  P_star <- P_star_0
  loglik <- 0.0
  d_diffuse <- NA_integer_
  in_diffuse <- (nunit > 0L)

  for (t in seq_len(n_T)) {
    y_t <- Y[, t]
    if (anyNA(y_t)) return(fail)
    v <- y_t - d_obs - as.numeric(ZZ %*% s)

    if (in_diffuse) {
      F_inf  <- .sym(ZZ %*% P_inf  %*% tZZ)
      F_star <- .sym(ZZ %*% P_star %*% tZZ + HHme)
      scale_star <- max(1, max(abs(F_star)))

      if (max(abs(F_inf)) < diffuse_tol * scale_star) {
        ## -- Case A -----------------------------------------------------------
        step <- .do_stat_step(s, P_star, v)
        if (is.null(step)) return(fail)
        loglik <- loglik + step$ll
        P_inf_new <- .sym(tcrossprod(TT %*% P_inf, TT))
        step_store[[t]] <- list(
          case = "A",
          s_prev = s, P_star_prev = P_star, P_inf_prev = P_inf,
          v = v, Fi = step$Fi, K = step$K, A = step$A, B = step$B,
          Fiv = step$Fiv
        )
        s      <- step$s
        P_star <- step$P
        P_inf  <- P_inf_new
      } else {
        ## -- Case B -----------------------------------------------------------
        rc     <- tryCatch(rcond(F_inf), error = function(e) 0)
        Fc_inf <- tryCatch(chol(F_inf), error = function(e) NULL)
        if (is.null(Fc_inf) || !is.finite(rc) || rc <= diffuse_tol) return(fail)

        F_inf_inv <- chol2inv(Fc_inf)
        ll_t      <- -0.5 * 2 * sum(log(diag(Fc_inf)))
        if (!is.finite(ll_t)) return(fail)
        loglik <- loglik + ll_t

        M_inf  <- TT %*% P_inf  %*% tZZ
        M_star <- TT %*% P_star %*% tZZ + SS
        K0     <- M_inf %*% F_inf_inv

        s_new <- as.numeric(TT %*% s) + as.numeric(K0 %*% v)
        P_inf_new  <- .sym(tcrossprod(TT %*% P_inf, TT) -
                             tcrossprod(K0 %*% F_inf, K0))
        P_star_new <- .sym(tcrossprod(TT %*% P_star, TT) + QQ -
                             K0 %*% t(M_star) - M_star %*% t(K0) +
                             tcrossprod(K0 %*% F_star, K0))

        step_store[[t]] <- list(
          case = "B",
          s_prev = s, P_star_prev = P_star, P_inf_prev = P_inf,
          v = v, F_inf = F_inf, F_inf_inv = F_inf_inv,
          F_star = F_star, M_inf = M_inf, M_star = M_star, K0 = K0
        )
        s      <- s_new
        P_inf  <- P_inf_new
        P_star <- P_star_new
      }

      if (max(abs(P_inf)) < conv_tol * max(1, max(abs(P_star)))) {
        d_diffuse  <- t
        in_diffuse <- FALSE
      } else if (t >= cap) {
        warning(".kf_loglik_adjoint_diffuse: P_inf did not converge within ",
                cap, " periods.")
        return(fail)
      }

    } else {
      ## -- Stationary tail ---------------------------------------------------
      step <- .do_stat_step(s, P_star, v)
      if (is.null(step)) return(fail)
      loglik <- loglik + step$ll
      step_store[[t]] <- list(
        case = "stat",
        s_prev = s, P_star_prev = P_star,
        v = v, Fi = step$Fi, K = step$K, A = step$A, B = step$B,
        Fiv = step$Fiv
      )
      s      <- step$s
      P_star <- step$P
    }
  }

  ## ---- Backward sweep ----------------------------------------------------- #
  bar_s      <- numeric(n_state)
  bar_P_star <- matrix(0, n_state, n_state)
  bar_P_inf  <- matrix(0, n_state, n_state)

  G_TT  <- matrix(0, n_state, n_state)
  G_RR  <- matrix(0, n_state, n_exo)
  G_ZZ  <- matrix(0, n_obs,   n_state)
  G_DD  <- matrix(0, n_obs,   n_exo)
  g_d   <- numeric(n_obs)
  G_Sig <- matrix(0, n_exo,   n_exo)

  for (t in n_T:1) {
    st <- step_store[[t]]
    if (is.null(st)) next

    if (st$case == "stat" || st$case == "A") {
      ## ---- Standard KF step adjoint (operates on bar_P_star) ----------------
      ## This is the dense adjoint from gradient-adjoint-kf.R, applied to P_star.
      s_prev <- st$s_prev
      P_prev <- st$P_star_prev
      v      <- st$v
      Fi     <- st$Fi
      K      <- st$K
      A      <- st$A
      B      <- st$B
      Fiv    <- st$Fiv

      ## Adjoint of P_t = sym(A P_prev A' + B Sig B')
      bar_P_star <- .sym(bar_P_star)
      bar_A      <- 2 * bar_P_star %*% A %*% P_prev
      bar_B      <- 2 * bar_P_star %*% B %*% Sigma_e
      bar_Pprev_AP <- t(A) %*% bar_P_star %*% A
      G_Sig <- G_Sig + t(B) %*% bar_P_star %*% B

      ## Adjoint of s_t = TT s_prev + K v
      G_TT       <- G_TT + outer(bar_s, s_prev)
      bar_K      <- outer(bar_s, v)
      bar_v      <- as.numeric(t(K) %*% bar_s)
      bar_s_prev <- as.numeric(t(TT) %*% bar_s)

      ## Adjoint of ll_t = -0.5*(log|F| + v'Fi v)
      bar_F <- -0.5 * (Fi - outer(Fiv, Fiv))
      bar_v <- bar_v - Fiv

      ## Adjoint of A = TT - K ZZ
      G_TT  <- G_TT  + bar_A
      G_ZZ  <- G_ZZ  - t(K) %*% bar_A
      bar_K <- bar_K - bar_A %*% tZZ

      ## Adjoint of B = RR - K DD
      G_RR  <- G_RR  + bar_B
      G_DD  <- G_DD  - t(K) %*% bar_B
      bar_K <- bar_K - bar_B %*% t(DD)

      ## Adjoint of K = (TT P_prev ZZ' + SS) Fi
      bar_F    <- bar_F - t(K) %*% bar_K %*% Fi
      bar_F    <- .sym(bar_F)
      bar_Mnum <- bar_K %*% Fi

      G_TT       <- G_TT + bar_Mnum %*% ZZ %*% P_prev
      G_ZZ       <- G_ZZ + t(bar_Mnum) %*% TT %*% P_prev
      bar_Pprev_Mnum <- t(TT) %*% bar_Mnum %*% ZZ
      G_RR  <- G_RR  + bar_Mnum %*% DD %*% Sigma_e
      G_DD  <- G_DD  + t(bar_Mnum) %*% RR %*% Sigma_e
      G_Sig <- G_Sig + t(RR) %*% bar_Mnum %*% DD

      ## Adjoint of F = sym(ZZ P_prev ZZ' + HH + me_diag)
      bar_Pprev_F <- t(ZZ) %*% bar_F %*% ZZ
      G_ZZ  <- G_ZZ  + 2 * bar_F %*% ZZ %*% P_prev
      G_DD  <- G_DD  + 2 * bar_F %*% DD %*% Sigma_e
      G_Sig <- G_Sig + t(DD) %*% bar_F %*% DD

      ## Adjoint of v = y - d - ZZ s_prev
      g_d        <- g_d - bar_v
      G_ZZ       <- G_ZZ - outer(bar_v, s_prev)
      bar_s_prev <- bar_s_prev - as.numeric(tZZ %*% bar_v)

      bar_P_star_new <- .sym(bar_Pprev_AP + bar_Pprev_Mnum + bar_Pprev_F)

      ## In Case A: also propagate bar_P_inf through P_inf <- TT P_inf TT'
      if (st$case == "A") {
        P_inf_prev    <- st$P_inf_prev
        ## Adjoint of P_inf <- sym(TT P_inf TT'): bar_P_inf_prev = TT' bar_P_inf TT
        bar_P_inf_sym <- .sym(bar_P_inf)
        ## G_TT from this: d/dTT tr(bar_P_inf TT P_inf_prev TT') = 2 tr(bar_P_inf TT P_inf_prev dTT')
        ##  => G_TT += 2 bar_P_inf TT P_inf_prev  (sym adjoint)
        G_TT      <- G_TT + 2 * bar_P_inf_sym %*% TT %*% P_inf_prev
        bar_P_inf <- t(TT) %*% bar_P_inf_sym %*% TT
      }

      bar_s      <- bar_s_prev
      bar_P_star <- bar_P_star_new

    } else if (st$case == "B") {
      ## ---- Case B adjoint ---------------------------------------------------
      ## Forward quantities stored:
      s_prev      <- st$s_prev
      P_star_prev <- st$P_star_prev
      P_inf_prev  <- st$P_inf_prev
      v           <- st$v
      F_inf       <- st$F_inf
      Fi          <- st$F_inf_inv   # F_inf^{-1}
      F_star      <- st$F_star
      M_inf       <- st$M_inf
      M_star      <- st$M_star
      K0          <- st$K0
      tK0         <- t(K0)

      ## We unpack the accumulated bar_P_star as bar_{P_star_t} (the new P_star)
      ## and bar_P_inf as bar_{P_inf_t}.
      bar_Pstar_t <- .sym(bar_P_star)  # adjoint of P_star AFTER this step
      bar_Pinf_t  <- .sym(bar_P_inf)   # adjoint of P_inf  AFTER this step

      ## ----- 1. Adjoint of P_star' = sym(TT P_star TT' + QQ
      ##                                  - K0 M_star' - M_star K0' + K0 F_star K0') ----

      ## d/d(TT P_star TT') = sym: contrib 2 bar_Pstar_t * TT * P_star_prev
      G_TT <- G_TT + 2 * bar_Pstar_t %*% TT %*% P_star_prev
      bar_Pstar_from_TT <- t(TT) %*% bar_Pstar_t %*% TT  # [n x n]

      ## d/dQQ = bar_Pstar_t  (QQ enters P_star directly)
      ## QQ = RR Sig RR' => G_Sig += RR' bar_Pstar_t RR; G_RR += 2 bar_Pstar_t RR Sig
      G_Sig <- G_Sig + t(RR) %*% bar_Pstar_t %*% RR
      G_RR  <- G_RR  + 2 * bar_Pstar_t %*% RR %*% Sigma_e

      ## d/d(-K0 M_star' - M_star K0') = -sym(K0 M_star'):
      ## tr(bar_Pstar_t d(-sym(K0 M_star')))
      ##   = -2 tr(bar_Pstar_t dK0 M_star') [sym adjoint]
      ##   = -2 tr(M_star' bar_Pstar_t dK0)  [cyclic]
      ##   => bar_K0 from this: -2 * bar_Pstar_t %*% M_star  [n x q]
      ## d/dM_star: -2 tr(K0' bar_Pstar_t dM_star)
      ##   => bar_M_star: -2 * K0' %*% bar_Pstar_t  [q x n] -> as n x q: -(2*bar_Pstar_t %*% K0)'
      ##   Actually bar_M_star has shape [n x q]; the formula:
      ##   tr(grad_{M_star}' dM_star) = -2 tr(K0' bar_Pstar_t dM_star)
      ##   => grad_{M_star} = -2 bar_Pstar_t %*% K0  [n x q] after transposing: NO.
      ##   tr(K0' bar_Pstar_t dM_star) where K0 (n x q), bar_P (n x n), dM_star (n x q):
      ##   = tr((K0' bar_Pstar_t) dM_star) where (K0' bar_P) is (q x n), dM* is (n x q)
      ##   = sum_{ij} [K0' bar_P]_{ij} dM*_{ji}  [trace of q x n times n x q]
      ##   = Frobenius inner product of (K0' bar_P)^T and dM*
      ##   => bar_M_star = -(K0' bar_Pstar_t)' * 2 = -2 bar_Pstar_t %*% K0  [n x q] ... wait:
      ##   (K0' bar_P)^T = bar_P' K0 = bar_P K0  [n x q] (since bar_P symmetric)
      ##   So bar_M_star = -2 bar_Pstar_t %*% K0  [n x q]  -- YES, this is (n x q).
      bar_K0    <- -2 * bar_Pstar_t %*% M_star               # [n x q]
      bar_Mstar <- -2 * bar_Pstar_t %*% K0                   # [n x q]

      ## d/d(K0 F_star K0') = sym:
      ## tr(bar_P d(K0 F_star K0')) = 2 tr(bar_P dK0 F_star K0') + tr(bar_P K0 dF_star K0')
      ##   = 2 tr(F_star K0' bar_P dK0) + tr(K0' bar_P K0 dF_star)
      ## => bar_K0 += 2 bar_Pstar_t %*% K0 %*% F_star  [n x q]
      ## => bar_F_star = tK0 %*% bar_Pstar_t %*% K0    [q x q]
      bar_K0     <- bar_K0 + 2 * bar_Pstar_t %*% K0 %*% F_star
      bar_Fstar  <- tK0 %*% bar_Pstar_t %*% K0              # [q x q]

      ## ----- 2. Adjoint of P_inf' = sym(TT P_inf TT' - K0 F_inf K0') ---------
      ## d/d(TT P_inf TT'): contrib 2 bar_Pinf_t * TT * P_inf_prev
      G_TT <- G_TT + 2 * bar_Pinf_t %*% TT %*% P_inf_prev
      bar_Pinf_from_TT <- t(TT) %*% bar_Pinf_t %*% TT       # [n x n]

      ## d/d(-K0 F_inf K0') = -sym(K0 F_inf K0'):
      ## => bar_K0 += -2 * bar_Pinf_t %*% K0 %*% F_inf  [n x q]
      ## => bar_F_inf = -tK0 %*% bar_Pinf_t %*% K0  [q x q]  (initialise here)
      bar_K0    <- bar_K0 - 2 * bar_Pinf_t %*% K0 %*% F_inf
      bar_F_inf <- -tK0 %*% bar_Pinf_t %*% K0               # [q x q]

      ## ----- 3. Adjoint of ll_t = -0.5 * log|F_inf| --------------------------
      ## d/dF_inf (-0.5 log|F_inf|) = -0.5 * F_inf^{-1}
      bar_F_inf <- bar_F_inf - 0.5 * Fi                      # [q x q]

      ## ----- 4. Adjoint of s' = TT s_prev + K0 v ----------------------------
      G_TT       <- G_TT + outer(bar_s, s_prev)
      bar_K0     <- bar_K0 + outer(bar_s, v)                 # [n x q]
      bar_v      <- as.numeric(tK0 %*% bar_s)
      bar_s_prev <- as.numeric(t(TT) %*% bar_s)

      ## ----- 5. Adjoint of v = y - d - ZZ s_prev ----------------------------
      g_d        <- g_d - bar_v
      G_ZZ       <- G_ZZ - outer(bar_v, s_prev)
      bar_s_prev <- bar_s_prev - as.numeric(tZZ %*% bar_v)

      ## ----- 6. Adjoint of K0 = M_inf F_inf^{-1} ----------------------------
      ## bar_K0 fully accumulated now. Propagate:
      ## bar_M_inf = bar_K0 %*% Fi  [n x q]
      ## bar_F_inf += -K0' bar_K0 Fi  [q x q]  (from F_inf^{-1} dependence)
      bar_Minf   <- bar_K0 %*% Fi                            # [n x q]
      bar_F_inf  <- bar_F_inf - tK0 %*% bar_K0 %*% Fi       # [q x q]
      bar_F_inf  <- .sym(bar_F_inf)

      ## ----- 7. Adjoint of M_star = TT P_star ZZ' + SS ----------------------
      ## bar_Mstar [n x q] from step 1.
      ## From TT P_star ZZ' (same pattern as Mnum in dense adjoint):
      G_TT  <- G_TT  + bar_Mstar %*% ZZ %*% P_star_prev
      G_ZZ  <- G_ZZ  + t(bar_Mstar) %*% TT %*% P_star_prev
      bar_Pstar_from_Mstar <- t(TT) %*% bar_Mstar %*% ZZ    # [n x n]
      ## From SS = RR Sig DD':
      G_RR  <- G_RR  + bar_Mstar %*% DD %*% Sigma_e
      G_DD  <- G_DD  + t(bar_Mstar) %*% RR %*% Sigma_e
      G_Sig <- G_Sig + t(RR) %*% bar_Mstar %*% DD

      ## ----- 8. Adjoint of M_inf = TT P_inf ZZ' -----------------------------
      G_TT  <- G_TT  + bar_Minf %*% ZZ %*% P_inf_prev
      G_ZZ  <- G_ZZ  + t(bar_Minf) %*% TT %*% P_inf_prev
      bar_Pinf_from_Minf <- t(TT) %*% bar_Minf %*% ZZ       # [n x n]

      ## ----- 9. Adjoint of F_star = sym(ZZ P_star ZZ' + HH + me_diag) -------
      bar_Fstar <- .sym(bar_Fstar)
      G_ZZ  <- G_ZZ  + 2 * bar_Fstar %*% ZZ %*% P_star_prev
      bar_Pstar_from_Fstar <- t(ZZ) %*% bar_Fstar %*% ZZ    # [n x n]
      G_DD  <- G_DD  + 2 * bar_Fstar %*% DD %*% Sigma_e
      G_Sig <- G_Sig + t(DD) %*% bar_Fstar %*% DD

      ## ----- 10. Adjoint of F_inf = sym(ZZ P_inf ZZ') -----------------------
      ## bar_F_inf already symmetrized.
      G_ZZ  <- G_ZZ  + 2 * bar_F_inf %*% ZZ %*% P_inf_prev
      bar_Pinf_from_Finf <- t(ZZ) %*% bar_F_inf %*% ZZ      # [n x n]

      ## ----- Collect bar_P_{t-1} ---------------------------------------------
      bar_P_inf  <- .sym(bar_Pinf_from_TT + bar_Pinf_from_Minf + bar_Pinf_from_Finf)
      bar_P_star <- .sym(bar_Pstar_from_TT + bar_Pstar_from_Mstar + bar_Pstar_from_Fstar)
      bar_s      <- bar_s_prev
    }
  }  # end backward loop

  ## At t=0: bar_P_star = bar_P_star_0, bar_P_inf = bar_P_inf_0, bar_s_0 unused.
  bar_P_inf_0  <- bar_P_inf
  bar_P_star_0 <- bar_P_star

  ## ============================================================
  ## Adjoint through P_star_0 (Lyapunov in Schur basis)
  ## and P_inf_0 (spectral projector)
  ## ============================================================

  ## Re-derive the Schur decomposition used in .kf_diffuse_P0 (SAME swap logic).
  schur <- Matrix::Schur(Matrix::Matrix(TT, sparse = FALSE))
  Ts    <- as.matrix(schur@T)
  U     <- as.matrix(schur@Q)

  ur_target <- 0L
  n <- n_state
  for (i in seq_len(n)) {
    if (abs(abs(diag(Ts)[i]) - 1) < ur_tol) {
      j <- i
      while (j > ur_target + 1L) {
        sw <- .swap_schur_11(Ts, U, j - 1L)
        Ts <- sw$T; U <- sw$Q
        j  <- j - 1L
      }
      ur_target <- ur_target + 1L
    }
  }
  ## Now U, Ts match exactly what .kf_diffuse_P0 used.

  stage1_ok <- TRUE   # forward pass + diffuse/stationary recursion adjoint done
  stage2_ok <- TRUE   # full analytic init adjoint below (FALSE only on failure)

  ## ---- Full analytic adjoint through the diffuse initialization ------------
  ## P_inf_0 = U J U' (J = blockdiag(I_u, 0)) and P_star_0 = U_s Pa_ss U_s' are
  ## BOTH invariant to within-block Schur rotations, so the only hidden variable
  ## is the inter-block coupling Omega_su, which solves the Sylvester equation
  ##   T_ss Omega_su - Omega_su T_uu = -(U' dTT U)_su .
  ## We assemble the exact init-gradient by running the analytic FORWARD
  ## directional derivative of (P_inf_0, P_star_0) on each basis direction and
  ## contracting with (bar_P_inf_0, bar_P_star_0). The forward derivative is
  ## exact (Sylvester + Lyapunov derivatives); assembling the Jacobian
  ## column-by-column is one-time (init only), O(n^2) small solves. This is the
  ## full Schur-rotation init adjoint (general non-triangular TT included).
  if (nunit > 0L) {
    idx_u <- seq_len(nunit)
    has_stable <- (nunit < n)
    idx_s <- if (has_stable) (nunit + 1L):n else integer(0)
    U_u  <- U[, idx_u, drop = FALSE]
    T_uu <- Ts[idx_u, idx_u, drop = FALSE]
    if (has_stable) {
      U_s   <- U[, idx_s, drop = FALSE]
      T_ss  <- Ts[idx_s, idx_s, drop = FALSE]
      QQ_ss <- t(U_s) %*% QQ %*% U_s
      Pa_ss <- tryCatch(.sym(.solve_lyapunov(T_ss, QQ_ss)),
                        error = function(e) NULL)
      if (is.null(Pa_ss) || !all(is.finite(Pa_ss))) stage2_ok <- FALSE
    }

    ## Sylvester solve  T_ss X - X T_uu = C  (small, one-time) via kron.
    solve_syl <- function(C) {
      ns <- nrow(T_ss); nu <- ncol(T_uu)
      A  <- kronecker(diag(nu), T_ss) - kronecker(t(T_uu), diag(ns))
      matrix(solve(A, as.vector(C)), ns, nu)
    }

    ## Forward directional derivative of (P_inf_0, P_star_0) given (dTT, dQQ).
    fwd <- function(dTT, dQQ) {
      if (has_stable && !is.null(dTT)) {
        Gsu   <- (t(U) %*% dTT %*% U)[idx_s, idx_u, drop = FALSE]
        Om_su <- solve_syl(-Gsu)                  # T_ss Om - Om T_uu = -G_su
      } else {
        Om_su <- matrix(0, length(idx_s), nunit)
      }
      dU_u   <- if (has_stable) U_s %*% Om_su else matrix(0, n, nunit)
      dP_inf <- .sym(dU_u %*% t(U_u) + U_u %*% t(dU_u))
      dP_star <- matrix(0, n, n)
      if (has_stable) {
        dU_s  <- -U_u %*% t(Om_su)
        dTm   <- if (is.null(dTT)) matrix(0, n, n) else dTT
        dQm   <- if (is.null(dQQ)) matrix(0, n, n) else dQQ
        dT_ss <- t(dU_s) %*% TT %*% U_s + t(U_s) %*% dTm %*% U_s +
                 t(U_s) %*% TT %*% dU_s
        dQ_ss <- t(dU_s) %*% QQ %*% U_s + t(U_s) %*% dQm %*% U_s +
                 t(U_s) %*% QQ %*% dU_s
        rhs   <- dT_ss %*% Pa_ss %*% t(T_ss) + T_ss %*% Pa_ss %*% t(dT_ss) + dQ_ss
        dPa   <- .sym(.solve_lyapunov(T_ss, .sym(rhs)))
        dP_star <- .sym(dU_s %*% Pa_ss %*% t(U_s) + U_s %*% dPa %*% t(U_s) +
                        U_s %*% Pa_ss %*% t(dU_s))
      }
      list(dP_inf = dP_inf, dP_star = dP_star)
    }

    if (isTRUE(stage2_ok)) {
      ## bar_TT_init: contract bars with the forward derivative per TT basis dir.
      bar_TT_init <- matrix(0, n, n)
      E <- matrix(0, n, n)
      for (a in seq_len(n)) for (b in seq_len(n)) {
        E[a, b] <- 1
        fd <- tryCatch(fwd(E, NULL), error = function(e) NULL)
        if (is.null(fd)) { stage2_ok <- FALSE; break }
        bar_TT_init[a, b] <- sum(bar_P_inf_0 * fd$dP_inf) +
                             sum(bar_P_star_0 * fd$dP_star)
        E[a, b] <- 0
      }
      if (isTRUE(stage2_ok)) G_TT <- G_TT + bar_TT_init

      ## bar_QQ_init: only P_star_0 depends on QQ.
      if (isTRUE(stage2_ok) && has_stable) {
        bar_QQ_init <- matrix(0, n, n)
        E <- matrix(0, n, n)
        for (a in seq_len(n)) for (b in seq_len(n)) {
          E[a, b] <- 1
          fd <- fwd(NULL, E)
          bar_QQ_init[a, b] <- sum(bar_P_star_0 * fd$dP_star)
          E[a, b] <- 0
        }
        bar_QQ_init <- .sym(bar_QQ_init)
        ## QQ = RR Sigma_e RR'  ->  G_Sig, G_RR (same map as the dense adjoint).
        G_Sig <- G_Sig + t(RR) %*% bar_QQ_init %*% RR
        G_RR  <- G_RR  + 2 * bar_QQ_init %*% RR %*% Sigma_e
      }
    }
  }

  if (!isTRUE(stage2_ok)) {
    has_dTT <- any(vapply(d_ss_list, function(p)
      !is.null(p) && !is.null(p[["dTT"]]), logical(1)))
    if (has_dTT)
      warning(".kf_loglik_adjoint_diffuse: diffuse init adjoint (Sylvester/",
              "Lyapunov) failed; TT gradient omits the init sensitivity.")
  }
  ## ============================================================
  ## Final contraction: grad[j] = <G_TT, dTT_j> + ... Frobenius
  ## ============================================================
  grad <- numeric(n_par)
  for (j in seq_len(n_par)) {
    dpar <- d_ss_list[[j]]
    if (is.null(dpar)) { grad[j] <- 0; next }

    gj <- 0
    dTT  <- dpar[["dTT"]];      if (!is.null(dTT))  gj <- gj + sum(G_TT  * dTT)
    dRR  <- dpar[["dRR"]];      if (!is.null(dRR))  gj <- gj + sum(G_RR  * dRR)
    dZZ  <- dpar[["dZZ"]];      if (!is.null(dZZ))  gj <- gj + sum(G_ZZ  * dZZ)
    dDD  <- dpar[["dDD"]];      if (!is.null(dDD))  gj <- gj + sum(G_DD  * dDD)
    dd_j <- dpar[["dd"]];       if (!is.null(dd_j)) gj <- gj + sum(g_d   * dd_j)
    dSig <- dpar[["dSigma_e"]]; if (!is.null(dSig)) gj <- gj + sum(G_Sig * dSig)

    grad[j] <- gj
  }

  list(
    loglik       = loglik,
    grad         = grad,
    stage1_ok    = stage1_ok,
    stage2_ok    = stage2_ok,
    bar_P_inf_0  = bar_P_inf_0,
    bar_P_star_0 = bar_P_star_0,
    nunit        = nunit,
    d_diffuse    = d_diffuse,
    min_eig_margin       = min_eig_margin,
    near_regime_boundary = near_regime_boundary
  )
}
