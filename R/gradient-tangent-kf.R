## R/gradient-tangent-kf.R
## --------------------------------------------------------------------------
## Analytic TANGENT (forward-sensitivity) Kalman filter.
##
## This is one of two building blocks for full implicit-differentiation
## likelihood gradients (Childers, Fernandez-Villaverde, Perla, Rackauckas &
## Wu 2022, NBER w30573): given the derivatives of the state-space matrices
## (TT, RR, ZZ, DD, d, Sigma_e) with respect to a structural parameter theta_j
## -- typically obtained from a separate solution-derivative ("implicit
## function theorem on the policy function") layer -- this file propagates
## the corresponding derivatives of the filtered state, covariance, and the
## per-period log-likelihood contributions through the SAME recursion used by
## .kf_step() in R/kalman-filter.R, accumulating the exact gradient of the
## Gaussian log-likelihood.
##
## .kf_loglik_score_sigma() (R/analytic-gradient.R) is the special case where
## only Sigma_e moves (TT, RR, ZZ, DD, d held fixed); .kf_loglik_tangent()
## below generalizes that recursion to arbitrary derivatives of every
## state-space matrix, including dd (the observation-constant derivative).
##
## Conventions mirror R/kalman-filter.R's .kf_step EXACTLY:
##   QQ = RR Sigma_e RR',  HH = DD Sigma_e DD',  SS = RR Sigma_e DD'
##   v_t = y_t - d - ZZ s            (s = predicted state at time t)
##   F   = ZZ P ZZ' + HH + me_diag   (symmetrised)
##   K   = (TT P ZZ' + SS) F^{-1}
##   s'  = TT s + K v
##   P'  = (TT - K ZZ) P (TT - K ZZ)' + (RR - K DD) Sigma_e (RR - K DD)'
##         (symmetrised; "Joseph form")
##   P0  = solve_lyapunov(TT, QQ),  s0 = 0
##
## Plain R, no C++: correctness first. An Rcpp port of the per-step tangent
## recursion can follow the kf_score_sigma_cpp pattern (R/analytic-gradient.R)
## once this reference implementation is validated against finite differences.
## --------------------------------------------------------------------------


## Extract a derivative piece from a per-parameter list, defaulting to "zero"
## (a matrix/vector of the same shape as `zero`) when `nm` is absent or NULL.
## Shared by the C++ dispatch shim (.kf_loglik_tangent) and the R reference
## loop's local `.dpiece`.
## @noRd
.dpiece2 <- function(dpar, nm, zero) {
  val <- dpar[[nm]]
  if (is.null(val)) zero else val
}

#' TRUE when the compiled tangent Kalman-filter recursion is available.
#' @noRd
.HAS_RCPP_KF_TANGENT <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("kf_tangent_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

#' Analytic tangent (forward-sensitivity) Kalman-filter log-likelihood gradient
#'
#' Runs the exact per-step Kalman filter (the same recursion as
#' \code{.kf_step} in \code{R/kalman-filter.R}) together with one
#' forward-sensitivity ("tangent") recursion per parameter, simultaneously
#' accumulating the Gaussian log-likelihood and its exact gradient with
#' respect to each parameter \code{theta_j} given the supplied derivatives of
#' the state-space matrices.
#'
#' This generalizes \code{.kf_loglik_score_sigma} (which only handles
#' parameters that move \code{Sigma_e}) to arbitrary derivatives of
#' \code{TT, RR, ZZ, DD, d, Sigma_e}; any of the six matrix derivatives may be
#' \code{NULL} (treated as the zero matrix/vector of the appropriate shape).
#'
#' @param Y observation matrix, \code{n_obs x n_T} (no missing values --
#'   \code{NA}s are not supported and trigger an error).
#' @param ss list with the base state-space matrices: \code{TT, RR, ZZ, DD}
#'   (matrices, shaped as in \code{kalman_filter}), \code{d} (named numeric
#'   vector of observation steady states, length \code{n_obs}), and
#'   \code{Sigma_e} (shock covariance, \code{n_exo x n_exo}).
#' @param d_ss_list list of length \code{n_par}; element \code{j} is itself a
#'   list with (any subset of) \code{dTT, dRR, dZZ, dDD, dd, dSigma_e} giving
#'   \code{d <matrix> / d theta_j}. Missing/\code{NULL} entries are treated as
#'   zero (no dependence of that matrix on \code{theta_j}).
#' @param me_variance scalar measurement-error variance added to the diagonal
#'   of \code{F} (as in \code{kalman_filter}); does not enter any gradient
#'   (it is parameter-independent here, mirroring \code{me_diag}).
#'
#' @return list with:
#'   \item{loglik}{scalar Gaussian log-likelihood (matches the "dare"/standard
#'     per-step Kalman filter to numerical precision).}
#'   \item{grad}{numeric vector of length \code{n_par}, the exact gradient of
#'     \code{loglik} with respect to each \code{theta_j}.}
#'
#' If the base filter's \code{F} matrix is not positive definite (Cholesky
#' failure) at any step, or a per-period log-likelihood contribution falls
#' below \code{.KF_LL_MIN}, this returns \code{list(loglik = -Inf, grad =
#' rep(NA_real_, n_par))} -- gracefully, with no error -- mirroring the
#' \code{.kf_step} guard in \code{R/kalman-filter.R}.
#'
#' @noRd
.kf_loglik_tangent <- function(Y, ss, d_ss_list, me_variance = 0,
                               me_extra = NULL, shock_scale = NULL) {

  TT <- ss$TT; RR <- ss$RR; ZZ <- ss$ZZ; DD <- ss$DD
  d  <- as.numeric(ss$d); Sigma_e <- ss$Sigma_e

  n_state <- nrow(TT)
  n_obs   <- nrow(ZZ)
  n_exo   <- ncol(RR)
  n_par   <- length(d_ss_list)

  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  if (anyNA(Y)) stop(".kf_loglik_tangent: Y must not contain missing values.")
  n_T <- ncol(Y)

  has_me_extra    <- !is.null(me_extra)
  ## Any TRUE observation noise (base me_variance or per-period me_extra):
  ## both feed F AND the Joseph covariance term (F3-D).
  has_me_true     <- has_me_extra || me_variance != 0
  has_shock_scale <- !is.null(shock_scale)

  ## -- Fast path: compiled per-step recursion (kf_tangent_cpp) -------------
  ## Densify d_ss_list into per-parameter arma::cube/matrix arguments (zero
  ## slices/columns for NULL entries) and dispatch to the C++ backend, which
  ## mirrors this R function line-for-line (see src/kf_tangent.cpp). The R
  ## loop below remains the reference fallback (see test-gradient-tangent-kf
  ## for the parity checks); bit-identical to 1e-12 / 1e-10.
  ## tv inputs (me_extra / shock_scale) are now handled by C++ -- sentinels
  ## (all-ones / all-zeros) are passed when the tv arg is absent so the
  ## constant path is arithmetically equivalent (IEEE 754 exact).
  ## Force-R escape hatch: set option dynhr.use_rcpp = FALSE.
  if (.HAS_RCPP_KF_TANGENT()) {
    zTT0 <- matrix(0, n_state, n_state)
    zRR0 <- matrix(0, n_state, n_exo)
    zZZ0 <- matrix(0, n_obs, n_state)
    zDD0 <- matrix(0, n_obs, n_exo)
    zSig0 <- matrix(0, n_exo, n_exo)

    dTT_cube <- array(0, c(n_state, n_state, n_par))
    dRR_cube <- array(0, c(n_state, n_exo, n_par))
    dZZ_cube <- array(0, c(n_obs, n_state, n_par))
    dDD_cube <- array(0, c(n_obs, n_exo, n_par))
    dd_mat   <- matrix(0, n_obs, n_par)
    dSigma_cube <- array(0, c(n_exo, n_exo, n_par))

    for (j in seq_len(n_par)) {
      dpar <- d_ss_list[[j]]
      if (is.null(dpar)) dpar <- list()
      dTT_cube[, , j]    <- .dpiece2(dpar, "dTT", zTT0)
      dRR_cube[, , j]    <- .dpiece2(dpar, "dRR", zRR0)
      dZZ_cube[, , j]    <- .dpiece2(dpar, "dZZ", zZZ0)
      dDD_cube[, , j]    <- .dpiece2(dpar, "dDD", zDD0)
      dd_mat[, j]        <- .dpiece2(dpar, "dd", numeric(n_obs))
      dSigma_cube[, , j] <- .dpiece2(dpar, "dSigma_e", zSig0)
    }

    ## Sentinels for tv args: all-ones shock_scale => Se_t = Sigma_e (exact);
    ## all-zeros me_extra => me_diag_t = me_diag (exact).
    shock_scale_cpp <- if (has_shock_scale) shock_scale else matrix(1.0, n_exo, n_T)
    me_extra_cpp    <- if (has_me_extra)    me_extra    else matrix(0.0, n_obs, n_T)

    out <- kf_tangent_cpp(Y, TT, RR, ZZ, DD, d, Sigma_e,
                          dTT_cube, dRR_cube, dZZ_cube, dDD_cube, dd_mat,
                          dSigma_cube, me_variance, .KF_LL_MIN,
                          shock_scale_cpp, me_extra_cpp)

    if (!isTRUE(out$ok)) {
      return(list(loglik = -Inf, grad = rep(NA_real_, n_par)))
    }
    return(list(loglik = out$loglik, grad = as.numeric(out$grad)))
  }

  tZZ <- t(ZZ)
  ## Baseline noise matrices (use Sigma_e unchanged -- P0/dP0 stay on baseline)
  QQ  <- tcrossprod(RR %*% Sigma_e, RR)
  HH  <- tcrossprod(DD %*% Sigma_e, DD)
  SS  <- RR %*% Sigma_e %*% t(DD)
  me_diag  <- me_variance * diag(n_obs)
  ll_const <- -0.5 * n_obs * log(2 * pi)

  fail <- list(loglik = -Inf, grad = rep(NA_real_, n_par))

  ## -- Helper: extract a derivative piece, defaulting to "zero" -----------
  ## Returns a matrix/vector of the right shape filled with zeros if `nm` is
  ## absent or NULL in `dpar`.
  .dpiece <- function(dpar, nm, zero) {
    val <- dpar[[nm]]
    if (is.null(val)) zero else val
  }

  zTT <- matrix(0, n_state, n_state)
  zRR <- matrix(0, n_state, n_exo)
  zZZ <- matrix(0, n_obs, n_state)
  zDD <- matrix(0, n_obs, n_exo)
  zd  <- numeric(n_obs)
  zSig <- matrix(0, n_exo, n_exo)

  ## -- Base initial conditions ---------------------------------------------
  P0 <- .solve_lyapunov(TT, QQ)
  if (!all(is.finite(P0))) return(fail)
  s0 <- numeric(n_state)

  ## -- Per-parameter pieces: dQQ, dHH, dSS, dP0, ds0 -------------------------
  dTT_l <- vector("list", n_par); dRR_l <- vector("list", n_par)
  dZZ_l <- vector("list", n_par); dDD_l <- vector("list", n_par)
  dd_l  <- vector("list", n_par); dSig_l <- vector("list", n_par)
  dQQ_l <- vector("list", n_par); dHH_l <- vector("list", n_par)
  dSS_l <- vector("list", n_par)
  dP_l  <- vector("list", n_par); ds_l  <- vector("list", n_par)

  for (j in seq_len(n_par)) {
    dpar <- d_ss_list[[j]]
    if (is.null(dpar)) dpar <- list()

    dTT  <- .dpiece(dpar, "dTT", zTT)
    dRR  <- .dpiece(dpar, "dRR", zRR)
    dZZ  <- .dpiece(dpar, "dZZ", zZZ)
    dDD  <- .dpiece(dpar, "dDD", zDD)
    dd_j <- .dpiece(dpar, "dd",  zd)
    dSig <- .dpiece(dpar, "dSigma_e", zSig)

    dTT_l[[j]] <- dTT;  dRR_l[[j]] <- dRR
    dZZ_l[[j]] <- dZZ;  dDD_l[[j]] <- dDD
    dd_l[[j]]  <- as.numeric(dd_j); dSig_l[[j]] <- dSig

    ## dQQ = dRR Sig RR' + RR dSig RR' + RR Sig dRR'   (symmetrise)
    dQQ <- tcrossprod(dRR %*% Sigma_e, RR) +
           tcrossprod(RR %*% dSig, RR) +
           tcrossprod(RR %*% Sigma_e, dRR)
    dQQ <- .sym(dQQ)

    ## dHH = dDD Sig DD' + DD dSig DD' + DD Sig dDD'   (symmetrise)
    dHH <- tcrossprod(dDD %*% Sigma_e, DD) +
           tcrossprod(DD %*% dSig, DD) +
           tcrossprod(DD %*% Sigma_e, dDD)
    dHH <- .sym(dHH)

    ## dSS = dRR Sig DD' + RR dSig DD' + RR Sig dDD'   (NOT symmetric in general)
    dSS <- dRR %*% Sigma_e %*% t(DD) +
           RR %*% dSig %*% t(DD) +
           RR %*% Sigma_e %*% t(dDD)

    dQQ_l[[j]] <- dQQ; dHH_l[[j]] <- dHH; dSS_l[[j]] <- dSS

    ## dP0 = solve_lyapunov(TT, RHS), RHS = dTT P0 TT' + TT P0 dTT' + dQQ
    rhs <- dTT %*% P0 %*% t(TT) + TT %*% P0 %*% t(dTT) + dQQ
    rhs <- .sym(rhs)
    dP0 <- .solve_lyapunov(TT, rhs)
    if (!all(is.finite(dP0))) return(fail)

    dP_l[[j]] <- dP0
    ds_l[[j]] <- numeric(n_state)
  }

  s <- s0
  P <- P0
  loglik <- 0
  grad   <- numeric(n_par)

  for (t in seq_len(n_T)) {
    ## -- Per-period tv substitutions (shock_scale / me_extra) --------------
    if (has_shock_scale) {
      sc_t  <- shock_scale[, t]
      Se_t  <- Sigma_e * outer(sc_t, sc_t)    # D_t Sigma_e D_t
      HH_t  <- tcrossprod(DD %*% Se_t, DD)
      SS_t  <- RR %*% Se_t %*% t(DD)
    } else {
      Se_t <- Sigma_e; HH_t <- HH; SS_t <- SS
    }
    me_x_t    <- if (has_me_extra) me_variance + me_extra[, t]
                 else rep(me_variance, n_obs)
    me_diag_t <- diag(me_x_t, n_obs)

    ## -- Base step (mirrors .kf_step exactly) ------------------------------
    PZ <- P %*% tZZ
    Ft <- ZZ %*% PZ + HH_t + me_diag_t
    Ft <- .sym(Ft)

    Fc <- tryCatch(chol(Ft), error = function(e) NULL)
    if (is.null(Fc)) return(fail)
    Fi  <- chol2inv(Fc)
    ldf <- 2 * sum(log(diag(Fc)))

    v  <- Y[, t] - d - as.numeric(ZZ %*% s)
    Fiv <- as.numeric(Fi %*% v)

    ll_t <- ll_const - 0.5 * (ldf + sum(v * Fiv))
    if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) return(fail)
    loglik <- loglik + ll_t

    K    <- (TT %*% PZ + SS_t) %*% Fi
    A    <- TT - K %*% ZZ          # TT - K ZZ
    B    <- RR - K %*% DD          # RR - K DD

    s_n <- as.numeric(TT %*% s) + as.numeric(K %*% v)
    P_n <- tcrossprod(A %*% P, A) + tcrossprod(B %*% Se_t, B)
    ## Joseph true-noise term for the FULL ME diagonal (mirrors
    ## kalman_filter's standard path): P' += K diag(me_x_t) K'. Kme is
    ## hoisted for the per-parameter tangent recursion below.
    Kme <- NULL
    if (has_me_true) {
      Kme    <- K %*% me_diag_t                    # K diag(me_x_t)
      P_n    <- P_n + tcrossprod(Kme, K)           # + K me_x K'
    }
    P_n <- .sym(P_n)

    ## -- Hoisted parameter-independent pieces for the tangent recursion ----
    TtP   <- TT %*% P            # TT P  (used for TT dP and dTT P)
    AP    <- A %*% P
    BSig  <- B %*% Se_t          # B Se_t (tv)

    for (j in seq_len(n_par)) {
      dTT  <- dTT_l[[j]]; dRR <- dRR_l[[j]]; dZZ <- dZZ_l[[j]]
      dDD  <- dDD_l[[j]]; dd_j <- dd_l[[j]]; dSig <- dSig_l[[j]]
      dP   <- dP_l[[j]];  ds  <- ds_l[[j]]

      ## Per-period derivatives of HH/SS: when shock_scale is active, substitute
      ## Se_t for Sigma_e; otherwise reuse the pre-loop constants (Landmine L1).
      if (has_shock_scale) {
        dSig_t <- dSig * outer(sc_t, sc_t)  # D_t dSig_j D_t
        dHH_t <- .sym(tcrossprod(dDD %*% Se_t, DD) +
                       tcrossprod(DD %*% dSig_t, DD) +
                       tcrossprod(DD %*% Se_t, dDD))
        dSS_t <- dRR %*% Se_t %*% t(DD) +
                 RR %*% dSig_t %*% t(DD) +
                 RR %*% Se_t %*% t(dDD)
      } else {
        dHH_t <- dHH_l[[j]]; dSS_t <- dSS_l[[j]]; dSig_t <- dSig
      }

      tdZZ <- t(dZZ)

      ## dv = -dd - dZZ s - ZZ ds
      dv <- -dd_j - as.numeric(dZZ %*% s) - as.numeric(ZZ %*% ds)

      ## dF = dZZ P ZZ' + ZZ dP ZZ' + ZZ P dZZ' + dHH_t   (symmetrise)
      ## Note: me_extra is theta-independent; its derivative is zero, so dF does
      ## NOT include any d(me_extra[,t])/d(theta_j) term.
      dF <- tcrossprod(dZZ %*% P, ZZ) +
            ZZ %*% tcrossprod(dP, ZZ) +
            ZZ %*% tcrossprod(P, dZZ) +
            dHH_t
      dF <- .sym(dF)

      ## dll_t = -0.5 * ( tr(Fi dF) - (Fi v)' dF (Fi v) ) - (Fi v)' dv
      Fi_dF_Fiv <- as.numeric(dF %*% Fiv)
      dll_t <- -0.5 * (sum(Fi * dF) - sum(Fiv * Fi_dF_Fiv)) - sum(Fiv * dv)
      grad[j] <- grad[j] + dll_t

      ## dK = (dTT P ZZ' + TT dP ZZ' + TT P dZZ' + dSS_t) Fi - K dF Fi
      dPZ <- dP %*% tZZ
      dK_num <- tcrossprod(dTT %*% P, ZZ) +
                TT %*% dPZ +
                TtP %*% tdZZ +
                dSS_t
      dK <- dK_num %*% Fi - K %*% (dF %*% Fi)

      ## ds' = dTT s + TT ds + dK v + K dv
      ds_n <- as.numeric(dTT %*% s) + as.numeric(TT %*% ds) +
              as.numeric(dK %*% v) + as.numeric(K %*% dv)

      ## dA = dTT - dK ZZ - K dZZ;  dB = dRR - dK DD - K dDD
      dA <- dTT - dK %*% ZZ - K %*% dZZ
      dB <- dRR - dK %*% DD - K %*% dDD

      ## dP' = dA P A' + A dP A' + A P dA' + dB Se_t B' + B dSig_t B' + B Se_t dB'
      ## (Se_t = Sigma_e when no shock_scale; BSig = B * Se_t already)
      dP_n <- tcrossprod(dA %*% P, A) +
              A %*% tcrossprod(dP, A) +
              AP %*% t(dA) +
              tcrossprod(dB %*% Se_t, B) +
              B %*% tcrossprod(dSig_t, B) +
              BSig %*% t(dB)
      ## Tangent of the ME Joseph term P' += K me_x K' (the ME diagonal is
      ## data, not differentiated): dP' += dK me_x K' + K me_x dK'.
      if (has_me_true)
        dP_n <- dP_n + tcrossprod(dK %*% me_diag_t, K) + Kme %*% t(dK)
      dP_n <- .sym(dP_n)

      ds_l[[j]] <- ds_n
      dP_l[[j]] <- dP_n
    }

    s <- s_n
    P <- P_n
  }

  list(loglik = loglik, grad = grad)
}
