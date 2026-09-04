## R/gradient-adjoint-ss.R
## --------------------------------------------------------------------------
## Steady-state-aware variant of the dense adjoint gradient.
##
## Exploits convergence of the Riccati recursion P_t -> P_inf (DARE fixed
## point) to cut storage from O(T*n^2) to O(t_conv*n^2 + T*(n+q)).
##
## In a stationary model, after t_conv steps the covariance P_t converges and
## {F_t, K_t, A_t, B_t} become CONSTANT (= P_inf, Fi_inf, K_inf, A_inf, B_inf).
## Only the data-dependent quantities {s_{t-1}, v_t} keep varying in the tail.
##
## Forward pass:
##   - For t <= t_conv: store per-step {s, P, v, Fi, K, A, B} exactly as the
##     dense adjoint does.
##   - Once ||P_new - P|| < ss_tol: lock {P_inf, Fi_inf, K_inf, A_inf, B_inf},
##     and for t > t_conv store ONLY {s_{t-1}, v_t} (size n+q per step).
##
## Backward sweep:
##   - IDENTICAL math to the dense adjoint for every step.
##   - For t > t_conv, P_t is replaced by P_inf and {Fi,K,A,B} by the locked
##     constants -- the same approximation the steady-state KF already makes
##     for the loglik, controlled by ss_tol.
##
## Storage: O(t_conv * n^2 + T * (n + q))  vs  O(T * n^2) for dense.
## Convergence detection mirrors kalman_standard_loop_cpp: t > 1 AND
##   max(abs(P_new - P)) < ss_tol.
##
## State-space convention: identical to gradient-adjoint-kf.R.
## --------------------------------------------------------------------------

#' Steady-state adjoint Kalman-filter log-likelihood gradient
#'
#' Drop-in replacement for \code{.kf_loglik_adjoint} that exploits
#' covariance convergence to cut n^2-matrix storage from O(T) to O(t_conv).
#'
#' @param Y          n_obs x n_T observation matrix (no NAs).
#' @param ss         list: TT, RR, ZZ, DD, d, Sigma_e.
#' @param d_ss_list  list of length n_par; element j has (any of) dTT, dRR,
#'                   dZZ, dDD, dd, dSigma_e.
#' @param me_variance scalar measurement-error variance (parameter-independent).
#' @param ss_tol    convergence threshold for max|P_new - P| (default 1e-12).
#'
#' @return list(loglik, grad, t_conv, n_T) -- same loglik/grad contract as
#'   .kf_loglik_adjoint, plus t_conv and n_T for diagnostics.
#' @noRd
.kf_loglik_adjoint_ss <- function(Y, ss, d_ss_list, me_variance = 0,
                                   ss_tol = 1e-12) {

  TT <- ss$TT; RR <- ss$RR; ZZ <- ss$ZZ; DD <- ss$DD
  d  <- as.numeric(ss$d); Sigma_e <- ss$Sigma_e

  n_state <- nrow(TT)
  n_obs   <- nrow(ZZ)
  n_exo   <- ncol(RR)
  n_par   <- length(d_ss_list)

  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  if (anyNA(Y)) stop(".kf_loglik_adjoint_ss: Y must not contain missing values.")
  n_T <- ncol(Y)

  tZZ     <- t(ZZ)
  QQ      <- tcrossprod(RR %*% Sigma_e, RR)
  HH      <- tcrossprod(DD %*% Sigma_e, DD)
  SS      <- RR %*% Sigma_e %*% t(DD)
  me_diag <- me_variance * diag(n_obs)
  ll_const <- -0.5 * n_obs * log(2 * pi)

  fail <- list(loglik = -Inf, grad = rep(NA_real_, n_par),
               t_conv = NA_integer_, n_T = n_T)

  ## -- Stationary initialisation (Lyapunov P_0) -----------------------------
  P0 <- .solve_lyapunov(TT, QQ)
  if (!all(is.finite(P0))) return(fail)

  ## -- Forward pass ----------------------------------------------------------
  ## For t in [1, t_conv]: store {s, P, v, Fi, K, A, B} (full n^2 matrices).
  ## For t in (t_conv, T]:  store only {s_{t-1} (n-vector), v_t (q-vector)}.
  ## Converged constants: P_inf, Fi_inf, K_inf, A_inf, B_inf.

  ## Transient storage (indices 1..t_conv once known)
  s_store  <- vector("list", n_T)   # s_{t-1} entering step t
  P_store  <- vector("list", n_T)   # P_{t-1} entering step t
  v_store  <- vector("list", n_T)   # innovation v_t
  Fi_store <- vector("list", n_T)   # F_t^{-1}
  K_store  <- vector("list", n_T)   # K_t
  A_store  <- vector("list", n_T)   # A_t = TT - K ZZ
  B_store  <- vector("list", n_T)   # B_t = RR - K DD

  ## Locked-tail storage (small vectors; allocated lazily once ss reached)
  sv_s <- vector("list", n_T)   # s_{t-1} for t in tail (n-vector)
  sv_v <- vector("list", n_T)   # v_t     for t in tail (q-vector)

  ## Converged constants (set when ss_reached)
  P_inf  <- NULL; Fi_inf <- NULL; K_inf <- NULL
  A_inf  <- NULL; B_inf  <- NULL; ldf_inf <- NULL

  s <- numeric(n_state)
  P <- P0
  loglik  <- 0.0
  t_conv  <- NA_integer_   # step at which lock is set
  ss_reached <- FALSE

  for (t in seq_len(n_T)) {

    if (!ss_reached) {
      ## ---------- Transient step: full Riccati + storage --------------------
      s_store[[t]] <- s
      P_store[[t]] <- P

      PZ  <- P %*% tZZ
      Ft  <- .sym(ZZ %*% PZ + HH + me_diag)

      Fc  <- tryCatch(chol(Ft), error = function(e) NULL)
      if (is.null(Fc)) return(fail)
      Fi  <- chol2inv(Fc)
      ldf <- 2 * sum(log(diag(Fc)))

      v    <- Y[, t] - d - as.numeric(ZZ %*% s)
      Fiv  <- as.numeric(Fi %*% v)

      ll_t <- ll_const - 0.5 * (ldf + sum(v * Fiv))
      if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) return(fail)
      loglik <- loglik + ll_t

      K <- (TT %*% PZ + SS) %*% Fi
      A <- TT - K %*% ZZ
      B <- RR - K %*% DD

      Fi_store[[t]] <- Fi
      v_store[[t]]  <- v
      K_store[[t]]  <- K
      A_store[[t]]  <- A
      B_store[[t]]  <- B

      s_new <- as.numeric(TT %*% s) + as.numeric(K %*% v)
      P_raw <- tcrossprod(A %*% P, A) + tcrossprod(B %*% Sigma_e, B)
      ## TRUE measurement-noise law (F3-D): P' += K me_diag K'.
      if (me_variance != 0) P_raw <- P_raw + me_variance * tcrossprod(K)
      P_new <- .sym(P_raw)

      ## Convergence check (require t > 1, mirror kalman_standard_loop_cpp)
      if (t > 1 && max(abs(P_new - P)) < ss_tol) {
        ss_reached <- TRUE
        t_conv  <- t
        P_inf   <- P_new
        Fi_inf  <- Fi
        K_inf   <- K
        A_inf   <- A
        B_inf   <- B
        ldf_inf <- ldf     # log-det at the lock step (constant for all tail)
      }

      s <- s_new
      P <- P_new

    } else {
      ## ---------- Locked-tail step: cheap recursion, minimal storage ---------
      sv_s[[t]] <- s           # s_{t-1}

      v   <- Y[, t] - d - as.numeric(ZZ %*% s)
      Fiv <- as.numeric(Fi_inf %*% v)
      ll_t <- ll_const - 0.5 * (ldf_inf + sum(v * Fiv))
      if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) return(fail)
      loglik <- loglik + ll_t

      sv_v[[t]] <- v
      s <- as.numeric(TT %*% s) + as.numeric(K_inf %*% v)
    }
  }

  ## If SS never reached (short series or slow convergence) treat all steps
  ## as transient -- t_conv = n_T, tail is empty.
  if (!ss_reached) t_conv <- n_T

  ## -- Backward sweep -------------------------------------------------------
  ## bar_s and bar_P carry adjoints from step t+1 (or zero at t=T).
  bar_s <- numeric(n_state)
  bar_P <- matrix(0, n_state, n_state)

  G_TT  <- matrix(0, n_state, n_state)
  G_RR  <- matrix(0, n_state, n_exo)
  G_ZZ  <- matrix(0, n_obs, n_state)
  G_DD  <- matrix(0, n_obs, n_exo)
  g_d   <- numeric(n_obs)
  G_Sig <- matrix(0, n_exo, n_exo)

  for (t in n_T:1) {

    ## Retrieve per-step quantities: use locked constants in tail.
    if (t > t_conv) {
      ## Tail step: get locked {Fi, K, A, B, P_prev = P_inf}
      s_prev <- sv_s[[t]]
      P_prev <- P_inf
      v      <- sv_v[[t]]
      Fi     <- Fi_inf
      K      <- K_inf
      A      <- A_inf
      B      <- B_inf
    } else {
      ## Transient step: per-step storage
      s_prev <- s_store[[t]]
      P_prev <- P_store[[t]]
      v      <- v_store[[t]]
      Fi     <- Fi_store[[t]]
      K      <- K_store[[t]]
      A      <- A_store[[t]]
      B      <- B_store[[t]]
    }

    Fiv <- as.numeric(Fi %*% v)

    ## ---- Step 1: adjoint of P_t = sym(A P_prev A' + B Sig B') -------------
    bar_P <- .sym(bar_P)

    bar_A <- 2 * bar_P %*% A %*% P_prev                   # [n x n]
    bar_B <- 2 * bar_P %*% B %*% Sigma_e                  # [n x p]
    bar_P_prev_from_AP <- t(A) %*% bar_P %*% A            # [n x n]
    G_Sig <- G_Sig + t(B) %*% bar_P %*% B                 # [p x p]

    ## ---- Step 2: adjoint of s_t = TT s_prev + K v --------------------------
    G_TT    <- G_TT + outer(bar_s, s_prev)                 # [n x n]
    bar_K   <- outer(bar_s, v)                             # [n x q]

    ## Adjoint of the ME Joseph term in P_t (P_t += K me_diag K'; me_diag is
    ## DATA, not differentiated): with bar_P symmetrized in Step 1,
    ## d tr(bar_P K me K') / dK = 2 bar_P K me.
    if (me_variance != 0)
      bar_K <- bar_K + 2 * me_variance * (bar_P %*% K)     # [n x q]
    bar_v   <- as.numeric(t(K) %*% bar_s)                 # [q]
    bar_s_prev <- as.numeric(t(TT) %*% bar_s)             # [n]

    ## ---- Step 3: adjoint of ll_t = -0.5*(log|F| + v' Fi v) -----------------
    bar_F <- -0.5 * (Fi - outer(Fiv, Fiv))                # [q x q]
    bar_v <- bar_v - Fiv                                   # [q]

    ## ---- Step 4: adjoint of A = TT - K ZZ ----------------------------------
    G_TT  <- G_TT  + bar_A
    G_ZZ  <- G_ZZ  - t(K) %*% bar_A                       # [q x n]
    bar_K <- bar_K - bar_A %*% tZZ                        # [n x q]

    ## ---- Step 5: adjoint of B = RR - K DD ----------------------------------
    G_RR  <- G_RR  + bar_B
    G_DD  <- G_DD  - t(K) %*% bar_B                       # [q x p]
    bar_K <- bar_K - bar_B %*% t(DD)                      # [n x q]

    ## ---- Step 6: adjoint of K = Mnum Fi (Mnum = TT P_prev ZZ' + SS) --------
    bar_F   <- bar_F - t(K) %*% bar_K %*% Fi              # [q x q]
    bar_F   <- .sym(bar_F)
    bar_Mnum <- bar_K %*% Fi                              # [n x q]

    G_TT  <- G_TT  + bar_Mnum %*% ZZ %*% P_prev           # [n x n]
    G_ZZ  <- G_ZZ  + t(bar_Mnum) %*% TT %*% P_prev        # [q x n]
    bar_P_prev_from_Mnum <- t(TT) %*% bar_Mnum %*% ZZ    # [n x n]

    G_RR  <- G_RR  + bar_Mnum %*% DD %*% Sigma_e          # [n x p]
    G_DD  <- G_DD  + t(bar_Mnum) %*% RR %*% Sigma_e       # [q x p]
    G_Sig <- G_Sig + t(RR) %*% bar_Mnum %*% DD            # [p x p]

    ## ---- Step 7: adjoint of F = sym(ZZ P_prev ZZ' + HH + me_diag) ----------
    bar_P_prev_from_F <- t(ZZ) %*% bar_F %*% ZZ           # [n x n]
    G_ZZ  <- G_ZZ  + 2 * bar_F %*% ZZ %*% P_prev         # [q x n]
    G_DD  <- G_DD  + 2 * bar_F %*% DD %*% Sigma_e        # [q x p]
    G_Sig <- G_Sig + t(DD) %*% bar_F %*% DD              # [p x p]

    ## ---- Step 8: adjoint of v_t = y_t - d - ZZ s_prev ----------------------
    g_d      <- g_d - bar_v
    G_ZZ     <- G_ZZ - outer(bar_v, s_prev)               # [q x n]
    bar_s_prev <- bar_s_prev - as.numeric(t(ZZ) %*% bar_v) # [n]

    ## ---- Step 9: collect bar_P_{t-1} ----------------------------------------
    bar_P_prev_new <- .sym(bar_P_prev_from_AP +
                             bar_P_prev_from_Mnum +
                             bar_P_prev_from_F)

    ## ---- Update carry variables for next (earlier) step ---------------------
    bar_s <- bar_s_prev
    bar_P <- bar_P_prev_new
  }

  ## ------------------------------------------------------------------
  ## Adjoint through P_0 = solve_lyapunov(TT, QQ)
  ## ------------------------------------------------------------------
  bar_P0 <- .sym(bar_P)

  if (n_state == 1) {
    denom  <- 1 - TT[1, 1]^2
    bar_QQ <- matrix(bar_P0[1, 1] / denom, 1, 1)
    G_TT[1, 1] <- G_TT[1, 1] + bar_P0[1, 1] * 2 * TT[1, 1] * P0[1, 1] / denom
  } else {
    bar_QQ <- tryCatch(.solve_lyapunov(t(TT), bar_P0), error = function(e) NULL)
    if (is.null(bar_QQ) || !all(is.finite(bar_QQ))) return(fail)
    bar_QQ <- .sym(bar_QQ)
    G_TT   <- G_TT + 2 * bar_QQ %*% TT %*% P0            # [n x n]
  }

  G_Sig <- G_Sig + t(RR) %*% bar_QQ %*% RR               # [p x p]
  G_RR  <- G_RR  + 2 * bar_QQ %*% RR %*% Sigma_e         # [n x p]

  ## ------------------------------------------------------------------
  ## Final contraction
  ## ------------------------------------------------------------------
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

  list(loglik = loglik, grad = grad, t_conv = t_conv, n_T = n_T)
}
