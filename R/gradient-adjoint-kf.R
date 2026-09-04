## R/gradient-adjoint-kf.R
## --------------------------------------------------------------------------
## Analytic ADJOINT (reverse-mode) Kalman filter gradient.
##
## Companion to R/gradient-tangent-kf.R.  Where the tangent filter runs one
## forward-sensitivity pass per parameter (O(n_par) cost), the adjoint runs
## one forward pass that stores intermediate quantities and one backward sweep
## that accumulates gradient matrices G_TT / G_RR / G_ZZ / G_DD / g_d /
## G_Sig, then contracts them against the per-parameter solution derivatives
## in d_ss_list.  Total cost is O(1) in n_par -- the backward sweep visits
## the same per-step arithmetic as one tangent pass, and the final Frobenius
## contractions are O(n_state^2 * n_par).
##
## State-space convention (identical to gradient-tangent-kf.R):
##   QQ = RR Sigma_e RR',  HH = DD Sigma_e DD',  SS = RR Sigma_e DD'
##   v_t = y_t - d - ZZ s_{t-1|t-1}     (s = previous filtered state)
##   F_t = sym(ZZ P_{t-1} ZZ' + HH + me_diag)
##   K_t = (TT P_{t-1} ZZ' + SS) F_t^{-1}
##   A_t = TT - K_t ZZ,   B_t = RR - K_t DD
##   s_t = TT s_{t-1} + K_t v_t
##   P_t = sym(A_t P_{t-1} A_t' + B_t Sigma_e B_t')
##   P_0 = solve_lyapunov(TT, QQ),   s_0 = 0
##
## Adjoint notation: bar_X = d loglik / d X.
## Adjoint of sym(X): bar_X = sym(bar_{sym(X)}).
##
## Key adjoint rules used throughout (derived from first principles):
##
##   P_t = sym(A P A' + B Sig B'):
##     bar_A += bar_P A P + P A' bar_P         [n x n, P sym]
##     bar_B += 2 * bar_P B Sig                [n x p, Sig sym]
##     bar_{P_prev} += A' bar_P A              [n x n]
##     bar_Sig += B' bar_P B                   [p x p]
##
##   s_t = TT s + K v:
##     G_TT  += outer(bar_s, s_prev)           [n x n]
##     bar_K += outer(bar_s, v)                [n x q]
##     bar_v += K' bar_s                       [q]
##     bar_{s_prev} += TT' bar_s               [n]
##
##   ll_t = -0.5*(log|F| + v' Fi v):
##     bar_F += -0.5*(Fi - outer(Fiv, Fiv))    [q x q, Fi=F^{-1}]
##     bar_v += -Fiv                           [q]
##
##   K = Mnum Fi  (Mnum = TT P ZZ' + SS):
##     bar_Mnum = bar_K Fi                     [n x q]
##     bar_F   += -K' bar_K Fi                 [q x q]
##     (bar_Mnum contributes to G_TT, G_ZZ, bar_P_prev via Mnum = TT P ZZ')
##
##   Mnum = TT P ZZ' + SS:
##     G_TT  += bar_Mnum ZZ P_prev             [n x n]   (from TT contribution)
##     G_ZZ  += bar_Mnum' TT P_prev            [q x n]   (from ZZ contribution)
##     bar_P_prev += TT' bar_Mnum ZZ           [n x n]   (from P contribution)
##
##   F = sym(ZZ P ZZ' + HH + me_diag):
##     bar_F symmetrised first; then:
##     G_ZZ  += 2 * bar_F ZZ P_prev            [q x n]
##     G_DD  += 2 * bar_F DD Sigma_e           [q x p]
##     G_Sig += DD' bar_F DD                   [p x p]  (from HH = DD Sig DD')
##     bar_P_prev += ZZ' bar_F ZZ              [n x n]
##
##   A = TT - K ZZ:
##     G_TT  += bar_A                          [n x n]
##     G_ZZ  -= K' bar_A                       [q x n]
##     bar_K -= bar_A ZZ'                      [n x q]
##
##   B = RR - K DD:
##     G_RR  += bar_B                          [n x p]
##     G_DD  -= K' bar_B                       [q x p]
##     bar_K -= bar_B DD'                      [n x q]
##
##   v_t = y - d - ZZ s_prev:
##     g_d   -= bar_v                          [q]
##     G_ZZ  -= outer(bar_v, s_prev)           [q x n]
##     bar_{s_prev} -= ZZ' bar_v              [n]
##
##   P_0 Lyapunov adjoint: solve(TT', bar_P_0) -> bar_QQ; then:
##     G_TT  += 2 * bar_QQ TT P_0             [n x n]
##     G_RR  += 2 * bar_QQ RR Sig             [n x p]  (from QQ = RR Sig RR')
##     G_Sig += RR' bar_QQ RR                 [p x p]
## --------------------------------------------------------------------------

#' TRUE when the compiled adjoint Kalman-filter recursion is available.
#' @noRd
.HAS_RCPP_KF_ADJOINT <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("kf_adjoint_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

#' Analytic adjoint (reverse-mode) Kalman-filter log-likelihood gradient
#'
#' Runs the exact Kalman filter (same recursion as \code{.kf_step} in
#' \code{R/kalman-filter.R}), stores the per-step quantities needed for the
#' backward sweep, then executes one reverse sweep accumulating the gradient
#' matrices (G_TT, G_RR, G_ZZ, G_DD, g_d, G_Sig) and contracts them against
#' \code{d_ss_list} to produce the full parameter gradient.
#'
#' Drop-in replacement for \code{.kf_loglik_tangent}: identical argument list
#' and return value.  The loglik is numerically identical; the gradient agrees
#' with the tangent to ~1e-10 absolute error (both are exact).
#'
#' @param Y          n_obs x n_T observation matrix (no NAs).
#' @param ss         list: TT, RR, ZZ, DD, d, Sigma_e.
#' @param d_ss_list  list of length n_par; element j has (any of) dTT, dRR,
#'                   dZZ, dDD, dd, dSigma_e.
#' @param me_variance scalar measurement-error variance (parameter-independent).
#'
#' @param return_bars when TRUE, additionally return the raw adjoint (bar)
#'   matrices wrt the state-space system as \code{bars = list(G_TT, G_RR,
#'   G_ZZ, G_DD, g_d, G_Sig)} -- the inputs \code{.solution_adjoint()}
#'   (Tier-18 A2) contracts against the analytic primitive derivatives.
#'   Served by both kernels: the compiled fast path exports the bars it
#'   accumulates internally (A2 phase 2); the R kernel remains the spec.
#'
#' @return list(loglik, grad) -- same contract as .kf_loglik_tangent.
#'   With \code{return_bars = TRUE}, also \code{bars}.
#' @noRd
.kf_loglik_adjoint <- function(Y, ss, d_ss_list, me_variance = 0,
                               me_extra = NULL, shock_scale = NULL,
                               return_bars = FALSE) {

  TT <- ss$TT; RR <- ss$RR; ZZ <- ss$ZZ; DD <- ss$DD
  d  <- as.numeric(ss$d); Sigma_e <- ss$Sigma_e

  n_state <- nrow(TT)
  n_obs   <- nrow(ZZ)
  n_exo   <- ncol(RR)
  n_par   <- length(d_ss_list)

  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  if (anyNA(Y)) stop(".kf_loglik_adjoint: Y must not contain missing values.")
  n_T <- ncol(Y)

  has_me_extra    <- !is.null(me_extra)
  ## Any TRUE observation noise at all (base me_variance or per-period
  ## me_extra): both feed F AND the Joseph covariance term (F3-D).
  has_me_true     <- has_me_extra || me_variance != 0
  has_shock_scale <- !is.null(shock_scale)

  ## -- Fast path: compiled adjoint recursion (kf_adjoint_cpp) ---------------
  ## tv inputs (me_extra / shock_scale) are now handled by C++ -- sentinels
  ## (all-ones / all-zeros) are passed when the tv arg is absent so the
  ## constant path is arithmetically equivalent (IEEE 754 exact).
  ## Force-R escape hatch: set option dynhr.use_rcpp = FALSE.
  if (.HAS_RCPP_KF_ADJOINT()) {
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

    ## Sentinels for tv args: all-ones shock_scale => Se_t = Sigma_e (exact);
    ## all-zeros me_extra => me_diag_t = me_diag (exact).
    shock_scale_cpp <- if (has_shock_scale) shock_scale else matrix(1.0, n_exo, n_T)
    me_extra_cpp    <- if (has_me_extra)    me_extra    else matrix(0.0, n_obs, n_T)

    out <- kf_adjoint_cpp(Y, TT, RR, ZZ, DD, d, Sigma_e,
                          dTT_cube, dRR_cube, dZZ_cube, dDD_cube, dd_mat,
                          dSigma_cube, me_variance, .KF_LL_MIN,
                          shock_scale_cpp, me_extra_cpp, return_bars)

    if (!isTRUE(out$ok)) {
      return(list(loglik = -Inf, grad = rep(NA_real_, n_par)))
    }
    res <- list(loglik = out$loglik, grad = as.numeric(out$grad))
    if (return_bars)
      res$bars <- list(G_TT  = out$bars$G_TT, G_RR = out$bars$G_RR,
                       G_ZZ  = out$bars$G_ZZ, G_DD = out$bars$G_DD,
                       g_d   = as.numeric(out$bars$g_d),
                       G_Sig = out$bars$G_Sig)
    return(res)
  }

  ## ---------------------------------------------------------------------------
  ## R reference implementation
  ## ---------------------------------------------------------------------------

  tZZ <- t(ZZ)
  QQ  <- tcrossprod(RR %*% Sigma_e, RR)
  HH  <- tcrossprod(DD %*% Sigma_e, DD)
  SS  <- RR %*% Sigma_e %*% t(DD)
  me_diag  <- me_variance * diag(n_obs)
  ll_const <- -0.5 * n_obs * log(2 * pi)

  fail <- list(loglik = -Inf, grad = rep(NA_real_, n_par))

  ## -- Stationary initialisation (Lyapunov P_0) -----------------------------
  P0 <- .solve_lyapunov(TT, QQ)
  if (!all(is.finite(P0)))
    ## The kron vec-solve's rcond gate fires on HIGHLY NON-NORMAL stable TT
    ## (e.g. Reiter-HANK transition matrices): rcond(I - TT %x% TT) underflows
    ## machine eps while the Lyapunov equation itself is well-posed. Doubling
    ## recovery (solve_lyapunov tries doubling first and still returns NaN for
    ## genuine unit/explosive roots, so the fail contract is preserved).
    P0 <- solve_lyapunov(TT, QQ)
  if (!all(is.finite(P0))) return(fail)

  ## -- Forward pass: run the filter, storing per-step quantities ------------
  ## Stored at step t: s_{t-1} (state entering t), P_{t-1} (covariance),
  ## v_t, Fi_t = F_t^{-1}, K_t, A_t = TT - K ZZ, B_t = RR - K DD.
  ## When shock_scale is active, Se_store[[t]] holds D_t Sigma_e D_t for use
  ## in the backward sweep.
  s_store  <- vector("list", n_T)
  P_store  <- vector("list", n_T)
  v_store  <- vector("list", n_T)
  Fi_store <- vector("list", n_T)
  K_store  <- vector("list", n_T)
  A_store  <- vector("list", n_T)
  B_store  <- vector("list", n_T)
  Se_store <- if (has_shock_scale) vector("list", n_T) else NULL

  s <- numeric(n_state)
  P <- P0
  loglik <- 0.0

  for (t in seq_len(n_T)) {
    s_store[[t]] <- s
    P_store[[t]] <- P

    ## Per-period tv substitutions
    if (has_shock_scale) {
      sc_t <- shock_scale[, t]
      Se_t <- Sigma_e * outer(sc_t, sc_t)
      HH_t <- tcrossprod(DD %*% Se_t, DD)
      SS_t <- RR %*% Se_t %*% t(DD)
      Se_store[[t]] <- Se_t
    } else {
      Se_t <- Sigma_e; HH_t <- HH; SS_t <- SS
    }
    ## Full ME diagonal for this period as a VECTOR: base me_variance plus
    ## this period's me_extra. Both are TRUE observation noise (F3-D), so both
    ## enter F AND the Joseph covariance update.
    me_vec_t <- if (has_me_extra) me_variance + me_extra[, t]
                else rep(me_variance, n_obs)
    me_diag_t <- diag(me_vec_t, n_obs)

    PZ <- P %*% tZZ
    Ft <- ZZ %*% PZ + HH_t + me_diag_t
    Ft <- .sym(Ft)

    Fc <- tryCatch(chol(Ft), error = function(e) NULL)
    if (is.null(Fc)) return(fail)
    Fi  <- chol2inv(Fc)
    ldf <- 2 * sum(log(diag(Fc)))

    v    <- Y[, t] - d - as.numeric(ZZ %*% s)
    Fiv  <- as.numeric(Fi %*% v)

    ll_t <- ll_const - 0.5 * (ldf + sum(v * Fiv))
    if (!is.finite(ll_t) || ll_t < .KF_LL_MIN) return(fail)
    loglik <- loglik + ll_t

    K  <- (TT %*% PZ + SS_t) %*% Fi
    A  <- TT - K %*% ZZ
    B  <- RR - K %*% DD

    Fi_store[[t]] <- Fi
    v_store[[t]]  <- v
    K_store[[t]]  <- K
    A_store[[t]]  <- A
    B_store[[t]]  <- B

    s <- as.numeric(TT %*% s) + as.numeric(K %*% v)
    P <- tcrossprod(A %*% P, A) + tcrossprod(B %*% Se_t, B)
    ## Joseph true-noise term for the FULL ME diagonal (mirrors
    ## kalman_filter's standard path): P' += K diag(me_vec_t) K'.
    if (has_me_true)
      P <- P + (K %*% diag(me_vec_t, n_obs)) %*% t(K)
    P <- .sym(P)
  }

  ## -- Backward sweep -------------------------------------------------------
  ## bar_s and bar_P carry the adjoint of the state (s_t, P_t) that was
  ## produced by step t and consumed by step t+1 (or implicitly by loglik).
  ## Both initialise to zero (loglik does not depend on (s_T, P_T) post-loop).
  ##
  ## Gradient matrices accumulate d loglik / d {TT, RR, ZZ, DD, d, Sigma_e}.
  bar_s   <- numeric(n_state)
  bar_P   <- matrix(0, n_state, n_state)

  G_TT  <- matrix(0, n_state, n_state)
  G_RR  <- matrix(0, n_state, n_exo)
  G_ZZ  <- matrix(0, n_obs, n_state)
  G_DD  <- matrix(0, n_obs, n_exo)
  g_d   <- numeric(n_obs)
  G_Sig <- matrix(0, n_exo, n_exo)

  for (t in n_T:1) {
    s_prev <- s_store[[t]]   # s_{t-1} entering step t
    P_prev <- P_store[[t]]   # P_{t-1} entering step t
    v   <- v_store[[t]]
    Fi  <- Fi_store[[t]]
    K   <- K_store[[t]]
    A   <- A_store[[t]]      # TT - K ZZ
    B   <- B_store[[t]]      # RR - K DD

    ## Retrieve per-period Se_t and sc_t for backward sweep substitutions.
    if (has_shock_scale) {
      sc_t <- shock_scale[, t]
      Se_t <- Se_store[[t]]
    } else {
      Se_t <- Sigma_e
    }

    Fiv <- as.numeric(Fi %*% v)

    ## ---- Step 1: adjoint of P_t = sym(A P_prev A' + B Se_t B') -------------
    ## bar_P enters as adjoint of P_t from the next step.
    ## Adjoint through sym(): bar_P <- sym(bar_P) (already symmetric in
    ## practice due to symmetry of all inputs, but symmetrize to be safe).
    bar_P <- .sym(bar_P)

    ## bar_A = 2 bar_P A P_prev: with W = bar_P symmetric,
    ## d tr(W A P A') = tr(W dA P A') + tr(W A P dA'); cycling each term and
    ## using W = W', P = P' gives identical contributions W A P, hence the 2.
    bar_A <- 2 * bar_P %*% A %*% P_prev                  # [n x n]

    ## bar_B = 2 * bar_P * B * Se_t   [n x p]  (tv: Se_t replaces Sigma_e)
    bar_B <- 2 * bar_P %*% B %*% Se_t                     # [n x p]

    ## bar_P_prev from A P A': d/dP tr(bar_P A P A') = tr(A' bar_P A dP) => += A' bar_P A
    bar_P_prev_from_AP <- t(A) %*% bar_P %*% A            # [n x n]

    ## G_Sig from B Se_t B': accumulate d loglik_t / d Sigma_e.
    ## With Se_t = D_t Sigma_e D_t: d Se_t_{kl} / d Sigma_e_{ij} = sc_t_i sc_t_j.
    ## So d loglik_t / d Sigma_e = outer(sc_t,sc_t) .* (B' bar_P B) when shock_scale active.
    ## When no shock_scale, Se_t = Sigma_e and this reduces to B' bar_P B.
    if (has_shock_scale) {
      G_Sig <- G_Sig + outer(sc_t, sc_t) * (t(B) %*% bar_P %*% B)   # [p x p]
    } else {
      G_Sig <- G_Sig + t(B) %*% bar_P %*% B                           # [p x p]
    }

    ## ---- Step 2: adjoint of s_t = TT s_prev + K v --------------------------
    ## G_TT += outer(bar_s, s_prev): d/dTT tr(bar_s' TT s_prev) => += outer(bar_s, s_prev)
    G_TT <- G_TT + outer(bar_s, s_prev)                   # [n x n]

    ## bar_K_from_s: d/dK tr(bar_s' K v) => += outer(bar_s, v)
    bar_K <- outer(bar_s, v)                               # [n x q]

    ## Adjoint of the ME Joseph term in P_t (P_t += K me K', with
    ## me = diag(me_variance + me_extra[, t]) being DATA, not differentiated):
    ## with bar_P symmetrized above, d tr(bar_P K me K') / dK = 2 bar_P K me.
    if (has_me_true) {
      me_vec_b <- if (has_me_extra) me_variance + me_extra[, t]
                  else rep(me_variance, n_obs)
      bar_K <- bar_K + 2 * bar_P %*% K %*% diag(me_vec_b, n_obs)
    }

    ## bar_v_from_s: d/dv tr(bar_s' K v) => K' bar_s
    bar_v <- as.numeric(t(K) %*% bar_s)                   # [q]

    ## bar_s_prev from TT s_prev: d/ds_prev tr(bar_s' TT s_prev) => TT' bar_s
    bar_s_prev <- as.numeric(t(TT) %*% bar_s)             # [n]

    ## ---- Step 3: adjoint of ll_t = -0.5*(log|F| + v' Fi v) -----------------
    ## bar_F += -0.5*(Fi - outer(Fiv, Fiv))
    bar_F <- -0.5 * (Fi - outer(Fiv, Fiv))                # [q x q]

    ## bar_v += -Fiv
    bar_v <- bar_v - Fiv                                   # [q]

    ## ---- Step 4: adjoint of A = TT - K ZZ ----------------------------------
    ## G_TT += bar_A
    G_TT <- G_TT + bar_A

    ## G_ZZ -= K' bar_A:  d/dZZ tr(bar_A' (-K ZZ)) = tr(-K' bar_A dZZ) => -= K' bar_A
    G_ZZ <- G_ZZ - t(K) %*% bar_A                         # [q x n]

    ## bar_K from A: d/dK tr(bar_A' (-K ZZ)) = tr(-ZZ bar_A' dK) => -= bar_A ZZ'
    bar_K <- bar_K - bar_A %*% tZZ                         # [n x q]

    ## ---- Step 5: adjoint of B = RR - K DD ----------------------------------
    ## G_RR += bar_B
    G_RR <- G_RR + bar_B

    ## G_DD -= K' bar_B
    G_DD <- G_DD - t(K) %*% bar_B                         # [q x p]

    ## bar_K from B: d/dK tr(bar_B' (-K DD)) => -= bar_B DD'
    bar_K <- bar_K - bar_B %*% t(DD)                      # [n x q]

    ## ---- Step 6: adjoint of K = Mnum Fi (Mnum = TT P_prev ZZ' + SS_t) --------
    ## bar_F through Fi = F^{-1}: => bar_F += -K' bar_K Fi   [Fi symmetric]
    bar_F <- bar_F - t(K) %*% bar_K %*% Fi                # [q x q]

    ## Symmetrize bar_F (it comes from symmetric F via sym())
    bar_F <- .sym(bar_F)

    ## bar_Mnum = bar_K Fi
    bar_Mnum <- bar_K %*% Fi                               # [n x q]

    ## From Mnum = TT P_prev ZZ' + SS_t:
    G_TT <- G_TT + bar_Mnum %*% ZZ %*% P_prev             # [n x n]
    G_ZZ <- G_ZZ + t(bar_Mnum) %*% TT %*% P_prev          # [q x n]
    bar_P_prev_from_Mnum <- t(TT) %*% bar_Mnum %*% ZZ     # [n x n]

    ## From SS_t = RR Se_t DD' in Mnum: use Se_t instead of Sigma_e for G_RR/G_DD.
    ## G_Sig from SS_t: outer(sc_t,sc_t) .* (RR' bar_Mnum DD) when tv.
    G_RR  <- G_RR  + bar_Mnum %*% DD %*% Se_t             # [n x p]
    G_DD  <- G_DD  + t(bar_Mnum) %*% RR %*% Se_t          # [q x p]
    if (has_shock_scale) {
      G_Sig <- G_Sig + outer(sc_t, sc_t) * (t(RR) %*% bar_Mnum %*% DD)   # [p x p]
    } else {
      G_Sig <- G_Sig + t(RR) %*% bar_Mnum %*% DD                           # [p x p]
    }

    ## ---- Step 7: adjoint of F = sym(ZZ P_prev ZZ' + HH_t + me_diag_t) ------
    ## bar_F already symmetrized above.
    ## From ZZ P ZZ': bar_P_prev += ZZ' bar_F ZZ
    bar_P_prev_from_F <- t(ZZ) %*% bar_F %*% ZZ            # [n x n]

    ## G_ZZ from ZZ P ZZ': += 2 bar_F ZZ P_prev
    G_ZZ <- G_ZZ + 2 * bar_F %*% ZZ %*% P_prev            # [q x n]

    ## From HH_t = DD Se_t DD': G_DD += 2 bar_F DD Se_t; G_Sig += outer(sc_t,sc_t) .* DD' bar_F DD
    G_DD  <- G_DD  + 2 * bar_F %*% DD %*% Se_t            # [q x p]
    if (has_shock_scale) {
      G_Sig <- G_Sig + outer(sc_t, sc_t) * (t(DD) %*% bar_F %*% DD)       # [p x p]
    } else {
      G_Sig <- G_Sig + t(DD) %*% bar_F %*% DD                               # [p x p]
    }

    ## ---- Step 8: adjoint of v_t = y_t - d - ZZ s_prev ----------------------
    ## g_d -= bar_v (d/d(d) = -1 for each observation equation)
    g_d  <- g_d  - bar_v

    ## G_ZZ from ZZ s_prev: d/dZZ tr(bar_v' (-ZZ s_prev)) = -tr(bar_v' dZZ s_prev)
    ##   = -tr(dZZ s_prev bar_v') [cyclic] => G_ZZ -= outer(bar_v, s_prev)
    G_ZZ <- G_ZZ - outer(bar_v, s_prev)                    # [q x n]

    ## bar_s_prev from ZZ: d/ds tr(bar_v' (-ZZ s)) => -= ZZ' bar_v
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
  ##
  ## Lyapunov equation: P_0 - TT P_0 TT' = QQ.
  ## The adjoint operator L*(bar_QQ) = bar_QQ - TT' bar_QQ TT is the
  ## same system with TT' replacing TT.  So bar_QQ = solve_lyapunov(TT', bar_P_0).
  ##
  ## After bar_QQ, propagate through QQ = RR Sigma_e RR' and through
  ## the TT dependence of P_0 on TT (via the Lyapunov RHS).
  ## ------------------------------------------------------------------
  bar_P0 <- .sym(bar_P)   # bar_P after the full backward loop = bar_P_0

  if (n_state == 1) {
    ## Scalar: P_0 = QQ / (1 - TT^2); d P_0 / dTT = 2 TT QQ / (1-TT^2)^2
    denom   <- 1 - TT[1, 1]^2
    bar_QQ  <- matrix(bar_P0[1, 1] / denom, 1, 1)
    ## G_TT from P_0(TT): bar_P_0 * dP_0/dTT = bar_P_0 * 2 TT P_0 / (1-TT^2)
    G_TT[1, 1] <- G_TT[1, 1] + bar_P0[1, 1] * 2 * TT[1, 1] * P0[1, 1] / denom
  } else {
    ## General: bar_QQ = solve_lyapunov(TT', bar_P_0)
    ## (I - TT'^T ⊗ TT'^T) vec(bar_QQ) = vec(bar_P_0),
    ## which is the Lyapunov equation P - TT' P TT = bar_P_0.
    bar_QQ <- tryCatch(.solve_lyapunov(t(TT), bar_P0), error = function(e) NULL)
    if (is.null(bar_QQ) || !all(is.finite(bar_QQ)))
      ## Same non-normal-TT recovery as the forward P0 solve above; doubling
      ## converges for any symmetric (possibly indefinite) bar_P0 when the
      ## spectral radius is < 1.
      bar_QQ <- tryCatch(solve_lyapunov(t(TT), bar_P0), error = function(e) NULL)
    if (is.null(bar_QQ) || !all(is.finite(bar_QQ))) return(fail)
    bar_QQ <- .sym(bar_QQ)

    ## G_TT from P_0's dependence on TT (via LHS of Lyapunov):
    ## dP_0 / dTT contracts with bar_P_0 as:
    ## <bar_P_0, d/dTT solve_lyapunov(TT, QQ)> = <bar_QQ, dTT P_0 TT' + TT P_0 dTT'>
    ## = 2 tr(bar_QQ dTT P_0 TT') [by symmetry + bar_QQ sym]
    ## = 2 tr(dTT P_0 TT' bar_QQ) [cyclic] => G_TT += 2 (P_0 TT' bar_QQ)^T = 2 bar_QQ TT P_0
    G_TT <- G_TT + 2 * bar_QQ %*% TT %*% P0              # [n x n]
  }

  ## G_Sig from QQ = RR Sigma_e RR': bar_Sig += RR' bar_QQ RR
  ## d/dSig tr(bar_QQ RR dSig RR') = tr(RR' bar_QQ RR dSig) => += RR' bar_QQ RR
  G_Sig <- G_Sig + t(RR) %*% bar_QQ %*% RR              # [p x p]

  ## G_RR from QQ = RR Sigma_e RR': 2 * bar_QQ * RR * Sigma_e
  ## d/dRR tr(bar_QQ (RR Sig RR')) = 2 tr(bar_QQ dRR Sig RR') [sym + bar_QQ sym]
  ## = 2 tr(dRR Sig RR' bar_QQ) [cyclic] => G_RR += 2 (Sig RR' bar_QQ)^T = 2 bar_QQ RR Sig
  G_RR <- G_RR + 2 * bar_QQ %*% RR %*% Sigma_e          # [n x p]

  ## ------------------------------------------------------------------
  ## Final contraction: grad[j] = <G_TT, dTT_j> + ... Frobenius inner products
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

  out <- list(loglik = loglik, grad = grad)
  if (return_bars)
    out$bars <- list(G_TT = G_TT, G_RR = G_RR, G_ZZ = G_ZZ, G_DD = G_DD,
                     g_d = g_d, G_Sig = G_Sig)
  out
}
