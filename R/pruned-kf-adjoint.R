## R/pruned-kf-adjoint.R
## --------------------------------------------------------------------------
## B3 (phase a): reverse-mode adjoint of .pruned_kf_correlated() with respect
## to ALL NINE matrix/vector inputs (Tlin, ZZ, d_y, c_drift, QQ, HH, SS, mu0,
## Sxi0). Does NOT chain to theta (that is phase b).
##
## .pruned_kf_correlated is defined in R/pruned-state-space.R:479 -- this
## file mirrors its forward recursion EXACTLY (same partial-NA subsetting,
## same symmetrization points) while stashing per-period intermediates, then
## runs a reverse sweep that accumulates dL/d(input) via standard adjoint
## (reverse-mode AD) rules for matrix products, inverses, and Cholesky
## log-determinants.
##
## Convention: gradients w.r.t. QQ/HH/Sxi0 (and all other matrix inputs) are
## returned as dL/dX_ij treating every entry of X as independent -- i.e. NOT
## folding in any symmetric-parameterization Jacobian. This is exactly the
## convention numDeriv::grad sees when perturbing single entries of the full
## matrix, which is how the validation test in
## tests/testthat/test-pruned-kf-adjoint.R checks it.
## --------------------------------------------------------------------------

## Adjoint cookbook used below (all in "row-major dL/dX_ij" convention, i.e.
## if L is scalar and Y = f(X), dbarX (the adjoint of X) has the same shape
## as X and dbarX_ij = dL/dX_ij):
##
##   Y = A %*% B                =>  dbarA = dbarY %*% t(B);  dbarB = t(A) %*% dbarY
##   Y = t(A)                   =>  dbarA = t(dbarY)
##   y = A %*% x (matrix*vector)=>  dbarA = outer(dbarY, x); dbarx = t(A) %*% dbarY
##   Y = A + B                  =>  dbarA += dbarY; dbarB += dbarY
##   s = sum(x * y) (dot)       =>  dbarx += s_bar * y; dbary += s_bar * x
##   Y = symmetrize(X) = (X+t(X))/2 => dbarX = (dbarY + t(dbarY)) / 2
##   For F_inv = solve(F), with F symmetric and downstream scalar loss:
##     dbarF (from usage of F_inv only) = - F_inv %*% dbar(F_inv) %*% F_inv
##     (then symmetrize since F_inv is symmetric-implied but we treat F_inv's
##      own adjoint using the plain matrix-inverse rule, valid regardless of
##      symmetry of the adjoint itself)
##   For logdet(F) = log(det(F)):  dbar(logdet) contributes s_bar * F_inv to
##     dbarF (standard: d/dF log|F| = F^{-T} = F_inv when F symmetric).
##
## We do NOT explicitly form F_chol adjoints -- log_det_F and F_inv are both
## re-expressed as direct functions of F (log(det(F)) and solve(F)) for
## adjoint purposes, which is mathematically equivalent to going through the
## Cholesky factor and is far simpler to implement correctly.

#' Reverse-mode adjoint of the pruned correlated-noise Kalman filter
#'
#' Computes the scalar log-likelihood exactly as \code{.pruned_kf_correlated}
#' does, while stashing per-period intermediates on a forward pass, then
#' runs a reverse sweep to accumulate \code{d(loglik)/d(input)} for all nine
#' matrix/vector inputs. Does NOT chain through to model parameters (theta);
#' that is a later phase.
#'
#' @param Y n_obs x T data matrix (may contain NA for missing observables).
#' @param Tlin d x d transition matrix.
#' @param ZZ n_obs x d observation matrix.
#' @param d_y n_obs observation intercept.
#' @param c_drift d state drift.
#' @param QQ d x d state noise covariance.
#' @param HH n_obs x n_obs observation noise covariance.
#' @param SS d x n_obs state/obs noise cross-covariance.
#' @param mu0 d initial state mean.
#' @param Sxi0 d x d initial state covariance.
#' @return list(loglik = scalar, grad = list(Tlin=, ZZ=, d_y=, c_drift=,
#'   QQ=, HH=, SS=, mu0=, Sxi0=)), each grad element matching the shape of
#'   its input, holding d(loglik)/d(input)_ij (all entries independent).
#' @noRd
.pruned_kf_correlated_adjoint <- function(Y, Tlin, ZZ, d_y, c_drift, QQ, HH,
                                          SS, mu0, Sxi0) {
  n_obs <- nrow(ZZ)
  n_T   <- ncol(Y)
  d_dim <- nrow(Tlin)

  ## ---- Forward pass: run the filter, stash per-period intermediates ------
  xi <- mu0
  P  <- Sxi0
  loglik <- 0

  ## Per-period record. `kind` is "missing" or "obs".
  recs <- vector("list", n_T)

  for (t in seq_len(n_T)) {
    y_t <- Y[, t]
    o <- which(!is.na(y_t))

    xi_in <- xi   # xi_{t|t-1} entering this period
    P_in  <- P    # P_{t|t-1} entering this period

    if (length(o) == 0L) {
      xi <- as.numeric(Tlin %*% xi_in + c_drift)
      P_pre <- Tlin %*% P_in %*% t(Tlin) + QQ
      P <- (P_pre + t(P_pre)) * 0.5
      recs[[t]] <- list(kind = "missing", xi_in = xi_in, P_in = P_in)
      next
    }

    n_t <- length(o)
    ll_const_t <- -0.5 * n_t * log(2 * pi)

    ZZ_o  <- ZZ[o, , drop = FALSE]
    d_y_o <- d_y[o]
    HH_o  <- HH[o, o, drop = FALSE]
    SS_o  <- SS[, o, drop = FALSE]

    v  <- y_t[o] - d_y_o - as.numeric(ZZ_o %*% xi_in)
    F_pre <- ZZ_o %*% P_in %*% t(ZZ_o) + HH_o
    F  <- (F_pre + t(F_pre)) * 0.5

    F_chol <- tryCatch(chol(F), error = function(e) NULL)
    if (is.null(F_chol))
      stop(".pruned_kf_correlated_adjoint: F not PD at t = ", t)
    log_det_F <- 2 * sum(log(diag(F_chol)))
    F_inv <- chol2inv(F_chol)

    ll_t <- ll_const_t - 0.5 * log_det_F - 0.5 * sum(v * (F_inv %*% v))
    if (!is.finite(ll_t))
      stop(".pruned_kf_correlated_adjoint: non-finite ll_t at t = ", t)
    loglik <- loglik + ll_t

    M <- Tlin %*% P_in %*% t(ZZ_o) + SS_o
    K <- M %*% F_inv
    xi <- as.numeric(Tlin %*% xi_in + c_drift + K %*% v)
    P_pre <- Tlin %*% P_in %*% t(Tlin) + QQ - K %*% t(M)
    P <- (P_pre + t(P_pre)) * 0.5

    recs[[t]] <- list(kind = "obs", o = o, xi_in = xi_in, P_in = P_in,
                       v = v, F = F, F_inv = F_inv, M = M, K = K,
                       ZZ_o = ZZ_o, SS_o = SS_o)
  }

  ## ---- Reverse sweep -------------------------------------------------------
  ## Global adjoints (accumulated across periods).
  g_Tlin    <- matrix(0, d_dim, d_dim)
  g_ZZ      <- matrix(0, n_obs, d_dim)
  g_d_y     <- numeric(n_obs)
  g_c_drift <- numeric(d_dim)
  g_QQ      <- matrix(0, d_dim, d_dim)
  g_HH      <- matrix(0, n_obs, n_obs)
  g_SS      <- matrix(0, d_dim, n_obs)

  ## Adjoints of the *carried* state (xi_{t|t-1}, P_{t|t-1}) flowing backward
  ## from period t+1 into period t. Start at zero: loglik does not depend on
  ## the terminal (post-loop) xi, P directly.
  bar_xi_out <- numeric(d_dim)
  bar_P_out  <- matrix(0, d_dim, d_dim)

  for (t in rev(seq_len(n_T))) {
    rec <- recs[[t]]

    if (rec$kind == "missing") {
      ## Forward: xi_out = Tlin %*% xi_in + c_drift
      ##          P_pre  = Tlin %*% P_in %*% t(Tlin) + QQ
      ##          P_out  = sym(P_pre)
      bar_P_pre <- (bar_P_out + t(bar_P_out)) * 0.5

      ## P_pre = Tlin P_in Tlin' + QQ
      g_QQ <- g_QQ + bar_P_pre
      ## d/dTlin of Tlin %*% P_in %*% t(Tlin): standard quadratic-form rule
      g_Tlin <- g_Tlin + bar_P_pre %*% Tlin %*% t(rec$P_in) +
                          t(bar_P_pre) %*% Tlin %*% rec$P_in
      bar_P_in <- t(Tlin) %*% bar_P_pre %*% Tlin

      ## xi_out = Tlin %*% xi_in + c_drift
      g_c_drift <- g_c_drift + bar_xi_out
      g_Tlin <- g_Tlin + outer(bar_xi_out, rec$xi_in)
      bar_xi_in <- t(Tlin) %*% bar_xi_out

      bar_xi_out <- as.numeric(bar_xi_in)
      bar_P_out  <- bar_P_in
      next
    }

    ## ---- observed period -----------------------------------------------
    o     <- rec$o
    xi_in <- rec$xi_in
    P_in  <- rec$P_in
    v     <- rec$v
    F     <- rec$F
    F_inv <- rec$F_inv
    M     <- rec$M
    K     <- rec$K
    ZZ_o  <- rec$ZZ_o
    SS_o  <- rec$SS_o

    ## == 1. Adjoints from ll_t = const -0.5*log_det_F -0.5*v'F_inv v ======
    ## d ll_t / d log_det_F = -0.5   =>  dbar F (via logdet) = -0.5 * F_inv
    bar_F <- -0.5 * F_inv

    ## d ll_t / d(v' F_inv v):
    Finv_v <- F_inv %*% v
    ## contribution to bar_v (quadratic form v'Av, A=F_inv, symmetric):
    bar_v <- -0.5 * 2 * as.numeric(Finv_v)   ## = -Finv_v
    ## Direct closed-form for q = v' F_inv v treated as a function of F
    ## alone (v is a FIXED vector here, not re-differentiated through
    ## F_inv): d(v' F^{-1} v)/dF = -F^{-1} v v' F^{-1}, so this quadratic
    ## term contributes -0.5 * (-F_inv v v' F_inv) = +0.5*outer(Finv_v,Finv_v)
    ## to bar_F directly. (NB: using bar_F_inv = -0.5*outer(v, v) -- i.e.
    ## treating q = sum_ij F_inv_ij v_i v_j linearly in F_inv -- and then
    ## folding via bar_F += -F_inv %*% bar_F_inv %*% F_inv gives the SAME
    ## answer; outer(Finv_v, Finv_v) is just F_inv %*% outer(v,v) %*% F_inv
    ## expanded. Do NOT use outer(Finv_v, Finv_v) as bar_F_inv itself --
    ## that double-applies F_inv and is wrong.)
    bar_F_inv <- -0.5 * outer(v, v)
    bar_F <- bar_F - t(F_inv) %*% bar_F_inv %*% t(F_inv)

    ## == 2. Adjoints from the "out" state (xi_out, P_out) via chain from
    ##       later periods (bar_xi_out, bar_P_out coming from t+1) =========
    bar_P_pre <- (bar_P_out + t(bar_P_out)) * 0.5   ## sym adjoint
    bar_xi_full <- bar_xi_out                        ## xi_out already a vector

    ## P_pre = Tlin P_in Tlin' + QQ - K M'
    g_QQ <- g_QQ + bar_P_pre
    g_Tlin <- g_Tlin + bar_P_pre %*% Tlin %*% t(P_in) +
                        t(bar_P_pre) %*% Tlin %*% P_in
    bar_P_in_from_Ppre <- t(Tlin) %*% bar_P_pre %*% Tlin

    ## -K M' term: Y = -K %*% t(M). For Y = A %*% t(B): dbarA = bar_Y %*% B,
    ## dbarB = t(bar_Y) %*% A.
    bar_K_from_P <- -bar_P_pre %*% M
    bar_M_from_P <- -t(bar_P_pre) %*% K

    ## xi_out = Tlin %*% xi_in + c_drift + K %*% v
    g_c_drift <- g_c_drift + bar_xi_full
    g_Tlin <- g_Tlin + outer(bar_xi_full, xi_in)
    bar_xi_in_from_xiout <- as.numeric(t(Tlin) %*% bar_xi_full)
    bar_K_from_xi <- outer(bar_xi_full, v)
    bar_v_from_xi <- as.numeric(t(K) %*% bar_xi_full)

    bar_K <- bar_K_from_P + bar_K_from_xi
    bar_M <- bar_M_from_P

    ## == 3. K = M %*% F_inv ================================================
    bar_M <- bar_M + bar_K %*% t(F_inv)
    bar_F_inv <- bar_F_inv + t(M) %*% bar_K
    ## fold this additional F_inv adjoint back into F
    bar_F <- bar_F - F_inv %*% (t(M) %*% bar_K) %*% F_inv

    ## == 4. M = Tlin %*% P_in %*% t(ZZ_o) + SS_o ===========================
    g_SS[, o] <- g_SS[, o] + bar_M
    ## Tlin %*% P_in %*% t(ZZ_o): treat as A %*% B %*% t(C), A=Tlin,B=P_in,C=ZZ_o
    AB <- Tlin %*% P_in                       # d x d
    g_Tlin <- g_Tlin + bar_M %*% ZZ_o %*% t(P_in)
    bar_P_in_from_M <- t(Tlin) %*% bar_M %*% ZZ_o
    g_ZZ[o, ] <- g_ZZ[o, ] + t(bar_M) %*% AB

    ## == 5. v = y_o - d_y_o - ZZ_o %*% xi_in ===============================
    bar_v_total <- bar_v + bar_v_from_xi
    g_d_y[o] <- g_d_y[o] - bar_v_total
    g_ZZ[o, ] <- g_ZZ[o, ] - outer(bar_v_total, xi_in)
    bar_xi_in_from_v <- as.numeric(-t(ZZ_o) %*% bar_v_total)

    ## == 6. F = sym(F_pre), F_pre = ZZ_o %*% P_in %*% t(ZZ_o) + HH_o =======
    bar_F_pre <- (bar_F + t(bar_F)) * 0.5
    g_HH[o, o] <- g_HH[o, o] + bar_F_pre
    g_ZZ[o, ] <- g_ZZ[o, ] + bar_F_pre %*% ZZ_o %*% t(P_in) +
                              t(bar_F_pre) %*% ZZ_o %*% P_in
    bar_P_in_from_F <- t(ZZ_o) %*% bar_F_pre %*% ZZ_o

    ## == combine all bar_P_in and bar_xi_in contributions ==================
    bar_P_in <- bar_P_in_from_Ppre + bar_P_in_from_M + bar_P_in_from_F
    bar_xi_in <- bar_xi_in_from_xiout + bar_xi_in_from_v

    bar_xi_out <- bar_xi_in
    bar_P_out  <- bar_P_in
  }

  ## After the loop, bar_xi_out / bar_P_out are the adjoints of xi_{1|0},
  ## P_{1|0}, i.e. exactly mu0 and Sxi0.
  g_mu0  <- as.numeric(bar_xi_out)
  g_Sxi0 <- bar_P_out

  list(
    loglik = loglik,
    grad = list(
      Tlin    = g_Tlin,
      ZZ      = g_ZZ,
      d_y     = g_d_y,
      c_drift = g_c_drift,
      QQ      = g_QQ,
      HH      = g_HH,
      SS      = g_SS,
      mu0     = g_mu0,
      Sxi0    = g_Sxi0
    )
  )
}
