## R/smoother-diffuse.R
## ---------------------------------------------------------------------------
## EXACT diffuse smoothing, by the sequential (Koopman-Durbin) treatment on the
## augmented state x_t = [s_{t-1}; eps_t].
##
## WHY AUGMENTED. dynhr's state space has correlated noise -- the shocks enter
## the observation equation directly, y_t = Z s_{t-1} + D eps_t -- and the
## textbook diffuse smoothing recursions are written for the uncorrelated case.
## On the augmented state that correlation disappears: eps_t IS a state
## component, Zb = [Z, D] observes it with no noise beyond the measurement
## error, and the recursions apply as written. The same trick the filter
## already uses (.kf_univariate_dispatch), for the same reason.
##
## It also makes the disturbance smoother free: E[eps_t | y_{1:T}] is just the
## eps block of the smoothed augmented state, so there is no second recursion
## to keep consistent with the first.
##
## WHY SEQUENTIAL. One observable at a time is what makes a diffuse phase and
## missing data compose: a missing observable is skipped, a component with no
## forecast variance is skipped, and a component with diffuse forecast variance
## takes the diffuse branch -- three independent decisions per (period,
## observable) rather than one decision per period about a matrix. The
## multivariate exact-diffuse recursion has to invert F_inf, which is singular
## exactly when some (not all) components of the period are diffuse.
##
## THE RECURSIONS ARE DERIVED, NOT QUOTED. Everything below is the kappa -> Inf
## limit of the ordinary recursion at P = P_star + kappa P_inf, expanded in
## 1/kappa:
##
##   M/F   = K0 + K1/kappa + O(kappa^-2),  K0 = M_inf/F_inf,
##                                         K1 = (M_star - K0 F_star)/F_inf
##   L     = I - (M/F) Z = L0 + L1/kappa,  L0 = I - K0 Z,  L1 = -K1 Z
##   r     = r0 + r1/kappa,   N = N0 + N1/kappa + N2/kappa^2
##
## Matching orders in r = Z'v/F + L'r+ and N = Z'Z/F + L'N+L gives the diffuse
## branch below; the O(1) terms of x = a + P r and V = P - P N P give
##   x_hat = a + P_star r0 + P_inf r1
##   V     = P_star - P_star N0 P_star - P_inf N1 P_star - P_star N1 P_inf
##                  - P_inf N2 P_inf
## and the kappa^2 / kappa^1 terms vanish because P_inf N0 = 0 and
## P_inf N1 P_inf = 0 along the recursion. Deriving it this way is also what
## fixes the convention: the diffuse log-likelihood term is
## -0.5 log F_inf per diffuse observable, with no log(2*pi) and no v^2/F --
## exactly the convention kalman_filter(lik_init = "diffuse",
## method = "univariate") already uses, so the two are directly comparable.
##
## Every L is a rank-one update (L0 = I - K0 Z with Z a single row), so the
## N recursion is done as rank-one updates -- O(nb^2) per observable rather
## than the O(nb^3) three matrix products the literal form would cost.
## ---------------------------------------------------------------------------


#' Exact diffuse sequential smoother on the augmented state
#'
#' @param Y   \code{n_obs x T}, already demeaned. \code{NA} = missing.
#' @param Zb  \code{n_obs x nb} augmented observation matrix \code{[Z, D]}.
#' @param Tb  \code{nb x nb} augmented transition \code{[[T, R], [0, 0]]}.
#' @param G   \code{n_state x nb} map \code{x_t -> s_t}, i.e. \code{[T, R]}.
#' @param Sigma_list length-T list of \code{n_shock x n_shock} shock
#'   covariances (already carrying any \code{shock_scale}).
#' @param a1,P_star1,P_inf1 initial augmented mean / proper covariance /
#'   diffuse covariance for \code{x_1 = [s_0; eps_1]}.
#' @param me  \code{n_obs x T} total measurement-error variance per cell.
#' @param s_idx,e_idx augmented-state index blocks.
#' @return List with the smoothed / filtered / predicted moments in the
#'   ORIGINAL state indexing, plus \code{loglik}, \code{d_diffuse} and the
#'   per-period count of skipped (exactly predictable) components.
#' @noRd
.smoother_diffuse_seq <- function(Y, Zb, Tb, G, Sigma_list, a1, P_star1,
                                  P_inf1, me, s_idx, e_idx,
                                  kalman_tol = 1e-10, diffuse_tol = 1e-10,
                                  conv_tol = 1e-8) {
  n_obs <- nrow(Y); n_T <- ncol(Y); nb <- ncol(Zb)
  n_s   <- length(s_idx); n_e <- length(e_idx)
  tTb   <- t(Tb)
  ## Strip dimnames: Y[i, t] on a named row is a NAMED scalar, and the name
  ## rides the innovation into the log-likelihood -- the same trap
  ## .kf_univariate_loop_R documents. Caught here by test-smoother-singular-F,
  ## where c(filter = ., smoother = .) came back as "smoother.y".
  dimnames(Y) <- NULL

  ## ---- Forward pass -------------------------------------------------------
  a_st <- vector("list", n_T); Ps_st <- vector("list", n_T)
  Pi_st <- vector("list", n_T)
  a_en <- vector("list", n_T); Pe_en <- vector("list", n_T)
  ## Per (period, observable): branch and the quantities the backward pass
  ## needs. 0 = skipped, 1 = diffuse, 2 = ordinary.
  br   <- matrix(0L, n_obs, n_T)
  vv   <- matrix(0,  n_obs, n_T)
  Finf <- matrix(0,  n_obs, n_T)
  Fst  <- matrix(0,  n_obs, n_T)
  K0s  <- vector("list", n_T)      # nb x n_obs
  K1s  <- vector("list", n_T)

  a <- a1; P_star <- P_star1; P_inf <- P_inf1
  diffuse <- max(abs(P_inf)) > 0
  d_diffuse <- NA_integer_
  loglik <- 0
  n_skipped <- integer(n_T)

  for (t in seq_len(n_T)) {
    a_st[[t]] <- a; Ps_st[[t]] <- P_star; Pi_st[[t]] <- P_inf
    K0t <- matrix(0, nb, n_obs); K1t <- matrix(0, nb, n_obs)

    for (i in seq_len(n_obs)) {
      y_i <- Y[i, t]
      if (!is.finite(y_i)) next                      # missing: not an event
      z  <- Zb[i, ]
      v  <- y_i - sum(z * a)
      Ms <- drop(P_star %*% z); Fs <- sum(z * Ms) + me[i, t]
      Mi <- if (diffuse) drop(P_inf %*% z) else numeric(nb)
      Fi <- if (diffuse) sum(z * Mi) else 0

      if (diffuse && Fi > diffuse_tol * max(1, Fs)) {
        K0 <- Mi / Fi
        K1 <- (Ms - K0 * Fs) / Fi
        a  <- a + K0 * v
        P_star <- P_star + tcrossprod(K0) * Fs -
          (tcrossprod(K0, Ms) + tcrossprod(Ms, K0))
        P_inf  <- P_inf - tcrossprod(Mi) / Fi
        loglik <- loglik - 0.5 * log(Fi)
        br[i, t] <- 1L; vv[i, t] <- v; Finf[i, t] <- Fi; Fst[i, t] <- Fs
        K0t[, i] <- K0; K1t[, i] <- K1
      } else if (Fs > kalman_tol) {
        K0 <- Ms / Fs
        a  <- a + K0 * v
        P_star <- P_star - tcrossprod(Ms) / Fs
        loglik <- loglik - 0.5 * (log(2 * pi) + log(Fs) + v * v / Fs)
        br[i, t] <- 2L; vv[i, t] <- v; Fst[i, t] <- Fs
        K0t[, i] <- K0
      } else {
        ## Zero forecast variance: the component is exactly predictable from
        ## what has already been processed and carries no information.
        n_skipped[t] <- n_skipped[t] + 1L
      }
    }
    K0s[[t]] <- K0t; K1s[[t]] <- K1t
    a_en[[t]] <- a; Pe_en[[t]] <- P_star

    ## Time update. The eps block is rebuilt from this period's Sigma_e, which
    ## is what makes shock_scale a per-period quantity here.
    a <- drop(Tb %*% a)
    P_star <- Tb %*% P_star %*% tTb
    P_star[e_idx, e_idx] <- P_star[e_idx, e_idx] +
      Sigma_list[[min(t + 1L, n_T)]]
    P_star <- (P_star + t(P_star)) * 0.5
    if (diffuse) {
      P_inf <- Tb %*% P_inf %*% tTb
      P_inf <- (P_inf + t(P_inf)) * 0.5
      if (max(abs(P_inf)) < conv_tol * max(1, max(abs(P_star)))) {
        diffuse <- FALSE
        d_diffuse <- t
      }
    }
  }

  ## ---- Backward pass ------------------------------------------------------
  ## Every L in this pass is a rank-one update of the identity, L = I - k z',
  ## so each L'(A)L is done in O(nb^2) instead of two nb x nb products. A is
  ## always symmetric here (N0, N1, N2 all are along the recursion).
  LAL <- function(A, k, z) {                    # (I - k z')' A (I - k z')
    Ak <- drop(A %*% k)
    A - outer(z, Ak) - outer(Ak, z) + outer(z, z) * sum(k * Ak)
  }
  LcAL <- function(A, k1, k2, z) {              # (-k1 z')' A (I - k2 z')
    Ak1 <- drop(A %*% k1)
    -outer(z, Ak1) + outer(z, z) * sum(k2 * Ak1)
  }

  r0 <- numeric(nb); r1 <- numeric(nb)
  N0 <- matrix(0, nb, nb); N1 <- matrix(0, nb, nb); N2 <- matrix(0, nb, nb)
  xhat <- vector("list", n_T); Vs <- vector("list", n_T)

  for (t in n_T:1L) {
    K0t <- K0s[[t]]; K1t <- K1s[[t]]
    for (i in n_obs:1L) {
      b <- br[i, t]
      if (b == 0L) next
      z <- Zb[i, ]; u <- K0t[, i]; zz <- outer(z, z)
      if (b == 2L) {
        ## Ordinary step: L = I - K0 z'. r1 / N1 / N2 ride the same L.
        Fs <- Fst[i, t]
        r0 <- z * (vv[i, t] / Fs) + (r0 - z * sum(u * r0))
        r1 <- r1 - z * sum(u * r1)
        N2 <- LAL(N2, u, z)
        N1 <- LAL(N1, u, z)
        N0 <- zz / Fs + LAL(N0, u, z)
      } else {
        ## Diffuse step: L0 = I - K0 z', L1 = -K1 z'.
        w  <- K1t[, i]
        Fi <- Finf[i, t]; Fs <- Fst[i, t]
        r0_in <- r0
        r1 <- z * (vv[i, t] / Fi) + (r1 - z * sum(u * r1)) - z * sum(w * r0_in)
        r0 <- r0_in - z * sum(u * r0_in)
        N0_in <- N0; N1_in <- N1; N2_in <- N2
        cross1 <- LcAL(N1_in, w, u, z)          # L1' N1 L0
        cross0 <- LcAL(N0_in, w, u, z)          # L1' N0 L0
        N2 <- -zz * (Fs / Fi^2) + LAL(N2_in, u, z) + cross1 + t(cross1) +
          zz * sum(w * drop(N0_in %*% w))       # L1' N0 L1
        N1 <- zz / Fi + LAL(N1_in, u, z) + cross0 + t(cross0)
        N0 <- LAL(N0_in, u, z)
      }
    }
    ## Smoothed moments of x_t = [s_{t-1}; eps_t].
    Ps <- Ps_st[[t]]; Pi <- Pi_st[[t]]
    xhat[[t]] <- a_st[[t]] + drop(Ps %*% r0) + drop(Pi %*% r1)
    PiN1Ps <- Pi %*% N1 %*% Ps
    V <- Ps - Ps %*% N0 %*% Ps - PiN1Ps - t(PiN1Ps) - Pi %*% N2 %*% Pi
    Vs[[t]] <- (V + t(V)) * 0.5
    ## Move to period t-1.
    r0 <- drop(tTb %*% r0); r1 <- drop(tTb %*% r1)
    N0 <- tTb %*% N0 %*% Tb; N1 <- tTb %*% N1 %*% Tb; N2 <- tTb %*% N2 %*% Tb
  }

  ## ---- Map back to the original state indexing ---------------------------
  ## s_t = G x_t with G = [T, R]: one linear map turns every augmented moment
  ## into the reported one, and the shocks are already a block of x.
  tG <- t(G)
  ## vapply collapses to a plain vector when the result length is 1, so build
  ## the T x n matrices explicitly rather than transposing -- a one-state model
  ## would otherwise come back 1 x T.
  rows <- function(lst, f, n)
    matrix(vapply(lst, f, numeric(n)), nrow = n_T, ncol = n, byrow = TRUE)
  s_sm <- rows(xhat, function(x) drop(G %*% x), n_s)
  e_sm <- rows(xhat, function(x) x[e_idx], n_e)
  s_fl <- rows(a_en, function(x) drop(G %*% x), n_s)
  s_pr <- rows(a_st, function(x) drop(G %*% x), n_s)
  arr <- function(lst, f) {
    out <- array(0, dim = c(n_s, n_s, n_T))
    for (t in seq_len(n_T)) out[, , t] <- f(lst[[t]])
    out
  }
  GVG <- function(V) G %*% V %*% tG
  list(loglik = loglik,
       smoothed_states = s_sm, smoothed_shocks = e_sm,
       filtered_states = s_fl, predicted_states = s_pr,
       filtered_cov  = arr(Pe_en, GVG),
       predicted_cov = arr(Ps_st, function(P) GVG(P)),
       smoothed_cov  = arr(Vs, GVG),
       smoothed_initial     = xhat[[1L]][s_idx],
       smoothed_initial_cov = Vs[[1L]][s_idx, s_idx, drop = FALSE],
       d_diffuse = d_diffuse, n_skipped = n_skipped,
       ## The sample ended with P_inf still alive: some diffuse direction was
       ## never identified by the data, so the smoothed state is not
       ## determined in it and the reported covariance understates that (it
       ## carries the proper part only). The caller is told; there is no
       ## honest number to return instead.
       diffuse_failed = diffuse)
}
