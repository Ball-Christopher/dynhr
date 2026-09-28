## R/sv-grid-filter.R
## --------------------------------------------------------------------------
## Deterministic likelihood ORACLES for the SV-on-shocks RB-PF (R/sv-rbpf.R).
##
## The RB-PF returns a particle estimate whose exponential is unbiased for the
## marginal likelihood p(y_{1:T}) = E_h[ p(y_{1:T} | h_{1:T}) ], where
## p(y | h) is the exact Kalman likelihood given the volatility path and h is
## the latent AR(1) log-variance of each SV shock (R/stochastic-volatility.R).
## Two independent computations of that same p(y) live here:
##
##  1. .sv_grid_loglik() -- Kitagawa's non-Gaussian filter (Kitagawa 1987;
##     re-examined in arXiv 2607.27576, Sec. 2.2 eqs. (9)-(12) and the
##     renormalisation (14)): the predictive / filtering densities of h_t are
##     represented on a fixed grid, the prediction is a discrete convolution
##     with the Gaussian AR(1) kernel, and the filter update multiplies by the
##     observation density at each node. log p(y_t | y_{1:t-1}) is the log of
##     the update's normalising constant; the sum is a smooth, deterministic
##     log-likelihood. Deviation from the paper: the kernel is POINT-evaluated
##     at the grid nodes (trapezoidal rule) rather than cell-averaged; for the
##     smooth Gaussian kernel the trapezoidal rule converges spectrally in the
##     node spacing, which is what makes this usable as an oracle. The kernel
##     columns are renormalised to unit mass (the paper's eq. (14) step), so the
##     truncation to mu +/- width * sd_stationary loses no probability.
##
##  2. .sv_chan_loglik() -- Chan (arXiv 2608.21619) tangent-twisted proposal,
##     Algorithm 1 (Sec. 3.4) with the closed forms of Sec. 3.1 (eq. (8)-(10))
##     and the "practical grid" of Sec. 3.2 (nodes m_t + w_j s_t, w on
##     [-4, 6], m = posterior mode, s = Gaussian-approximation sd). Its
##     importance weights are bounded by the attained constant C_1
##     (Theorem 1, Corollary 1), so mean(w) is an UNBIASED likelihood
##     estimator with relative variance <= (1/p - 1)/M -- a second, stochastic
##     but exact-in-expectation oracle, and log C_1 >= log p(y) is a
##     deterministic upper bound.
##
## EXACTNESS CLASS -- read this before using either oracle.
## Both methods need the observations to be conditionally independent given
## the CURRENT volatility state: p(y_t | h_{1:t}, y_{1:t-1}) = p(y_t | h_t,
## y_{1:t-1}). In an SV-DSGE that holds only when the Kalman filtering
## distribution of the Gaussian block is the SAME for every volatility path,
## i.e. the post-update (s_{t+1}, P_{t+1}) do not depend on h_t. That is the
## "state revealed by the observables" class: e.g. AR(1)/VAR(1) observables
## without measurement error, where x_t = y_t is known after the update and
## P_{t+1} = 0 (kf_step's cancellation snap makes it an exact zero). Outside
## it -- a latent Gaussian state, or measurement error -- the conditional
## Gaussian covariance carries the whole h-path, a grid over h_t with ONE
## Kalman filter per node is a moment-collapsing APPROXIMATION (not an
## oracle), and Chan's Markov/concavity assumptions fail. .sv_revealed_path()
## therefore VERIFIES the class at every period (kf_step at probe volatility
## values spanning +/- 3 stationary sd per SV shock must give identical
## (s', P')) and aborts with class "dynhr_sv_oracle_not_exact" otherwise.
## Within the class, period 1 is still exact: its (s, P) = (0, P0) prior is
## shared by every path.
##
## Chan additionally requires a SCALAR state (exactly one SV shock) and a
## concave observation log-density l_t(h). With one observable,
## F_t(h) = a_t + b e^h, and (derived in the comments of
## .sv_chan_obs_logdens) l_t is concave on all of R iff a_t = 0 or
## v_t^2 <= a_t. In the revealed class a_t = 0 for t >= 2 but a_1 = ZZ P0 ZZ'
## > 0, so the FIRST innovation must satisfy v_1^2 <= a_1; otherwise the
## method does not apply and .sv_chan_loglik() aborts
## ("dynhr_sv_chan_not_logconcave") rather than return a number whose bound
## (Theorem 1) no longer holds.
##
## Scope: internal (noRd). These are TEST ORACLES for a narrow exact class,
## not general-purpose estimators; exporting them would advertise a
## likelihood that silently cannot exist for the latent-state DSGEs users
## actually estimate. tests/testthat/test-fix-0925-sv-oracles.R drives them.
## --------------------------------------------------------------------------


## ---- shared: exact-class Kalman path ---------------------------------------

#' Kalman path shared by every volatility path (exactness-class check)
#'
#' Runs the prediction-form Kalman recursion of \code{kf_step()} and, at each
#' period, verifies that the post-update moments do not depend on the
#' volatility state (probe scales at mu + c(-3, 0, 3) stationary sd per SV
#' shock). Returns the per-period innovation and the volatility-FREE part of
#' the forecast covariance, from which the per-node observation density is
#' \eqn{N(v_t; 0, A_t + DD\, Se(h)\, DD')}.
#'
#' @inheritParams .sv_rbpf_loglik
#' @param tol Relative tolerance on the node-independence of (s', P').
#' @return \code{list(v = n_obs x T, A = list of T n_obs x n_obs)}.
#' @noRd
.sv_revealed_path <- function(Y, TT, ZZ, RR, DD, Sigma_e, d, P0, sv_idx, hyper,
                              n_exo, me_diag = NULL, tol = 1e-8) {
  n_obs <- nrow(Y); n_T <- ncol(Y); n_s <- nrow(TT); n_sv <- length(sv_idx)
  mu <- hyper[, "mu"]; rho <- hyper[, "rho"]; seta <- hyper[, "sigma_eta"]
  sd0 <- seta / sqrt(1 - rho^2)

  ## Probe volatility states: every combination of mu + {-3, 0, 3} sd0.
  probe_h <- as.matrix(expand.grid(lapply(seq_len(n_sv), function(k)
    mu[k] + c(0, -3, 3) * sd0[k])))
  probe_scale <- lapply(seq_len(nrow(probe_h)), function(i) {
    sc <- rep(1, n_exo); sc[sv_idx] <- exp(probe_h[i, ] / 2); sc })

  me_mat <- if (is.null(me_diag)) 0 else diag(rep_len(me_diag, n_obs),
                                              nrow = n_obs)
  s <- numeric(n_s); P <- P0
  v_out <- matrix(0, n_obs, n_T)
  A_out <- vector("list", n_T)
  for (t in seq_len(n_T)) {
    v <- Y[, t] - as.numeric(ZZ %*% s)
    if (!is.null(d)) v <- v - d
    A <- ZZ %*% P %*% t(ZZ) + me_mat
    v_out[, t] <- v
    A_out[[t]] <- (A + t(A)) * 0.5

    steps <- lapply(probe_scale, function(sc)
      kf_step(s, P, Y[, t], TT, ZZ, RR, DD, Sigma_e, scale = sc, d = d,
              me_diag = me_diag))
    if (any(vapply(steps, is.null, logical(1))))
      .dynhr_abort(sprintf(paste0(
        ".sv_revealed_path: non-PD forecast covariance at a probe volatility ",
        "state in period %d."), t), class = "dynhr_sv_oracle_not_exact")
    s_ref <- steps[[1L]]$s; P_ref <- steps[[1L]]$P
    s_scl <- 1 + max(abs(s_ref), abs(s))
    P_scl <- 1 + max(abs(P_ref), abs(P))
    dev <- max(vapply(steps, function(st)
      max(max(abs(st$s - s_ref)) / s_scl, max(abs(st$P - P_ref)) / P_scl),
      numeric(1)))
    if (dev > tol)
      .dynhr_abort(sprintf(paste0(
        "SV likelihood oracle: the Kalman filtering distribution depends on the ",
        "volatility state at period %d (relative deviation %.2e across probe ",
        "states). The grid / Chan oracles are exact only when the observables ",
        "reveal the Gaussian state (no latent state, no measurement error); ",
        "for this model a per-node Kalman filter would be an approximation."),
        t, dev), class = "dynhr_sv_oracle_not_exact")
    s <- s_ref; P <- P_ref
  }
  list(v = v_out, A = A_out)
}


## ---- 1. Kitagawa non-Gaussian grid filter ----------------------------------

#' Log observation density at every grid node (vectorised over nodes)
#'
#' @param v Innovation, length n_obs.
#' @param A Volatility-free forecast covariance, n_obs x n_obs.
#' @param DD,Sigma_e Shock loading and baseline covariance.
#' @param scale n_exo x N matrix of shock scales, one column per node.
#' @return Length-N vector of \eqn{\log N(v; 0, A + DD Se DD')}; \code{-Inf}
#'   where the forecast covariance is not positive definite.
#' @noRd
.sv_grid_obs_loglik <- function(v, A, DD, Sigma_e, scale) {
  n_obs <- length(v); n_exo <- ncol(DD); N <- ncol(scale)
  ## F_flat[, i] = vec(A) + sum_{k,l} Sigma_kl s_k s_l vec(D_k D_l').
  F_flat <- matrix(as.numeric(A), n_obs * n_obs, N)
  for (k in seq_len(n_exo)) for (l in seq_len(n_exo)) {
    if (Sigma_e[k, l] == 0) next
    M_kl <- as.numeric(tcrossprod(DD[, k], DD[, l])) * Sigma_e[k, l]
    F_flat <- F_flat + outer(M_kl, scale[k, ] * scale[l, ])
  }
  c0 <- -0.5 * n_obs * log(2 * pi)
  if (n_obs == 1L) {
    f <- F_flat[1L, ]
    out <- ifelse(f > 0, c0 - 0.5 * (log(f) + v^2 / f), -Inf)
  } else if (n_obs == 2L) {
    f11 <- F_flat[1L, ]; f21 <- F_flat[2L, ]; f22 <- F_flat[4L, ]
    det <- f11 * f22 - f21^2
    maha <- (f22 * v[1]^2 - 2 * f21 * v[1] * v[2] + f11 * v[2]^2) / det
    out <- ifelse(f11 > 0 & det > 0, c0 - 0.5 * (log(det) + maha), -Inf)
  } else {
    out <- vapply(seq_len(N), function(i) {
      Fi <- matrix(F_flat[, i], n_obs, n_obs)
      ev <- eigen((Fi + t(Fi)) * 0.5, symmetric = TRUE)
      if (!all(is.finite(ev$values)) || min(ev$values) <= 0) return(-Inf)
      c0 - 0.5 * (sum(log(ev$values)) +
                    sum(drop(crossprod(ev$vectors, v))^2 / ev$values))
    }, numeric(1))
  }
  out[!is.finite(out)] <- -Inf
  out
}


#' Deterministic grid-filter log-likelihood for SV-on-shocks (oracle)
#'
#' Kitagawa's non-Gaussian filter over the latent log-variances (one or two
#' SV shocks), Rao-Blackwellised over the Gaussian DSGE block, which in the
#' exactness class is a single shared Kalman path (see the file header).
#' Targets exactly the quantity whose exponential the RB-PF estimates without
#' bias: same AR(1) law for h (stationary start), same P0, same
#' \code{kf_step()} conditional likelihood.
#'
#' @inheritParams .sv_rbpf_loglik
#' @param n_grid Grid points per SV shock (recycled to \code{length(sv_idx)}).
#' @param width  Half-width of the grid in stationary sd of h.
#' @return Scalar log-likelihood (\code{-Inf} if every node is infeasible).
#' @noRd
.sv_grid_loglik <- function(Y, TT, ZZ, RR, DD, Sigma_e, d, P0, sv_idx, hyper,
                            n_exo, n_grid = 200L, width = 8, me_diag = NULL) {
  n_sv <- length(sv_idx)
  if (n_sv < 1L || n_sv > 2L)
    .dynhr_abort(".sv_grid_loglik: supports 1 or 2 SV shocks (got ", n_sv,
                 "); the grid cost grows as n_grid^n_sv.")
  n_grid <- rep_len(as.integer(n_grid), n_sv)
  if (any(n_grid < 3L))
    .dynhr_abort(".sv_grid_loglik: 'n_grid' must be >= 3.")
  if (!is.matrix(Y)) Y <- matrix(Y, nrow = 1L)

  path <- .sv_revealed_path(Y, TT, ZZ, RR, DD, Sigma_e, d, P0, sv_idx, hyper,
                            n_exo, me_diag = me_diag)

  mu <- hyper[, "mu"]; rho <- hyper[, "rho"]; seta <- hyper[, "sigma_eta"]
  sd0 <- seta / sqrt(1 - rho^2)
  grids <- lapply(seq_len(n_sv), function(k)
    mu[k] + sd0[k] * seq(-width, width, length.out = n_grid[k]))
  ## Transition kernel K[i, j] = P(h_t = x_i | h_{t-1} = x_j), point-evaluated
  ## (trapezoidal rule) and column-normalised to unit mass (eq. (14)).
  kern <- lapply(seq_len(n_sv), function(k) {
    x <- grids[[k]]
    K <- stats::dnorm(outer(x, mu[k] + rho[k] * (x - mu[k]), "-"),
                      sd = seta[k])
    sweep(K, 2L, colSums(K), "/")
  })
  init <- lapply(seq_len(n_sv), function(k) {
    w <- stats::dnorm(grids[[k]], mu[k], sd0[k]); w / sum(w) })

  ## Node shock scales, nodes enumerated with the FIRST SV shock fastest
  ## (column-major, matching a G1 x G2 mass matrix).
  node_h <- as.matrix(expand.grid(grids))
  scale <- matrix(1, n_exo, nrow(node_h))
  for (k in seq_len(n_sv)) scale[sv_idx[k], ] <- exp(node_h[, k] / 2)

  ## Filtered masses at t = 0: the stationary law (the RB-PF draws h_0 there).
  filt <- if (n_sv == 1L) init[[1L]] else outer(init[[1L]], init[[2L]])
  loglik <- 0
  for (t in seq_len(ncol(Y))) {
    pred <- if (n_sv == 1L) drop(kern[[1L]] %*% filt) else
      kern[[1L]] %*% filt %*% t(kern[[2L]])
    ll_node <- .sv_grid_obs_loglik(path$v[, t], path$A[[t]], DD, Sigma_e, scale)
    lw <- log(as.numeric(pred)) + ll_node
    m <- max(lw)
    if (!is.finite(m)) return(-Inf)
    w <- exp(lw - m)
    loglik <- loglik + m + log(sum(w))
    filt <- w / sum(w)
    if (n_sv == 2L) filt <- matrix(filt, n_grid[1L], n_grid[2L])
  }
  loglik
}


## ---- 2. Chan (2608.21619) tangent-twisted unbiased likelihood --------------

#' Stable log(1 - exp(x)) for x <= 0
#' @noRd
.sv_chan_log1mexp <- function(x) {
  ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x)))
}

#' log(Phi(hi) - Phi(lo)) evaluated in the tail that keeps precision
#' @noRd
.sv_chan_log_pdiff <- function(lo, hi) {
  up <- lo > 0
  out <- numeric(length(lo))
  if (any(up)) {                                  # both in the upper tail
    ls_lo <- stats::pnorm(lo[up], lower.tail = FALSE, log.p = TRUE)
    ls_hi <- stats::pnorm(hi[up], lower.tail = FALSE, log.p = TRUE)
    out[up] <- ls_lo + .sv_chan_log1mexp(ls_hi - ls_lo)
  }
  if (any(!up)) {
    lp_lo <- stats::pnorm(lo[!up], log.p = TRUE)
    lp_hi <- stats::pnorm(hi[!up], log.p = TRUE)
    out[!up] <- lp_hi + .sv_chan_log1mexp(lp_lo - lp_hi)
  }
  out[!(hi > lo)] <- -Inf
  out
}

#' Standard-normal draws truncated to (lo, hi], by inversion in log space
#' @noRd
.sv_chan_rtrunc <- function(lo, hi) {
  n <- length(lo); U <- stats::runif(n); x <- numeric(n)
  up <- lo > 0
  if (any(up)) {
    ls_lo <- stats::pnorm(lo[up], lower.tail = FALSE, log.p = TRUE)
    ls_hi <- stats::pnorm(hi[up], lower.tail = FALSE, log.p = TRUE)
    lp <- ls_lo + log1p(-U[up] * (1 - exp(ls_hi - ls_lo)))
    x[up] <- stats::qnorm(lp, lower.tail = FALSE, log.p = TRUE)
  }
  if (any(!up)) {
    lp_lo <- stats::pnorm(lo[!up], log.p = TRUE)
    lp_hi <- stats::pnorm(hi[!up], log.p = TRUE)
    r <- exp(lp_lo - lp_hi)
    lp <- lp_hi + log(r + U[!up] * (1 - r))
    x[!up] <- stats::qnorm(lp, log.p = TRUE)
  }
  pmin(pmax(x, lo), hi)
}

#' Observation log-density l_t(u), l_t'(u), l_t''(u) for one observable
#'
#' F(u) = a + b e^u. With a = 0: l = -0.5 log(2 pi b) - u/2 - (v^2/b) e^-u / 2
#' (Chan's leading SV example). With a > 0, write c = b/a, kappa = v^2/a,
#' q = plogis(u + log c): l = -0.5 log(2 pi a) - kappa/2 - 0.5 log(1 + c e^u)
#' + kappa q / 2, l' = q (kappa (1 - q) - 1) / 2,
#' l'' = q (1 - q) (kappa (1 - 2q) - 1) / 2, which is <= 0 for every u iff
#' kappa <= 1 (the q -> 0 limit is the binding one).
#' @noRd
.sv_chan_obs_logdens <- function(u, a, b, v) {
  if (a == 0) {
    r <- v^2 / b; e <- r * exp(-u)
    return(list(l = -0.5 * log(2 * pi * b) - 0.5 * u - 0.5 * e,
                d1 = -0.5 + 0.5 * e, d2 = -0.5 * e))
  }
  kap <- v^2 / a; z <- u + log(b / a)
  q <- stats::plogis(z)
  sp <- ifelse(z > 0, z + log1p(exp(-z)), log1p(exp(z)))
  list(l = -0.5 * log(2 * pi * a) - 0.5 * kap - 0.5 * sp + 0.5 * kap * q,
       d1 = 0.5 * q * (kap * (1 - q) - 1),
       d2 = 0.5 * q * (1 - q) * (kap * (1 - 2 * q) - 1))
}

#' Lower envelope of tangent lines (Chan Lemma 3 + the merge/prune rule)
#'
#' @param v,f,fp Nodes, f at nodes, f' at nodes.
#' @return list(g = slopes (decreasing), al = intercepts, z = breakpoints
#'   with z[1] = -Inf, z[J + 1] = Inf).
#' @noRd
.sv_chan_envelope <- function(v, f, fp) {
  g <- fp; al <- f - fp * v
  o <- order(-g, al); g <- g[o]; al <- al[o]
  ## Merge (numerically) equal slopes, keeping the lowest intercept: a lower
  ## tangent is still >= f everywhere, so psi >= f is preserved.
  keep <- c(TRUE, diff(g) < -1e-12 * (1 + abs(g[-1L])))
  if (!all(keep)) {
    grp <- cumsum(keep)
    al <- as.numeric(tapply(al, grp, min)); g <- g[keep]
  }
  ## Lower envelope of lines ordered by decreasing slope (active left->right).
  xint <- function(i, j) (al[j] - al[i]) / (g[i] - g[j])
  st <- 1L
  for (j in seq_along(g)[-1L]) {
    while (length(st) >= 2L &&
           xint(st[length(st) - 1L], st[length(st)]) >= xint(st[length(st)], j))
      st <- st[-length(st)]
    st <- c(st, j)
  }
  g <- g[st]; al <- al[st]
  J <- length(g)
  z <- c(-Inf, if (J > 1L) (al[-1L] - al[-J]) / (g[-J] - g[-1L]), Inf)
  list(g = g, al = al, z = z)
}

#' Twisted-transition log weights (Chan eq. (8)) at conditioning means m
#'
#' @param env Envelope from .sv_chan_envelope().
#' @param m   Vector of transition means a_t(u).
#' @param om  Transition variance.
#' @return list(lw = length(m) x J log weights, b = component means,
#'   lo/hi = standardised interval bounds, logC = row log-sums).
#' @noRd
.sv_chan_weights <- function(env, m, om) {
  J <- length(env$g); n <- length(m); so <- sqrt(om)
  b  <- outer(m, env$g * om, "+")
  zl <- matrix(env$z[-(J + 1L)], n, J, byrow = TRUE)
  zh <- matrix(env$z[-1L],       n, J, byrow = TRUE)
  lo <- (zl - b) / so; hi <- (zh - b) / so
  lw <- matrix(env$al + 0.5 * env$g^2 * om, n, J, byrow = TRUE) +
    outer(m, env$g) +
    matrix(.sv_chan_log_pdiff(as.numeric(lo), as.numeric(hi)), n, J)
  mx <- apply(lw, 1L, max)
  logC <- mx + log(rowSums(exp(lw - mx)))
  list(lw = lw, b = b, lo = lo, hi = hi, logC = logC)
}

#' Chan tangent-twisted unbiased likelihood estimator for SV-on-shocks
#'
#' Scalar-state (one SV shock), one-observable models in the exactness class
#' (file header). Builds the tangent twists by Chan's backward recursion,
#' draws \code{n_draws} independent paths from the twisted proposal and
#' returns the log of the mean importance weight -- an unbiased estimator of
#' the likelihood (Corollary 1) -- together with the dominating constant.
#'
#' @inheritParams .sv_rbpf_loglik
#' @param n_draws Number of proposal paths M.
#' @param n_nodes Nodes per date G (default \code{max(10, ceiling(4 sqrt(T)))},
#'   Chan's G proportional to sqrt(T)).
#' @param node_range The fixed grid w on which nodes m_t + w s_t sit.
#' @return list(loglik = log mean weight, log_C1 = log dominating constant
#'   (>= the true log-likelihood), accept_rate = Rao-Blackwellised estimate of
#'   the rejection-sampling acceptance probability, max_resid = largest
#'   sum_t d_t(h_t) seen (<= 0 in exact arithmetic), log_w = the M log
#'   weights).
#' @noRd
.sv_chan_loglik <- function(Y, TT, ZZ, RR, DD, Sigma_e, d, P0, sv_idx, hyper,
                            n_exo, n_draws = 1000L, n_nodes = NULL,
                            node_range = c(-4, 6), me_diag = NULL) {
  if (!is.matrix(Y)) Y <- matrix(Y, nrow = 1L)
  if (length(sv_idx) != 1L)
    .dynhr_abort(".sv_chan_loglik: Chan's exact method needs a SCALAR state ",
                 "(exactly one SV shock); got ", length(sv_idx), ".",
                 class = "dynhr_sv_chan_not_applicable")
  if (nrow(Y) != 1L)
    .dynhr_abort(".sv_chan_loglik: implemented for one observable (got ",
                 nrow(Y), ").", class = "dynhr_sv_chan_not_applicable")
  k <- sv_idx
  if (any(Sigma_e[k, -k] != 0))
    .dynhr_abort(".sv_chan_loglik: the SV shock must be uncorrelated with the ",
                 "other shocks (a cross term e^{h/2} breaks the a + b e^h form).",
                 class = "dynhr_sv_chan_not_applicable")
  path <- .sv_revealed_path(Y, TT, ZZ, RR, DD, Sigma_e, d, P0, sv_idx, hyper,
                            n_exo, me_diag = me_diag)
  n_T <- ncol(Y)
  mu <- hyper[1L, "mu"]; rho <- hyper[1L, "rho"]; om <- hyper[1L, "sigma_eta"]^2
  om1 <- om / (1 - rho^2)

  ## F_t(u) = a_t + b e^u with b from the SV shock, a_t from P_t, ME and the
  ## non-SV shocks.
  b <- DD[1L, k]^2 * Sigma_e[k, k]
  if (!(b > 0))
    .dynhr_abort(".sv_chan_loglik: the SV shock does not load on the ",
                 "observable.", class = "dynhr_sv_chan_not_applicable")
  a <- vapply(seq_len(n_T), function(t) {
    Dn <- DD[1L, -k, drop = FALSE]
    path$A[[t]][1L, 1L] +
      if (ncol(Dn)) drop(Dn %*% Sigma_e[-k, -k, drop = FALSE] %*% t(Dn)) else 0
  }, numeric(1))
  v <- path$v[1L, ]
  if (any(a < 0))
    .dynhr_abort(".sv_chan_loglik: negative volatility-free variance.",
                 class = "dynhr_sv_chan_not_applicable")
  bad <- which(a > 0 & v^2 > a)
  if (length(bad))
    .dynhr_abort(sprintf(paste0(
      ".sv_chan_loglik: observation log-density is not concave in h at ",
      "period(s) %s (innovation^2 > volatility-free variance), violating ",
      "Chan's Assumption 1(ii); the tangent bound would not hold."),
      paste(bad, collapse = ", ")), class = "dynhr_sv_chan_not_logconcave")
  obs <- function(t, u) .sv_chan_obs_logdens(u, a[t], b, v[t])
  amean <- function(u) mu + rho * (u - mu)

  ## ---- Step 1: posterior mode m and Gaussian-approximation sd s ----------
  Q <- matrix(0, n_T, n_T)
  if (n_T == 1L) Q[1L, 1L] <- 1 / om1 else {
    diag(Q) <- c(1, rep(1 + rho^2, n_T - 2L), 1) / om
    for (t in seq_len(n_T - 1L)) Q[t, t + 1L] <- Q[t + 1L, t] <- -rho / om
  }
  obj <- function(h) sum(vapply(seq_len(n_T), function(t) obs(t, h[t])$l,
                                numeric(1))) -
    0.5 * drop(crossprod(h - mu, Q %*% (h - mu)))
  h <- rep(mu, n_T); f_h <- obj(h)
  for (it in seq_len(200L)) {
    od <- lapply(seq_len(n_T), function(t) obs(t, h[t]))
    gr <- vapply(od, `[[`, numeric(1), "d1") - drop(Q %*% (h - mu))
    H <- diag(vapply(od, `[[`, numeric(1), "d2"), n_T) - Q
    step <- -solve(H, gr)
    lam <- 1
    repeat {
      h_new <- h + lam * step; f_new <- obj(h_new)
      if (is.finite(f_new) && f_new >= f_h - 1e-12 * abs(f_h)) break
      lam <- lam / 2
      if (lam < 1e-8) break
    }
    h <- h_new; f_h <- f_new
    if (max(abs(lam * step)) < 1e-10) break
  }
  od <- lapply(seq_len(n_T), function(t) obs(t, h[t]))
  H <- diag(vapply(od, `[[`, numeric(1), "d2"), n_T) - Q
  s_mode <- sqrt(pmax(diag(solve(-H)), .Machine$double.eps))

  ## ---- Step 2: backward pass (tangent twists) ------------------------------
  if (is.null(n_nodes)) n_nodes <- max(10L, ceiling(4 * sqrt(n_T)))
  wgrid <- seq(node_range[1L], node_range[2L], length.out = n_nodes)
  nodes <- lapply(seq_len(n_T), function(t) h[t] + wgrid * s_mode[t])
  env <- vector("list", n_T)
  ot <- obs(n_T, nodes[[n_T]]); f_t <- ot$l; fp_t <- ot$d1
  for (t in rev(seq_len(n_T))) {
    env[[t]] <- .sv_chan_envelope(nodes[[t]], f_t, fp_t)
    if (t == 1L) break
    u <- nodes[[t - 1L]]
    cw <- .sv_chan_weights(env[[t]], amean(u), om)
    pw <- exp(cw$lw - cw$logC)
    ## Truncated-normal means mu_tj = b + sqrt(om) (phi(lo) - phi(hi)) / Z.
    lz <- cw$lw - matrix(env[[t]]$al + 0.5 * env[[t]]$g^2 * om, length(u),
                         length(env[[t]]$g), byrow = TRUE) -
      outer(amean(u), env[[t]]$g)
    tn <- cw$b + sqrt(om) * (exp(stats::dnorm(cw$lo, log = TRUE) - lz) -
                               exp(stats::dnorm(cw$hi, log = TRUE) - lz))
    tn[!is.finite(tn)] <- cw$b[!is.finite(tn)]
    Mt <- rowSums(pw * tn)
    op <- obs(t - 1L, u)
    f_t  <- op$l + cw$logC
    fp_t <- op$d1 + rho * (Mt - amean(u)) / om
  }
  c1 <- .sv_chan_weights(env[[1L]], mu, om1)
  log_C1 <- c1$logC

  ## ---- Step 3: forward proposals and log weights --------------------------
  M <- as.integer(n_draws)
  psi <- function(e, x) apply(outer(x, e$g) +
                                matrix(e$al, length(x), length(e$g),
                                       byrow = TRUE), 1L, min)
  draw <- function(cw) {
    J <- ncol(cw$lw); n <- nrow(cw$lw)
    pw <- exp(cw$lw - cw$logC)
    cp <- pw %*% upper.tri(diag(J), diag = TRUE)     # row-wise cumsum
    j <- pmin(1L + rowSums(cp < stats::runif(n)), J)
    idx <- cbind(seq_len(n), j)
    list(x = cw$b[idx] + sqrt(cw$om) * .sv_chan_rtrunc(cw$lo[idx], cw$hi[idx]))
  }
  c1r <- .sv_chan_weights(env[[1L]], rep(mu, M), om1); c1r$om <- om1
  hp <- draw(c1r)$x
  L <- numeric(M)
  if (n_T >= 2L) for (t in 2:n_T) {
    cw <- .sv_chan_weights(env[[t]], amean(hp), om); cw$om <- om
    L <- L + obs(t - 1L, hp)$l + cw$logC - psi(env[[t - 1L]], hp)
    hp <- draw(cw)$x
  }
  L <- L + obs(n_T, hp)$l - psi(env[[n_T]], hp)

  log_w <- log_C1 + L
  mx <- max(log_w)
  list(loglik = mx + log(mean(exp(log_w - mx))),
       log_C1 = log_C1,
       accept_rate = mean(exp(pmin(L, 0))),
       max_resid = max(L),
       log_w = log_w)
}
