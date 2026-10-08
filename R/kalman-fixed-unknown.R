## R/kalman-fixed-unknown.R
## ---------------------------------------------------------------------------
## lik_init = "fixed_unknown": the augmented ("fixed unknown initial state")
## Kalman filter of de Jong (1991) / Rosenberg (1973) -- IRIS's default
## initialisation (`initMeanUnit = 'optimal'`; IRIS 2018 @model/mykalman.m,
## +kalman/init.m, ped.m, correct.m, oolik.m).
##
## Lagged timing, as everywhere in kalman-filter.R:
##   s_t = T s_{t-1} + R e_t,     y_t = d + Z s_{t-1} + D e_t (+ u_t).
## The unit-root subspace of T gets a FIXED UNKNOWN level instead of a diffuse
## prior:
##   s_0 = a0 + A_0 delta,  Var(s_0) = P_0 (zero along the unit directions)
## The filter runs ONCE with that proper covariance, carrying
## A_t = d s_{t|t} / d delta:
##   A_t = T A_{t-1} - K_t Z_o A_{t-1}   (observed rows Z_o)
##   A_t = T A_{t-1}                     (all-missing period)
## and accumulates S = sum M' F^-1 M, s = sum M' F^-1 v with M = Z_o A_{t-1}
## (the innovation's response to delta is -M). Then
##   delta = pinv(S) s                (GLS; unidentified directions get 0)
##   s_{t|t}   += A_t delta,   s_{t|t-1} += T A_{t-1} delta
##   loglik     = ll0 + 0.5 delta' s  (concentrated over delta)
## The concentrated loglik carries NO -0.5 log det S term: it is IRIS's
## loglik(..., 'relative=', false). (IRIS's default relative = true also
## profiles a common scale factor out of the covariances; not reproduced.)
##
## COORDINATES (0.9.4.64). "Zero variance along the unit directions, Lyapunov
## on the ORTHOGONAL complement" is not invariant to the choice of state
## coordinates: the orthogonal complement depends on them, and a random
## component along the unit directions cannot be absorbed by the fixed delta.
## (The diffuse prior and the smoothed states ARE invariant.) IRIS builds the
## init on ITS state vector xb: every transition variable that is not
## forward-looking -- static ones (identities such as l = n + p + er)
## included -- plus the lags; measurement variables are not in it. dynhr's
## state vector omits static variables, so the init is computed in IRIS's
## coordinates X and mapped to dynhr's states s = C x (C a 0/1 selector,
## dynhr's states being a subset of X):
##   T_X = ghx[X, ] C,  R_X = ghu[X, ],  (U_Xu, P_X*) = .kf_diffuse_P0(T_X, R_X Q R_X')
##   A_0 = C U_Xu,      P_0 = C P_X* C'
## (.fu_init_coordinates / .fu_init). On the potential-output model (11 dynhr
## states, 13 IRIS states) this moved the filtered states from 4e-2 to 4e-12
## of IRIS and the loglik from 0.1 to 1e-9 nats.
##
## SINGULAR F. An observable that is an exact identity of others (data
## consistent) makes F singular. The multivariate recursion then hands over
## to a SEQUENTIAL (univariate, Koopman-Durbin) form of the same augmented
## filter on x_t = [s_{t-1}; e_t]: one component at a time, a component whose
## conditional variance is <= kalman_tol is DROPPED (the univariate filter's
## rule), and A_t rides every update. A dropped component is counted
## informative when its CORRECTED innovation v - M delta is not negligible
## (.kf_informative_skip on that). With no component dropped the sequential
## recursion is the multivariate one to round-off.
##
## The smoothed states of this filter are those of the exact diffuse smoother
## (the diffuse prior is the flat-prior limit of the GLS estimate of delta),
## which is how kalman_smoother(lik_init = "fixed_unknown") is implemented.
## ---------------------------------------------------------------------------

## Pseudo-inverse of a symmetric PSD matrix (eigen; MASS::ginv's relative
## tolerance). Directions with no information get a zero estimate, as IRIS's
## pinv does.
.fu_pinv_sym <- function(S, tol = sqrt(.Machine$double.eps)) {
  k <- nrow(S)
  if (k == 0L) return(list(inv = matrix(0, 0, 0), rank = 0L))
  e <- eigen((S + t(S)) / 2, symmetric = TRUE)
  ev <- e$values
  keep <- ev > tol * max(ev[1L], 0)
  if (!any(keep)) return(list(inv = matrix(0, k, k), rank = 0L))
  V <- e$vectors[, keep, drop = FALSE]
  list(inv = V %*% (t(V) / ev[keep]), rank = sum(keep))
}

## Prefixes of the auxiliary variables that carry a LEAD (R/aux-expansion-
## monolith.R, R/parse-mod.R): forward-looking by construction.
.FU_LEAD_AUX_PREFIX <- c("AUX_LEAD_", "AUX_EXO_LEAD_", "AUX_EXPECT_LEAD_")

## Which endogenous variables are forward-looking (appear with a lead)?
## From the lead-lag incidence (its t+1 row), else the variable
## classification; NULL when the model carries neither.
.fu_forward_looking <- function(model, endo_names) {
  lli <- model$lead_lag_incidence
  if (is.matrix(lli) && ncol(lli) == length(endo_names)) {
    if (!is.null(colnames(lli))) {
      ix <- match(endo_names, colnames(lli))
      if (anyNA(ix)) return(NULL)
      lli <- lli[, ix, drop = FALSE]
    }
    rn <- rownames(lli)
    lead_row <- if (!is.null(rn) && "t+1" %in% rn) match("t+1", rn)
                else if (nrow(lli) == 3L) 3L
                else NA_integer_
    if (!is.na(lead_row)) return(as.numeric(lli[lead_row, ]) > 0)
  }
  vc <- model$variable_classification
  if (is.list(vc) && !is.data.frame(vc) &&
      any(c("predetermined", "forward", "mixed", "static") %in% names(vc)))
    return(endo_names %in% c(vc$forward, vc$mixed))
  NULL
}

## IRIS's state coordinates for the fixed-unknown init: every endogenous
## variable that is not forward-looking (no lead; lead auxiliaries excluded)
## and is not a measurement variable, static ones INCLUDED, lag auxiliaries
## included; dynhr's own states are always in. A measurement variable is an
## observable (`obs_vars`, or the model's declared varobs) that is not itself
## a dynhr state -- an observable that IS a state (no separate measurement
## equation) is a transition variable in IRIS's sense and stays. Returns the
## names in declaration order, or NULL when the model does not say which
## variables are forward-looking (the init is then built on dynhr's states).
## @noRd
.fu_init_coordinates <- function(model, endo_names, state_names, obs_vars) {
  if (is.null(model) || !length(endo_names)) return(NULL)
  fwd <- .fu_forward_looking(model, endo_names)
  if (is.null(fwd)) return(NULL)
  lead_aux <- Reduce(`|`, lapply(.FU_LEAD_AUX_PREFIX,
                                 function(p) startsWith(endo_names, p)))
  meas <- setdiff(unique(c(obs_vars, model$varobs, model$varobs_names)),
                  state_names)
  keep <- (!fwd & !lead_aux & !(endo_names %in% meas)) |
    endo_names %in% state_names
  endo_names[keep]
}

## The fixed-unknown initialisation: the unit-root block A_0 (n_state x
## n_unit) and the proper covariance P_0 on the complement, built in the
## coordinates `coords` (see the header) when they can be, otherwise on
## dynhr's states. `state_names` are the filter's states in ghx-column order
## (the order of TT's rows and columns). The IRIS coordinates are used only
## when the selection reproduces TT exactly -- a transformed or augmented
## state space falls back to the state coordinates, and says so in
## `$coordinates`.
## @return list(U_unit, P_star, nunit, coordinates)
## @noRd
.fu_init <- function(TT, RR, Sigma_e, state_names, coords = NULL,
                     ghx = NULL, ghu = NULL, endo_names = NULL) {
  ok_X <- !is.null(coords) && !is.null(ghx) && !is.null(ghu) &&
    !is.null(endo_names) && !is.null(state_names) &&
    nrow(ghx) == length(endo_names) && nrow(ghu) == length(endo_names) &&
    ncol(ghx) == nrow(TT) && ncol(ghu) == ncol(RR) &&
    length(state_names) == nrow(TT) &&
    all(coords %in% endo_names) && all(state_names %in% coords)
  if (ok_X) {
    is_ <- match(state_names, endo_names)
    sc  <- max(1, abs(TT))
    ok_X <- max(abs(ghx[is_, , drop = FALSE] - TT)) <= 1e-12 * sc &&
      max(abs(ghu[is_, , drop = FALSE] - RR)) <= 1e-12 * max(1, abs(RR))
  }
  if (!ok_X) {
    dp <- .kf_diffuse_P0(TT, tcrossprod(RR %*% Sigma_e, RR))
    return(list(U_unit = dp$U_unit, P_star = dp$P_star, nunit = dp$nunit,
                coordinates = state_names))
  }
  ix   <- match(coords, endo_names)
  Csel <- matrix(0, length(state_names), length(coords))
  Csel[cbind(seq_along(state_names), match(state_names, coords))] <- 1
  TX <- ghx[ix, , drop = FALSE] %*% Csel
  RX <- ghu[ix, , drop = FALSE]
  dpx <- .kf_diffuse_P0(TX, tcrossprod(RX %*% Sigma_e, RX))
  list(U_unit = Csel %*% dpx$U_unit,
       P_star = .sym(Csel %*% dpx$P_star %*% t(Csel)),
       nunit = dpx$nunit, coordinates = coords)
}

## The fixed-unknown filter on DEVIATION data.
##
## @param Y        n_obs x T, deviations from the observation intercept (and
##                 any trend / deterministic path already taken off); NA =
##                 missing.
## @param TT,RR,ZZ,DD  lagged-form system matrices.
## @param Sigma_e  baseline shock covariance.
## @param a0       length n_state, the KNOWN part of s_0 (deviations).
## @param me_mat   n_obs x T measurement-error variances (base + extra).
## @param shock_scale n_exo x T std scale factors, or NULL.
## @param init     .fu_init(): the unit-root block (U_unit) and the proper
##                 part (P_star).
## @param kalman_tol drop threshold of the sequential fallback (the
##                 univariate filter's).
## @return list(ok, loglik, loglik0, ll_contrib, filtered, predicted, s0,
##   final_state, final_cov, delta, init_unit_estimate, n_unit, rank_S,
##   n_obs_used, recursion, dropped, dropped_informative). `filtered` /
##   `predicted` are n_state x T and CORRECTED; `recursion` is
##   "multivariate", or "sequential" when a singular F handed over to the
##   component-by-component form; `ok = FALSE` (loglik = -Inf) on a
##   non-finite step.
## @noRd
.kf_fixed_unknown_core <- function(Y, TT, RR, ZZ, DD, Sigma_e, a0, me_mat,
                                   shock_scale = NULL, init,
                                   return_filtered = TRUE,
                                   kalman_tol = 1e-10) {
  run <- .fu_multivariate(Y, TT, RR, ZZ, DD, Sigma_e, a0, me_mat,
                          shock_scale, init, return_filtered)
  if (isTRUE(run$singular))
    run <- .fu_sequential(Y, TT, RR, ZZ, DD, Sigma_e, a0, me_mat,
                          shock_scale, init, return_filtered, kalman_tol)
  if (!isTRUE(run$ok)) return(list(ok = FALSE, loglik = -Inf))
  .fu_finish(run, TT, a0, init, kalman_tol, return_filtered)
}

## Multivariate per-step recursion. Returns the uncorrected pieces, or
## list(singular = TRUE) at the first singular innovation covariance.
.fu_multivariate <- function(Y, TT, RR, ZZ, DD, Sigma_e, a0, me_mat,
                             shock_scale, init, return_filtered) {
  n   <- nrow(TT)
  n_T <- ncol(Y)
  Uu  <- init$U_unit
  k   <- ncol(Uu)

  a <- as.numeric(a0)
  P <- init$P_star
  A <- Uu
  c0 <- numeric(n_T)                       # -0.5 (n_t log 2 pi + log det F_t + v'F^-1 v)
  b_list <- vector("list", n_T)            # M' F^-1 v
  C_list <- vector("list", n_T)            # M' F^-1 M
  S_mat <- matrix(0, k, k); s_vec <- numeric(k)
  filt <- if (return_filtered) matrix(0, n, n_T) else NULL
  Af   <- if (return_filtered) array(0, c(n, k, n_T)) else NULL
  n_used <- 0L

  for (t in seq_len(n_T)) {
    Q_t <- if (is.null(shock_scale)) Sigma_e
           else Sigma_e * outer(shock_scale[, t], shock_scale[, t])
    TA <- TT %*% A
    o <- which(is.finite(Y[, t]))
    if (!length(o)) {
      a <- as.numeric(TT %*% a)
      P <- tcrossprod(TT %*% P, TT) + tcrossprod(RR %*% Q_t, RR)
      P <- (P + t(P)) * 0.5
      A <- TA
    } else {
      Z <- ZZ[o, , drop = FALSE]; D <- DD[o, , drop = FALSE]
      h <- me_mat[o, t]
      v <- Y[o, t] - as.numeric(Z %*% a)
      Ft <- Z %*% P %*% t(Z) + tcrossprod(D %*% Q_t, D) +
        diag(h, nrow = length(o))
      Ft <- (Ft + t(Ft)) * 0.5
      if (!all(is.finite(Ft))) return(list(ok = FALSE))
      ## Pivoted Cholesky does not stop() on a singular matrix: the rank says
      ## so, and .kf_F_singular() applies the package's scale-invariant rule.
      Fc <- suppressWarnings(chol(Ft, pivot = TRUE))
      if (attr(Fc, "rank") < length(o)) return(list(singular = TRUE))
      Fi_p <- chol2inv(Fc)
      if (.kf_F_singular(Fc, Ft, Fi = Fi_p)) return(list(singular = TRUE))
      ip <- order(attr(Fc, "pivot"))
      Fi <- Fi_p[ip, ip, drop = FALSE]
      ldf <- 2 * sum(log(diag(Fc)))
      G  <- TT %*% P %*% t(Z) + RR %*% Q_t %*% t(D)
      K  <- G %*% Fi
      M  <- Z %*% A
      Fv <- as.numeric(Fi %*% v)
      FM <- Fi %*% M
      c0[t] <- -0.5 * (length(o) * log(2 * pi) + ldf + sum(v * Fv))
      b_list[[t]] <- as.numeric(crossprod(M, Fv))
      C_list[[t]] <- crossprod(M, FM)
      S_mat <- S_mat + C_list[[t]]
      s_vec <- s_vec + b_list[[t]]
      n_used <- n_used + length(o)
      a <- as.numeric(TT %*% a) + as.numeric(K %*% v)
      A <- TA - K %*% M
      ## Joseph form, as .kf_step: exact for any gain, symmetric by construction.
      TmKZ <- TT - K %*% Z; RmKD <- RR - K %*% D
      P <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Q_t, RmKD)
      if (any(h != 0)) P <- P + K %*% (h * t(K))
      P <- (P + t(P)) * 0.5
    }
    if (return_filtered) { filt[, t] <- a; Af[, , t] <- A }
  }
  list(ok = TRUE, recursion = "multivariate", c0 = c0, b_list = b_list,
       C_list = C_list, S = S_mat, s = s_vec, filt = filt, Af = Af,
       a = a, A = A, P = P, n_used = n_used, skipped = NULL)
}

## Sequential (univariate) form of the same augmented recursion, on
## x_t = [s_{t-1}; e_t] with y_t = [Z D] x_t and x_{t+1} = [[T, R], [0, 0]] x_t
## + [0; e_{t+1}] (the .kf_univariate_loop_R layout). A component whose
## conditional variance is <= kalman_tol is dropped and recorded; whether it
## was informative is decided on its CORRECTED innovation in .fu_finish.
.fu_sequential <- function(Y, TT, RR, ZZ, DD, Sigma_e, a0, me_mat,
                           shock_scale, init, return_filtered, kalman_tol) {
  n    <- nrow(TT)
  n_T  <- ncol(Y)
  n_e  <- ncol(RR)
  nb   <- n + n_e
  s_ix <- seq_len(n); e_ix <- n + seq_len(n_e)
  Uu   <- init$U_unit
  k    <- ncol(Uu)
  log2pi <- log(2 * pi)
  Q_at <- function(t) if (is.null(shock_scale)) Sigma_e
                      else Sigma_e * outer(shock_scale[, t], shock_scale[, t])
  Zb <- cbind(ZZ, DD)
  Tb <- rbind(cbind(TT, RR), matrix(0, n_e, nb))

  a <- c(as.numeric(a0), numeric(n_e))
  P <- matrix(0, nb, nb)
  P[s_ix, s_ix] <- init$P_star
  P[e_ix, e_ix] <- Q_at(1L)
  A <- rbind(Uu, matrix(0, n_e, k))
  c0 <- numeric(n_T)
  b_list <- vector("list", n_T)
  C_list <- vector("list", n_T)
  S_mat <- matrix(0, k, k); s_vec <- numeric(k)
  filt <- if (return_filtered) matrix(0, n, n_T) else NULL
  Af   <- if (return_filtered) array(0, c(n, k, n_T)) else NULL
  n_used <- 0L
  sk_t <- integer(0); sk_v <- numeric(0); sk_y <- numeric(0)
  sk_M <- matrix(0, k, 0)

  for (t in seq_len(n_T)) {
    ll_t <- 0
    b_t <- numeric(k); C_t <- matrix(0, k, k)
    for (i in which(is.finite(Y[, t]))) {
      Zi <- Zb[i, ]
      y_i <- Y[i, t]
      v  <- y_i - sum(Zi * a)
      Mi <- as.numeric(crossprod(A, Zi))
      Ki <- as.numeric(P %*% Zi)
      Fi <- sum(Zi * Ki) + me_mat[i, t]
      if (!is.finite(Fi) || !is.finite(v)) return(list(ok = FALSE))
      if (Fi > kalman_tol) {
        ll_t <- ll_t - 0.5 * (log2pi + log(Fi) + v * v / Fi)
        b_t  <- b_t + Mi * (v / Fi)
        C_t  <- C_t + tcrossprod(Mi) / Fi
        a <- a + Ki * (v / Fi)
        A <- A - tcrossprod(Ki, Mi) / Fi
        P <- P - tcrossprod(Ki) / Fi
        n_used <- n_used + 1L
      } else {
        sk_t <- c(sk_t, t); sk_v <- c(sk_v, v); sk_y <- c(sk_y, y_i)
        sk_M <- cbind(sk_M, Mi)
      }
    }
    c0[t] <- ll_t
    b_list[[t]] <- b_t; C_list[[t]] <- C_t
    S_mat <- S_mat + C_t; s_vec <- s_vec + b_t
    a <- as.numeric(Tb %*% a)
    A <- Tb %*% A
    P <- Tb %*% P %*% t(Tb)
    P[e_ix, e_ix] <- Q_at(min(t + 1L, n_T))
    P <- (P + t(P)) * 0.5
    if (return_filtered) {
      filt[, t] <- a[s_ix]
      Af[, , t] <- A[s_ix, , drop = FALSE]
    }
  }
  list(ok = TRUE, recursion = "sequential", c0 = c0, b_list = b_list,
       C_list = C_list, S = S_mat, s = s_vec, filt = filt, Af = Af,
       a = a[s_ix], A = A[s_ix, , drop = FALSE], P = P[s_ix, s_ix, drop = FALSE],
       n_used = n_used,
       skipped = list(t = sk_t, v = sk_v, y = sk_y, M = sk_M))
}

## GLS estimate of delta, the concentrated loglik and the corrected paths.
.fu_finish <- function(run, TT, a0, init, kalman_tol, return_filtered) {
  n_T <- length(run$c0)
  n   <- nrow(TT)
  k   <- ncol(init$U_unit)
  pv    <- .fu_pinv_sym(run$S)
  delta <- as.numeric(pv$inv %*% run$s)
  ll0   <- sum(run$c0)
  loglik <- ll0 + 0.5 * sum(delta * run$s)
  if (!is.finite(loglik)) return(list(ok = FALSE, loglik = -Inf))

  ## Contributions of the CORRECTED innovations v_t - M_t delta. They sum to
  ## the concentrated loglik exactly: with delta = S^+ s (s in range(S)),
  ## sum (v - M delta)' F^-1 (v - M delta) = sum v'F^-1 v - delta's.
  ll_contrib <- run$c0
  for (t in seq_len(n_T)) if (!is.null(run$b_list[[t]]))
    ll_contrib[t] <- run$c0[t] + sum(delta * run$b_list[[t]]) -
      0.5 * sum(delta * (run$C_list[[t]] %*% delta))

  ## Dropped components (sequential form only): informative when the
  ## corrected innovation is not negligible (.kf_informative_skip).
  dropped <- dropped_inf <- integer(n_T)
  sk <- run$skipped
  if (length(sk$t)) {
    vc <- sk$v - as.numeric(crossprod(sk$M, delta))
    for (j in seq_along(sk$t)) {
      dropped[sk$t[j]] <- dropped[sk$t[j]] + 1L
      if (.kf_informative_skip(vc[j], sk$y[j], kalman_tol))
        dropped_inf[sk$t[j]] <- dropped_inf[sk$t[j]] + 1L
    }
  }

  unit_est <- as.numeric(init$U_unit %*% delta)
  s0 <- as.numeric(a0) + unit_est
  filt <- pred <- NULL
  if (return_filtered) {
    filt <- run$filt
    if (k > 0L) for (t in seq_len(n_T))
      filt[, t] <- filt[, t] + matrix(run$Af[, , t], n, k) %*% delta
    ## s_{t|t-1} = T s_{t-1|t-1} (E e_t = 0), corrected path included.
    pred <- TT %*% cbind(s0, filt[, -n_T, drop = FALSE])
    dimnames(pred) <- NULL
  }
  list(ok = TRUE, loglik = loglik, loglik0 = ll0, ll_contrib = ll_contrib,
       filtered = filt, predicted = pred, s0 = s0,
       final_state = run$a + as.numeric(run$A %*% delta),
       ## MSE of s_{T|T} with delta estimated: P_T + A_T S^+ A_T' (de Jong
       ## 1991) -- the exact diffuse filter's P_T once delta is identified.
       final_cov = .sym(run$P + run$A %*% pv$inv %*% t(run$A)),
       delta = delta, init_unit_estimate = unit_est, n_unit = k,
       rank_S = pv$rank, n_obs_used = run$n_used,
       recursion = run$recursion,
       dropped = dropped, dropped_informative = dropped_inf)
}
