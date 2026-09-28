## R/obc-filter.R
## --------------------------------------------------------------------------
## OBC piecewise-linear Kalman filter (Dynare's OccBin PKF).
##
## Provides:
##   kalman_filter_obc_pkf() -- piecewise-linear Kalman filter (exported):
##                              Giovannini, Pfeiffer & Ratto (2021), the
##                              algorithm of Dynare 7.1's OccBin likelihood
##                              (occbin/+occbin/kalman_update_engine.m and
##                              kalman_update_algo_1.m)
##   kalman_filter_obc()     -- Kalman filter along a GIVEN regime path
##   pkf_extract_shock()     -- eps_{t|t} from KF quantities (GPR eq. 6)
##   pkf_backward_one_step() -- s_{t-1|t} (GPR eq. 7, one step)
##   pkf_check_binding()     -- slack-rule regime check at t only (no longer
##                              used by any filter: the PKF since W49, the
##                              particle filters since W50)
##
## STATE SPACE (dynhr's lag-1 timing).  In period t the model follows the
## TIME-VARYING piecewise-linear rule of the regime sequence expected at t
## (.obc_pwl_rules(), R/obc-binding.R):
##   y_t = ghx_t s_{t-1} + ghu_t eps_t + c_t
##   s_t = TT_t s_{t-1} + RR_t eps_t + c_state_t,   TT_t = ghx_t[state, ] ...
##   obs_t = ZZ_t s_{t-1} + DD_t eps_t + d + c_obs_t + me
## The rule of period t depends on the WHOLE expected sequence (a period that
## binds with k more binding periods ahead has a different rule from the last
## period of a spell), so the filter carries, per period, the regime sequence
## expected at t over a check-ahead horizon (Dynare: 200 periods).
##
## W49 (2026-09-25) replaced the pre-0.9.3.93 filters, which applied ONE
## policy per regime (next period assumed slack) and checked the regime of
## period t only against the SLACK rule -- wrong for every spell of 2+
## periods and for any period that anticipates a later binding period.
##
## kf_store completeness invariant (kept from 5fd8b64): for every t,
##   kf_store$L[[t]] and kf_store$RR[[t]] are matrices; NULL entries are
##   stored with single-bracket list(NULL) so the list keeps its length.
## --------------------------------------------------------------------------


# =============================================================================
# Kalman filter along a given regime path
# =============================================================================

#' Kalman filter for an OBC model along a GIVEN regime path
#'
#' Evaluates the Gaussian likelihood of the piecewise-linear state space
#' along \code{regime_path}.  The rule of period t is the time-varying
#' OccBin rule of the regime sequence expected at t:
#' \itemize{
#'   \item when \code{regime_cache} carries the per-period rules of a
#'     \code{kalman_filter_obc_pkf()} run on the same regime path (as
#'     \code{obc_guess_verify()} leaves it), those rules -- the regime
#'     sequences the filter EXPECTED at each t -- are used, and the
#'     log-likelihood is the PKF's (except in a period with several regime
#'     solutions, whose densities the PKF sums);
#'   \item otherwise agents are taken to foresee the given path: the
#'     expected sequence at t is \code{regime_path[t:T]} (slack after T).
#' }
#' Before W49 (0.9.3.93) every binding period used the one-period policy of
#' \code{obc_ensure_policy()} (next period slack).
#'
#' @param Y             Observation matrix (n_obs x T); columns = time periods
#' @param dr_slack      Slack-regime DecisionRules
#' @param regime_cache  R environment seeded by obc_ensure_policy() (and, for
#'                      the PKF rules, filled by kalman_filter_obc_pkf())
#' @param model         dynhr_mod
#' @param params        Named numeric parameter vector
#' @param obs_vars      Character vector of observed variable names
#' @param regime_path   Integer vector (length T): bitfield regime per period
#' @param me_variance   Scalar measurement error variance added to F
#' @param return_filtered Logical: return the n_state x T filtered states
#' @return List with $loglik, $filtered_states (or NULL), $n_obs, $n_T
#' @noRd
kalman_filter_obc <- function(Y, dr_slack, regime_cache,
                               model, params, obs_vars, regime_path,
                               me_variance = 1e-8, return_filtered = FALSE) {
  endo    <- dr_slack$endo_names
  n_obs   <- length(obs_vars)
  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    .dynhr_abort("kalman_filter_obc: observed variables not in model: ",
                 paste(obs_vars[is.na(obs_idx)], collapse = ", "),
                 class = "dynhr_error_obc_filter")

  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)
  regime_path <- as.integer(regime_path)
  if (length(regime_path) != n_T)
    .dynhr_abort(sprintf(paste0(
      "kalman_filter_obc: regime_path length (%d) != number of time ",
      "periods (%d)."), length(regime_path), n_T),
      class = "dynhr_error_obc_filter")

  rules <- .obc_path_rules(regime_cache, regime_path)

  Sigma_e <- .get_shock_cov(model, dr_slack$exo_names, params)
  pk <- .obc_pkf_prep(NULL, dr_slack, obs_idx, Sigma_e, me_variance, 1L)
  d  <- .obc_pkf_obs_ss(dr_slack, obs_vars)

  s <- numeric(pk$n_state)
  P <- .obc_pkf_P0(pk)
  loglik   <- 0
  filtered <- if (return_filtered) matrix(0, pk$n_state, n_T) else NULL
  for (t in seq_len(n_T)) {
    up <- .obc_pkf_kf_step(pk, .obc_pkf_mats(pk, rules[[t]]), s, P, Y[, t], d)
    if (!up$ok)
      return(list(loglik = -Inf, filtered_states = NULL, n_obs = n_obs,
                  n_T = n_T))
    loglik <- loglik + up$ll
    s <- up$s_new
    P <- up$P_new
    if (return_filtered) filtered[, t] <- s
  }
  if (return_filtered) rownames(filtered) <- endo[dr_slack$state_idx]
  list(loglik = loglik, filtered_states = filtered, n_obs = n_obs, n_T = n_T)
}

#' Per-period rules along a regime path
#'
#' The PKF's rules when \code{regime_cache} holds those of a
#' kalman_filter_obc_pkf() run with the same regime path (the regime
#' sequences the filter EXPECTED in each period); otherwise perfect foresight
#' of the given path: the rule of period t is that of the expected sequence
#' \code{regime_path[t:T]} (slack after T).
#' @return List (length T): NULL (slack rule) or list(ghx, ghu, c)
#' @noRd
.obc_path_rules <- function(regime_cache, regime_path) {
  regime_path <- as.integer(regime_path)
  n_T <- length(regime_path)
  if (exists(".pkf_rules", envir = regime_cache, inherits = FALSE)) {
    pr <- get(".pkf_rules", envir = regime_cache, inherits = FALSE)
    if (identical(as.integer(pr$regime_path), regime_path))
      return(pr$rules)
  }
  rules <- vector("list", n_T)
  if (all(regime_path == 0L)) return(rules)
  ctx <- .obc_pwl_cache_context(regime_cache)
  for (t in seq_len(n_T)) {
    if (all(regime_path[t:n_T] == 0L)) break
    rules[t] <- list(.obc_pwl_rules(ctx, regime_path[t:n_T])[[1L]])
  }
  rules
}

# =============================================================================
# PKF helpers — Giovannini, Pfeiffer, Ratto (2021) §2.2
#
# Eqs. 6-7 (shock extraction and one-step backward smoother) as standalone
# helpers; .obc_pkf_kf_step() computes the same two quantities inline.
# pkf_check_binding() (period t against the SLACK rule) is no longer used by
# kalman_filter_obc_pkf(), which solves the whole expected regime sequence
# (W49), nor by the particle filters in R/obc-ppf.R, which solve it per
# particle with .obc_pkf_solve()'s algorithm (W50).
#
# The "inversion" label: eq. 6 inverts the innovation to recover the
# underlying shock eps_{t|t}, which is then used (together with the
# backward-smoothed state) to verify the OBC regime — rather than checking
# the slack prediction with zero shocks as obc_should_bind() does.
# =============================================================================

#' Extract the period-t shock estimate from Kalman filter quantities
#'
#' Implements GPR (2021) eq. 6:
#'   eps_{t|t} = Sigma_e * t(DD) * F_inv * v
#'
#' With dynhr's correlated-noise timing (y_t = ZZ*s_{t-1} + DD*eps_t),
#' DD encodes H*R in paper notation, so Cov(eps_t, v_t) = Sigma_e * t(DD).
#' Handles partial-observation subsetting: pass DD[obs_ok, , drop=FALSE].
#'
#' @param Sigma_e n_exo x n_exo shock covariance matrix
#' @param DD      n_obs_t x n_exo observation-shock loading (ghu[obs_ok, ])
#' @param F_inv   n_obs_t x n_obs_t inverse innovation covariance
#' @param v       length-n_obs_t innovation vector
#' @return length-n_exo vector eps_{t|t}
#' @noRd
pkf_extract_shock <- function(Sigma_e, DD, F_inv, v) {
  drop(Sigma_e %*% t(DD) %*% (F_inv %*% v))
}


#' One-step backward smoother: update s_{t-1|t-1} to s_{t-1|t}
#'
#' Implements the first step of the backward recursion in GPR (2021) eq. 7,
#' initialised with r_{t+1} = 0 so only period-t innovations contribute:
#'
#'   r_t       = t(ZZ) * F_inv * v          (r_{t+1} = 0 term vanishes)
#'   s_{t-1|t} = s_{t-1|t-1} + P_{t-1|t-1} * r_t
#'
#' Equivalently: s_{t-1|t} = s + Cov(s_{t-1}, v_t) * F^{-1} * v
#' where Cov(s_{t-1}, v_t) = P * t(ZZ) (the cross-covariance under the
#' lag-1 observation convention y_t = ZZ*s_{t-1} + DD*eps_t).
#'
#' @param s_prev  length-n_state filtered state s_{t-1|t-1}
#' @param P_prev  n_state x n_state filtered covariance P_{t-1|t-1}
#' @param ZZ      n_obs_t x n_state observation-state loading (ghx[obs_ok, ])
#' @param F_inv   n_obs_t x n_obs_t inverse innovation covariance
#' @param v       length-n_obs_t innovation vector
#' @return length-n_state backward-smoothed state s_{t-1|t}
#' @noRd
pkf_backward_one_step <- function(s_prev, P_prev, ZZ, F_inv, v) {
  drop(s_prev + P_prev %*% (t(ZZ) %*% (F_inv %*% v)))
}


#' Check which OBC constraints bind given backward-smoothed state and shock
#'
#' Predicts the constrained variable at t using the SLACK policy applied to
#' (s_backward, eps_hat):
#'
#'   var_t = ghx_slack[var_idx, ] * s_backward + ghu_slack[var_idx, ] * eps_hat
#'
#' This is the shock-aware replacement for obc_should_bind(), which uses
#' the lagged filtered state with zero shock instead.  Because eps_hat is
#' extracted from the actual observations (via pkf_extract_shock), the regime
#' check correctly reflects the driving force behind any constraint violation.
#'
#' @param s_backward length-n_state backward-smoothed state s_{t-1|t}
#' @param eps_hat    length-n_exo extracted shock eps_{t|t}
#' @param specs      OBC spec list from obc_parse_tags
#' @param dr_slack   Slack-regime DecisionRules
#' @return logical vector (length k): TRUE where constraint j binds
#' @noRd
pkf_check_binding <- function(s_backward, eps_hat, specs, dr_slack) {
  # The MCP bound is a LEVEL (package convention, 0.9.4 ledger A6) while
  # ghx/ghu produce DEVIATIONS.  Convert the bound once, through the SAME
  # shared helper every other regime/binding check uses
  # (obc_solve_binding, obc-boehl, obc-regime, obc-lcp), so the
  # particle/Kalman filter path cannot desynchronise from the solvers:
  #   var_t >= b   <=>   (var_t - var_ss) >= b - var_ss
  # .obc_bound_dev() also handles a `dynhr_steady` $ys, a NULL/short $ys and
  # a non-finite steady-state entry, all of which the old inline
  # `var_dev + dr_slack$ys[s$var_idx]` silently got wrong.
  bnd <- .obc_bound_dev(specs, dr_slack)
  vapply(seq_along(specs), function(j) {
    s <- specs[[j]]
    var_dev <- sum(dr_slack$ghx[s$var_idx, ] * s_backward) +
               sum(dr_slack$ghu[s$var_idx, ] * eps_hat)
    if (s$op == ">") var_dev < bnd[[j]] else var_dev > bnd[[j]]
  }, logical(1))
}



# =============================================================================
# Piecewise-linear Kalman filter: building blocks
# =============================================================================

#' Filter context: slack matrices, engine context, slack-tail check matrix
#'
#' @param ctx      .obc_pwl_context() (NULL: no regime solves, e.g. for
#'                 kalman_filter_obc along a given path)
#' @param dr_slack Slack-regime DecisionRules
#' @param obs_idx  Observable positions in the endogenous vector
#' @param Sigma_e  Shock covariance
#' @param me       Measurement-error variance (scalar)
#' @param horizon  Check-ahead horizon of the regime solves
#' @noRd
.obc_pkf_prep <- function(ctx, dr_slack, obs_idx, Sigma_e, me, horizon) {
  ## Every OBC Kalman/particle filter builds its measurement noise here as
  ## me * I: a per-observable vector is refused with a classed error rather
  ## than recycled (kalman_filter() is the path that takes H = diag(me)).
  me <- .kf_me_variance(me, character(length(obs_idx)),
                        "OBC filter (kalman_filter_obc / kalman_filter_obc_pkf / ppf_likelihood)",
                        allow_vector = FALSE)
  si <- dr_slack$state_idx
  pk <- list(ctx = ctx, si = si, obs_idx = as.integer(obs_idx),
             n_state = length(si), n_exo = ncol(dr_slack$ghu),
             n_endo = nrow(dr_slack$ghx), Sigma_e = Sigma_e, me = me,
             H = as.integer(horizon), ghx = dr_slack$ghx, ghu = dr_slack$ghu)
  pk$slack <- .obc_pkf_mats(pk, NULL)
  ## Slack tail: rows of ghx[var, ] %*% TT_s^(k-2), k = 2..H (period-major,
  ## spec-minor), so the constrained variables of an all-slack path started
  ## at s_1 are one matrix product (the common case of the regime solve).
  if (!is.null(ctx) && ctx$n_spec > 0L && pk$H > 1L) {
    TTs <- pk$slack$TT
    G   <- dr_slack$ghx[ctx$var, , drop = FALSE]
    M   <- matrix(0, ctx$n_spec * (pk$H - 1L), pk$n_state)
    for (k in seq_len(pk$H - 1L)) {
      M[(k - 1L) * ctx$n_spec + seq_len(ctx$n_spec), ] <- G
      G <- G %*% TTs
    }
    pk$tail <- M
  }
  pk
}

#' State-space matrices of one period rule (NULL rule = slack)
#' @noRd
.obc_pkf_mats <- function(pk, rule) {
  if (is.null(rule)) {
    ghx <- pk$ghx; ghu <- pk$ghu; cc <- numeric(pk$n_endo)
  } else {
    ghx <- rule$ghx; ghu <- rule$ghu; cc <- as.numeric(rule$c)
  }
  si <- pk$si; oi <- pk$obs_idx
  list(ghx = ghx, ghu = ghu, c = cc,
       TT = ghx[si, , drop = FALSE], RR = ghu[si, , drop = FALSE],
       ZZ = ghx[oi, , drop = FALSE], DD = ghu[oi, , drop = FALSE],
       cs = cc[si], co = cc[oi])
}

#' Steady-state level of the observables (deviation-form filter offset)
#' @noRd
.obc_pkf_obs_ss <- function(dr_slack, obs_vars) {
  ys <- if (inherits(dr_slack$ys, "dynhr_steady")) dr_slack$ys$values
        else dr_slack$ys
  d <- as.numeric(ys[obs_vars])
  d[!is.finite(d)] <- 0
  d
}

#' Initial state covariance: the slack regime's stationary covariance
#' (Dynare lik_init = 1)
#' @noRd
.obc_pkf_P0 <- function(pk) {
  m <- pk$slack
  P <- solve_lyapunov(m$TT, tcrossprod(m$RR %*% pk$Sigma_e, m$RR))
  if (any(!is.finite(P))) P <- diag(1e6, pk$n_state)
  P
}

#' One Kalman step of period t under the rule matrices m
#'
#' Prediction with the period-t rule, update with the period-t data (missing
#' entries dropped), the shock estimate eps_{t|t} = Sigma_e DD' F^-1 v and the
#' one-step smoothed lagged state s_{t-1|t} = s_sm + P ZZ' F^-1 v that
#' Dynare's kalman_update_algo_1 feeds to the regime solve.  s_sm is the
#' previous period's Kalman-updated state; it differs from s (the state the
#' prediction starts from) only after a period whose regime iteration cycled,
#' where Dynare carries the solver's path forward as the state but smooths
#' from the Kalman update (kalman_update_algo_1: a = out.piecewise vs
#' alphahat = a1 + P1 r).  Accordingly s_new, the state carried forward, is
#' the solver's path TT s_{t-1|t} + RR eps_{t|t} + c_state (Dynare's
#' out.piecewise), and s_kf the Kalman update TT s + K v + c_state; the two
#' coincide when s = s_sm.
#' @return list(ok, ll, s_new, s_kf, P_new, eps, s_back, v, F_inv, L, ZZ, DD)
#' @noRd
.obc_pkf_kf_step <- function(pk, m, s, P, y, d, s_sm = s) {
  Sig <- pk$Sigma_e
  v_full <- y - drop(m$ZZ %*% s) - d - m$co
  ok_obs <- which(!is.na(v_full))
  if (length(ok_obs) == 0L) {
    P_new <- tcrossprod(m$TT %*% P, m$TT) + tcrossprod(m$RR %*% Sig, m$RR)
    return(list(ok = TRUE, ll = 0, s_new = drop(m$TT %*% s_sm) + m$cs,
                s_kf = drop(m$TT %*% s) + m$cs,
                P_new = (P_new + t(P_new)) * 0.5, eps = numeric(pk$n_exo),
                s_back = s_sm, v = NULL, F_inv = NULL, L = m$TT, K = NULL,
                obs_ok = integer(0), ZZ = NULL,
                DD = NULL))
  }
  v  <- v_full[ok_obs]
  ZZ <- m$ZZ[ok_obs, , drop = FALSE]
  DD <- m$DD[ok_obs, , drop = FALSE]
  n  <- length(ok_obs)
  Fm <- ZZ %*% P %*% t(ZZ) + DD %*% Sig %*% t(DD) + pk$me * diag(n)
  Fm <- (Fm + t(Fm)) * 0.5
  ## Dynare rejects a rank-deficient F (error 326): a failed candidate.
  sc <- sqrt(pmax(diag(Fm), 0))
  if (any(!is.finite(Fm)) || any(sc <= 0) ||
      rcond(Fm / tcrossprod(sc)) < 1e-15)
    return(list(ok = FALSE))
  Fc    <- chol(Fm)
  F_inv <- chol2inv(Fc)
  Fiv   <- drop(F_inv %*% v)
  ll    <- -0.5 * (n * log(2 * pi) + 2 * sum(log(diag(Fc))) + sum(v * Fiv))
  K     <- (m$TT %*% P %*% t(ZZ) + m$RR %*% Sig %*% t(DD)) %*% F_inv
  L     <- m$TT - K %*% ZZ
  RmKD  <- m$RR - K %*% DD
  P_new <- tcrossprod(L %*% P, L) + tcrossprod(RmKD %*% Sig, RmKD) +
           pk$me * tcrossprod(K)
  eps    <- drop(Sig %*% t(DD) %*% Fiv)
  s_back <- drop(s_sm + P %*% (t(ZZ) %*% Fiv))
  list(ok = TRUE, ll = ll,
       s_new = drop(m$TT %*% s_back) + drop(m$RR %*% eps) + m$cs,
       s_kf  = drop(m$TT %*% s) + drop(K %*% v) + m$cs,
       P_new = (P_new + t(P_new)) * 0.5,
       eps = eps, s_back = s_back,
       v = v, F_inv = F_inv, L = L, K = K, obs_ok = ok_obs, ZZ = ZZ,
       DD = DD)
}

#' Regime sequence without trailing slack periods (integer(0) = all slack)
#' @noRd
.obc_pkf_trim <- function(x) {
  x <- as.integer(x)
  nz <- which(x != 0L)
  if (length(nz) == 0L) integer(0) else x[seq_len(max(nz))]
}

#' Period-1 rule of an expected regime sequence (NULL = slack); FALSE when
#' the sequence has a singular system
#' @noRd
.obc_pkf_rule <- function(pk, seq) {
  if (length(seq) == 0L) return(NULL)
  r <- .obc_pwl_rules(pk$ctx, seq, strict = FALSE)
  if (is.null(r)) return(FALSE)
  r[[1L]]
}

#' OccBin solve of the regime sequence expected at t
#'
#' Shock eps in period 1 (period t), none after; initial state s0
#' (= s_{t-1|t}); horizon pk$H, doubled while the last period binds.  An
#' all-slack initial guess is first checked with the precomputed slack tail
#' (one matrix product) before the full guess-and-verify.  maxit = 30 is
#' Dynare's options_.occbin.likelihood.maxit.
#' @param tol relative round-off band of the regime decisions
#'   (.obc_gap_binds(), R/obc-binding.R)
#' @return list(converged, seq) with seq trimmed
#' @noRd
.obc_pkf_solve <- function(pk, eps, s0, init = NULL, maxit = 30L,
                           tol = 1e-8) {
  ctx <- pk$ctx
  if (length(init) == 0L && !is.null(pk$tail)) {
    y1 <- drop(pk$ghx %*% s0) + drop(pk$ghu %*% eps)
    x  <- c(y1[ctx$var], drop(pk$tail %*% y1[pk$si]))
    bx <- .obc_gap_binds(ctx$sgn * (x - ctx$bnd), x, ctx$bnd, tol)
    if (!any(bx))
      return(list(converged = TRUE, seq = integer(0)))
  }
  H <- max(pk$H, length(init) + 1L)
  repeat {
    e <- matrix(0, pk$n_exo, H)
    e[, 1L] <- eps
    ini <- if (length(init) > 0L) c(init, integer(H - length(init))) else NULL
    r <- .obc_pwl_solve(ctx, e, s0, ini, max_iter = maxit, tol = tol,
                        strict = FALSE)
    if (!isTRUE(r$converged)) return(list(converged = FALSE, seq = NULL))
    if (r$regime_path[H] == 0L) break
    if (H >= 16L * pk$H)                      # binding at the horizon's end
      return(list(converged = FALSE, seq = NULL))
    init <- .obc_pkf_trim(r$regime_path)
    H <- 2L * H
  }
  list(converged = TRUE, seq = .obc_pkf_trim(r$regime_path))
}

#' One PKF update from a starting regime guess (Dynare kalman_update_algo_1)
#'
#' Kalman step with the period-t rule of the guess seq0; OccBin solve of the
#' expected regime sequence from (eps_{t|t}, s_{t-1|t}); repeat with the
#' solved sequence until it reproduces the sequence whose rule produced the
#' update (at most max_iter re-updates).  A cycle is resolved as in Dynare
#' (filter.periodic_solution = true): the cycle member with the highest
#' likelihood is imposed without verification (Dynare re-solves with
#' maxit = 1, which returns the guessed regime's path), its rule gives the
#' likelihood, and the state is its path from the latest update's shock and
#' smoothed state.
#' @param guess TRUE: seq0 is a forced guess (Dynare guess_regime; the first
#'   solve starts from it); FALSE: the first solve starts all-slack
#' @return list(ok, seq, up, rule) -- up the accepted Kalman step
#' @noRd
.obc_pkf_algo1 <- function(pk, s, P, y, d, seq0, guess, max_iter,
                           s_sm = s) {
  fail <- list(ok = FALSE)
  rule <- .obc_pkf_rule(pk, seq0)
  if (identical(rule, FALSE)) return(fail)
  up <- .obc_pkf_kf_step(pk, .obc_pkf_mats(pk, rule), s, P, y, d, s_sm)
  if (!up$ok) return(fail)
  sol <- .obc_pkf_solve(pk, up$eps, up$s_back, if (guess) seq0 else NULL)
  if (!sol$converged)
    sol <- .obc_pkf_solve(pk, up$eps, up$s_back, if (guess) NULL else seq0)
  if (!sol$converged) return(fail)
  hist    <- list(seq0)
  hist_ll <- up$ll
  cur <- seq0
  new <- sol$seq
  niter <- 1L
  while (!identical(new, cur) && niter <= max_iter) {
    niter <- niter + 1L
    rule  <- .obc_pkf_rule(pk, new)
    if (identical(rule, FALSE)) return(fail)
    up <- .obc_pkf_kf_step(pk, .obc_pkf_mats(pk, rule), s, P, y, d, s_sm)
    if (!up$ok) return(fail)
    hist[[niter]]  <- new
    hist_ll[niter] <- up$ll
    sol <- .obc_pkf_solve(pk, up$eps, up$s_back, new)
    if (!sol$converged) sol <- .obc_pkf_solve(pk, up$eps, up$s_back, NULL)
    if (!sol$converged) return(fail)
    cur <- new
    new <- sol$seq
    per <- which(vapply(hist[seq_len(niter - 1L)], identical, logical(1),
                        new))
    if (length(per) > 0L) {
      idx  <- seq.int(per[1L], niter)
      best <- hist[[idx[which.max(hist_ll[idx])]]]
      rule <- .obc_pkf_rule(pk, best)
      if (identical(rule, FALSE)) return(fail)
      m    <- .obc_pkf_mats(pk, rule)
      ## Dynare's state is the solver's path: the accepted rule applied to
      ## the shock and smoothed state of the update BEFORE this re-update
      ## (only at a non-periodic fixed point do the two coincide).
      s_path <- drop(m$TT %*% up$s_back) + drop(m$RR %*% up$eps) + m$cs
      up <- .obc_pkf_kf_step(pk, m, s, P, y, d, s_sm)
      if (!up$ok) return(fail)
      up$s_new <- s_path
      cur <- new <- best
      break
    }
  }
  if (!identical(new, cur)) return(fail)            # Dynare error 331
  list(ok = TRUE, seq = new, up = up, rule = rule)
}

#' Brute-force regime guesses of Dynare's kalman_update_engine
#'
#' Groups of guesses (tried in order until a group yields a solution):
#' one constraint -- for k = 1..5, binding in periods 5..4+4k ("in
#' expectation"), then 1..4k ("now"); two constraints -- Dynare's loops over
#' the other constraint's pattern (slack; binding now 1..4; binding in
#' expectation 5..8), the duration k, the constraint and now/expected.  With
#' more than two constraints the two-constraint scheme is applied to each
#' constraint with every other constraint following the same pattern.
#' @return list of groups, each a list of trimmed bitfield sequences
#' @noRd
.obc_pkf_brute_guesses <- function(n_spec) {
  pat <- function(kind, k) {
    if (kind == 0L) integer(0)
    else if (kind == 1L) c(integer(4L), rep(1L, 4L * k))  # expected: 5..4+4k
    else rep(1L, 4L * k)                                   # now: 1..4k
  }
  comb <- function(parts) {                 # per-spec 0/1 vectors -> bitfield
    n <- max(c(0L, lengths(parts)))
    out <- integer(n)
    for (j in seq_along(parts)) {
      b <- parts[[j]]
      if (length(b)) out[seq_along(b)] <- out[seq_along(b)] + b * 2L^(j - 1L)
    }
    .obc_pkf_trim(out)
  }
  groups <- list()
  if (n_spec == 1L) {
    for (k in 1:5) groups[[k]] <- list(pat(1L, k), pat(2L, k))
    return(groups)
  }
  for (jk in 0:1) for (k in 1:5) {
    g <- list()
    for (jr in seq_len(n_spec)) for (kk in 1:2) {
      others <- if (jk == 0L) list(0L) else list(c(2L, 1L), c(1L, 1L))
      for (o in others) {
        parts <- vector("list", n_spec)
        for (j in seq_len(n_spec))
          parts[[j]] <- if (j == jr) pat(kk, k) else pat(o[1L], if (jk) jk else 1L)
        g[[length(g) + 1L]] <- comb(parts)
      }
    }
    groups[[length(groups) + 1L]] <- g
  }
  groups
}

#' One period of the PKF (Dynare kalman_update_engine, multivariate)
#'
#' Candidate 1 starts from the all-slack guess; candidate 2 from the regime
#' sequence expected in t-1 for t onwards (when it binds somewhere).  When a
#' candidate fails, Dynare's brute-force guesses are tried.  The state
#' carried forward is the highest-likelihood candidate's; when candidates
#' converge to DIFFERENT regime sequences the period's likelihood is the sum
#' of their densities (Dynare: likx = -2 log sum exp(-lik/2)).
#' @return list(ok, ll, seq, up, rule, n_solutions)
#' @noRd
.obc_pkf_period <- function(pk, s, P, y, d, prev_seq, max_iter, s_sm = s) {
  cands <- list()
  add <- function(cands, c) {
    for (x in cands) if (identical(x$seq, c$seq)) return(cands)
    c(cands, list(c))
  }
  best <- NULL
  c1 <- .obc_pkf_algo1(pk, s, P, y, d, integer(0), FALSE, max_iter, s_sm)
  info0 <- !c1$ok
  if (c1$ok) { cands <- list(c1); best <- c1 }
  info1 <- info0
  if (length(prev_seq) > 0L) {
    c2 <- .obc_pkf_algo1(pk, s, P, y, d, prev_seq, FALSE, max_iter, s_sm)
    info1 <- !c2$ok
    if (c2$ok) {
      cands <- add(cands, c2)
      if (is.null(best) || c2$up$ll > best$up$ll) best <- c2
    }
  }
  if (info0 || info1) {
    found <- NULL
    for (g in .obc_pkf_brute_guesses(pk$ctx$n_spec)) {
      for (gs in g) {
        cb <- .obc_pkf_algo1(pk, s, P, y, d, gs, TRUE, max_iter, s_sm)
        if (cb$ok) {
          found <- cb
          if (is.null(best) || !identical(cb$seq, best$seq)) break
        }
      }
      if (!is.null(found)) break
    }
    if (!is.null(found)) {
      cands <- add(cands, found)
      if (is.null(best) || found$up$ll > best$up$ll) best <- found
    }
  }
  if (is.null(best)) return(list(ok = FALSE))
  lls <- vapply(cands, function(x) x$up$ll, numeric(1))
  ll  <- if (length(lls) > 1L) max(lls) + log(sum(exp(lls - max(lls))))
         else best$up$ll
  list(ok = TRUE, ll = ll, seq = best$seq, up = best$up, rule = best$rule,
       n_solutions = length(lls))
}


# =============================================================================
# Piecewise-linear Kalman filter (Dynare OccBin PKF)
# =============================================================================

#' Piecewise-linear Kalman filter for OBC models (Dynare's OccBin PKF)
#'
#' Implements the piecewise-linear Kalman filter of Giovannini, Pfeiffer and
#' Ratto (2021) as Dynare 7.1 runs it for \code{occbin_setup;
#' estimation(...)} (\code{options_.occbin.likelihood.status = true}; files
#' \code{+occbin/kalman_update_engine.m} and \code{kalman_update_algo_1.m}).
#' In every period t:
#' \enumerate{
#'   \item a guess of the regime sequence expected at t (over a check-ahead
#'     horizon) gives the time-varying OccBin rule of period t;
#'   \item a Kalman step with that rule gives the likelihood contribution,
#'     the shock estimate \eqn{\epsilon_{t|t}} and the one-step smoothed
#'     lagged state \eqn{s_{t-1|t}};
#'   \item the OccBin guess-and-verify solve from \eqn{s_{t-1|t}} with shock
#'     \eqn{\epsilon_{t|t}} (and no further shocks) gives the expected regime
#'     sequence; steps 1-3 repeat until it reproduces the guess
#'     (\code{max_inner} re-updates, Dynare's
#'     \code{likelihood.max_number_of_iterations = 10}).
#' }
#' The update is started from the all-slack guess and, when it binds
#' somewhere, from the sequence expected in t-1; if one of them fails the
#' brute-force guesses of Dynare's engine are tried.  Different converged
#' regime sequences contribute the SUM of their densities; the state follows
#' the most likely one.  A period in which no guess converges makes the
#' log-likelihood \code{-Inf} (Dynare rejects the draw); \code{$failed_period}
#' records it.
#'
#' Changed in 0.9.3.93 (W49): the previous filter used ONE policy per
#' regime (the next period assumed slack) and checked only period t against
#' the slack rule, which is wrong for spells of two or more periods.
#' \code{regime_path_init} now plays the role of Dynare's
#' \code{likelihood.init_regime_history}: the expected regime sequence for
#' the first period (it used to warm-start every period, which made the
#' likelihood depend on the warm start).
#'
#' @param data            n_obs x T observation matrix
#' @param dr_slack        Slack-regime DecisionRules
#' @param regime_cache    R environment seeded by obc_ensure_policy(0, ...);
#'                        the filter stores its per-period rules there
#'                        (\code{.pkf_rules}) for kalman_filter_obc() and
#'                        historical_decomposition_obc()
#' @param sys             System matrices (from extract_system_matrices_fast)
#' @param model           dynhr_mod
#' @param params          Named numeric parameter vector
#' @param obs_vars        Character vector of observed variable names
#' @param specs           OBC spec list from obc_parse_tags
#' @param obs_idx         Integer vector: observable positions in endo vector
#' @param regime_path_init NULL, or an integer bitfield vector: the regime
#'                        sequence expected in period 1 (an extra starting
#'                        guess for period 1 only)
#' @param me_variance     Measurement error variance added to every
#'                        observable (default 1e-8)
#' @param max_inner       Maximum re-updates per period and starting guess
#'                        (default 10, Dynare's default)
#' @param horizon         Check-ahead horizon of the regime solves (default
#'                        200, Dynare's \code{check_ahead_periods}); doubled
#'                        while a solution binds in its last period
#' @param return_filtered Logical: return the n_state x T filtered states
#' @param return_shocks   Logical: return the n_exo x T shocks
#'                        \eqn{\epsilon_{t|t}}
#' @param return_store    Logical: return the per-period quantities needed by
#'                        pkf_smoother_obc()
#' @param return_P_last   Logical: return \eqn{P_{T|T}} in \code{$P_last}
#' @return List with:
#'   \item{loglik}{total log-likelihood (\code{-Inf} on failure)}
#'   \item{loglik_t}{per-period log-likelihood contributions (NA after a
#'     failed period)}
#'   \item{regime_path}{integer vector (length T): regime in period t}
#'   \item{regime_expected}{list (length T): the regime sequence expected in
#'     period t, trimmed after its last binding period (\code{integer(0)}:
#'     all slack)}
#'   \item{n_solutions}{integer vector: distinct regime solutions per period}
#'   \item{failed_period}{NA, or the period in which no guess converged}
#'   \item{filtered_states, filtered_shocks, kf_store, P_last}{as requested}
#'   \item{n_obs, n_T}{dimensions}
#'
#' @references
#'   Giovannini, M., Pfeiffer, P. and Ratto, M. (2021). Efficient and robust
#'     inference of models with occasionally binding constraints. JRC Working
#'     Papers in Economics and Finance 2021/03.
#' @export
kalman_filter_obc_pkf <- function(data, dr_slack, regime_cache, sys,
                                   model, params, obs_vars, specs,
                                   obs_idx          = NULL,
                                   regime_path_init = NULL,
                                   me_variance      = 1e-8,
                                   max_inner        = 10L,
                                   horizon          = 200L,
                                   return_filtered  = FALSE,
                                   return_shocks    = FALSE,
                                   return_store     = FALSE,
                                   return_P_last    = FALSE) {
  endo  <- dr_slack$endo_names
  exo   <- dr_slack$exo_names
  n_obs <- length(obs_vars)
  if (is.null(obs_idx)) obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    .dynhr_abort("kalman_filter_obc_pkf: observed variables not in model: ",
                 paste(obs_vars[is.na(obs_idx)], collapse = ", "),
                 class = "dynhr_error_obc_filter")
  if (!exists(".pwl_src", envir = regime_cache, inherits = FALSE))
    obc_ensure_policy(0L, regime_cache, sys, dr_slack, specs, obs_idx)

  Sigma_e <- .get_shock_cov(model, exo, params)
  pk <- .obc_pkf_prep(.obc_pwl_context(sys, dr_slack, specs), dr_slack,
                      obs_idx, Sigma_e, me_variance, horizon)
  d  <- .obc_pkf_obs_ss(dr_slack, obs_vars)
  n_state <- pk$n_state
  n_exo   <- pk$n_exo

  if (is.null(dim(data))) data <- matrix(data, nrow = n_obs)
  if (nrow(data) != n_obs) data <- t(data)
  n_T <- ncol(data)

  s     <- numeric(n_state)
  s_sm  <- s
  s_lin <- s        # plain Kalman recursion along the accepted rules (store)
  P     <- .obc_pkf_P0(pk)
  prev <- if (length(regime_path_init) > 0L)
            .obc_pkf_trim(regime_path_init) else integer(0)

  loglik    <- 0
  ll_t      <- rep(NA_real_, n_T)
  regime    <- integer(n_T)
  expected  <- vector("list", n_T)
  rules     <- vector("list", n_T)
  n_sol     <- integer(n_T)
  filtered  <- if (return_filtered) matrix(0, n_state, n_T) else NULL
  shocks    <- if (return_shocks)   matrix(0, n_exo,   n_T) else NULL
  kf_store  <- if (return_store) list(
    n_T = n_T, n_state = n_state, n_exo = n_exo, Sigma_e = Sigma_e,
    endo_state_names = endo[dr_slack$state_idx], exo_names = exo,
    v = vector("list", n_T), F_inv = vector("list", n_T),
    L = vector("list", n_T), P_in = vector("list", n_T),
    s_in = matrix(0, n_state, n_T),
    TT = vector("list", n_T), RR = vector("list", n_T),
    ZZ = vector("list", n_T), DD = vector("list", n_T),
    rules = NULL) else NULL

  for (t in seq_len(n_T)) {
    res <- .obc_pkf_period(pk, s, P, data[, t], d, prev, max_inner, s_sm)
    if (!res$ok)
      return(list(loglik = -Inf, loglik_t = ll_t, regime_path = regime,
                  regime_expected = expected, n_solutions = n_sol,
                  failed_period = t, filtered_states = NULL,
                  filtered_shocks = NULL, kf_store = NULL, P_last = NULL,
                  n_obs = n_obs, n_T = n_T))
    up <- res$up
    loglik      <- loglik + res$ll
    ll_t[t]     <- res$ll
    regime[t]   <- if (length(res$seq)) res$seq[1L] else 0L
    expected[t] <- list(res$seq)
    rules[t]    <- list(res$rule)
    n_sol[t]    <- res$n_solutions
    if (return_store) {
      ## The smoother gets the plain Kalman recursion along the accepted
      ## rules (s_lin, its innovations); it equals the filter's state except
      ## after a period whose regime iteration cycled (see .obc_pkf_kf_step).
      m <- .obc_pkf_mats(pk, res$rule)
      v_lin <- NULL
      s_lin_new <- drop(m$TT %*% s_lin) + m$cs
      if (length(up$obs_ok)) {
        v_lin <- (data[, t] - drop(m$ZZ %*% s_lin) - d - m$co)[up$obs_ok]
        s_lin_new <- s_lin_new + drop(up$K %*% v_lin)
      }
      kf_store$P_in[[t]]  <- P
      kf_store$s_in[, t]  <- s_lin
      kf_store$L[[t]]     <- up$L
      kf_store$TT[[t]]    <- m$TT
      kf_store$RR[[t]]    <- m$RR
      kf_store$v[t]       <- list(v_lin)
      kf_store$F_inv[t]   <- list(up$F_inv)
      kf_store$ZZ[t]      <- list(up$ZZ)
      kf_store$DD[t]      <- list(up$DD)
      s_lin <- s_lin_new
    }
    s    <- up$s_new
    s_sm <- up$s_kf
    P    <- up$P_new
    prev <- if (length(res$seq) > 1L) .obc_pkf_trim(res$seq[-1L]) else integer(0)
    if (return_filtered) filtered[, t] <- s
    if (return_shocks)   shocks[, t]   <- up$eps
  }

  ## Per-period rules for kalman_filter_obc() / historical_decomposition_obc()
  assign(".pkf_rules", list(regime_path = regime, rules = rules,
                            obs_idx = as.integer(obs_idx)),
         envir = regime_cache)
  if (return_store) kf_store$rules <- rules
  if (return_filtered) rownames(filtered) <- endo[dr_slack$state_idx]
  if (return_shocks)   rownames(shocks)   <- exo

  list(
    loglik          = loglik,
    loglik_t        = ll_t,
    regime_path     = regime,
    regime_expected = expected,
    n_solutions     = n_sol,
    failed_period   = NA_integer_,
    filtered_states = filtered,
    filtered_shocks = shocks,
    kf_store        = kf_store,
    P_last          = if (return_P_last) P else NULL,
    n_obs = n_obs, n_T = n_T
  )
}
