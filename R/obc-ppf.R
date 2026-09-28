## R/obc-ppf.R
## --------------------------------------------------------------------------
## OBC Piecewise Particle Filter (PPF): bootstrap and conditionally-optimal
## (COPF) proposals for OccBin piecewise-linear models.
##
## TRANSITION (W50, 2026-09-25).  A particle carries s_{t-1}.  In period t it
## draws eps_t and follows the OccBin rule of the regime SEQUENCE expected
## from ITS OWN (s_{t-1}, eps_t): the guess-and-verify solve of Dynare's
## OccBin (shock eps_t in period t, no shock after, check-ahead horizon 200,
## slack after it) gives the expected sequence, and the period-t rule of that
## sequence is the time-varying rule of .obc_pwl_rules() (R/obc-binding.R):
##   y_t = ghx_t s_{t-1} + ghu_t eps_t + c_t,   s_t = y_t[state],
##   obs_t = y_t[obs] + d + me.
## This is the per-period step of the piecewise-linear Kalman filter
## (kalman_filter_obc_pkf(), R/obc-filter.R) applied to a POINT state and a
## drawn shock instead of the filtered mean: the solve starts all-slack and,
## when that does not converge, from the sequence the particle expected in
## t-1 for t onwards (the PKF's two starting guesses).  A particle whose
## solve does not converge has no model solution and gets weight zero.
##
## Vectorisation: particles are grouped by their current regime-sequence
## guess; each group's rules are built once and the guess is verified for
## the whole group with matrix products (the slack continuation after the
## last guessed binding period is one product with the precomputed slack
## tail of .obc_pkf_prep()).  An all-slack solution -- the common case -- is
## one product for the whole cloud.
##
## Before W50 each particle checked period t only, against the SLACK rule
## (pkf_check_binding()), and a binding period used the one-period policy of
## obc_ensure_policy() (next period slack): exact for one-period spells only.
##
## References:
##   Aruoba, Cuba-Borda, Higa-Flores, Schorfheide & Villalvazo (2021, RED)
##   Guerrieri & Iacoviello (2015, JME); Dynare 7.1 +occbin
##
## Provides:
##   ppf_likelihood()              -- per-evaluation PF loglik
##   make_log_posterior_obc_ppf()  -- closure factory (mirrors PKF version)
## --------------------------------------------------------------------------


## ============================================================================
## Regime-sequence engine for particle clouds
## ============================================================================

#' Key of a trimmed regime sequence ("" = all slack)
#' @noRd
.ppf_seq_key <- function(seq) {
  if (length(seq) == 0L) "" else paste(seq, collapse = ",")
}

#' Regime sequence of a key (integer(0) = all slack)
#' @noRd
.ppf_key_seq <- function(key) {
  if (is.na(key) || !nzchar(key)) integer(0)
  else as.integer(strsplit(key, ",", fixed = TRUE)[[1L]])
}

#' Particle-filter engine: PKF context, eq-row matrices, rule caches
#'
#' @param sys, dr_slack, specs, obs_idx  model pieces
#' @param Sigma_e     shock covariance
#' @param me_variance measurement-error variance
#' @param horizon     check-ahead horizon of the regime solves (Dynare: 200)
#' @noRd
.ppf_engine <- function(sys, dr_slack, specs, obs_idx, Sigma_e, me_variance,
                        horizon = 200L) {
  ctx <- .obc_pwl_context(sys, dr_slack, specs)
  pk  <- .obc_pkf_prep(ctx, dr_slack, obs_idx, Sigma_e, me_variance, horizon)
  eq  <- ctx$eq
  list(pk = pk, ctx = ctx,
       Fm_e = ctx$Fm[eq, , drop = FALSE], F0_e = ctx$F0[eq, , drop = FALSE],
       Fp_e = ctx$Fp[eq, , drop = FALSE], Fe_e = ctx$Fe[eq, , drop = FALSE],
       w = 2L^(seq_len(ctx$n_spec) - 1L),
       rules = new.env(parent = emptyenv(), hash = TRUE),
       mats  = new.env(parent = emptyenv(), hash = TRUE))
}

#' Time-varying rules of a regime sequence (cached; NULL = singular system)
#' @noRd
.ppf_rules <- function(eng, key) {
  ek <- paste0("s", key)                 # environment names must be non-empty
  if (exists(ek, envir = eng$rules, inherits = FALSE))
    return(get(ek, envir = eng$rules, inherits = FALSE))
  r <- if (!nzchar(key)) list()
       else .obc_pwl_rules(eng$ctx, .ppf_key_seq(key), strict = FALSE)
  assign(ek, r, envir = eng$rules)
  r
}

#' State-space matrices of the period rule of an expected sequence (cached)
#'
#' A singular sequence falls back to the slack matrices (only used as a
#' COPF proposal guess; a solved sequence is never singular).
#' @noRd
.ppf_mats <- function(eng, key) {
  ek <- paste0("s", key)
  if (exists(ek, envir = eng$mats, inherits = FALSE))
    return(get(ek, envir = eng$mats, inherits = FALSE))
  r <- .ppf_rules(eng, key)
  m <- if (length(r) == 0L) eng$pk$slack else .obc_pkf_mats(eng$pk, r[[1L]])
  assign(ek, m, envir = eng$mats)
  m
}

#' Expected sequences one period later: drop the first period, trim
#' @noRd
.ppf_shift_keys <- function(keys) {
  keys[is.na(keys)] <- ""
  u <- unique(keys)
  su <- vapply(u, function(k) {
    s <- .ppf_key_seq(k)
    if (length(s) <= 1L) "" else .ppf_seq_key(.obc_pkf_trim(s[-1L]))
  }, character(1), USE.NAMES = FALSE)
  su[match(keys, u)]
}

#' Verify one regime-sequence guess for a group of particles
#'
#' One pass of the OccBin guess-and-verify (.obc_pwl_solve() iteration) over
#' the check-ahead horizon pk$H for every column of (S, E): the guess's rules
#' in periods 1..L (L = last guessed binding period), the slack continuation
#' L+1..H through the precomputed slack tail, the complementarity check of
#' .obc_pwl_check() (a slack period binds when the constrained variable
#' violates its bound; a binding period is released when the multiplier of
#' its relaxed equation has the wrong sign).
#' @return list(ok = FALSE) for a singular guess, else list(ok, conv, edge,
#'   newkey, Y1): conv = the check reproduces the guess; edge = converged
#'   but binding in period H (the PKF solve then doubles the horizon);
#'   newkey = the checked sequence of the non-converged columns; Y1 =
#'   n_endo x n period-1 values under the guess's rule
#' @noRd
.ppf_eval_guess <- function(eng, key, S, E, tol) {
  pk  <- eng$pk
  ctx <- eng$ctx
  rules <- .ppf_rules(eng, key)
  if (is.null(rules)) return(list(ok = FALSE))
  sq <- .ppf_key_seq(key)
  L  <- length(sq)
  H  <- pk$H
  n  <- ncol(S)
  ns <- ctx$n_spec
  Lx <- max(L, 1L)
  bits <- matrix(0L, H, n)
  Sp <- S
  Y1 <- NULL
  for (t in seq_len(Lx)) {
    ru <- if (t <= L) rules[[t]] else NULL
    if (is.null(ru)) {
      Yt <- ctx$ghx %*% Sp
      if (t == 1L) Yt <- Yt + ctx$ghu %*% E
    } else {
      Yt <- ru$ghx %*% Sp + ru$c
      if (t == 1L) Yt <- Yt + ru$ghu %*% E
    }
    St  <- Yt[ctx$si, , drop = FALSE]
    Xt  <- Yt[ctx$var, , drop = FALSE]
    gap <- ctx$sgn * (Xt - ctx$bnd)
    Bn  <- .obc_gap_binds(gap, Xt, ctx$bnd, tol)
    rb  <- if (t <= L) obc_regime_flags(sq[t], ns) else logical(ns)
    if (any(rb)) {
      ru1 <- if (t < L) rules[[t + 1L]] else NULL
      Ey1 <- if (is.null(ru1)) ctx$ghx %*% St else ru1$ghx %*% St + ru1$c
      Fr  <- eng$Fm_e %*% Sp + eng$F0_e %*% Yt + eng$Fp_e %*% Ey1
      Mg  <- abs(eng$Fm_e) %*% abs(Sp) + abs(eng$F0_e) %*% abs(Yt) +
             abs(eng$Fp_e) %*% abs(Ey1)
      if (t == 1L) {
        Fr <- Fr + eng$Fe_e %*% E
        Mg <- Mg + abs(eng$Fe_e) %*% abs(E)
      }
      Bn[rb, ] <- .obc_mult_keeps((ctx$sgn * Fr)[rb, , drop = FALSE],
                                  Mg[rb, , drop = FALSE], tol)
    }
    bits[t, ] <- as.integer(colSums(Bn * eng$w))
    if (t == 1L) Y1 <- Yt
    Sp <- St
  }
  if (H > Lx && !is.null(pk$tail)) {
    ## slack continuation: period Lx + k is tail block k applied to s_{Lx}
    np <- H - Lx
    X  <- pk$tail[seq_len(np * ns), , drop = FALSE] %*% Sp
    Bt <- .obc_gap_binds(rep(ctx$sgn, np) * (X - rep(ctx$bnd, np)), X,
                         rep(ctx$bnd, np), tol)
    tb <- if (ns == 1L) Bt
          else rowsum(Bt * rep(eng$w, np), rep(seq_len(np), each = ns),
                      reorder = FALSE)
    bits[Lx + seq_len(np), ] <- as.integer(tb)
  }
  gb <- c(sq, integer(H - L))
  conv <- colSums(bits != gb) == 0L
  edge <- conv & bits[H, ] != 0L
  newkey <- rep(NA_character_, n)
  for (j in which(!conv)) {
    nz <- which(bits[, j] != 0L)
    newkey[j] <- if (length(nz) == 0L) ""
                 else paste(bits[seq_len(max(nz)), j], collapse = ",")
  }
  list(ok = TRUE, conv = conv, edge = edge, newkey = newkey, Y1 = Y1)
}

#' OccBin guess-and-verify for a particle cloud (vectorised .obc_pkf_solve)
#'
#' Same iterates as .obc_pkf_solve(pk, E[, i], S[, i], init) for every
#' particle i: at most maxit guesses, a repeated guess (cycle) or a singular
#' guess is a failure.  Particles whose solution binds in the last period of
#' the horizon (or whose initial guess is longer than it) are handed to
#' .obc_pkf_solve() itself, which doubles the horizon.
#' @param init_keys character, per particle: the starting guess ("" slack)
#' @param tol relative round-off band of the regime decisions
#'   (.obc_gap_binds(), R/obc-binding.R)
#' @return list(key = solved sequence per particle (NA: no solution),
#'   Y1 = n_endo x N period-t values under the solved rule)
#' @noRd
.ppf_solve_batch <- function(eng, S, E, init_keys, maxit = 30L, tol = 1e-8) {
  pk  <- eng$pk
  N   <- ncol(S)
  key <- rep(NA_character_, N)
  Y1  <- matrix(NA_real_, pk$n_endo, N)
  if (eng$ctx$n_spec == 0L) {
    return(list(key = rep("", N),
                Y1 = eng$ctx$ghx %*% S + eng$ctx$ghu %*% E))
  }
  init_len <- ifelse(nzchar(init_keys),
                     nchar(init_keys) - nchar(gsub(",", "", init_keys,
                                                   fixed = TRUE)) + 1L, 0L)
  undecided <- init_len >= pk$H
  cur    <- init_keys
  seen   <- vector("list", N)
  active <- which(!undecided)
  for (iter in seq_len(maxit)) {
    if (length(active) == 0L) break
    nxt <- integer(0)
    grp <- split(active, cur[active])
    for (g in seq_along(grp)) {
      k   <- names(grp)[g]
      idx <- grp[[g]]
      ev  <- .ppf_eval_guess(eng, k, S[, idx, drop = FALSE],
                             E[, idx, drop = FALSE], tol)
      if (!ev$ok) next                       # singular guess: no solution
      fin <- ev$conv & !ev$edge
      if (any(fin)) {
        key[idx[fin]]   <- k
        Y1[, idx[fin]]  <- ev$Y1[, fin, drop = FALSE]
      }
      undecided[idx[ev$edge]] <- TRUE
      for (j in which(!ev$conv)) {
        i <- idx[j]
        seen[[i]] <- c(seen[[i]], k)
        if (ev$newkey[j] %in% seen[[i]]) next          # cycle: no solution
        cur[i] <- ev$newkey[j]
        nxt <- c(nxt, i)
      }
    }
    active <- nxt
  }
  for (i in which(undecided)) {
    sol <- .obc_pkf_solve(pk, E[, i], S[, i], .ppf_key_seq(init_keys[i]),
                          maxit = maxit, tol = tol)
    if (!isTRUE(sol$converged)) next
    k <- .ppf_seq_key(sol$seq)
    r <- .ppf_rules(eng, k)
    if (is.null(r)) next
    ru <- if (length(r)) r[[1L]] else NULL
    key[i]  <- k
    Y1[, i] <- if (is.null(ru)) drop(pk$ghx %*% S[, i]) + drop(pk$ghu %*% E[, i])
               else drop(ru$ghx %*% S[, i]) + drop(ru$ghu %*% E[, i]) + ru$c
  }
  list(key = key, Y1 = Y1)
}

#' Particle transition: expected regime sequence and period-t values
#'
#' The solve starts all-slack; a particle whose solve does not converge is
#' re-solved from the sequence it expected in t-1 for t onwards (prev_keys),
#' as the PKF's second starting guess.  NA key: no solution.
#' @param S n_state x N states s_{t-1}; E n_exo x N shocks eps_t
#' @param prev_keys character N: expected sequences carried from t-1
#'   (already shifted to start at t)
#' @return list(key, Y1 = n_endo x N)
#' @noRd
.ppf_transition <- function(eng, S, E, prev_keys = NULL) {
  N  <- ncol(S)
  tr <- .ppf_solve_batch(eng, S, E, rep("", N))
  if (!is.null(prev_keys)) {
    retry <- which(is.na(tr$key) & !is.na(prev_keys) & nzchar(prev_keys))
    if (length(retry) > 0L) {
      r2 <- .ppf_solve_batch(eng, S[, retry, drop = FALSE],
                             E[, retry, drop = FALSE], prev_keys[retry])
      tr$key[retry]  <- r2$key
      tr$Y1[, retry] <- r2$Y1
    }
  }
  tr
}

#' Gaussian measurement log-density of every particle (NA observations
#' dropped; the missing pattern is common to all particles)
#' @param pred n_obs x N predicted observables (levels)
#' @noRd
.ppf_log_meas <- function(y_t, pred, me_variance) {
  ok <- which(!is.na(y_t))
  n_ok <- length(ok)
  if (n_ok == 0L) return(numeric(ncol(pred)))
  innov <- y_t[ok] - pred[ok, , drop = FALSE]
  -0.5 * n_ok * log(2 * pi) - 0.5 * n_ok * log(me_variance) -
    0.5 * colSums(innov^2) / me_variance
}

#' Resample and propagate a weighted cloud
#' @return list(particles, log_lik_contrib, keys) or NULL (all weights zero)
#' @noRd
.ppf_resample <- function(log_w, Y1, keys, si) {
  N <- length(log_w)
  log_lik_contrib <- .smc_log_sum_exp(log_w) - log(N)
  if (!is.finite(log_lik_contrib)) return(NULL)
  w_norm <- exp(log_w - max(log_w))
  w_norm <- w_norm / sum(w_norm)
  idx <- .smc_systematic_resample(w_norm, N)
  list(particles = Y1[si, idx, drop = FALSE],
       log_lik_contrib = log_lik_contrib, keys = keys[idx])
}


## ============================================================================
## Per-period bootstrap PPF step (single period, vectorised over particles)
## ============================================================================

#' Bootstrap PPF: process one time period
#'
#' Each particle draws eps_t^i ~ N(0, Sigma_e), solves the regime sequence
#' expected from (s_{t-1}^i, eps_t^i) and follows its period-t rule.
#' Particle weight = Gaussian p(y_t | s_{t-1}^i, eps_t^i).  Accumulates
#' log_lik_contrib BEFORE resampling; resamples at the end of the period.
#'
#' @param particles n_state x N matrix of particles (s_{t-1}^i)
#' @param y_t       length-n_obs observation vector (may contain NAs)
#' @param L_e       n_exo x n_exo lower-triangular Cholesky of Sigma_e
#' @param eng       .ppf_engine()
#' @param d_obs     length-n_obs observable steady-state level
#' @param me_variance Scalar measurement error variance (must be > 0)
#' @param prev_keys character N: each particle's expected sequence from t-1,
#'   shifted to start at t (NULL: all slack)
#' @return list(particles = n_state x N updated, log_lik_contrib = scalar,
#'   keys = the resampled particles' expected sequences at t,
#'   n_failed = particles without a regime solution)
#' @noRd
.ppf_run_period <- function(particles, y_t, L_e, eng, d_obs, me_variance,
                            prev_keys = NULL) {
  N     <- ncol(particles)
  n_exo <- ncol(L_e)
  shocks <- L_e %*% matrix(rnorm(n_exo * N), nrow = n_exo)   # n_exo x N

  tr    <- .ppf_transition(eng, particles, shocks, prev_keys)
  fail  <- is.na(tr$key)
  log_w <- .ppf_log_meas(y_t, tr$Y1[eng$pk$obs_idx, , drop = FALSE] + d_obs,
                         me_variance)
  log_w[fail] <- -Inf

  rs <- .ppf_resample(log_w, tr$Y1, tr$key, eng$pk$si)
  if (is.null(rs)) return(list(particles = particles, log_lik_contrib = -Inf,
                               keys = tr$key, n_failed = sum(fail)))
  rs$n_failed <- sum(fail)
  rs
}


## ============================================================================
## Per-period COPF step (conditionally-optimal proposal)
## ============================================================================

#' Cholesky factor of a symmetric matrix, NULL unless clearly positive
#' definite (explicit eigenvalue test instead of catching a chol() error)
#' @noRd
.ppf_chol <- function(A) {
  if (any(!is.finite(A))) return(NULL)
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) <= max(abs(ev)) * nrow(A) * .Machine$double.eps) return(NULL)
  chol(A)
}

#' COPF quantities of one guessed rule for the observed subset
#' @return list(Omega, L_Omega, F_inv, log_det_F) with NULL entries when the
#'   corresponding matrix is not numerically positive definite
#' @noRd
.copf_quantities <- function(DD_ok, Sigma_e, Sigma_e_inv, me_variance) {
  n_ok <- nrow(DD_ok)
  out  <- list(Omega = NULL, L_Omega = NULL, F_inv = NULL, log_det_F = NULL)
  ch_Oi <- .ppf_chol(crossprod(DD_ok) / me_variance + Sigma_e_inv)
  if (!is.null(ch_Oi)) {
    Om   <- chol2inv(ch_Oi)
    ch_O <- .ppf_chol(Om)
    if (!is.null(ch_O)) {
      out$Omega   <- Om
      out$L_Omega <- t(ch_O)
    }
  }
  ch_F <- .ppf_chol(DD_ok %*% Sigma_e %*% t(DD_ok) + me_variance * diag(n_ok))
  if (!is.null(ch_F)) {
    out$F_inv     <- chol2inv(ch_F)
    out$log_det_F <- 2 * sum(log(diag(ch_F)))
  }
  out
}

#' COPF: process one time period using the conditionally-optimal Gaussian proposal
#'
#' Proposal.  For a guessed expected regime sequence g with period-t rule
#' (ZZ_g, DD_g, c_g), the posterior of eps_t given (s_{t-1}^i, y_t) in the
#' linear-Gaussian model of that rule is N(mu_g^i, Omega_g).  Each particle's
#' guess starts from guess_keys (the ancestor's expectation shifted to t, or
#' all slack) and is iterated to a fixed point ON THE POSTERIOR MEAN, as the
#' PKF does for the filtered mean: g <- the sequence solved from
#' (s_{t-1}^i, mu_g^i), at most copf_iter times.  The guess is a deterministic
#' function of (s_{t-1}^i, its carried expectation, y_t), so
#' q = N(mu_g, Omega_g) is a valid proposal; without the iteration every
#' particle entering a spell proposed from the slack rule and mismatched.
#'
#' Weight.  eps_t^i ~ q; the particle's transition is solved from
#' (s_{t-1}^i, eps_t^i) exactly as in the bootstrap filter
#' (.ppf_transition()).  When the solved sequence is the guess the weight is
#' the marginal N(v_g^i; 0, F_g) (the Gaussian identity); otherwise the
#' generally-valid p(y_t | s, eps, solved rule) p(eps) / q(eps).  A guess whose
#' Omega_g is degenerate draws from the prior with the bootstrap weight.
#'
#' Accumulates log_lik_contrib BEFORE resampling (same as bootstrap PPF).
#'
#' @param particles     n_state x N matrix of particles (s_{t-1}^i)
#' @param y_t           length-n_obs observation vector (may contain NAs)
#' @param L_e           n_exo x n_exo lower-triangular Cholesky of Sigma_e
#' @param Sigma_e       n_exo x n_exo shock covariance
#' @param Sigma_e_inv   n_exo x n_exo inverse of Sigma_e (pre-computed once)
#' @param eng           .ppf_engine()
#' @param d_obs         length-n_obs observable steady-state level
#' @param me_variance   Scalar measurement error variance (must be > 0)
#' @param prev_keys     character N: expected sequences carried from t-1,
#'   shifted to t (NULL: all slack).  Always used by the transition.
#' @param guess_keys    character N: starting proposal guesses (NULL: all
#'   slack)
#' @param U_copf        CPM: list(z_copf = n_exo x N, z_fallback = n_exo x N);
#'   NULL draws fresh from the RNG
#' @param copf_iter     maximum guess updates on the posterior mean
#' @return list(particles, log_lik_contrib, n_fallback = proposal
#'   mismatches + degenerate-Omega prior draws, keys = resampled expected
#'   sequences at t, n_failed, z_copf_used, z_fallback_used)
#' @noRd
.copf_run_period <- function(particles, y_t, L_e, Sigma_e, Sigma_e_inv,
                             eng, d_obs, me_variance,
                             prev_keys = NULL, guess_keys = NULL,
                             U_copf = NULL, copf_iter = 10L) {
  N     <- ncol(particles)
  n_exo <- ncol(L_e)
  pk    <- eng$pk

  ## Both z_copf and z_fallback are always drawn/supplied (CPM: fixed U
  ## structure); only the per-particle SELECTION is data-dependent.
  z_copf     <- if (!is.null(U_copf)) U_copf$z_copf     else matrix(rnorm(n_exo * N), nrow = n_exo)
  z_fallback <- if (!is.null(U_copf)) U_copf$z_fallback else matrix(rnorm(n_exo * N), nrow = n_exo)

  if (is.null(prev_keys)) prev_keys <- rep("", N)
  g <- if (is.null(guess_keys)) rep("", N) else guess_keys
  g[is.na(g)] <- ""
  ok   <- which(!is.na(y_t))
  n_ok <- length(ok)

  ## Proposal of a guess (per period: the missing pattern is common)
  qenv <- new.env(parent = emptyenv(), hash = TRUE)
  getq <- function(k) {
    ek <- paste0("s", k)
    if (exists(ek, envir = qenv, inherits = FALSE))
      return(get(ek, envir = qenv, inherits = FALSE))
    m <- .ppf_mats(eng, k)
    q <- .copf_quantities(m$DD[ok, , drop = FALSE], Sigma_e, Sigma_e_inv,
                          me_variance)
    q$DD_ok  <- m$DD[ok, , drop = FALSE]
    q$ZZ_ok  <- m$ZZ[ok, , drop = FALSE]
    q$off_ok <- d_obs[ok] + m$co[ok]
    assign(ek, q, envir = qenv)
    q
  }
  post <- function(q, idx) {
    V <- y_t[ok] - q$ZZ_ok %*% particles[, idx, drop = FALSE] - q$off_ok
    list(V = V, mu = q$Omega %*% (t(q$DD_ok) %*% V / me_variance))
  }

  ## Guess iteration on the posterior mean
  if (n_ok > 0L) {
    todo <- seq_len(N)
    for (it in seq_len(copf_iter)) {
      mu  <- matrix(0, n_exo, length(todo))
      use <- logical(length(todo))
      grp <- split(seq_along(todo), g[todo])
      for (j in seq_along(grp)) {
        q <- getq(names(grp)[j])
        if (is.null(q$L_Omega)) next
        pos <- grp[[j]]
        mu[, pos] <- post(q, todo[pos])$mu
        use[pos]  <- TRUE
      }
      if (!any(use)) break
      sub <- todo[use]
      tr0 <- .ppf_transition(eng, particles[, sub, drop = FALSE],
                             mu[, use, drop = FALSE], prev_keys[sub])
      new <- tr0$key
      new[is.na(new)] <- g[sub][is.na(new)]
      ch  <- new != g[sub]
      g[sub] <- new
      todo <- sub[ch]
      if (length(todo) == 0L) break
    }
  }

  shocks <- matrix(0, n_exo, N)
  ## per particle: 1 = COPF draw, 2 = prior draw (all missing / degenerate)
  mode   <- integer(N)
  lw_acc <- rep(NA_real_, N)          # N(v_g; 0, F_g) when the guess verifies
  lq     <- rep(NA_real_, N)          # log q(eps) up to the 2*pi constant
  grp    <- split(seq_len(N), g)
  for (j in seq_along(grp)) {
    idx <- grp[[j]]
    if (n_ok == 0L) {
      shocks[, idx] <- L_e %*% z_copf[, idx, drop = FALSE]
      mode[idx] <- 2L
      next
    }
    q <- getq(names(grp)[j])
    if (is.null(q$L_Omega)) {
      shocks[, idx] <- L_e %*% z_fallback[, idx, drop = FALSE]
      mode[idx] <- 2L
      next
    }
    pm <- post(q, idx)
    ep <- pm$mu + q$L_Omega %*% z_copf[, idx, drop = FALSE]
    shocks[, idx] <- ep
    mode[idx] <- 1L
    u_q <- forwardsolve(q$L_Omega, ep - pm$mu)
    lq[idx] <- -sum(log(diag(q$L_Omega))) - 0.5 * colSums(u_q^2)
    lw_acc[idx] <- if (is.null(q$F_inv)) -Inf else
      -0.5 * n_ok * log(2 * pi) - 0.5 * q$log_det_F -
        0.5 * colSums(pm$V * (q$F_inv %*% pm$V))
  }

  tr   <- .ppf_transition(eng, particles, shocks, prev_keys)
  fail <- is.na(tr$key)
  lp_y <- .ppf_log_meas(y_t, tr$Y1[pk$obs_idx, , drop = FALSE] + d_obs,
                        me_variance)
  log_w <- lp_y                                     # prior draws: bootstrap
  acc  <- mode == 1L & !fail & tr$key == g
  mis  <- mode == 1L & !fail & !acc
  log_w[acc] <- lw_acc[acc]
  if (any(mis)) {
    ## p(y | s, eps, solved rule) p(eps) / q(eps); the 2*pi normalisers of
    ## p(eps) and q(eps) cancel
    u_p <- forwardsolve(L_e, shocks[, mis, drop = FALSE])
    log_p_eps <- -sum(log(diag(L_e))) - 0.5 * colSums(u_p^2)
    log_w[mis] <- lp_y[mis] + log_p_eps - lq[mis]
  }
  log_w[fail] <- -Inf
  n_fallback <- sum(mis) + sum(mode == 2L & n_ok > 0L)

  rs <- .ppf_resample(log_w, tr$Y1, tr$key, pk$si)
  if (is.null(rs))
    rs <- list(particles = particles, log_lik_contrib = -Inf, keys = tr$key)
  rs$n_fallback      <- n_fallback
  rs$n_failed        <- sum(fail)
  rs$z_copf_used     <- z_copf
  rs$z_fallback_used <- z_fallback
  rs
}

## ============================================================================
## Main PPF likelihood function
## ============================================================================

#' Piecewise Particle Filter log-likelihood for OBC models
#'
#' Evaluates the marginal log-likelihood log p(Y | theta) with a particle
#' filter over the state of the OccBin piecewise-linear model.
#'
#' At each period t, N particles carry s_{t-1}^i; each particle:
#'   1. Draws eps_t^i from the prior N(0, Sigma_e) (bootstrap) or from the
#'      conditionally-optimal Gaussian proposal of a guessed rule (COPF)
#'   2. Solves the regime SEQUENCE expected from (s_{t-1}^i, eps_t^i)
#'      (OccBin guess-and-verify, check-ahead \code{horizon}) and follows the
#'      period-t time-varying rule of that sequence
#'   3. Computes weight w_t^i = p(y_t | s_{t-1}^i, eps_t^i) (times
#'      p(eps)/q(eps) for the COPF)
#'   4. Accumulates log_lik += log(mean w_t^i) BEFORE resampling
#'   5. Systematically resamples and propagates the state
#' A period with no observation propagates the cloud (weights equal).
#'
#' Changed in W50 (2026-09-25): step 2 used to check period t only against
#' the slack rule and apply the one-period binding policy (next period
#' slack), which is wrong for spells of two or more periods.
#'
#' @param Y             n_obs x T observation matrix
#' @param dr_slack      Slack-regime DecisionRules
#' @param regime_cache  Regime cache seeded by obc_ensure_policy(0, ...)
#' @param sys           System matrices
#' @param model         dynhr_mod
#' @param params        Named numeric parameter vector
#' @param obs_vars      Character vector of observed variable names
#' @param specs         OBC spec list from obc_parse_tags
#' @param obs_idx       Integer vector of observable indices
#' @param N             Number of particles (default 1000)
#' @param me_variance   Measurement error variance (must be > 0; default 1e-4)
#' @param proposal      "bootstrap" or "copf"
#' @param regime_guess  COPF proposal guess: "ancestor" (the resampled
#'   ancestor's expected sequence, shifted to t) or "slack"
#' @param horizon       Check-ahead horizon of the regime solves (default 200,
#'   Dynare's; doubled while a solution binds in its last period)
#' @param return_particles Logical (default FALSE).  When TRUE, append the
#'   terminal \code{n_state x N} particle cloud as \code{$particles} in the
#'   return list.  Has no effect on the loglik path (byte-identical when FALSE).
#' @param seed          Integer RNG seed for reproducibility (NULL = no seed)
#' @return List with:
#'   $loglik       scalar log-likelihood
#'   $n_obs        integer
#'   $n_T          integer
#'   $n_failed     particle-periods without a regime solution (weight zero)
#'   $particles    n_state x N terminal particle cloud (only when return_particles = TRUE)
#' @noRd
ppf_likelihood <- function(Y, dr_slack, regime_cache, sys,
                            model, params, obs_vars, specs,
                            obs_idx          = NULL,
                            N                = 1000L,
                            me_variance      = 1e-4,
                            proposal         = c("bootstrap", "copf"),
                            regime_guess     = c("ancestor", "slack"),
                            U_copf_list      = NULL,  # CPM: list of per-period U_copf structures
                            horizon          = 200L,
                            return_particles  = FALSE,
                            seed             = NULL) {

  proposal     <- match.arg(proposal)
  regime_guess <- match.arg(regime_guess)
  if (!is.null(seed)) {
    ## Local seed: deterministic in theta but leaves the CALLER's RNG stream
    ## untouched (an outer sampler must not replay the same proposals; see
    ## .with_local_seed in tpf-likelihood.R and NEWS 0.9.2.0003).
    .ge_ <- globalenv()
    .had_ <- exists(".Random.seed", envir = .ge_, inherits = FALSE)
    .old_ <- if (.had_) get(".Random.seed", envir = .ge_, inherits = FALSE) else NULL
    on.exit({
      if (.had_) assign(".Random.seed", .old_, envir = .ge_)
      else if (exists(".Random.seed", envir = .ge_, inherits = FALSE))
        rm(list = ".Random.seed", envir = .ge_)
    }, add = TRUE)
    set.seed(seed)
  }

  endo    <- dr_slack$endo_names
  exo     <- dr_slack$exo_names
  n_state <- length(dr_slack$state_idx)
  n_exo   <- length(exo)
  n_obs   <- length(obs_vars)

  if (is.null(obs_idx)) obs_idx <- match(obs_vars, endo)
  ## The measurement density below is N(y; ., me_variance I): a vector would
  ## be silently recycled into it, so refuse it (classed) instead.
  me_variance <- .kf_me_variance(me_variance, obs_vars, "ppf_likelihood",
                                 allow_vector = FALSE)

  ## Data orientation: n_obs x T
  if (is.null(dim(Y))) Y <- matrix(Y, nrow = n_obs)
  if (nrow(Y) != n_obs) Y <- t(Y)
  n_T <- ncol(Y)

  ## Shock covariance + Cholesky factor
  Sigma_e <- .get_shock_cov(model, exo, params)
  L_e     <- t(chol(Sigma_e))   # lower-triangular L s.t. L %*% t(L) = Sigma_e

  ## chol(Sigma_e) succeeded just above (L_e), so no failure branch here.
  Sigma_e_inv <- NULL
  if (proposal == "copf") Sigma_e_inv <- chol2inv(t(L_e))

  ## Regime-sequence engine (theta-dependent: built per call)
  if (!exists(".pwl_src", envir = regime_cache, inherits = FALSE))
    obc_ensure_policy(0L, regime_cache, sys, dr_slack, specs, obs_idx)
  eng   <- .ppf_engine(sys, dr_slack, specs, obs_idx, Sigma_e, me_variance,
                       horizon)
  d_obs <- .obc_pkf_obs_ss(dr_slack, obs_vars)

  ## ---- Initialise particles from slack-policy stationary distribution -----
  pol_s <- eng$pk$slack
  ## Both fallbacks are in the model's units (W79): P_0 scales as Sigma_e.
  ## They used to be diag(1e-6) (no stationary P_0) and a factor diag(1e-3)
  ## (chol failure) -- absolute, behind tryCatch: with every shock std and
  ## the data x 1e-4 a PSD-singular P_0 (an exact linear dependence among
  ## the states) started the cloud thousands of state stds wide, 71 nats
  ## off the rescale identity loglik(c) + N log(c) = const.
  ##  * solve_lyapunov() signals a non-stationary slack policy with NaN (it
  ##    does not throw): start from the one-period covariance QQ_s instead.
  ##  * chol() only where P_0 is positive definite by its eigenvalues
  ##    (lambda_min > 20 n^2.5 eps lambda_max is sufficient for chol to run
  ##    to completion: Higham 2002 Thm 10.7 with van der Sluis'
  ##    kappa(D P_0 D) <= n kappa(P_0)); a PSD-singular (or round-off
  ##    negative) P_0 gets its symmetric square root, negative part clipped.
  QQ_s  <- tcrossprod(pol_s$RR %*% Sigma_e, pol_s$RR)
  P_0   <- solve_lyapunov(pol_s$TT, QQ_s)
  if (!all(is.finite(P_0))) P_0 <- QQ_s
  ev_P  <- eigen(P_0, symmetric = TRUE)
  L_P   <- if (min(ev_P$values) >
               20 * n_state^2.5 * .Machine$double.eps * max(abs(ev_P$values)))
             t(chol(P_0))
           else
             ev_P$vectors %*% (sqrt(pmax(ev_P$values, 0)) * t(ev_P$vectors))
  particles <- L_P %*% matrix(rnorm(n_state * N), nrow = n_state)

  ## ---- Main filter loop ---------------------------------------------------
  loglik           <- 0
  total_fallback   <- 0L
  total_failed     <- 0L
  keys             <- rep("", N)   # expected sequences (all slack before t=1)
  U_copf_realized  <- vector("list", n_T)  # CPM: collect used z_copf/z_fallback per period

  for (t in seq_len(n_T)) {
    y_t  <- Y[, t]
    prev <- .ppf_shift_keys(keys)

    if (proposal == "bootstrap") {
      res <- .ppf_run_period(
        particles   = particles,
        y_t         = y_t,
        L_e         = L_e,
        eng         = eng,
        d_obs       = d_obs,
        me_variance = me_variance,
        prev_keys   = prev
      )
    } else {
      ## The proposal guess only: the transition always carries `prev`.
      ## CPM: thread per-period U_copf if supplied.
      U_copf_t <- if (!is.null(U_copf_list)) U_copf_list[[t]] else NULL
      res <- .copf_run_period(
        particles   = particles,
        y_t         = y_t,
        L_e         = L_e,
        Sigma_e     = Sigma_e,
        Sigma_e_inv = Sigma_e_inv,
        eng         = eng,
        d_obs       = d_obs,
        me_variance = me_variance,
        prev_keys   = prev,
        guess_keys  = if (regime_guess == "ancestor") prev else NULL,
        U_copf      = U_copf_t
      )
      total_fallback <- total_fallback + res$n_fallback
      ## Record used z_copf/z_fallback matrices (CPM: deterministic U structure)
      U_copf_realized[[t]] <- list(z_copf = res$z_copf_used, z_fallback = res$z_fallback_used)
    }

    loglik       <- loglik + res$log_lik_contrib
    total_failed <- total_failed + res$n_failed
    particles    <- res$particles
    keys         <- res$keys

    if (!is.finite(loglik))
      return(list(loglik = -Inf, n_obs = n_obs, n_T = n_T,
                  n_failed = total_failed))
  }

  out <- list(loglik = loglik, n_obs = n_obs, n_T = n_T,
              n_failed = total_failed)
  if (proposal == "copf") {
    out$n_fallback        <- total_fallback
    out$U_copf_realized   <- U_copf_realized  # CPM: used z_copf/z_fallback per period
  }
  ## return_particles = TRUE: append terminal n_state x N cloud (for ctx terminal-state
  ## dispatch; default FALSE preserves byte-identical loglik-only path).
  if (isTRUE(return_particles)) out$particles <- particles
  out
}


## ============================================================================
## Log-posterior factory
## ============================================================================

#' Create a PPF-based log-posterior evaluator for OBC models
#'
#' Particle filter variant of \code{make_log_posterior_obc_pkf}.
#' Each evaluation builds a FRESH regime_cache (theta-dependent) and runs
#' the piecewise particle filter with N particles.  Each particle follows the
#' OccBin piecewise-linear solution from its own state: in every period the
#' regime sequence expected from the particle's state and drawn shock is
#' solved by guess-and-verify (as \code{kalman_filter_obc_pkf} does for the
#' filtered state) and the period's time-varying rule is applied.
#'
#' Changed in W50 (2026-09-25): the particles used to check the current
#' period only, against the slack rule, and to apply the one-period binding
#' policy (next period slack), which is wrong for spells of two or more
#' periods; and a period with no observation was skipped instead of
#' propagated.
#'
#' me_variance > 0 is required (hard stop): the bootstrap weights degenerate
#' when me_variance = 0 with n_obs < n_exo.
#'
#' @param model       dynhr_mod
#' @param data        Observation matrix (T x n_obs or n_obs x T)
#' @param prior_spec  Prior specification data.frame
#' @param obs_vars    Character vector of observed variable names
#' @param compiled    dynhr_compiled
#' @param specs       OBC spec list (from obc_parse_tags), or NULL to parse
#' @param me_variance Measurement error variance (must be > 0; default 1e-4)
#' @param N           Number of particles per evaluation (default 1000)
#' @param proposal    "bootstrap" (prior proposal) or "copf" (conditionally
#'   optimal proposal with regime verification and general-ratio fallback)
#' @param regime_guess For proposal = "copf": the starting guess of each
#'   particle's expected regime sequence, which is then iterated on the mean
#'   of the proposal: "ancestor" (default) starts from the sequence the
#'   resampled ancestor expected (all slack at t = 1); "slack" starts all
#'   slack. Pure variance reduction -- weights are valid for any guess.
#' @param seed        Integer RNG seed (NULL = not fixed; each call differs)
#' @param power       Power-posterior (generalised-Bayes) tempering exponent
#'   \eqn{\zeta}: \code{$logpost} becomes
#'   \eqn{\log p(\theta) + \zeta \cdot \log L(\theta)} while \code{$loglik}
#'   keeps the RAW (untempered) particle-filter estimate -- so PMMH's
#'   unbiasedness argument and any marginal-likelihood use of \code{$loglik}
#'   are unaffected. \code{NULL} (default) resolves the \code{power_posterior}
#'   package option ONCE, at factory time.
#' @return Function(theta) -> list(logpost, loglik, logprior)
#' @export
make_log_posterior_obc_ppf <- function(model, data, prior_spec, obs_vars,
                                        compiled, specs = NULL,
                                        me_variance  = 1e-4,
                                        N            = 1000L,
                                        proposal     = c("bootstrap", "copf"),
                                        regime_guess = c("ancestor", "slack"),
                                        seed         = NULL,
                                        power        = NULL) {

  ## Force promises (closure-capture safety; mirrors PKF and TPF factories)
  force(prior_spec); force(me_variance); force(N); force(seed)
  ## Resolve zeta ONCE here, not per draw (see .resolve_power_posterior).
  power <- .resolve_power_posterior(power, "make_log_posterior_obc_ppf")
  proposal     <- match.arg(proposal)
  regime_guess <- match.arg(regime_guess)

  ## A per-observable vector: refused with a classed error (the PPF/COPF
  ## weights implement H = me I only); an all-equal vector is the scalar.
  if (length(me_variance) > 1L)
    me_variance <- .kf_me_variance(me_variance, obs_vars,
                                   "make_log_posterior_obc_ppf",
                                   allow_vector = FALSE)
  ## Hard stop: me_variance = 0 degenerates weights for both bootstrap and COPF
  if (!is.numeric(me_variance) || length(me_variance) != 1L ||
      !is.finite(me_variance) || me_variance <= 0) {
    stop(
      "make_log_posterior_obc_ppf: 'me_variance' must be a finite positive scalar.\n",
      "Bootstrap PF weights p(y_t | ...) degenerate at me_variance = 0 when ",
      "n_obs < n_exo. Recommended: me_variance >= 1e-6."
    )
  }

  if (is.null(specs)) specs <- obc_parse_tags(model)

  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence)) {
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence
  }
  endo    <- model$var_names
  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx)))
    stop("obs_vars contains names not found in model$var_names: ",
         paste(obs_vars[is.na(obs_idx)], collapse = ", "))

  if (is.data.frame(data)) data <- as.matrix(data)
  Y <- if (nrow(data) == length(obs_vars)) data else t(data)

  ## Inner evaluator over the shared closure builder (R/posterior-closure.R):
  ## cold steady-state solve, always-eigen() stationarity guard, no system
  ## priors -- the OBC bootstrap/COPF particle filter is the only branch-
  ## specific part. `pass_dots = TRUE` returns the raw `function(theta, ...)`;
  ## the seeded wrapper below restores the public `function(theta)` signature.
  eval_one <- .make_posterior_closure(
    model, data, prior_spec, obs_vars, compiled,
    stationarity = "eigen",
    loglik_fn = function(sol, params, ss, theta, me_floor_check, ...) {
      ## FRESH regime_cache per draw (theta-dependent matrices)
      regime_cache <- new.env(parent = emptyenv(), hash = TRUE)
      obc_ensure_policy(0L, regime_cache, sol$sys, sol$dr, specs, obs_idx)

      pf <- ppf_likelihood(
        Y, sol$dr, regime_cache, sol$sys,
        model, params, obs_vars, specs,
        obs_idx      = obs_idx,
        N            = N,
        me_variance  = me_variance,
        proposal     = proposal,
        regime_guess = regime_guess,
        seed         = NULL   # the closure below already set the local seed
      )
      if (is.null(pf) || !is.finite(pf$loglik)) return(NULL)
      list(loglik = pf$loglik)
    },
    power          = power,
    warm_start     = FALSE,
    needs_me_floor = FALSE,
    pass_dots      = TRUE)

  ## ---- Closure: evaluated at each parameter draw -------------------------
  function(theta) {
    ## Fixed seed per theta for reproducibility (same seed = same loglik).
    ## LOCAL: deterministic in theta but leaves the CALLER's RNG stream
    ## untouched (an outer sampler must not replay the same proposals; see
    ## .with_local_seed in tpf-likelihood.R and NEWS 0.9.2.0003).
    .with_local_seed(seed, eval_one(theta))
  }
}
