## R/sampler-nuts-pooled.R
## --------------------------------------------------------------------------
## Multi-chain POOLED NUTS warmup (Lao 2026, arXiv:2607.23788, Sec. 3/5/7).
##
## M chains are warmed up in lockstep. Every warmup iteration each chain makes
## one NUTS transition (.nuts_transition(), sampler-nuts.R) at ONE shared step
## size; the M acceptance statistics are mean-pooled into ONE dual-averaging
## update ("replicated measurements for one controller time", Lao Sec. 7). At
## each mass-matrix-window endpoint of the dimension-derived schedule
## (.lowrank_warmup_windows(n_warmup, d, n_chains = M): 1-step init, first
## window ceil(8 (k_cap + 1) / M) per chain, 1.5x growth, final 15% step size
## only) the M per-chain buffers are POOLED with chain labels and ONE common
## metric is estimated from the per-chain-centred (within-chain, Lao eq. 5 W)
## draws -- and scores, for the Fisher estimators:
##   metric = "lowrank"     .lowrank_estimate(..., chain = labels): the
##                          diagonal -> low-rank-plus-diagonal promotion of the
##                          controller, now with pooled support N = M n, so
##                          both the rank cap floor(N / 8) - 1 and the BBP
##                          edge (1 + sqrt(d / N))^2 use the pooled count;
##   metric = "diagonal"    diag(W): pooled within-chain variance (the paper's
##                          pooled Welford-diagonal baseline);
##   metric = "fisher_diag" sqrt(var_W(x) / var_W(g)) (Seyboldt et al. 2026,
##                          Thm 2.2) on the pooled within-chain draws/scores.
## The step size is then re-found (per chain, geometric mean) and the dual
## averager reset, as in dynhr_nuts(). At the end of warmup the step size
## (dual-averaged) and the metric are FROZEN and identical for every chain;
## each chain then samples independently from its own RNG stream with
## dynhr_nuts(n_warmup = 1, step_size = eps, adapt_mass = FALSE, ...), so the
## post-warmup draws are ordinary fixed-kernel NUTS chains.
##
## Evidence recorded at each endpoint (Lao Sec. 3, 5.1):
##   * within/between split: Brooks-Gelman (1998) multivariate PSRF
##     (n-1)/n + (M+1)/M lambda_max(W^{-1} B/n) of the window;
##   * held-out score-position linearity: R^2 of the per-chain-centred scores
##     regressed on the per-chain-centred draws, every 5th pooled row held out
##     (exactly 1 for a Gaussian target).
## Persistent disagreement -- in the last two slow windows BOTH the MPSRF
## exceeds `pooled_control$disagreement` (default 1.5) AND lambda_max exceeds
## `pooled_control$excess` (default 4) times its iid-null scale
## (sqrt(d) + sqrt(M - 1))^2 / ((M - 1) n) (the Marchenko-Pastur edge of the
## chain-mean scatter; a fixed MPSRF cut alone fires on well-mixed short
## windows: 1.1-1.25 was measured on a unimodal d = 10 Gaussian) -- keeps the within-chain metric and advises a
## population / tempering sampler; poor linearity (R^2 below
## `pooled_control$linearity` in the last window) advises reparameterisation.
## Both are advisory (classed warnings + $advice); neither changes the kernel.
##
## Scope vs the paper: the controller's within-chain (W) route is implemented;
## the between-means low-rank route and the no-demotion latch are not (each
## window is re-estimated memorylessly, as dynhr_nuts(metric = "lowrank")),
## and the fixed diagnostic thresholds are dynhr's own (the paper does not
## state its thresholds).
## --------------------------------------------------------------------------


#' Validate and merge a named control list
#' @noRd
.nuts_pooled_ctrl <- function(user, defaults, what) {
  if (is.null(user)) return(defaults)
  if (!is.list(user) || (length(user) > 0L && is.null(names(user))) ||
      !all(names(user) %in% names(defaults))) {
    .dynhr_abort(paste0("`", what, "` must be a named list with elements in: ",
                        paste(names(defaults), collapse = ", "), "."),
                 class = "dynhr_error_invalid_argument")
  }
  defaults[names(user)] <- user
  defaults
}


#' Sampler-space log density / gradient closures (as built inside dynhr_nuts())
#' @noRd
.nuts_sampler_target <- function(log_post_fn, par_names, transform = NULL,
                                 grad_fn = NULL, grad_method = "forward") {
  target_fn <- if (!is.null(transform)) {
    make_transformed_logpost(log_post_fn, transform, include_jacobian = TRUE)
  } else {
    log_post_fn
  }
  lp_scalar <- function(theta) {
    names(theta) <- par_names
    res <- target_fn(theta)
    val <- if (is.list(res)) res$logpost else res
    if (!is.finite(val)) -1e300 else val
  }
  grad <- if (is.null(grad_fn)) {
    function(theta) .hmc_gradient(lp_scalar, theta, method = grad_method)
  } else if (!is.null(transform)) {
    make_transformed_grad(grad_fn, transform)
  } else {
    grad_fn
  }
  list(lp = lp_scalar, grad = grad)
}


#' Pool an n x d x M window array into an (n M) x d matrix + chain labels
#' @noRd
.nuts_pooled_stack <- function(arr, rows) {
  M <- dim(arr)[3L]
  d <- dim(arr)[2L]
  n <- length(rows)
  X <- do.call(rbind, lapply(seq_len(M), function(c)
    matrix(arr[rows, , c], nrow = n, ncol = d)))
  list(X = X, chain = rep(seq_len(M), each = n), n = n, M = M)
}


#' Subtract per-chain column means
#' @noRd
.nuts_pooled_centre <- function(X, chain) {
  means <- rowsum(X, chain) / as.numeric(table(chain))
  X - means[chain, , drop = FALSE]
}


#' One common metric from a pooled window (within-chain estimators)
#'
#' @return list(M_diag, M_inv_diag, M_inv, chol_M, rank, threshold, est)
#' @noRd
.nuts_pooled_metric <- function(X, G, chain, metric, lr_ctrl) {
  d   <- ncol(X)
  N   <- nrow(X)
  dof <- N - length(unique(chain))
  if (identical(metric, "lowrank")) {
    est <- .lowrank_estimate(X, G, chain = chain,
                             cutoff   = lr_ctrl$cutoff,
                             max_rank = lr_ctrl$max_rank,
                             gamma    = lr_ctrl$gamma)
    return(list(M_diag = rep(1, d), M_inv_diag = NULL,
                M_inv = est$metric, chol_M = est$metric,
                rank = est$rank, threshold = est$threshold, est = est))
  }
  xc    <- .nuts_pooled_centre(X, chain)
  var_x <- colSums(xc^2) / dof
  if (identical(metric, "fisher_diag")) {
    gc    <- .nuts_pooled_centre(G, chain)
    var_g <- colSums(gc^2) / dof
    ok_x  <- is.finite(var_x) & var_x >= 1e-12
    ok_g  <- is.finite(var_g) & var_g >= 1e-12
    vars  <- sqrt(var_x / var_g)
    bad   <- !(ok_x & ok_g) | !is.finite(vars)
    vars[bad] <- ifelse(ok_x[bad], var_x[bad], 1)
  } else {
    vars <- var_x
  }
  vars[!is.finite(vars) | vars < 1e-12] <- 1
  ## Stan's rule: the INVERSE mass is the (within-chain) posterior variance.
  list(M_diag = 1 / vars, M_inv_diag = vars, M_inv = NULL, chol_M = NULL,
       rank = NA_integer_, threshold = NA_real_, est = NULL)
}


#' Within/between and score-position evidence of one pooled window
#'
#' @return list(mpsrf, between_excess (lambda_1 over its iid-null scale),
#'   linearity_r2)
#' @noRd
.nuts_pooled_evidence <- function(X, G, chain, n, M) {
  d <- ncol(X)
  mpsrf <- NA_real_
  if (M >= 2L && n >= 2L) {
    xbar <- rowsum(X, chain) / n
    xc   <- X - xbar[chain, , drop = FALSE]
    W    <- crossprod(xc) / (M * (n - 1L))
    Bn   <- stats::cov(xbar)                    # B / n (Lao eq. 5)
    ew   <- eigen((W + t(W)) / 2, symmetric = TRUE)
    ev   <- pmax(ew$values, max(ew$values) * 1e-12, .Machine$double.xmin)
    W_mh <- ew$vectors %*% ((1 / sqrt(ev)) * t(ew$vectors))
    A    <- W_mh %*% Bn %*% W_mh
    lam1 <- max(eigen((A + t(A)) / 2, symmetric = TRUE, only.values = TRUE)$values)
    mpsrf <- (n - 1) / n + (M + 1) / M * lam1
    ## iid-null scale of lambda_1: B/n ~ W_{M-1}(Sigma)/((M-1) n), whose
    ## whitened top eigenvalue sits near the Marchenko-Pastur edge
    ## (sqrt(d) + sqrt(M-1))^2 / ((M-1) n).
    null_lam <- (sqrt(d) + sqrt(M - 1))^2 / ((M - 1) * n)
    excess   <- lam1 / null_lam
  } else {
    excess <- NA_real_
  }
  r2 <- NA_real_
  if (!is.null(G)) {
    xc <- .nuts_pooled_centre(X, chain)
    gc <- .nuts_pooled_centre(G, chain)
    N  <- nrow(X)
    test  <- seq.int(5L, N, by = 5L)
    train <- setdiff(seq_len(N), test)
    if (length(test) >= 2L && length(train) > d + 1L) {
      B <- qr.coef(qr(xc[train, , drop = FALSE]), gc[train, , drop = FALSE])
      B[!is.finite(B)] <- 0
      res <- gc[test, , drop = FALSE] - xc[test, , drop = FALSE] %*% B
      gt  <- gc[test, , drop = FALSE]
      sst <- colSums(sweep(gt, 2L, colMeans(gt))^2)
      ok  <- is.finite(sst) & sst > 0
      if (any(ok)) r2 <- mean(1 - colSums(res^2)[ok] / sst[ok])
    }
  }
  list(mpsrf = mpsrf, between_excess = excess, linearity_r2 = r2)
}


#' Pooled multi-chain NUTS warmup (lockstep; shared step size and metric)
#'
#' @param log_post_fn Log posterior (theta space), returning a list with
#'   \code{$logpost} or a scalar.
#' @param starts List of M named THETA-space starting vectors.
#' @param rng_states List of M \code{.Random.seed} vectors: chain m's stream.
#' @param n_warmup Warmup iterations per chain (row 1 = the start, as in
#'   \code{dynhr_nuts()}).
#' @param metric \code{"lowrank"}, \code{"diagonal"} or \code{"fisher_diag"}.
#' @param mass_diag Optional initial MASS diagonal (e.g. 1 / diag(V_mode)).
#' @return list(states, lp, rng, step_size, M_diag, M_inv_diag, M_inv, chol_M,
#'   metric, history, advice, warm_theta (n_warmup x d x M, theta space),
#'   warm_lp (n_warmup x M), warm_depth, warm_div, n_grad (per chain),
#'   starts, rng_init, par_names, n_warmup). The caller's RNG state is
#'   restored on exit.
#' @noRd
.nuts_pooled_warmup <- function(log_post_fn, starts, rng_states,
                                n_warmup      = 1000L,
                                metric        = c("lowrank", "diagonal", "fisher_diag"),
                                grad_fn       = NULL,
                                grad_method   = "forward",
                                max_treedepth = 8L,
                                target_accept = 0.80,
                                delta_max     = 1000,
                                mass_diag     = NULL,
                                transform     = NULL,
                                lowrank_control = NULL,
                                pooled_control  = NULL,
                                verbose       = FALSE) {
  metric <- match.arg(metric)
  n_warmup <- as.integer(n_warmup)
  M <- length(starts)
  if (M < 1L || length(rng_states) != M || n_warmup < 1L) {
    .dynhr_abort("`starts` and `rng_states` must be non-empty lists of equal ",
                 "length and `n_warmup` >= 1.",
                 class = "dynhr_error_invalid_argument")
  }
  lr_ctrl <- .nuts_pooled_ctrl(lowrank_control,
                               list(cutoff = 2, max_rank = NULL, gamma = 1e-5),
                               "lowrank_control")
  pc <- .nuts_pooled_ctrl(pooled_control,
                          list(disagreement = 1.5, excess = 4, linearity = 0.5),
                          "pooled_control")
  par_names <- names(starts[[1L]])
  d <- length(par_names)
  tg <- .nuts_sampler_target(log_post_fn, par_names, transform, grad_fn, grad_method)
  trace_lp <- function(state, lp) {
    if (!is.null(transform)) lp - transform$log_jacobian(state) else lp
  }
  to_theta <- function(state) {
    if (!is.null(transform)) transform$to_constrained(state) else state
  }

  ## Per-chain RNG streams are swapped in and out of the global state; the
  ## caller's state is put back on exit.
  genv    <- globalenv()
  rng_old <- get0(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (is.null(rng_old)) {
      if (exists(".Random.seed", envir = genv, inherits = FALSE))
        rm(".Random.seed", envir = genv)
    } else {
      assign(".Random.seed", rng_old, envir = genv)
    }
  }, add = TRUE)
  rng <- rng_states

  states <- lapply(starts, function(s) {
    st <- if (!is.null(transform)) transform$to_unconstrained(s) else s
    st <- as.numeric(st)
    names(st) <- par_names
    st
  })
  ## Fused value + gradient (W92, .hmc_fused_target()): checked against
  ## log_post_fn at chain 1's start; then every chain's start value and every
  ## leaf come from the fused function, gradients carried with the states.
  fz <- .hmc_fused_target(log_post_fn, grad_fn, states[[1L]], par_names,
                          transform = transform, verbose = verbose,
                          sampler = "NUTS pooled warmup")
  vg_fn <- if (is.null(fz)) NULL else fz$vg
  g_c   <- vector("list", M)
  if (is.null(vg_fn)) {
    lp_c <- vapply(states, tg$lp, numeric(1))
    ## W94: the separate-call path carries the gradients with the states too
    ## (taken once here; the step-size search and every transition reuse them)
    g_c  <- lapply(states, tg$grad)
  } else {
    ev   <- c(list(list(lp = fz$lp0, grad = fz$g0)),
              lapply(states[-1L], vg_fn))
    lp_c <- vapply(ev, function(e) e$lp, numeric(1))
    g_c  <- lapply(ev, function(e) e$grad)
  }

  if (!is.null(mass_diag) && length(mass_diag) == d) {
    M_diag     <- pmax(as.numeric(mass_diag), 1e-12)
    M_inv_diag <- 1 / M_diag
  } else {
    M_diag     <- rep(1, d)
    M_inv_diag <- rep(1, d)
  }
  M_inv  <- NULL
  chol_M <- NULL

  windows <- .lowrank_warmup_windows(n_warmup, d, n_chains = M)
  slow    <- windows[windows$type == "slow", , drop = FALSE]
  use_scores <- metric %in% c("lowrank", "fisher_diag")

  warm_state <- array(NA_real_, c(n_warmup, d, M))
  warm_theta <- array(NA_real_, c(n_warmup, d, M),
                      dimnames = list(NULL, par_names, NULL))
  warm_grad  <- if (use_scores) array(NA_real_, c(n_warmup, d, M)) else NULL
  warm_lp    <- matrix(NA_real_, n_warmup, M)
  warm_depth <- matrix(0L, n_warmup, M)
  warm_div   <- matrix(FALSE, n_warmup, M)
  n_grad     <- integer(M)
  for (c in seq_len(M)) {
    warm_state[1L, , c] <- states[[c]]
    warm_theta[1L, , c] <- to_theta(states[[c]])
    warm_lp[1L, c]      <- trace_lp(states[[c]], lp_c[c])
  }

  ## Shared initial step size: each chain's Hoffman-Gelman search (own RNG
  ## stream), combined by the geometric mean.
  find_eps <- function() {
    e <- numeric(M)
    for (c in seq_len(M)) {
      assign(".Random.seed", rng[[c]], envir = genv)
      e[c] <- .hmc_find_stepsize(states[[c]], tg$lp, tg$grad, M_inv_diag, M_diag,
                                 M_inv = M_inv, chol_M = chol_M,
                                 vg_fn = vg_fn, lp0 = lp_c[c], g0 = g_c[[c]])
      rng[[c]] <- get(".Random.seed", envir = genv)
    }
    list(eps = exp(mean(log(e))), rng = rng)
  }
  fe    <- find_eps()
  rng   <- fe$rng
  eps0  <- fe$eps
  mu      <- log(10 * eps0)
  eps_bar <- 1
  H_bar   <- 0
  gamma_da <- 0.05
  t0_da    <- 10
  kappa_da <- 0.75
  eps_m    <- eps0
  da_m     <- 0L

  history <- list()
  win_i   <- 1L
  n_slow  <- nrow(slow)
  if (verbose) .dynhr_inform(sprintf(
    "NUTS pooled warmup: %d chains, metric %s, %d slow windows, initial step_size = %.4e",
    M, metric, n_slow, eps0))

  for (it in seq_len(n_warmup)[-1L]) {
    in_slow <- win_i <= n_slow && it >= slow$start[win_i]
    alpha <- numeric(M)
    for (c in seq_len(M)) {
      assign(".Random.seed", rng[[c]], envir = genv)
      tr <- .nuts_transition(states[[c]], lp_c[c], eps_m, tg$lp, tg$grad,
                             M_diag, M_inv_diag, M_inv, chol_M,
                             max_treedepth, delta_max,
                             g_curr = g_c[[c]], vg_fn = vg_fn)
      th <- tr$theta
      names(th) <- par_names
      states[[c]] <- th
      lp_c[c]     <- if (is.null(vg_fn)) tg$lp(th) else tr$lp
      g_c[[c]]    <- tr$g   # carried with the state (W92 fused, W94 both)
      warm_state[it, , c] <- th
      warm_theta[it, , c] <- to_theta(th)
      warm_lp[it, c]      <- trace_lp(th, lp_c[c])
      warm_depth[it, c]   <- tr$depth
      warm_div[it, c]     <- tr$divergent
      n_grad[c]           <- n_grad[c] + tr$n_leaves
      alpha[c] <- if (tr$n_alpha > 0) tr$alpha_sum / tr$n_alpha else 0
      if (use_scores && in_slow) {
        warm_grad[it, , c] <- g_c[[c]]   # carried with the state: no extra call
      }
      rng[[c]] <- get(".Random.seed", envir = genv)
    }

    ## ONE dual-averaging update per controller time: the M acceptance
    ## statistics at the shared step size are mean-pooled (Lao Sec. 7).
    da_m <- da_m + 1L
    alpha_bar <- mean(alpha)
    w <- 1 / (da_m + t0_da)
    H_bar <- (1 - w) * H_bar + w * (target_accept - alpha_bar)
    log_eps_m <- mu - (sqrt(da_m) / gamma_da) * H_bar
    eps_m <- exp(log_eps_m)
    m_kappa <- da_m^(-kappa_da)
    eps_bar <- exp(m_kappa * log_eps_m + (1 - m_kappa) * log(eps_bar))

    ## Mass-matrix-window endpoint: pool the M buffers, estimate ONE metric,
    ## record the evidence, re-find the step size, reset dual averaging.
    if (win_i <= n_slow && it == slow$end[win_i]) {
      rows <- seq.int(slow$start[win_i], slow$end[win_i])
      px <- .nuts_pooled_stack(warm_state, rows)
      pg <- if (use_scores) .nuts_pooled_stack(warm_grad, rows)$X else NULL
      if (nrow(px$X) - M >= 2L) {
        pm <- .nuts_pooled_metric(px$X, pg, px$chain, metric, lr_ctrl)
        M_diag     <- pm$M_diag
        M_inv_diag <- pm$M_inv_diag
        M_inv      <- pm$M_inv
        chol_M     <- pm$chol_M
        ev <- .nuts_pooled_evidence(px$X, pg, px$chain, px$n, M)
        fe   <- find_eps()
        rng  <- fe$rng
        eps0 <- fe$eps
        mu      <- log(10 * eps0)
        eps_bar <- 1
        H_bar   <- 0
        eps_m   <- eps0
        da_m    <- 0L
        history[[length(history) + 1L]] <- data.frame(
          start = rows[1L], end = rows[length(rows)], n = px$n, N = nrow(px$X),
          rank = pm$rank, threshold = pm$threshold, mpsrf = ev$mpsrf,
          between_excess = ev$between_excess,
          linearity_r2 = ev$linearity_r2, step_size = eps0)
        if (verbose) .dynhr_inform(sprintf(
          "NUTS pooled: window [%d,%d] (N = %d) -> %s%s, MPSRF %.3f, reset step_size = %.4e",
          rows[1L], rows[length(rows)], nrow(px$X), metric,
          if (is.na(pm$rank)) "" else sprintf(" rank %d", pm$rank),
          ev$mpsrf, eps0))
      }
      win_i <- win_i + 1L
    }
  }
  ## Freeze: the dual-averaged step size (dynhr_nuts() convention).
  step_size <- if (n_warmup >= 2L) eps_bar else eps0

  history <- if (length(history)) do.call(rbind, history) else NULL
  mp   <- if (is.null(history)) numeric(0) else utils::tail(history$mpsrf, 2L)
  bx   <- if (is.null(history)) numeric(0) else utils::tail(history$between_excess, 2L)
  persistent <- length(mp) > 0L && all(is.finite(mp)) && all(is.finite(bx)) &&
    all(mp > pc$disagreement) && all(bx > pc$excess)
  lin  <- if (is.null(history)) NA_real_ else utils::tail(history$linearity_r2, 1L)
  reparam <- length(lin) == 1L && is.finite(lin) && lin < pc$linearity
  advice <- list(
    persistent_disagreement = persistent,
    handoff        = if (persistent) "population" else "none",
    reparameterize = reparam,
    mpsrf          = if (length(mp)) mp[length(mp)] else NA_real_,
    linearity_r2   = lin)
  if (persistent) {
    .dynhr_warn(sprintf(paste0(
      "Pooled NUTS warmup: persistent within/between-chain disagreement ",
      "(multivariate PSRF %.3g > %g in the last warmup windows). One constant ",
      "metric does not describe the regions the chains occupy; the within-chain ",
      "metric was kept. Consider a population or tempering sampler (SMC, DIME) ",
      "for regional exploration."), advice$mpsrf, pc$disagreement),
      class = "dynhr_warning_pooled_disagreement")
  }
  if (reparam) {
    .dynhr_warn(sprintf(paste0(
      "Pooled NUTS warmup: poor held-out score-position linearity (R^2 = %.3g < %g): ",
      "a constant linear preconditioner fits the local geometry badly; consider ",
      "reparameterising (funnels, scale coupling)."), lin, pc$linearity),
      class = "dynhr_warning_pooled_nonlinear")
  }

  list(states = states, lp = lp_c, rng = rng, step_size = step_size,
       M_diag = M_diag, M_inv_diag = M_inv_diag, M_inv = M_inv, chol_M = chol_M,
       metric = metric, history = history, advice = advice,
       warm_theta = warm_theta, warm_lp = warm_lp, warm_depth = warm_depth,
       warm_div = warm_div, n_grad = n_grad, starts = starts,
       rng_init = rng_states, par_names = par_names, n_warmup = n_warmup)
}


#' Post-warmup sampling of one chain at the frozen pooled adaptation
#'
#' Restores chain \code{ch}'s RNG stream and runs \code{dynhr_nuts()} from the
#' chain's final warmup state with the frozen step size and metric
#' (\code{n_warmup = 1}: no adaptation, every retained row a fixed-kernel
#' transition). The warmup rows are prepended to \code{full_chain} /
#' \code{logpost_trace} and warmup gradients added to \code{n_grad_evals}, so
#' the result has the usual \code{dynhr_nuts()} shape.
#' @noRd
.nuts_pooled_sample_chain <- function(log_post_fn, warm, ch, n_draws,
                                      grad_fn = NULL, grad_method = "forward",
                                      transform = NULL, max_treedepth = 8L,
                                      target_accept = 0.80, delta_max = 1000,
                                      progressor = NULL) {
  par_names <- warm$par_names
  assign(".Random.seed", warm$rng[[ch]], envir = globalenv())
  st  <- warm$states[[ch]]
  th0 <- if (!is.null(transform)) transform$to_constrained(st) else st
  th0 <- as.numeric(th0)
  names(th0) <- par_names
  is_lr <- inherits(warm$M_inv, "dynhr_lowrank_metric")
  res <- dynhr_nuts(log_post_fn, th0,
                    n_draws       = n_draws,
                    n_warmup      = 1L,
                    step_size     = warm$step_size,
                    max_treedepth = max_treedepth,
                    target_accept = target_accept,
                    adapt_mass    = FALSE,
                    grad_fn       = grad_fn,
                    grad_method   = grad_method,
                    delta_max     = delta_max,
                    mass_diag     = if (is_lr) NULL else warm$M_diag,
                    M_inv         = warm$M_inv,
                    chol_M        = warm$chol_M,
                    verbose       = FALSE,
                    progressor    = progressor,
                    chain_id      = ch,
                    transform     = transform)
  n_warmup <- warm$n_warmup
  wt <- matrix(warm$warm_theta[, , ch], nrow = n_warmup,
               dimnames = list(NULL, par_names))
  if (!is.null(res$chain)) res$full_chain <- rbind(wt, res$chain)
  res$logpost_trace <- c(warm$warm_lp[, ch], res$post_logpost)
  res$n_burn        <- as.integer(n_warmup)
  res$n_grad_evals  <- res$n_grad_evals + warm$n_grad[ch]
  res$pooled <- list(n_chains  = length(warm$states),
                     metric    = warm$metric,
                     step_size = warm$step_size,
                     history   = warm$history,
                     advice    = warm$advice,
                     n_warmup_divergent = sum(warm$warm_div[, ch]))
  if (is_lr) {
    rk <- warm$history$rank
    res$lowrank <- list(
      rank     = warm$M_inv$rank,
      promoted = any(!is.na(rk) & rk >= 1L),
      metric   = warm$M_inv,
      history  = warm$history)
  }
  res
}


#' Per-chain summary table (run_nuts_mirai's chain_stats shape)
#' @noRd
.nuts_pooled_chain_stats <- function(chains, elapsed_min = NULL) {
  do.call(rbind, lapply(seq_along(chains), function(ch) {
    r <- chains[[ch]]
    ## accept_rate: the chain's mean NUTS acceptance statistic
    ## (dynhr_nuts()$acceptance_rate == $accept_stat since 0.9.3.127)
    data.frame(chain          = ch,
               accept_rate    = r$acceptance_rate %||% NA_real_,
               final_logpost  = utils::tail(r$post_logpost, 1L),
               n_divergent    = r$n_divergent %||% NA_integer_,
               mean_treedepth = r$mean_treedepth %||% NA_real_,
               elapsed_min    = if (is.null(elapsed_min)) NA_real_ else elapsed_min[ch],
               stringsAsFactors = FALSE)
  }))
}


#' Multi-chain NUTS with POOLED warmup adaptation (Lao 2026)
#'
#' Runs \code{n_chains} NUTS chains whose warmup is pooled: one shared step
#' size (dual averaging on the mean acceptance statistic of all chains) and
#' ONE common metric estimated at each window endpoint of Lao's
#' dimension-derived schedule from all chains' per-chain-centred warmup draws
#' (and scores for \code{"lowrank"} / \code{"fisher_diag"}). Step size and
#' metric are frozen at the end of warmup and are identical for every chain;
#' post-warmup each chain is an ordinary fixed-kernel NUTS chain on its own
#' RNG stream. See the header of R/sampler-nuts-pooled.R for the estimators,
#' the evidence diagnostics and the scope relative to the paper.
#'
#' Warmup runs the chains in lockstep in THIS process (serially over chains);
#' the post-warmup chains are independent and run serially here, or in
#' parallel via \code{run_nuts_mirai(adapt = "pooled")}.
#'
#' @param log_post_fn Log posterior (theta space).
#' @param theta_init Named start vector (all chains), an \code{n_chains x d}
#'   matrix with column names (one row per chain), or a list of named vectors.
#' @param n_chains Number of chains (ignored when \code{theta_init} is a
#'   matrix or list).
#' @param n_draws,n_warmup Post-warmup draws / warmup iterations per chain.
#' @param metric \code{"lowrank"} (default: the controller's diagonal ->
#'   low-rank-plus-diagonal promotion), \code{"diagonal"} or
#'   \code{"fisher_diag"}.
#' @param seeds Optional integer vector (one per chain) seeding each chain's
#'   RNG stream; NULL draws them from the caller's stream (so
#'   \code{set.seed()} before the call reproduces the run). The caller's
#'   stream is otherwise left untouched.
#' @param lowrank_control As in \code{dynhr_nuts()}.
#' @param pooled_control Optional list: \code{disagreement} (MPSRF above which,
#'   together with \code{excess} (default 4: lambda_max over its iid-null
#'   scale), the last two windows count as persistent disagreement, default
#'   1.5) and
#'   \code{linearity} (held-out score-position R^2 below which
#'   reparameterisation is advised, default 0.5).
#' @return list(chains (per-chain \code{dynhr_nuts()}-shaped results, each
#'   with a \code{$pooled} element), chain_stats, adaptation (step_size,
#'   metric, inverse_mass, history, advice, seeds)).
#' @noRd
dynhr_nuts_pooled <- function(log_post_fn, theta_init,
                              n_chains      = 4L,
                              n_draws       = 1000L,
                              n_warmup      = 1000L,
                              metric        = c("lowrank", "diagonal", "fisher_diag"),
                              grad_fn       = NULL,
                              grad_method   = "forward",
                              max_treedepth = 8L,
                              target_accept = 0.80,
                              delta_max     = 1000,
                              mass_diag     = NULL,
                              seeds         = NULL,
                              transform     = NULL,
                              lowrank_control = NULL,
                              pooled_control  = NULL,
                              verbose       = FALSE) {
  stopifnot(is.function(log_post_fn))
  metric <- match.arg(metric)
  starts <- if (is.list(theta_init)) {
    theta_init
  } else if (is.matrix(theta_init)) {
    lapply(seq_len(nrow(theta_init)), function(i) theta_init[i, ])
  } else {
    rep(list(theta_init), as.integer(n_chains))
  }
  M <- length(starts)
  if (is.null(names(starts[[1L]]))) {
    .dynhr_abort("`theta_init` must carry parameter names.",
                 class = "dynhr_error_invalid_argument")
  }
  if (is.null(seeds)) {
    seeds <- sample.int(.Machine$integer.max, M)
  } else if (length(seeds) != M) {
    .dynhr_abort(sprintf("`seeds` must have one entry per chain (%d).", M),
                 class = "dynhr_error_invalid_argument")
  }
  genv    <- globalenv()
  rng_old <- get0(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (is.null(rng_old)) {
      if (exists(".Random.seed", envir = genv, inherits = FALSE))
        rm(".Random.seed", envir = genv)
    } else {
      assign(".Random.seed", rng_old, envir = genv)
    }
  }, add = TRUE)
  rng_init <- lapply(seeds, function(s) {
    set.seed(s)
    get(".Random.seed", envir = genv)
  })

  warm <- .nuts_pooled_warmup(log_post_fn, starts, rng_init,
                              n_warmup = n_warmup, metric = metric,
                              grad_fn = grad_fn, grad_method = grad_method,
                              max_treedepth = max_treedepth,
                              target_accept = target_accept,
                              delta_max = delta_max, mass_diag = mass_diag,
                              transform = transform,
                              lowrank_control = lowrank_control,
                              pooled_control = pooled_control,
                              verbose = verbose)
  chains <- lapply(seq_len(M), function(ch)
    .nuts_pooled_sample_chain(log_post_fn, warm, ch, n_draws,
                              grad_fn = grad_fn, grad_method = grad_method,
                              transform = transform,
                              max_treedepth = max_treedepth,
                              target_accept = target_accept,
                              delta_max = delta_max))
  list(chains      = chains,
       chain_stats = .nuts_pooled_chain_stats(chains),
       adaptation  = .nuts_pooled_summary(warm, seeds))
}


#' Frozen-adaptation summary shared by the serial and mirai paths
#' @noRd
.nuts_pooled_summary <- function(warm, seeds = NULL) {
  list(step_size    = warm$step_size,
       metric       = warm$metric,
       inverse_mass = if (inherits(warm$M_inv, "dynhr_lowrank_metric"))
                        .lowrank_dense_inv(warm$M_inv) else warm$M_inv_diag,
       history      = warm$history,
       advice       = warm$advice,
       seeds        = seeds)
}


# ============================================================================
# mirai path: run_nuts_mirai(adapt = "pooled")
# ============================================================================

#' Chain start + RNG stream exactly as run_nuts_mirai's independent chain task
#'
#' \code{RNGkind("Mersenne-Twister", "Inversion", "Rejection")},
#' \code{set.seed(seed_base + ch)}; chain 1 starts at \code{theta_mode},
#' chains 2..N at \code{theta_mode + 0.5 chol(Sigma_prop)' z} (eta space
#' under \code{transform}; clamped to the prior box otherwise), falling back to
#' the mode when the log posterior there is not finite. Returns the start and
#' the RNG state after the dispersion draw, so the pooled chains start from
#' the same points and streams as the independent ones.
#' @noRd
.nuts_pooled_chain_start <- function(ch, log_post_fn, theta_mode, Sigma_prop,
                                     prior_spec, transform, seed_base) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed_base + ch)
  if (ch == 1L) {
    th0 <- theta_mode
  } else if (!is.null(transform)) {
    eta_mode <- transform$to_unconstrained(theta_mode)
    L    <- t(chol(Sigma_prop))
    z    <- rnorm(length(theta_mode))
    eta0 <- eta_mode + 0.5 * as.numeric(L %*% z)
    names(eta0) <- names(theta_mode)
    th0 <- transform$to_constrained(eta0)
  } else {
    L   <- t(chol(Sigma_prop))
    z   <- rnorm(length(theta_mode))
    th0 <- theta_mode + 0.5 * as.numeric(L %*% z)
    names(th0) <- names(theta_mode)
    for (i in seq_along(th0)) {
      th0[i] <- max(th0[i], prior_spec$lower[i] + 1e-6)
      th0[i] <- min(th0[i], prior_spec$upper[i] - 1e-6)
    }
  }
  lp0 <- log_post_fn(th0)
  lp0 <- if (is.list(lp0)) lp0$logpost else lp0
  if (!is.finite(lp0)) th0 <- theta_mode
  list(theta = th0, rng = get(".Random.seed", envir = globalenv()))
}


#' Pooled warmup as ONE daemon task (run_nuts_mirai(adapt = "pooled"))
#' @noRd
.nuts_pooled_warmup_task <- function(log_post_fn, grad_fn, theta_mode,
                                     Sigma_prop, prior_spec, transform,
                                     n_chains, seed_base, n_warmup, metric,
                                     mass_diag, max_treedepth, target_accept) {
  st <- lapply(seq_len(n_chains), function(ch)
    .nuts_pooled_chain_start(ch, log_post_fn, theta_mode, Sigma_prop,
                             prior_spec, transform, seed_base))
  .nuts_pooled_warmup(log_post_fn,
                      starts     = lapply(st, `[[`, "theta"),
                      rng_states = lapply(st, `[[`, "rng"),
                      n_warmup = n_warmup, metric = metric, grad_fn = grad_fn,
                      max_treedepth = max_treedepth,
                      target_accept = target_accept,
                      mass_diag = mass_diag, transform = transform)
}


#' Daemon-side gradient closure for the pooled tasks (mirrors chain_task)
#' @noRd
.nuts_pooled_worker_grad <- function(analytic_grad, prior_spec, obs_names,
                                     me_variance, me_extra, shock_scale,
                                     grad_method, likelihood, freq_band,
                                     transform, fuse = TRUE,
                                     system_priors = NULL) {
  if (!isTRUE(analytic_grad)) return(NULL)
  wm <- get0(".worker_model", envir = globalenv(), inherits = FALSE)
  wc <- get0(".worker_cm",    envir = globalenv(), inherits = FALSE)
  wy <- get0(".worker_Y",     envir = globalenv(), inherits = FALSE)
  if (is.null(wm) || is.null(wc) || is.null(wy)) return(NULL)
  ## lik_init: the init .worker_lp was built with (.mirai_pool_init).
  li <- get0(".worker_lik_init", envir = globalenv(), inherits = FALSE) %||%
    "auto"
  g <- make_posterior_grad(wm, wy, prior_spec, obs_names, wc,
                           me_variance = me_variance, me_extra = me_extra,
                           shock_scale = shock_scale, grad_method = grad_method,
                           likelihood = likelihood, freq_band = freq_band,
                           lik_init = li, system_priors = system_priors)
  ## THETA-space: .nuts_sampler_target() (warmup) and dynhr_nuts() (the
  ## frozen chains) both apply the eta chain rule when `transform` is set.
  ## W86: wrapping here as well applied it twice.
  ## W94: `system_priors` is the context's, the one .worker_lp carries
  ## (.mirai_pool_init), so this is the gradient -- and the fused value -- of
  ## the sampled target. (W92 dropped the fused value instead, and the
  ## gradient omitted the system prior.) `fuse = FALSE` still drops it.
  if (!isTRUE(fuse)) attr(g, "logpost_grad") <- NULL
  g
}


#' run_nuts_mirai(adapt = "pooled"): pooled warmup on one daemon, then the
#' frozen chains in parallel (one task per chain)
#'
#' The daemon pool must already be provisioned (\code{.worker_lp} bound).
#' Daemons run the INSTALLED dynhr, whose namespace must contain this file.
#' @return list(raw (per-chain task results, run_nuts_mirai's shape),
#'   adaptation)
#' @noRd
.nuts_pooled_mirai <- function(n_chains, n_draws, n_warmup, metric,
                               theta_mode, Sigma_prop, prior_spec, obs_names,
                               transform, seed_base, mass_diag,
                               max_treedepth, target_accept,
                               analytic_grad, grad_method, me_variance,
                               me_extra, shock_scale, likelihood, freq_band,
                               fuse = TRUE, system_priors = NULL) {
  grad_args <- list(analytic_grad = analytic_grad, prior_spec = prior_spec,
                    obs_names = obs_names, me_variance = me_variance,
                    me_extra = me_extra, shock_scale = shock_scale,
                    grad_method = grad_method, likelihood = likelihood,
                    freq_band = freq_band, transform = transform,
                    fuse = fuse, system_priors = system_priors)
  warm_h <- mirai::mirai({
    lp <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
    wg <- utils::getFromNamespace(".nuts_pooled_worker_grad", "dynhr")
    wt <- utils::getFromNamespace(".nuts_pooled_warmup_task", "dynhr")
    wt(lp, do.call(wg, grad_args), theta_mode, Sigma_prop, prior_spec,
       transform, n_chains, seed_base, n_warmup, metric, mass_diag,
       max_treedepth, target_accept)
  }, .args = list(grad_args = grad_args, theta_mode = theta_mode,
                  Sigma_prop = Sigma_prop, prior_spec = prior_spec,
                  transform = transform, n_chains = n_chains,
                  seed_base = seed_base, n_warmup = n_warmup, metric = metric,
                  mass_diag = mass_diag, max_treedepth = max_treedepth,
                  target_accept = target_accept))
  warm <- warm_h[]
  if (inherits(warm, "miraiError") || inherits(warm, "errorValue")) {
    .dynhr_abort("Pooled NUTS warmup failed on the daemon: ", as.character(warm),
                 class = "dynhr_error_sampler")
  }
  .dynhr_cat(sprintf("  Pooled warmup done (step_size %.4e, metric %s%s)\n",
                     warm$step_size, warm$metric,
                     if (inherits(warm$M_inv, "dynhr_lowrank_metric"))
                       sprintf(" rank %d", warm$M_inv$rank) else ""))

  sample_task <- function(ch) {
    lp <- get0(".worker_lp", envir = globalenv(), inherits = FALSE)
    wg <- utils::getFromNamespace(".nuts_pooled_worker_grad", "dynhr")
    sc <- utils::getFromNamespace(".nuts_pooled_sample_chain", "dynhr")
    t0 <- proc.time()[["elapsed"]]
    res <- sc(lp, warm, ch, n_draws, grad_fn = do.call(wg, grad_args),
              transform = grad_args$transform, max_treedepth = max_treedepth,
              target_accept = target_accept)
    list(chain_id = ch, result = res,
         elapsed_min = (proc.time()[["elapsed"]] - t0) / 60)
  }
  raw <- mirai::mirai_map(seq_len(n_chains), sample_task,
                          warm = warm, grad_args = grad_args, n_draws = n_draws,
                          max_treedepth = max_treedepth,
                          target_accept = target_accept)[]
  list(raw = raw, adaptation = .nuts_pooled_summary(warm))
}
