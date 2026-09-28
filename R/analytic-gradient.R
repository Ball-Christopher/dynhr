## R/analytic-gradient.R
## --------------------------------------------------------------------------
## Analytic gradient of the log-posterior for HMC / NUTS.
##
## The log-posterior is  log p(theta|Y) = loglik(theta) + logprior(theta).
##
##   * logprior: fully analytic score (.dlog_prior), per distribution.
##   * loglik:   the first-order decision rule is certainty-equivalent, so the
##               SHOCK-STD parameters (stderr -> sig_*) leave the state-space
##               matrices TT,RR,ZZ,DD unchanged and move only the shock
##               covariance Sigma_e = diag(sig^2). Their loglik gradient is
##               therefore analytic via a Kalman SCORE recursion that reuses a
##               SINGLE model solve (.kf_loglik_score_sigma).
##
##               The remaining (non-sigma) parameters DO change the decision
##               rule; the package has no analytic d(dr)/d(theta), so those
##               entries fall back to a relative-step central difference that
##               perturbs ONLY that parameter.
##
## make_posterior_grad() assembles these into a grad_fn(theta) suitable for the
## `grad_fn` argument of dynhr_nuts() / dynhr_hmc().
## --------------------------------------------------------------------------


# ============================================================================
# Analytic log-prior score  d/dx log p(x)
# ============================================================================

#' Derivative of a single parameter's log-prior density w.r.t. its value.
#' Mirrors the switch in log_prior() (R/prior-density.R) term for term.
#' @return scalar d/dx log-density (0 outside support / for flat priors).
#' @noRd
.dlog_prior_density1 <- function(x, dist, p1, p2, p3 = NA_real_,
                                 p4 = NA_real_) {
  ## Dynare p3/p4 (brief 23 A3): generalised beta on [p3, p4]; p3 is a SHIFT
  ## for gamma and both inverse gammas. Mirrors .lp_dist1() exactly.
  s <- if (is.na(p3)) 0 else p3
  switch(.normalize_dist(dist),
    "inv_gamma" =, "inv_gamma1" = {
      y <- x - s
      if (y <= 0) return(0)
      ## A12: share the ONE (alpha, theta) mapping, including the sd = Inf
      ## limit -- `(alpha - 1) * (p2^2 + p1^2)` was 0 * Inf = NaN there.
      ps <- .ig1_params(p1 - s, p2)
      # log f = log2 + a*log(theta) - lgamma(a) - (2a+1) log y - theta/y^2
      -(2 * ps$alpha + 1) / y + 2 * ps$theta / y^3
    },
    "inv_gamma2" = {
      y <- x - s
      if (y <= 0) return(0)
      ## A12: no `shape <= 2` special case. That branch fired exactly when
      ## sd = Inf and returned the gradient of the IMPROPER -log(x) density,
      ## which no longer matches log_prior().
      ps <- .ig2_params(p1 - s, p2)
      # log f = a*log(b) - lgamma(a) - (a+1) log y - b/y
      -(ps$shape + 1) / y + ps$scale / y^2
    },
    "beta" = {
      a <- if (is.na(p3)) 0 else p3
      b <- if (is.na(p4)) 1 else p4
      if (x <= a || x >= b) return(0)
      sh <- .beta_shapes(p1, p2, a, b)
      if (!is.finite(sh[1]) || !is.finite(sh[2]) || sh[1] <= 0 || sh[2] <= 0)
        return(0)
      len <- b - a
      y <- (x - a) / len
      ((sh[1] - 1) / y - (sh[2] - 1) / (1 - y)) / len
    },
    "gamma" = {
      y <- x - s
      if (y <= 0) return(0)
      m <- p1 - s
      shape <- (m / p2)^2
      rate  <- m / p2^2
      (shape - 1) / y - rate
    },
    "normal" = -(x - p1) / p2^2,
    "uniform" = 0,
    ## Mirror log_prior_density(): an unrecognised distribution must fail loud,
    ## not silently contribute a zero prior gradient (which would pair with the
    ## silent flat prior to corrupt gradient-based sampling/mode-finding).
    stop(".dlog_prior_density1: unsupported prior distribution \"", dist, "\". ",
         "Supported: beta, gamma, normal, inv_gamma (= inv_gamma1), inv_gamma2, uniform.",
         call. = FALSE)
  )
}

#' Analytic gradient of the total log-prior over a named parameter vector.
#' Returns 0 for any parameter out of bounds (the caller already rejects those).
#' @noRd
.dlog_prior <- function(theta, prior_spec) {
  g <- numeric(length(theta))
  names(g) <- names(theta)
  has_p3 <- "p3" %in% names(prior_spec)
  has_p4 <- "p4" %in% names(prior_spec)
  for (i in seq_len(nrow(prior_spec))) {
    nm <- prior_spec$name[i]
    if (!(nm %in% names(theta))) next
    g[nm] <- .dlog_prior_density1(theta[[nm]], prior_spec$distribution[i],
                                  prior_spec$p1[i], prior_spec$p2[i],
                                  if (has_p3) prior_spec$p3[i] else NA_real_,
                                  if (has_p4) prior_spec$p4[i] else NA_real_)
  }
  g
}

#' Conservative set of parameter names that can move .get_shock_cov().
#'
#' .get_shock_cov(model, exo, params) reads params only through: a parameter
#' named exactly like a shock (estimated stderr, Priority 0), a
#' comma-containing correlation key (estimated corr), and the identifiers in
#' the shocks block's expression/value strings (stderr/variance/corr/cov and
#' their _expr columns), evaluated against params. Any other parameter leaves
#' Sigma_e bit-identical, so its dSigma_e FD is exactly zero. Identifiers are
#' extracted by a regex over EVERY character column of the shocks tables -- a
#' superset of the true symbols (function names, shock names), which is safe:
#' extra names only cost an FD that returns zero.
#' @noRd
.sigma_e_param_deps <- function(model, exo_names, par_names) {
  strs <- character(0)
  for (tab in list(model$shocks$variances, model$shocks$correlations)) {
    if (!is.data.frame(tab) || nrow(tab) == 0L) next
    for (cn in names(tab)) {
      if (is.character(tab[[cn]]) || is.factor(tab[[cn]]))
        strs <- c(strs, as.character(tab[[cn]]))
    }
  }
  strs <- strs[!is.na(strs) & nzchar(strs)]
  toks <- unlist(regmatches(strs, gregexpr("[A-Za-z_.][A-Za-z0-9_.]*", strs)),
                 use.names = FALSE)
  par_names[par_names %in% exo_names |
            grepl(",", par_names, fixed = TRUE) |
            par_names %in% toks]
}

#' Validate every prior_spec distribution name at BUILD time (Item F3).
#'
#' \code{.dlog_prior_density1} is only invoked lazily, inside the per-draw
#' gradient closure returned by \code{\link{make_posterior_grad}}; a typo or
#' unsupported distribution name in \code{prior_spec} would otherwise not
#' \code{stop()} until the first (or a later, chain-dependent) gradient
#' evaluation, crashing an MCMC/NUTS chain mid-run instead of at closure
#' construction. This calls the same dispatcher once per row (the numeric
#' value passed does not matter -- the switch on \code{distribution} fires
#' regardless) so an unsupported name fails loud at
#' \code{make_posterior_grad()} time, naming the offending parameter.
#' @noRd
.validate_prior_spec_dist <- function(prior_spec) {
  for (i in seq_len(nrow(prior_spec))) {
    tryCatch(
      .dlog_prior_density1(prior_spec$p1[i], prior_spec$distribution[i],
                           prior_spec$p1[i], prior_spec$p2[i]),
      error = function(e) {
        stop("make_posterior_grad(): prior_spec entry \"", prior_spec$name[i],
             "\" has an unsupported distribution \"", prior_spec$distribution[i],
             "\" (caught at build time, not mid-chain). Supported: beta, gamma, ",
             "normal, inv_gamma (= inv_gamma1), inv_gamma2, uniform.",
             call. = FALSE)
      })
  }
  invisible(TRUE)
}


# ============================================================================
# Analytic Kalman score for the shock-std (sigma) parameters
# ============================================================================

#' TRUE when the compiled Kalman-score recursion is available.
#' @noRd
.HAS_RCPP_KF_SCORE <- function() {
  if (!isTRUE(getOption("dynhr.use_rcpp", TRUE))) return(FALSE)
  exists("kf_score_sigma_cpp", envir = asNamespace("dynhr"),
         inherits = FALSE, mode = "function")
}

#' Kalman log-likelihood SCORE w.r.t. parameters that move only Sigma_e.
#'
#' Differentiates the exact per-step Kalman filter (the "dare"/"standard"
#' recursion, which agree to < 1e-10 nats) w.r.t. a set of parameters whose only
#' effect on the state space is through the shock covariance Sigma_e (true of
#' shock-std parameters at first order, by certainty equivalence). The decision
#' rule -- hence TT, RR, ZZ, DD -- is held fixed and reused.
#'
#' For each such parameter k the caller supplies dSigma_k = d Sigma_e / d theta_k.
#' The derivative recursion propagates ds_k, dP_k alongside s, P and accumulates
#'   d ll_t = -0.5 ( tr(F^-1 dF) + 2 dv' F^-1 v - v' F^-1 dF F^-1 v ).
#'
#' @param Y observation matrix (n_obs x T) already oriented; d subtracted inside.
#' @param TT,RR,ZZ,DD state-space matrices from the decision rule.
#' @param Sigma_e shock covariance; d the observation constant (dr$ys[obs]).
#' @param me_diag measurement-error diagonal (me_variance * I).
#' @param dSigma_list named list of dSigma_e/dtheta_k (n_exo x n_exo each).
#' @return list(loglik, score) where score is a named numeric vector.
#' @noRd
.kf_loglik_score_sigma <- function(Y, TT, RR, ZZ, DD, Sigma_e, d, me_diag,
                                   dSigma_list) {
  n_state <- nrow(TT); n_obs <- nrow(ZZ); n_T <- ncol(Y)
  tZZ <- t(ZZ)
  QQ  <- tcrossprod(RR %*% Sigma_e, RR)
  HH  <- tcrossprod(DD %*% Sigma_e, DD)
  SS  <- RR %*% Sigma_e %*% t(DD)
  ll_const <- -0.5 * n_obs * log(2 * pi)

  K <- length(dSigma_list); knm <- names(dSigma_list)
  # Per-parameter derivative pieces of Q, H, SS and the Lyapunov-initialised dP.
  # These (and the Lyapunov solves) are cheap, one-time, and stay in R so the
  # C++ kernel needs no Lyapunov solver.
  dH <- dS <- vector("list", K)
  dsk <- rep(list(numeric(n_state)), K)
  dQ_vecs <- matrix(0, n_state * n_state, K)
  for (k in seq_len(K)) {
    dSe <- dSigma_list[[k]]
    dH[[k]] <- tcrossprod(DD %*% dSe, DD)
    dS[[k]] <- RR %*% dSe %*% t(DD)
    dQ_vecs[, k] <- as.vector(tcrossprod(RR %*% dSe, RR))
  }

  # P_1 = lyap(TT, QQ) and each dP_1 = lyap(TT, dQ_k) solve the SAME discrete
  # Lyapunov operator (I - TT (x) TT). For small state vectors, factor it ONCE
  # (kron is n^2 x n^2) and solve all 1 + K right-hand sides together -- much
  # cheaper than 1 + K separate iterative solves. Fall back to the iterative
  # solver when the state is large enough that the kron factorisation would
  # dominate.
  ## On a (near-)unit-root TT, M = I - kron(TT, TT) is singular and the
  ## stationary Lyapunov P0 does not exist -- mirror .solve_lyapunov()'s
  ## rcond guard (R/backend-monolith.R) and fail gracefully (NA score,
  ## -Inf loglik) instead of erroring out of solve().
  if (n_state <= 20L) {
    M <- diag(n_state * n_state) - kronecker(TT, TT)
    if (rcond(M) < .Machine$double.eps) {
      return(list(loglik = -Inf, score = setNames(rep(NA_real_, K), knm)))
    }
    Xs  <- solve(M, cbind(as.vector(QQ), dQ_vecs))
    P   <- matrix(Xs[, 1], n_state, n_state)
    dPk <- lapply(seq_len(K), function(k) matrix(Xs[, k + 1L], n_state, n_state))
  } else {
    P   <- solve_lyapunov(TT, QQ)
    dPk <- lapply(seq_len(K),
                  function(k) solve_lyapunov(TT, matrix(dQ_vecs[, k],
                                                        n_state, n_state)))
  }
  if (!all(is.finite(P)) || any(vapply(dPk, function(x) !all(is.finite(x)), logical(1)))) {
    return(list(loglik = -Inf, score = setNames(rep(NA_real_, K), knm)))
  }
  dll <- numeric(K)

  s <- numeric(n_state)
  loglik <- 0
  Yd <- Y - d

  ## Fast path: compiled per-step recursion (the hot loop).
  if (.HAS_RCPP_KF_SCORE()) {
    ## ss_tol = the forward filter's (kalman_filter(ss_tol = .LYAP_TOL)), so
    ## the score recursion locks at the same period as the objective
    out <- kf_score_sigma_cpp(Yd, TT, RR, ZZ, DD, Sigma_e, HH, SS, me_diag, P,
                              dSigma_list, dH, dS, dPk, ss_tol = .LYAP_TOL)
    ## On ok = FALSE the C++ loop breaks out mid-recursion and `score` holds
    ## whatever partial (e.g. all-zero) accumulator it had at the failing
    ## step, not a meaningful gradient. Report NA, matching every other KF
    ## gradient kernel's failure contract (.kf_loglik_adjoint et al.), so
    ## callers can't mistake a zero-by-construction score for a real one.
    if (!isTRUE(out$ok)) {
      return(list(loglik = -Inf, score = setNames(rep(NA_real_, K), knm)))
    }
    sc <- as.numeric(out$score); names(sc) <- knm
    return(list(loglik = out$loglik, score = sc))
  }

  ## R fallback (bit-equivalent reference; see test-kf-score-parity).
  ## me_diag is TRUE observation noise (F3-D): it enters F AND the Joseph
  ## covariance update.
  has_me_true <- any(me_diag != 0)
  for (t in seq_len(n_T)) {
    v   <- Yd[, t] - as.numeric(ZZ %*% s)
    PZ  <- P %*% tZZ
    F   <- ZZ %*% PZ + HH + me_diag
    F   <- (F + t(F)) * 0.5
    ## Same failure contract as the compiled kernel (ok = FALSE): a non-PD
    ## forecast covariance degrades to -Inf/NA instead of throwing from the
    ## reference path (this chol() was the last unguarded call of the
    ## b827f41 crash class).
    Fc  <- tryCatch(chol(F), error = function(e) NULL)
    if (is.null(Fc)) {
      return(list(loglik = -Inf, score = setNames(rep(NA_real_, K), knm)))
    }
    Fi <- chol2inv(Fc)
    loglik <- loglik + ll_const - 0.5 * (2 * sum(log(diag(Fc))) + sum(v * (Fi %*% v)))
    Kg   <- (TT %*% PZ + SS) %*% Fi
    TmKZ <- TT - Kg %*% ZZ
    RmKD <- RR - Kg %*% DD
    Fiv  <- Fi %*% v

    for (k in seq_len(K)) {
      dv    <- -as.numeric(ZZ %*% dsk[[k]])
      dPZ   <- dPk[[k]] %*% tZZ
      dF    <- ZZ %*% dPZ + dH[[k]]
      Fi_dF <- Fi %*% dF
      dll[k] <- dll[k] - 0.5 * (sum(diag(Fi_dF)) +
                                2 * sum(dv * Fiv) -
                                sum(Fiv * (dF %*% Fiv)))
      dKg   <- (TT %*% dPZ + dS[[k]]) %*% Fi - Kg %*% (dF %*% Fi)
      dTmKZ <- -dKg %*% ZZ
      dRmKD <- -dKg %*% DD
      dsk[[k]] <- as.numeric(TT %*% dsk[[k]] + dKg %*% v + Kg %*% dv)
      dPn <- dTmKZ %*% tcrossprod(P, TmKZ) + TmKZ %*% tcrossprod(dPk[[k]], TmKZ) +
             TmKZ %*% tcrossprod(P, dTmKZ) +
             dRmKD %*% tcrossprod(Sigma_e, RmKD) +
             RmKD %*% tcrossprod(dSigma_list[[k]], RmKD) +
             RmKD %*% tcrossprod(Sigma_e, dRmKD)
      ## Tangent of the ME Joseph term P' += K me K' (me is data):
      ## dP' += dK me K' + K me dK'.
      if (has_me_true)
        dPn <- dPn + tcrossprod(dKg %*% me_diag, Kg) +
               tcrossprod(Kg %*% me_diag, dKg)
      dPk[[k]] <- (dPn + t(dPn)) * 0.5
    }

    s <- as.numeric(TT %*% s + Kg %*% v)
    Pn <- tcrossprod(TmKZ %*% P, TmKZ) + tcrossprod(RmKD %*% Sigma_e, RmKD)
    ## TRUE measurement-noise law (F3-D): P' += K me K'.
    if (has_me_true) Pn <- Pn + tcrossprod(Kg %*% me_diag, Kg)
    P  <- (Pn + t(Pn)) * 0.5
  }

  names(dll) <- knm
  list(loglik = loglik, score = dll)
}


# ============================================================================
# Assembled analytic/semi-analytic gradient of the log-posterior
# ============================================================================

#' Build an analytic-gradient closure for the log-posterior
#'
#' Returns \code{function(theta) -> named numeric} giving d/dtheta log p(theta|Y),
#' suitable for the \code{grad_fn} argument of \code{\link{nuts}} /
#' \code{dynhr_hmc}. Composition:
#'   * log-prior: analytic score for every parameter;
#'   * log-lik:   analytic Kalman score (reusing one model solve) for the
#'     parameters whose perturbation leaves the decision rule unchanged
#'     (shock-std params, by certainty equivalence -- auto-detected once at
#'     construction); a relative-step central difference for the rest.
#'
#' @param model,data,prior_spec,obs_vars,compiled,me_variance as in
#'   \code{\link{make_posterior}}.
#' @param theta_ref reference parameter vector for the one-time certainty-
#'   equivalence classification (default: prior means).
#' @param verbose print the analytic/numerical parameter split.
#' @param grad_method \code{"auto"} (default), \code{"hybrid"},
#'   \code{"implicit"}, \code{"adjoint"} or \code{"adjoint_solution"}.
#'   \code{"auto"} (also the default of the \code{grad_method} package option
#'   -- see \code{\link{dynhr_set_options}} -- so mode finding, the Hessian,
#'   SBC and the samplers all take it since 0.9.3.127) picks, once at build
#'   time, the fastest
#'   EXACT method supported for the likelihood: \code{"adjoint_solution"} for
#'   the Gaussian likelihood and for the cumulant likelihood when
#'   \code{cumulant_orders} contains 3 but not 4 (the orders its reverse-mode
#'   path covers; e.g. rbc2shock, orders 1:3: 22 ms against 1921 ms for
#'   \code{"implicit"}), and \code{"implicit"} for the Whittle and pruned
#'   likelihoods and the other cumulant orders. (Cumulant orders within
#'   \code{1:2} -- no skewness or kurtosis -- are matched on the first-order
#'   rule, and every \code{grad_method} then takes the same exact analytic
#'   mean/covariance gradient: \code{solution_derivatives} plus a derivative
#'   Lyapunov equation; this used to be finite differences of the forward.)
#'   On the Gaussian likelihood it picks
#'   \code{"hybrid"} instead when no analytic Kalman kernel covers the data at
#'   \code{theta_ref} (missing observations together with \code{me_extra} or
#'   a unit root; a unit root together with \code{me_extra} or
#'   \code{shock_scale}), since every draw would then take the FD-hybrid
#'   fallback anyway, always under \code{lik_init = "kappa"} (see
#'   \code{lik_init}), and when the filter drops observation components at
#'   \code{theta_ref} (a singular innovation covariance; see Details). The
#'   rule is fitted to a benchmark (median of 20
#'   interleaved calls; 2 to 68 parameters, 1 to 37 states):
#'   \code{"adjoint_solution"} was the fastest, or tied within timing noise,
#'   of the exact methods on every model -- e.g. Smets-Wouters (2007),
#'   36 parameters: 16 ms against 32 (\code{"adjoint"}), 91
#'   (\code{"implicit"}) and 97 (\code{"hybrid"}); a 68-parameter
#'   small-open-economy model: 45 ms against 120, 391 and 316 -- so no size
#'   threshold is applied. Every returned closure carries the method actually
#'   used as \code{attr(g, "grad_method")} (and the requested one as
#'   \code{attr(g, "grad_method_requested")}).
#'
#'   \code{"hybrid"} is the behaviour described above: an analytic Kalman
#'   score for the shock-std (sigma) parameters and a relative-step central
#'   difference for the rest.
#'
#'   \code{"implicit"} additionally differentiates the decision rule itself
#'   via implicit differentiation of the perturbation fixed point
#'   (\code{solution_derivatives}, Childers, Fernandez-Villaverde,
#'   Perla, Rackauckas & Wu 2022, NBER w30573) and propagates those
#'   solution derivatives through the Kalman filter with a single tangent
#'   (forward-sensitivity) recursion (\code{.kf_loglik_tangent}) covering ALL
#'   parameters at once -- one model solve and one filter pass per gradient
#'   evaluation, instead of one solve+filter per non-sigma parameter.
#'   Per-parameter, this falls back to the \code{"hybrid"} relative-step FD
#'   whenever \code{solution_derivatives()} reports \code{ok = FALSE} for that
#'   parameter, or the tangent-filter gradient for it is non-finite. On models
#'   with a unit root in the state-transition matrix (the same
#'   \code{|eig(TT)| > 1 - 1e-6} test used by \code{kalman_filter(...,
#'   lik_init = "auto")}), the tangent filter's stationary-Lyapunov
#'   initialization does not apply; \code{"implicit"} then warns ONCE (at the
#'   first gradient evaluation) and falls back to \code{"hybrid"} entirely for
#'   the lifetime of the returned closure.
#'
#'   \code{"adjoint"} is identical to \code{"implicit"} except the tangent
#'   recursion is replaced by the reverse-mode adjoint Kalman filter
#'   (\code{.kf_loglik_adjoint}): one forward pass storing per-step filter
#'   quantities and one backward sweep, so the filter-side cost is O(1) in
#'   the number of parameters instead of O(n_par). The two agree to ~1e-12;
#'   prefer \code{"adjoint"} for models with many estimated parameters.
#'
#'   \code{"adjoint_solution"} (Tier 18 A2) extends \code{"adjoint"} with
#'   reverse mode through the perturbation solve as well: the adjoint Kalman
#'   filter exports its bar matrices wrt (TT, RR, ZZ, DD, d), and
#'   \code{.solution_adjoint()} turns them into the structural-parameter
#'   gradient with TWO transposed solves total plus one Frobenius contraction
#'   per parameter -- no per-parameter generalized-Sylvester solve at all
#'   (\code{"implicit"}/\code{"adjoint"} pay one backsolve + RHS assembly per
#'   structural parameter via \code{solution_derivatives}). The Sigma_e
#'   channel (estimated shock stds, stderr expressions) is still routed
#'   through the filter adjoint's \code{G_Sig}. Agrees with \code{"adjoint"}
#'   to ~1e-10. Draws needing the missing-data or exact-diffuse kernels (which
#'   do not export bars) fall back to the \code{"adjoint"} construction for
#'   that draw; whittle/pruned likelihoods treat it as \code{"implicit"}. For
#'   the cumulant likelihood it reverses orders 1-3 on the order-2 rule (the
#'   moment cotangent is taken back through the third-cumulant tensor
#'   Lyapunov, then one first-order and one order-2 solution adjoint);
#'   requests including order 4 take the \code{"implicit"} path, and orders
#'   within \code{1:2} the exact first-order mean/covariance gradient. It is what \code{"auto"} selects for the Gaussian
#'   likelihood and for cumulant orders containing 3 but not 4 (see above).
#' @param me_extra Optional n_obs x T matrix of per-period additive
#'   measurement-error variances (filter_tunes); the gradient is of the same
#'   tuned likelihood \code{make_log_posterior} evaluates.
#' @param shock_scale Optional n_exo x T matrix of known per-period shock-std
#'   scale factors (heteroskedastic_shocks).
#' @param likelihood \code{"gaussian"} (default), \code{"whittle"}, or
#'   \code{"cumulant"}; selects the likelihood whose gradient is built.
#' @param freq_band Numeric \code{c(lo, hi)} in radians; Whittle band
#'   restriction (ignored for the Gaussian likelihood).
#' @param cumulant_orders,cumulant_weight Cumulant orders to match (default
#'   \code{1:4}) and the weighting scheme; used only when
#'   \code{likelihood = "cumulant"} and must match the forward
#'   \code{make_log_posterior_cumulant} settings.
#' @param debias Logical; for \code{likelihood = "whittle"}, differentiate
#'   the debiased Whittle loglik (expected-periodogram form, Sykulski et al.
#'   2019 — the default of the estimation entry points). Must match the
#'   \code{debias} setting of the posterior being sampled. Ignored for the
#'   Gaussian likelihood.
#' @param pruned_order Integer perturbation order for the pruned state-space
#'   likelihood path (default \code{2L}); currently \code{2L} or \code{3L}.
#' @param power Power-posterior (generalised-Bayes) exponent zeta in (0, 1].
#'   \code{NULL} (default) resolves the \code{power_posterior} option ONCE at
#'   build time, exactly as the \code{make_log_posterior*} factories do, so the
#'   gradient is of the SAME tempered target
#'   \code{log p(theta) + zeta * loglik(theta) + log p_sys(theta)}: the
#'   likelihood score is multiplied by zeta, the prior score is not.
#' @param system_priors Optional \code{system_prior_spec} (as passed to
#'   \code{make_log_posterior}). Its log-density is part of the target, so its
#'   gradient is added: for the Gaussian and Whittle likelihoods by a central
#'   difference of the system-prior density alone (one first-order solve per
#'   side, no filter pass); for the cumulant and pruned likelihoods (whose
#'   system prior sees a lifted decision rule) by a central difference of the
#'   system-prior term the forward posterior closure itself evaluates.
#' @param lik_init Kalman filter \code{P0} initialization of the Gaussian
#'   likelihood being differentiated: \code{"auto"} (default),
#'   \code{"stationary"}, \code{"diffuse"} or \code{"kappa"}, exactly as in
#'   \code{\link{make_log_posterior}} -- it must match the \code{lik_init} of
#'   the posterior being sampled or optimized, because the initializations
#'   are different likelihoods (on a root within \code{1e-6} of the unit
#'   circle \code{"diffuse"} puts a diffuse prior on it where \code{"auto"}
#'   keeps the stationary one). The analytic kernels follow the init in force
#'   per draw: the stationary kernels, or the exact-diffuse adjoint where the
#'   init in force is diffuse; under \code{"stationary"} a draw with a unit
#'   root (whose log-posterior is \code{-Inf}) gets the prior score only.
#'   \code{"kappa"} (the big-\code{kappa} approximation) has no analytic
#'   kernel: \code{"implicit"}, \code{"adjoint"} and
#'   \code{"adjoint_solution"} refuse it (an error of class
#'   \code{dynhr_error_grad_lik_init}), \code{"auto"} resolves to
#'   \code{"hybrid"}, and \code{"hybrid"} then differentiates every parameter
#'   by finite differences of the kappa-initialized likelihood. Ignored by
#'   the Whittle, cumulant and pruned likelihoods (as by their posteriors).
#' @return \code{function(theta)} returning the gradient vector, suitable for the
#'   \code{grad_fn} argument of \code{\link{nuts}}, named and in
#'   \code{prior_spec$name} order. \code{theta} is read as the log-posterior
#'   closures read it: by name in any order, an unnamed \code{theta} of length
#'   \code{nrow(prior_spec)} in \code{prior_spec$name} order, anything else an
#'   error of class \code{dynhr_error_theta_names}. For the Gaussian
#'   likelihood it carries \code{attr(, "logpost_grad")}, a
#'   \code{function(theta)} returning \code{list(logpost, grad)} from ONE
#'   gradient pass: \code{grad} is the gradient above and \code{logpost} the
#'   log-posterior of \code{\link{make_log_posterior}} with the same
#'   arguments (\code{power}, \code{system_priors}, \code{lik_init},
#'   \code{me_extra}, \code{shock_scale}; no infeasibility penalty), assembled
#'   from the log-likelihood the pass already computed and \code{-Inf} exactly
#'   where that posterior is. The gradient samplers (\code{\link{nuts}}, HMC,
#'   ChEES, MALA) use it -- one evaluation per new position instead of a
#'   gradient plus a separate log-posterior call, the gradient carried along
#'   the trajectory -- once it agrees with their own \code{log_post_fn} at the
#'   start point to \code{1e-10} relative; otherwise they keep separate calls.
#'
#' @details
#' The shock-std parameters get an exact analytic Kalman score (validated to
#' machine precision); the remaining structural parameters use a relative-step
#' central difference. On models with a fast (C++) Kalman filter and a cheap
#' solve, the pure-R score recursion can be slower than the default numerical
#' forward-difference gradient -- it is provided for exactness (no finite-
#' difference truncation error) and as the basis for a future C++ port. On large
#' models where the solve dominates, it avoids the per-parameter re-solves and
#' wins. Pass it explicitly: \code{nuts(lp, theta0, grad_fn = make_posterior_grad(...))}.
#'
#' \code{grad_method = "implicit"} is the full implicit-differentiation
#' gradient of Childers et al. (2022): it amortizes the cost of differentiating
#' the decision rule across all non-sigma parameters (one shared Sylvester/QR
#' factorization, see \code{solution_derivatives}) and obtains the
#' likelihood gradient for every parameter from a single tangent Kalman-filter
#' pass (\code{.kf_loglik_tangent}), avoiding the per-parameter re-solve and
#' re-filter of \code{"hybrid"}.
#'
#' \strong{Singular innovation covariance.} When the Gaussian filter at
#' \code{theta_ref} falls back to the univariate filter (Dynare's rule:
#' \code{rcond(F) < kalman_tol} and either a diagonal entry of \code{F} or the
#' \code{rcond} of its correlation form below \code{kalman_tol = 1e-10}) and
#' drops observation components (conditional variance below the absolute
#' \code{kalman_tol}; reported as \code{$diagnostics$n_dropped} by
#' \code{\link{kalman_filter}}) while the analytic kernels still evaluate a
#' finite, different (undropped) likelihood -- typically an observable whose
#' variance, or whose variance conditional on the other observables, is
#' below \code{1e-10} while the others are of ordinary scale -- the dropped
#' set can change with the parameters and the log-posterior is then
#' discontinuous. (An EXACT
#' singularity, e.g. an observable that is an exact combination of others
#' without measurement error, is not affected: the kernels fail on it and
#' every draw takes finite differences of the posterior.) The analytic
#' kernels would differentiate the wrong function, so \code{"implicit"},
#' \code{"adjoint"} and \code{"adjoint_solution"} refuse such a model (class
#' \code{dynhr_error_grad_singular_F}); \code{"auto"} resolves to
#' \code{"hybrid"}, which then warns (class
#' \code{dynhr_warning_grad_singular_F}) and differentiates every parameter by
#' finite differences of the log-posterior. On every draw \code{"hybrid"}
#' also uses its analytic shock-std score only where the score recursion's
#' log-likelihood equals the posterior's own.
#'
#' @references Childers, D., Fernandez-Villaverde, J., Perla, J., Rackauckas,
#'   C., & Wu, P. (2022). \emph{Differentiable State-Space Models and
#'   Hamiltonian Monte Carlo Estimation}. NBER Working Paper No. 30573.
#' @seealso \code{\link{make_posterior}}, \code{\link{nuts}},
#'   \code{solution_derivatives}
#' @export
make_posterior_grad <- function(model, data, prior_spec, obs_vars, compiled,
                                me_variance = 0, theta_ref = NULL,
                                verbose = FALSE,
                                grad_method = c("auto", "hybrid", "implicit",
                                                "adjoint",
                                                "adjoint_solution"),
                                me_extra = NULL, shock_scale = NULL,
                                likelihood = "gaussian",
                                freq_band = c(0, pi),
                                cumulant_orders = 1:4,
                                cumulant_weight = "identity",
                                debias = TRUE,
                                pruned_order = 2L,
                                power = NULL,
                                system_priors = NULL,
                                lik_init = c("auto", "stationary",
                                             "diffuse", "kappa")) {
  ## The analytic / adjoint gradients build the measurement intercept from
  ## dr$ys and ignore observation_trends; callers catch this and fall back to
  ## finite differences of the (trend-aware) Gaussian posterior.
  .refuse_obs_trends(model, "make_posterior_grad()")
  grad_method <- match.arg(grad_method)
  lik_init    <- match.arg(lik_init)
  grad_method_requested <- grad_method   # "auto" is resolved after the base solve
  .validate_prior_spec_dist(prior_spec)  # F3: fail loud at BUILD time, not mid-chain
  ## Brief 23 A6: the gradient must be of the SAME target the posterior
  ## closures evaluate, logpost = logprior + zeta * loglik + log p_sys. Resolve
  ## zeta once, here, exactly as every make_log_posterior* factory does.
  power <- .resolve_power_posterior(power, "make_posterior_grad")
  if (!is.null(system_priors) && length(system_priors) == 0L)
    system_priors <- NULL
  sys_cache <- cache_system_structure(compiled)
  exo       <- model$varexo_names
  par_names <- prior_spec$name
  np        <- length(par_names)
  use_whittle  <- identical(likelihood, "whittle")
  use_cumulant <- identical(likelihood, "cumulant")
  use_pruned   <- identical(likelihood, "pruned")
  ## Fail loud on likelihoods with no gradient path, rather than silently
  ## returning a Gaussian gradient for them (pskf/student_t/tpf/ppf/copf have
  ## no analytic/FD gradient here; the sampler gate .ctx_allows_analytic_gradient
  ## already excludes them, so this only guards direct/user calls and typos).
  ## The analytic/adjoint score recursions implement H = me_variance * I
  ## only: a per-observable vector is refused (classed) at build time rather
  ## than differentiated as if it were a scalar. Samplers never get here with
  ## one -- .ctx_allows_analytic_gradient() routes them to FD of the posterior,
  ## which takes H = diag(me).
  me_variance <- .kf_me_variance(me_variance, obs_vars, "make_posterior_grad",
                                 allow_vector = FALSE)
  if (!likelihood %in% c("gaussian", "whittle", "cumulant", "pruned"))
    stop("make_posterior_grad(): no gradient path for likelihood '", likelihood,
         "'. Supported: gaussian, whittle, cumulant, pruned. ",
         "(pskf/student_t/tpf/ppf/copf have no gradient path.)")
  lp_fn <- if (use_whittle) {
    make_log_posterior_whittle(model, data, prior_spec, obs_vars, compiled,
                               me_variance = me_variance,
                               freq_band = freq_band,
                               debias = debias, power = power)
  } else if (use_cumulant) {
    make_log_posterior_cumulant(model, data, prior_spec, obs_vars, compiled,
                                me_variance = me_variance,
                                cumulant_orders = cumulant_orders,
                                cumulant_weight = cumulant_weight,
                                power = power)
  } else if (use_pruned) {
    ## Pruned-SS Gaussian KF on the AFVRR augmented state. Order 2: analytic
    ## adjoint-chain gradient (R/pruned-grad-chain.R), FD fallback per
    ## parameter for anything the chain does not cover. Order 3: the fold
    ## chain needs a derivative-Lyapunov pass through the order-3 pruned
    ## system (R/pruned-state-space-order3.R) that is OUT OF SCOPE here
    ## (deferred); make_posterior_grad falls back to numerical FD for ALL
    ## parameters at order 3, same as before this change.
    if (identical(pruned_order, 3L) || identical(pruned_order, 3)) {
      make_log_posterior_pruned3(model, data, prior_spec, obs_vars, compiled,
                                 me_variance = me_variance, power = power)
    } else {
      make_log_posterior_pruned(model, data, prior_spec, obs_vars, compiled,
                                me_variance = me_variance, power = power)
    }
  } else {
    make_log_posterior(model, data, prior_spec, obs_vars, compiled,
                       me_variance = me_variance,
                       me_extra = me_extra, shock_scale = shock_scale,
                       power = power, lik_init = lik_init)
  }
  ## lik_init only changes the Gaussian likelihood (the other posteriors
  ## ignore it); keep it only there so the kappa gate below cannot fire for
  ## a likelihood it does not apply to.
  use_kappa <- identical(likelihood, "gaussian") && identical(lik_init, "kappa")
  ## The big-kappa P0 (.build_P0) has no analytic Kalman kernel: the tangent /
  ## adjoint kernels start from the Lyapunov P0 or the exact-diffuse
  ## (P_inf, P_star) pair, which are different likelihoods (a kappa-dependent
  ## level AND, through the kappa P0's stable block, a different slope).
  ## Refuse rather than return the gradient of another posterior.
  if (use_kappa &&
      grad_method %in% c("implicit", "adjoint", "adjoint_solution"))
    .dynhr_abort("make_posterior_grad(): grad_method = \"", grad_method,
                 "\" has no analytic kernel for lik_init = \"kappa\" (the ",
                 "big-kappa P0 approximation). Use grad_method = \"hybrid\" ",
                 "or \"auto\" (finite differences of the kappa-initialized ",
                 "likelihood), or lik_init = \"diffuse\" (exact diffuse ",
                 "initialization, which has an exact adjoint).",
                 class = "dynhr_error_grad_lik_init")
  ## Time-varying inputs present? .kf_loglik_score_sigma (the analytic sigma
  ## score used by the "hybrid" closure) is NOT tv-aware: its loglik would be
  ## of a different likelihood than lp_fn, and mixing the two inside the FD
  ## step produces garbage gradients. Under tv, hybrid treats sigma params by
  ## consistent FD instead (tangent/adjoint paths are fully tv-aware).
  has_tv <- !is.null(me_extra) ||
            (!is.null(shock_scale) && !all(shock_scale == 1))

  ## Solve the model at a parameter vector, returning the decision rule (or NULL).
  .solve_dr <- function(theta) .grad_solve_dr(model, compiled, sys_cache, theta)

  ## --- Fused log-posterior (W92) ---------------------------------------------
  ## The Gaussian closures below already hold the log-likelihood at theta: the
  ## analytic kernel's (tang$loglik) or the objective's own (lp_fn(theta), the
  ## FD base of "hybrid" and of every fallback). They note it here, so that
  ## attr(grad, "logpost_grad") -- see .fused_of() at the end -- returns
  ## list(logpost, grad) from ONE pass instead of a gradient call plus a
  ## separate log-posterior call. Written only for the top-level theta (the FD
  ## perturbations go through .fd_loglik_grad1(), which does not note).
  fuse_env <- new.env(parent = emptyenv())
  .fuse_note <- function(loglik) fuse_env$loglik <- loglik

  ## --- Tempering + system prior (brief 23 A6) --------------------------------
  ## Every branch below builds a closure returning d/dtheta [logprior + loglik]
  ## (the untempered score; lp_fn is used only for its RAW $loglik). The
  ## target the posterior closures evaluate is
  ##   logpost = logprior + zeta * loglik + log p_sys,
  ## so .finish_grad() rescales the likelihood part by zeta and adds the
  ## system-prior score. zeta = 1 with no system prior returns the branch
  ## closure itself, untouched (bit-identical to the untempered gradient).
  ##
  ## System-prior density at theta. Gaussian/Whittle closures evaluate it on
  ## the first-order rule and .get_shock_cov() -- exactly what .solve_dr()
  ## gives, so no filter pass is needed. Cumulant/pruned closures evaluate it
  ## on a LIFTED (order-2/3) rule with their own Sigma_e; rather than
  ## re-implement those solves, take the system-prior term out of the forward
  ## closure itself: at zeta = 1, logpost - loglik - log_prior(theta) is the
  ## system-prior density in both the "lp" and the "extra" closure modes.
  .sp_logdens <- NULL
  if (!is.null(system_priors)) {
    if (use_cumulant || use_pruned) {
      lp_sp_fn <- if (use_cumulant) {
        make_log_posterior_cumulant(model, data, prior_spec, obs_vars, compiled,
                                    me_variance = me_variance,
                                    cumulant_orders = cumulant_orders,
                                    cumulant_weight = cumulant_weight,
                                    system_priors = system_priors, power = 1)
      } else if (identical(pruned_order, 3L) || identical(pruned_order, 3)) {
        make_log_posterior_pruned3(model, data, prior_spec, obs_vars, compiled,
                                   me_variance = me_variance,
                                   system_priors = system_priors, power = 1)
      } else {
        make_log_posterior_pruned(model, data, prior_spec, obs_vars, compiled,
                                  me_variance = me_variance,
                                  system_priors = system_priors, power = 1)
      }
      .sp_logdens <- function(theta) {
        r <- lp_sp_fn(theta)
        if (!is.finite(r$logpost) || !is.finite(r$loglik)) return(-Inf)
        r$logpost - r$loglik - log_prior(theta, prior_spec)
      }
    } else {
      .sp_logdens <- function(theta) {
        sol <- .solve_dr(theta)
        if (is.null(sol)) return(-Inf)
        .eval_system_priors(
          system_priors,
          list(theta   = theta,
               model   = model,
               dr      = sol$dr,
               Sigma_e = .get_shock_cov(model, exo, sol$params),
               params  = sol$params))
      }
    }
  }

  ## Central difference of the system-prior log-density, per parameter. A
  ## hard restriction is flat inside its region, so a side that falls outside
  ## (-Inf) takes the one-sided difference; both sides outside -> 0 (the draw
  ## itself is rejected by the posterior there, so its gradient is moot).
  .sp_grad <- function(theta) {
    gs <- setNames(numeric(np), par_names)
    f0 <- .sp_logdens(theta)
    fuse_env$sp <- f0   # the fused value's system-prior term at this theta

    for (nm in par_names) {
      h  <- 1e-5 * max(abs(theta[[nm]]), 1e-3)
      tp <- theta; tp[nm] <- tp[nm] + h
      tm <- theta; tm[nm] <- tm[nm] - h
      fp <- .sp_logdens(tp); fm <- .sp_logdens(tm)
      gs[nm] <- if (is.finite(fp) && is.finite(fm)) (fp - fm) / (2 * h)
                else if (is.finite(fp) && is.finite(f0)) (fp - f0) / h
                else if (is.finite(fm) && is.finite(f0)) (f0 - fm) / h
                else 0
    }
    gs
  }

  ## Every returned closure records the method that actually built it (the
  ## resolution of "auto", see .resolve_auto_grad_method()) and what the caller
  ## asked for, so a sampler / run record can report the gradient it used.
  .stamp_method <- function(gfn) {
    attr(gfn, "grad_method") <- grad_method
    attr(gfn, "grad_method_requested") <- grad_method_requested
    gfn
  }
  .finish_grad <- function(gfn) {
    if (power == 1 && is.null(system_priors)) return(.stamp_method(gfn))
    wrapped <- function(theta) {
      theta <- .theta_by_name(theta, par_names, "make_posterior_grad")
      g  <- gfn(theta)
      gp <- .dlog_prior(theta, prior_spec)
      ## g = gp + score_loglik  ->  gp + zeta * score_loglik
      g  <- gp + power * (g - gp)
      if (!is.null(system_priors)) g <- g + .sp_grad(theta)
      g
    }
    ## Carry the branch closure's diagnostics (kernel_stats) through.
    for (a in setdiff(names(attributes(gfn)), "srcref"))
      attr(wrapped, a) <- attr(gfn, a)
    .stamp_method(wrapped)
  }

  ## W92: the fused companion of a finished Gaussian gradient closure `gfn`,
  ## function(theta) -> list(logpost, grad): grad is gfn(theta) itself, and
  ## logpost is make_log_posterior's target assembled exactly as its closure
  ## assembles it (R/posterior-closure.R, system_prior_mode "lp"):
  ##   power * loglik + (log_prior + log p_sys),
  ## -Inf whenever the prior, the loglik or the system prior is -Inf. The
  ## loglik is the one the gradient pass noted (fuse_env): the analytic
  ## kernel's, or lp_fn's own on every "hybrid" / fallback draw. A draw on
  ## which the pass noted none (a failed solve, a "stationary" reject: the
  ## gradient is the prior score there) takes lp_fn's loglik, so the value is
  ## -Inf exactly where the objective's is. Attached as
  ## attr(grad, "logpost_grad"); the samplers use it only after checking it
  ## against their own log-posterior at the start point
  ## (.hmc_fused_target(), R/sampler-hmc.R).
  .fused_of <- function(gfn) {
    force(gfn)
    function(theta) {
      theta <- .theta_by_name(theta, par_names, "make_posterior_grad")
      fuse_env$loglik <- NULL
      fuse_env$sp     <- NULL
      g   <- gfn(theta)
      lpp <- log_prior(theta, prior_spec)
      if (!is.finite(lpp)) return(list(logpost = -Inf, grad = g))
      ll <- fuse_env$loglik
      if (is.null(ll)) ll <- lp_fn(theta)$loglik
      if (!is.finite(ll)) return(list(logpost = -Inf, grad = g))
      if (!is.null(system_priors)) {
        sp <- fuse_env$sp
        if (is.null(sp)) sp <- .sp_logdens(theta)
        if (!is.finite(sp)) return(list(logpost = -Inf, grad = g))
        lpp <- lpp + sp
      }
      list(logpost = power * ll + lpp, grad = g)
    }
  }
  .attach_fused <- function(gfn) {
    attr(gfn, "logpost_grad") <- .fused_of(gfn)
    gfn
  }

  ## --- One-time classification: which params leave the decision rule fixed? ---
  ## theta_ref is read positionally below (th[k], is_sigma[k]), so it is put
  ## into prior order by name first -- like every theta the closures take
  ## (.theta_by_name(), R/posterior.R).
  if (is.null(theta_ref)) theta_ref <- setNames(prior_spec$mean, par_names)
  else theta_ref <- .theta_by_name(theta_ref, par_names,
                                   "make_posterior_grad", "theta_ref")
  base <- .solve_dr(theta_ref)
  is_sigma <- rep(FALSE, np); names(is_sigma) <- par_names
  if (!is.null(base)) {
    ## The sigma-like paths below hold the WHOLE measurement system fixed except
    ## Sigma_e: TT/RR/ZZ/DD (ghx/ghu) AND the intercept d = ys[obs_vars]. A
    ## parameter that moves only the steady state of an observable (an
    ## estimated measurement constant, `y_obs = y + mu`) leaves ghx/ghu
    ## unchanged, so a ghx/ghu-only test classified it sigma-like and every
    ## grad_method silently dropped its likelihood score (returned the prior
    ## score alone). Include d in the invariance test.
    g0 <- list(ghx = base$dr$ghx, ghu = base$dr$ghu,
               d = as.numeric(base$dr$ys[obs_vars]))
    for (k in seq_len(np)) {
      th <- theta_ref; h <- 1e-5 * max(abs(th[k]), 1e-3)
      th[k] <- th[k] + h
      d2 <- .solve_dr(th)
      if (!is.null(d2)) {
        dchg <- max(abs(d2$dr$ghx - g0$ghx), abs(d2$dr$ghu - g0$ghu))
        d_chg <- max(abs(as.numeric(d2$dr$ys[obs_vars]) - g0$d) /
                       pmax(1, abs(g0$d)))
        ## decision rule AND observable steady state invariant -> sigma-like
        is_sigma[k] <- dchg < 1e-10 && d_chg < 1e-10
      }
    }
  }
  sig_names <- par_names[is_sigma]
  num_names <- par_names[!is_sigma]
  if (verbose)
    .dynhr_cat(sprintf("  Analytic gradient: %d analytic (Kalman score: %s), %d numerical (%s)\n",
                length(sig_names), paste(sig_names, collapse = ","),
                length(num_names), paste(num_names, collapse = ",")))

  ## Singular innovation covariance at theta_ref (W68): does the Gaussian
  ## forward DROP observation components there (kalman_filter's singular-F
  ## fallback to the univariate filter, which skips every component whose
  ## conditional variance is below the ABSOLUTE kalman_tol)? The dropped set
  ## can change with theta, so the objective can be discontinuous, and the
  ## analytic kernels -- dense multivariate recursions -- differentiate the
  ## undropped likelihood, a different function. No gradient method can
  ## repair that; the builder says so and uses FD of the objective. (W68 met
  ## this on art_zlb_mcp, whose well-conditioned but ~1e-9-scale F the old
  ## absolute pivot cut called singular; since W74 the fallback follows
  ## Dynare's scale-invariant rule (.kf_F_singular) and art_zlb_mcp stays
  ## multivariate. What remains is an F that Dynare also sends to the
  ## univariate filter: a component variance below kalman_tol beside
  ## components of ordinary scale.)
  sF <- if (identical(likelihood, "gaussian"))
    .grad_singular_F_mismatch(model, data, obs_vars, base, me_variance,
                              lik_init, me_extra, shock_scale)
  else NULL
  singular_F <- !is.null(sF)
  if (singular_F) {
    msg_sF <- paste0(
      "the Gaussian likelihood at theta_ref drops ", sF$n_dropped,
      " observation component(s) whose conditional variance is below ",
      "kalman_tol = ", format(.KF_ZERO_VAR_TOL), " (a numerically singular ",
      "innovation covariance: the filter falls back to the univariate ",
      "filter, which skips them). The dropped set changes with the ",
      "parameters, so the log-posterior is discontinuous there, and the ",
      "analytic Kalman kernels evaluate the undropped likelihood instead -- ",
      "a different function (", format(sF$ll_kernel, digits = 6), " vs ",
      format(sF$ll_forward, digits = 6), " at theta_ref). Rescale the ",
      "observables (e.g. to percent), add measurement error, or use a ",
      "gradient-free sampler.")
    if (grad_method %in% c("implicit", "adjoint", "adjoint_solution"))
      .dynhr_abort("make_posterior_grad(): grad_method = \"", grad_method,
                   "\" refused: ", msg_sF,
                   class = "dynhr_error_grad_singular_F")
    .dynhr_warn("make_posterior_grad(): ", msg_sF, " The returned gradient ",
                "is finite differences of the log-posterior for every ",
                "parameter (valid only between switches of the dropped set).",
                class = "dynhr_warning_grad_singular_F")
  }

  ## grad_method = "auto": resolve ONCE, here, from the likelihood, the data
  ## and the decision rule at theta_ref (see .resolve_auto_grad_method()).
  if (identical(grad_method, "auto")) {
    grad_method <- .resolve_auto_grad_method(
      likelihood, base = base, model = model, data = data,
      me_extra = me_extra, shock_scale = shock_scale,
      cumulant_orders = cumulant_orders, lik_init = lik_init,
      singular_F = singular_F)
    if (verbose)
      .dynhr_cat(sprintf("  Gradient method: auto -> %s\n", grad_method))
  }

  ## State-space matrices from a decision rule (mirrors kalman_filter()).
  obs_in_endo <- function(dr) match(obs_vars, dr$endo_names)
  me_diag <- me_variance * diag(length(obs_vars))

  ## dSigma_e/dtheta_nm for sigma-like params (and reused as the certainty-
  ## equivalence dSigma_e block under "implicit"). Parameters that cannot
  ## reach .get_shock_cov() at all get an exact zero (W59). See
  ## .sigma_e_param_deps() for the conservative dependency set.
  ##
  ## W91: the derivative is EXACT (.shock_cov_deriv_eval: forward-mode walk of
  ## .get_shock_cov's priority chain with stats::D() on the shocks-block
  ## expressions), evaluated once per params vector for every dependent
  ## parameter and memoised, so the per-parameter call sites below read it
  ## without recomputation. The central FD (two .get_shock_cov() calls per
  ## parameter; NZSIM: 20 x 2 calls, ~4 ms of a 26 ms gradient) remains only
  ## where no exact derivative exists: a model whose shocks-block expressions
  ## D() cannot differentiate exactly (plan NULL), or a draw where an entry is
  ## non-finite (sqrt of a zero variance).
  sigma_e_deps <- .sigma_e_param_deps(model, exo, par_names)
  zero_dSigma  <- matrix(0, length(exo), length(exo), dimnames = list(exo, exo))
  .dSigma_e_fd1 <- function(theta, params, nm) {
    h <- 1e-6 * max(abs(theta[[nm]]), 1e-3)
    tp <- theta; tm <- theta; tp[nm] <- tp[nm] + h; tm[nm] <- tm[nm] - h
    pp <- params; pm <- params
    pp[nm] <- tp[[nm]]; pm[nm] <- tm[[nm]]
    (.get_shock_cov(model, exo, pp) - .get_shock_cov(model, exo, pm)) / (2 * h)
  }
  sigma_plan <- if (length(sigma_e_deps)) {
    params_ref <- if (!is.null(base)) base$params
                  else .apply_theta_to_params(model, theta_ref)
    .shock_cov_deriv_plan(model, exo, params_ref, sigma_e_deps)
  } else NULL
  dS_memo <- new.env(parent = emptyenv())
  .dSigma_e_d <- function(theta, params, nm) {
    if (!(nm %in% sigma_e_deps)) return(zero_dSigma)
    if (!identical(dS_memo$params, params)) {
      ## A non-finite Sigma_e (an overflowing std, e.g. a sampler trial point
      ## at s = 4e154 -> s^2 = Inf) has no exact derivative, and the walk's
      ## drift guard would compare Inf - Inf (NA) and ERROR, crashing the
      ## sampler instead of registering a divergence (W92). Such a draw takes
      ## the per-parameter FD below, as a non-finite entry already did.
      S_ref <- .get_shock_cov(model, exo, params)
      dS_memo$val <- if (all(is.finite(S_ref)))
        .shock_cov_deriv_eval(sigma_plan, model, exo, params, S_ref) else NULL
      dS_memo$params <- params
    }
    d <- dS_memo$val[[nm]]
    if (is.null(d)) .dSigma_e_fd1(theta, params, nm) else d
  }

  ## Per-parameter relative-step central/forward FD of the loglik, holding the
  ## decision rule fixed at `theta` except for parameter `nm` (the "hybrid"
  ## fallback used both by grad_method = "hybrid" for num_names, and by
  ## grad_method = "implicit" for any parameter whose solution-derivative or
  ## tangent-filter result is unusable).
  .fd_loglik_grad1 <- function(theta, nm, base_ll) {
    h <- 1e-5 * max(abs(theta[[nm]]), 1e-3)
    tp <- theta; tp[nm] <- tp[nm] + h
    llp <- lp_fn(tp)$loglik
    if (is.finite(llp) && is.finite(base_ll)) (llp - base_ll) / h else NA_real_
  }

  ## --- "cumulant" gradient closure ------------------------------------------
  ## The cumulant likelihood is differentiated by cumulant_loglik_grad(), which
  ## needs the SAME decision rule the forward make_log_posterior_cumulant solves.
  ## The forward solves order 2 ONLY when an order-3/4 cumulant is requested
  ## (ghxx/ghss feed skewness/kurtosis); otherwise the bare first-order rule.
  ## We mirror that solve-order choice exactly below (order-1 core -> stationarity
  ## guard -> conditional order-2 solve), then take the implicit (analytic,
  ## tensor-Lyapunov) gradient. Validated against numDeriv of the forward logpost
  ## for BOTH order-1:2 and order-1:4 (test-gradient-exactness-audit.R). Any
  ## non-finite analytic entry takes the exact FD-of-forward fallback, so the
  ## gradient is always consistent with lp_fn.
  if (use_cumulant) {
    cumulant_grad_fn <- function(theta) {
      theta <- .theta_by_name(theta, par_names, "make_posterior_grad")
      g <- .dlog_prior(theta, prior_spec)        # analytic prior score
      params <- .apply_theta_to_params(model, theta)
      ss <- solve_steady_state(model, compiled, params, verbose = FALSE)
      if (is.null(ss) || !isTRUE(ss$converged)) return(g)   # infeasible: prior-only
      params <- ss$params %||% params                       # consistent p_c (Tier 13)
      sys <- extract_system_matrices_fast(sys_cache, ss$ss, params)
      d1  <- .solve_from_system(sys, model, compiled, ss$ss, params, FALSE)
      if (is.null(d1) || !isTRUE(d1$bk_satisfied)) return(g)
      ghx_state <- d1$ghx[d1$state_idx, , drop = FALSE]
      if (max(Mod(eigen(ghx_state, symmetric = FALSE, only.values = TRUE)$values)) >= 1) return(g)
      ## Mirror make_log_posterior_cumulant's solve-order choice EXACTLY: it
      ## solves order 2 only when an order-3/4 cumulant is requested (ghxx/ghss
      ## feed the skewness/kurtosis terms), else the bare first-order rule. The
      ## gradient must differentiate the SAME function the forward evaluates;
      ## unconditionally solving order 2 here made the order-1:2 gradient the
      ## derivative of the order-2 posterior (wrong sign on some params).
      solve_order <- if (any(cumulant_orders >= 3L)) 2L else 1L
      if (solve_order >= 2L) {
        Sigma_e <- .get_shock_cov(model, exo, params)
        dr_use <- tryCatch(
          solve_perturbation_order2(model, compiled, ss$ss, params, dr1 = d1,
                                    Sigma_e = Sigma_e, h = 1e-4, verbose = FALSE),
          error = function(e) .dynhr_reraise_bug(e, NULL))
        if (is.null(dr_use)) dr_use <- d1   # order-2 failed: orders 3-4 drop out
      } else {
        dr_use <- d1                        # first-order: matches the forward
      }
      ## grad_method == "adjoint_solution": reverse-mode (O(1) in P) path for
      ## orders 1-3 via .cumulant_loglik_grad_adjoint; any NA entry (e.g. order
      ## 4 requested, or a not-ok adjoint block) FD-fallbacks below, so the
      ## returned gradient stays exactly consistent with lp_fn. All other
      ## grad_method values keep the (unchanged) implicit path.
      ## Orders within 1:2 (first-order rule): the exact mean/covariance
      ## gradient (.cumulant_loglik_grad_order1) serves EVERY grad_method --
      ## the order-2-rule paths below need a DecisionRules2 and used to
      ## decline wholesale here, leaving per-parameter FD of the forward.
      gll <- if (solve_order == 1L) {
        dS_list <- lapply(par_names, function(nm) .dSigma_e_d(theta, params, nm))
        names(dS_list) <- par_names
        .cumulant_loglik_grad_order1(model, compiled, dr_use, params,
                                     par_names, obs_vars, data,
                                     orders = cumulant_orders,
                                     me_variance = me_variance,
                                     dSigma_e_list = dS_list)
      } else if (grad_method == "adjoint_solution") {
        g_adj <- tryCatch(
          .cumulant_loglik_grad_adjoint(model, compiled, dr_use, params,
                                        par_names, obs_vars, data,
                                        orders = cumulant_orders,
                                        me_variance = me_variance),
          error = function(e)
            .dynhr_reraise_bug(e, setNames(rep(NA_real_, np), par_names)))
        ## If the reverse path declined wholesale (all NA — e.g. order 4 or a
        ## non-DecisionRules2), fall back to the implicit path rather than a
        ## full per-param FD, matching the accuracy of the other methods.
        if (all(is.na(g_adj))) {
          tryCatch(
            cumulant_loglik_grad(model, compiled, dr_use, params, par_names,
                                 obs_vars, data, orders = cumulant_orders,
                                 me_variance = me_variance, deriv = "implicit"),
            error = function(e)
              .dynhr_reraise_bug(e, setNames(rep(NA_real_, np), par_names)))
        } else g_adj
      } else {
        tryCatch(
          cumulant_loglik_grad(model, compiled, dr_use, params, par_names,
                               obs_vars, data, orders = cumulant_orders,
                               me_variance = me_variance, deriv = "implicit"),
          error = function(e)
            .dynhr_reraise_bug(e, setNames(rep(NA_real_, np), par_names)))
      }
      base_ll <- NULL
      for (nm in par_names) {
        gj <- gll[[nm]]
        if (!is.null(gj) && is.finite(gj)) {
          g[nm] <- g[nm] + gj
        } else {
          if (is.null(base_ll)) base_ll <- lp_fn(theta)$loglik
          d1fd <- .fd_loglik_grad1(theta, nm, base_ll)
          if (is.finite(d1fd)) g[nm] <- g[nm] + d1fd
        }
      }
      g
    }
    return(.finish_grad(cumulant_grad_fn))
  }

  ## --- "pruned" gradient closure ---------------------------------------------
  ## Order 2: analytic adjoint-chain gradient (phase a: R/pruned-kf-adjoint.R,
  ## phase b: R/pruned-grad-chain.R). One order-2 solve +
  ## solution_derivatives_order2() call covers ALL structural parameters via a
  ## shared factorization (Childers et al. efficiency point, same as the
  ## Gaussian "implicit" path); each parameter also needs a dSigma_e (the
  ## exact .dSigma_e_d, same as the Gaussian paths). Any parameter whose
  ## solution_derivatives_order2() block comes back not-ok, or whose chained gradient is non-finite, falls back to
  ## .fd_loglik_grad1() against lp_fn -- so the returned gradient is always a
  ## mix of "exact chain" and "exact FD-of-forward", never silently wrong.
  ##
  ## Order 3: the FULL fold-chain derivative (analytic d(ghxxx)/dtheta etc.
  ## plus a derivative-Lyapunov pass through the order-3 augmented system) is
  ## OUT OF SCOPE (see R/pruned-grad-chain-order3.R's file header / the D1
  ## completion report's scope section). D1 instead implements a
  ## SEMI-ANALYTIC middle path (R/pruned-grad-chain-order3.R,
  ## .pgo3_grad_chain): ONE .pruned_kf_correlated_adjoint() filter pass
  ## supplies d(loglik)/d(9 SSM inputs), and central FD of the ASSEMBLY ONLY
  ## (order-3 solve + fold + stationary moments, NO filter pass) supplies
  ## d(inputs)/dtheta_j; the two are contracted via a Frobenius inner
  ## product per parameter. This is exact up to the assembly-FD truncation
  ## error (validated against FD-of-forward to 1e-4 relative, see
  ## test-pruned-grad-order3.R), and is cheaper than FD-of-forward whenever
  ## the filter's O(T) loop dominates the (fixed-cost) assembly. Any
  ## parameter whose base assembly or per-parameter assembly-FD fails falls
  ## back to exact FD-of-forward via lp_fn, so the returned gradient is
  ## always a mix of "semi-analytic" and "exact FD-of-forward", never
  ## silently wrong.
  if (use_pruned) {
    is_pruned3 <- identical(pruned_order, 3L) || identical(pruned_order, 3)
    pruned_grad_fn <- function(theta) {
      theta <- .theta_by_name(theta, par_names, "make_posterior_grad")
      g <- .dlog_prior(theta, prior_spec)        # analytic prior score (all params)

      if (is_pruned3) {
        Y3 <- if (is.null(dim(data))) matrix(data, nrow = length(obs_vars)) else data
        if (nrow(Y3) != length(obs_vars)) Y3 <- t(Y3)

        chain_res3 <- tryCatch(
          .pgo3_grad_chain(theta, model, compiled, par_names, Y3, obs_vars,
                           me_variance = me_variance),
          error = function(e) .dynhr_reraise_bug(e, NULL))

        base_ll <- if (!is.null(chain_res3)) chain_res3$loglik else
          tryCatch(lp_fn(theta)$loglik, error = function(e) .dynhr_reraise_bug(e, -Inf))
        if (!is.finite(base_ll)) return(g)

        fd_names3 <- character(0)
        for (nm in par_names) {
          gj <- if (!is.null(chain_res3)) chain_res3$grad[[nm]] else NA_real_
          if (!is.null(gj) && is.finite(gj)) {
            g[nm] <- g[nm] + gj
          } else {
            fd_names3 <- c(fd_names3, nm)
          }
        }
        for (nm in fd_names3) {
          d1 <- .fd_loglik_grad1(theta, nm, base_ll)
          if (is.finite(d1)) g[nm] <- g[nm] + d1
        }
        return(g)
      }

      params <- .apply_theta_to_params(model, theta)
      ss_result <- tryCatch(
        solve_steady_state(model, compiled, params, verbose = FALSE),
        error = function(e) .dynhr_reraise_bug(e, NULL))
      if (is.null(ss_result) || !isTRUE(ss_result$converged)) return(g)
      params <- ss_result$params %||% params

      sys <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)
      dr1 <- tryCatch(.solve_from_system(sys, model, compiled, ss_result$ss,
                                         params, FALSE),
                      error = function(e) .dynhr_reraise_bug(e, NULL))
      if (is.null(dr1) || !isTRUE(dr1$bk_satisfied)) return(g)

      dr2 <- tryCatch(
        solve_perturbation_order2(model, compiled, ss_result$ss, params, dr1,
                                  verbose = FALSE),
        error = function(e) .dynhr_reraise_bug(e, NULL))
      if (is.null(dr2)) return(g)

      pss <- tryCatch(pruned_state_space(dr2, model, params),
                      error = function(e) .dynhr_reraise_bug(e, NULL))
      if (is.null(pss)) return(g)

      Y <- if (is.null(dim(data))) matrix(data, nrow = length(obs_vars)) else data
      if (nrow(Y) != length(obs_vars)) Y <- t(Y)

      sd2 <- tryCatch(
        solution_derivatives_order2(model, compiled, dr2, params, par_names),
        error = function(e) .dynhr_reraise_bug(e, NULL))
      if (is.null(sd2)) {
        base_ll <- tryCatch(pruned_ss_loglik(pss, Y, obs_vars, me_variance = me_variance),
                            error = function(e) .dynhr_reraise_bug(e, -Inf))
        if (!is.finite(base_ll)) return(g)
        for (nm in par_names) {
          d1 <- .fd_loglik_grad1(theta, nm, base_ll)
          if (is.finite(d1)) g[nm] <- g[nm] + d1
        }
        return(g)
      }

      dSigma_e_list <- setNames(vector("list", length(par_names)), par_names)
      for (nm in par_names) dSigma_e_list[[nm]] <- .dSigma_e_d(theta, params, nm)

      chain_res <- tryCatch(
        .pruned_ss_loglik_grad_chain(pss, Y, obs_vars, sd2, dSigma_e_list,
                                     me_variance = me_variance,
                                     model = model, compiled = compiled,
                                     ss = ss_result$ss, dr1 = dr1, params = params),
        error = function(e) .dynhr_reraise_bug(e, NULL))

      base_ll <- if (!is.null(chain_res)) chain_res$loglik else
        tryCatch(pruned_ss_loglik(pss, Y, obs_vars, me_variance = me_variance),
                 error = function(e) .dynhr_reraise_bug(e, -Inf))
      if (!is.finite(base_ll)) return(g)

      fd_names <- character(0)
      for (nm in par_names) {
        gj <- if (!is.null(chain_res)) chain_res$grad[[nm]] else NA_real_
        if (!is.null(gj) && is.finite(gj)) {
          g[nm] <- g[nm] + gj
        } else {
          fd_names <- c(fd_names, nm)
        }
      }
      if (length(fd_names) > 0) {
        for (nm in fd_names) {
          d1 <- .fd_loglik_grad1(theta, nm, base_ll)
          if (is.finite(d1)) g[nm] <- g[nm] + d1
        }
      }
      g
    }
    return(.finish_grad(pruned_grad_fn))
  }

  ## --- "hybrid" gradient closure --------------------------------------------
  hybrid_grad_fn <- function(theta) {
    ## Every closure here reads theta BY NAME into prior order (brief 31 A2):
    ## the old `names(theta) <- par_names` relabelled a PERMUTED named theta
    ## positionally, silently scrambling it. Unnamed = prior order.
    theta <- .theta_by_name(theta, par_names, "make_posterior_grad")
    g <- .dlog_prior(theta, prior_spec)        # analytic prior score (all params)

    sol <- .solve_dr(theta)
    if (is.null(sol)) return(g)                # infeasible: prior-only (rare on-path)
    dr <- sol$dr; params <- sol$params
    si <- dr$state_idx; oi <- obs_in_endo(dr)
    TT <- dr$ghx[si, , drop = FALSE]; RR <- dr$ghu[si, , drop = FALSE]
    ZZ <- dr$ghx[oi, , drop = FALSE]; DD <- dr$ghu[oi, , drop = FALSE]
    dvec <- dr$ys[obs_vars]
    Sigma_e <- .get_shock_cov(model, exo, params)
    Y <- if (is.null(dim(data))) matrix(data, nrow = length(obs_vars)) else data
    if (nrow(Y) != length(obs_vars)) Y <- t(Y)

    ## The P0 initialization the objective runs at this draw (lik_init).
    ## "reject": lik_init = "stationary" on a unit root, where the posterior
    ## is -Inf -- no likelihood score exists, return the prior score.
    ## .kf_loglik_score_sigma starts from the stationary Lyapunov P0, so it is
    ## the score of the objective only when that is the init in force; under
    ## a diffuse (lik_init "diffuse" on a near-unit root, or a true unit
    ## root) or kappa init the sigma params take FD of lp_fn instead, which
    ## carries the same lik_init. (A near-unit root under "diffuse" used to
    ## get the finite -- but stationary -- score here.)
    init_now <- .grad_init_in_force(TT, RR, Sigma_e, lik_init)
    if (identical(init_now, "reject")) return(g)
    sigma_fd <- has_tv || singular_F || !identical(init_now, "stationary")

    ## Analytic Kalman score for the sigma-like params. The FD base value is
    ## the OBJECTIVE's own loglik (lp_fn), never the score recursion's: the
    ## two are different likelihoods whenever the forward filter drops
    ## components (the singular-F univariate fallback; W68, art_zlb_mcp:
    ## 1241 vs 3297), and (lp_fn(theta + h) - <other likelihood>) / h was the
    ## 1e8-relative "gradient". The comparison also decides whether the
    ## analytic sigma score is the objective's: if not, FD for those too.
    base_ll <- NULL
    if (length(sig_names) > 0 && !sigma_fd) {
      dS_list <- lapply(sig_names, function(nm) .dSigma_e_d(theta, params, nm))
      names(dS_list) <- sig_names
      sc <- .kf_loglik_score_sigma(Y, TT, RR, ZZ, DD, Sigma_e, dvec, me_diag, dS_list)
      base_ll <- lp_fn(theta)$loglik
      .fuse_note(base_ll)
      sc_ok <- is.finite(sc$loglik) && is.finite(base_ll) &&
        abs(sc$loglik - base_ll) <= 1e-6 * max(1, abs(base_ll))
      finite_sc <- is.finite(sc$score) & sc_ok
      g[sig_names[finite_sc]] <- g[sig_names[finite_sc]] + sc$score[finite_sc]
      ## Sigma params whose analytic score came back non-finite (e.g. a
      ## (near-)unit-root TT makes the stationary Lyapunov P0 undefined) or is
      ## not the objective's: relative-step FD for just those parameters.
      if (!all(finite_sc)) {
        for (nm in sig_names[!finite_sc]) {
          d1 <- .fd_loglik_grad1(theta, nm, base_ll)
          if (is.finite(d1)) g[nm] <- g[nm] + d1
        }
      }
    }

    ## Numerical loglik score for the structural (dr-changing) params -- and,
    ## under tv inputs, also for the sigma params (the analytic score above is
    ## skipped; FD against lp_fn is self-consistent because lp_fn carries the
    ## same me_extra/shock_scale).
    fd_names <- if (sigma_fd) c(sig_names, num_names) else num_names
    if (length(fd_names) > 0) {
      if (is.null(base_ll)) {
        base_ll <- lp_fn(theta)$loglik
        .fuse_note(base_ll)
      }
      for (nm in fd_names) {
        d1 <- .fd_loglik_grad1(theta, nm, base_ll)
        if (is.finite(d1)) g[nm] <- g[nm] + d1
      }
    }
    g
  }

  ## --- Whittle gradient closure (all grad_method values) -------------------
  ## For the Whittle likelihood the tangent/adjoint Kalman filter is replaced
  ## by .whittle_loglik_grad.  Three paths:
  ##
  ##  * "hybrid":  sigma-movers -> .whittle_loglik_grad with dSigma_e only
  ##               num_names    -> FD via lp_fn (same as Gaussian hybrid)
  ##
  ##  * "implicit"/"adjoint":
  ##               all params   -> d_ss_list built from solution_derivatives
  ##                               (sigma-movers get dSigma_e AND zero dTT/...)
  ##               FD fallback  -> lp_fn for any param where sd_res$ok==FALSE
  ##
  ## Under whittle, has_tv is always FALSE (whittle does not support me_extra/
  ## shock_scale; the ctx_allows_analytic_gradient guard ensures this).
  if (use_whittle) {
    ## Pre-compute the periodogram once (data is fixed for this grad builder).
    wdata <- if (is.null(dim(data))) matrix(data, ncol = 1L) else data
    if (nrow(wdata) < ncol(wdata)) wdata <- t(wdata)
    obs_col_idx <- match(obs_vars, colnames(wdata))
    if (anyNA(obs_col_idx)) obs_col_idx <- seq_along(obs_vars)
    Y_raw_w    <- wdata[, obs_col_idx, drop = FALSE]
    col_means_w <- colMeans(Y_raw_w, na.rm = TRUE)
    Y_dem_w    <- sweep(Y_raw_w, 2, col_means_w, "-")
    w_pdgm     <- .whittle_periodogram(Y_dem_w)

    whittle_grad_fn <- function(theta) {
      theta <- .theta_by_name(theta, par_names, "make_posterior_grad")
      g <- .dlog_prior(theta, prior_spec)

      sol <- .solve_dr(theta)
      if (is.null(sol)) return(g)
      dr <- sol$dr; params <- sol$params
      si <- dr$state_idx; oi <- obs_in_endo(dr)
      Sigma_e <- .get_shock_cov(model, exo, params)
      ## Construct dsge_ss at the Whittle-gradient boundary, then unpack to
      ## locals for .whittle_loglik_grad (bare-matrix hot path; no S3 dispatch).
      ss_wg <- new_dsge_ss(
        T_mat   = dr$ghx[si, , drop = FALSE],
        R_mat   = dr$ghu[si, , drop = FALSE],
        Z_mat   = dr$ghx[oi, , drop = FALSE],
        D_mat   = dr$ghu[oi, , drop = FALSE],
        Sigma_e = Sigma_e,
        timing  = "lagged"
      )
      TT <- ss_wg$T_mat; RR <- ss_wg$R_mat
      ZZ <- ss_wg$Z_mat; DD <- ss_wg$D_mat

      ## Build d_ss_list
      d_ss_list <- vector("list", np)
      names(d_ss_list) <- par_names

      ## Sigma-movers: dSigma_e from FD; TT/RR/ZZ/DD derivatives are zero
      ## (certainty equivalence). .whittle_loglik_grad fills zeros automatically
      ## when dTT/dRR/dZZ/dDD are absent from the slot list.
      for (nm in sig_names)
        d_ss_list[[nm]] <- list(dSigma_e = .dSigma_e_d(theta, params, nm))

      ## Solution-movers
      fd_names <- character(0)
      if (grad_method == "hybrid" || length(num_names) == 0L) {
        ## Hybrid: set all num_names to NULL -> fall through to FD below
        for (nm in num_names) d_ss_list[nm] <- list(NULL)
        fd_names <- num_names
      } else {
        ## implicit/adjoint: use solution_derivatives for all num_names
        sd_res <- tryCatch(
          solution_derivatives(model, compiled, dr, params,
                               param_names = num_names, obs_vars = obs_vars),
          error = function(e) .dynhr_reraise_bug(e, NULL)
        )
        for (nm in num_names) {
          d <- if (!is.null(sd_res)) sd_res$derivs[[nm]] else NULL
          if (!is.null(d) && isTRUE(d$ok)) {
            d_ss_list[[nm]] <- list(
              dTT = d$dTT, dRR = d$dRR, dZZ = d$dZZ, dDD = d$dDD,
              dSigma_e = .dSigma_e_d(theta, params, nm)
            )
          } else {
            d_ss_list[nm] <- list(NULL)
            fd_names <- c(fd_names, nm)
          }
        }
      }

      ## Debiased gradient: precompute EI_list, then dEI per active param.
      ## Requires P0 and K_mat; recomputed here at the same cost as the
      ## loglik-side ctau call (solve_lyapunov + T matrix products).
      EI_list_grad   <- NULL
      dEI_per_param  <- NULL
      use_debias_grad <- isTRUE(debias)
      if (use_debias_grad) {
        ## A NUMERICAL failure (e.g. the Lyapunov solve at a near-unit-root
        ## draw) drops to the non-debiased gradient; a programming error is
        ## re-raised -- the 0.9.3.84 index-shift bug hid here for a release.
        debias_ok <- tryCatch({
          T_len_w <- nrow(Y_dem_w)
          c_arr_g <- .whittle_compute_ctau(TT, RR, ZZ, DD, Sigma_e,
                                            T_len_w, me_variance)
          EI_list_grad <- .whittle_compute_EI(c_arr_g, w_pdgm$omega, T_len_w)
          ## P0 and K_mat for dctau
          Q_rr  <- RR %*% Sigma_e %*% t(RR)
          P0_g  <- solve_lyapunov(TT, Q_rr)
          K_g   <- TT %*% P0_g %*% t(ZZ) + RR %*% Sigma_e %*% t(DD)
          ## Per-param dEI
          dEI_per_param <- vector("list", np)
          names(dEI_per_param) <- par_names
          n_state_w <- nrow(TT); n_exo_w <- ncol(RR); n_obs_w <- nrow(ZZ)
          zero_TT_w <- matrix(0, n_state_w, n_state_w)
          zero_RR_w <- matrix(0, n_state_w, n_exo_w)
          zero_ZZ_w <- matrix(0, n_obs_w, n_state_w)
          zero_DD_w <- matrix(0, n_obs_w, n_exo_w)
          zero_Se_w <- matrix(0, n_exo_w, n_exo_w)
          for (k_nm in par_names) {
            d_k <- d_ss_list[[k_nm]]
            if (is.null(d_k)) {
              dEI_per_param[k_nm] <- list(NULL)  ## FD fallback: keep the slot ([[<- NULL would DELETE it and shift every later index)
              next
            }
            dTT_k <- d_k$dTT       %||% zero_TT_w
            dRR_k <- d_k$dRR       %||% zero_RR_w
            dZZ_k <- d_k$dZZ       %||% zero_ZZ_w
            dDD_k <- d_k$dDD       %||% zero_DD_w
            dSe_k <- d_k$dSigma_e  %||% zero_Se_w
            dcArr_k <- .whittle_compute_dctau(TT, RR, ZZ, DD, Sigma_e,
                                               P0_g, K_g,
                                               dTT_k, dRR_k, dZZ_k, dDD_k,
                                               dSe_k, T_len_w)
            dEI_per_param[[k_nm]] <- .whittle_compute_dEI(
              dcArr_k, w_pdgm$omega, T_len_w)
          }
          TRUE
        }, error = function(e) .dynhr_reraise_bug(e, FALSE))
        if (!isTRUE(debias_ok)) {
          EI_list_grad  <- NULL
          dEI_per_param <- NULL
          use_debias_grad <- FALSE
        }
      }

      ## Analytic Whittle gradient (NA entries -> FD)
      wg <- .whittle_loglik_grad(w_pdgm, TT, RR, ZZ, DD, Sigma_e,
                                  d_ss_list, freq_band = freq_band,
                                  me_variance = me_variance,
                                  debias = use_debias_grad,
                                  EI_list = EI_list_grad,
                                  dEI_per_param = dEI_per_param)
      ## Accumulate finite analytic entries
      finite_wg <- is.finite(wg)
      g[names(wg)[finite_wg]] <- g[names(wg)[finite_wg]] + wg[finite_wg]

      ## FD fallback for NULL-slot params + non-finite analytic entries
      fd_names <- c(fd_names, names(wg)[!is.na(wg) & !is.finite(wg)])
      all_fd   <- c(fd_names, par_names[na_slots_at_call(d_ss_list)])

      if (length(all_fd) > 0L) {
        base_ll <- lp_fn(theta)$loglik
        for (nm in unique(all_fd)) {
          d1 <- .fd_loglik_grad1(theta, nm, base_ll)
          if (is.finite(d1)) g[nm] <- g[nm] + d1
        }
      }
      g
    }
    ## Helper: which d_ss_list slots are NULL (NA in grad -> FD needed)?
    na_slots_at_call <- function(dssl) {
      vapply(dssl, is.null, logical(1))
    }
    return(.finish_grad(whittle_grad_fn))
  }

  if (grad_method == "hybrid")
    return(.attach_fused(.finish_grad(hybrid_grad_fn)))

  ## --- "implicit" gradient closure ------------------------------------------
  ## One model solve + one tangent Kalman-filter pass covers ALL parameters:
  ##   * sig_names  -> d_ss_list entry with ONLY dSigma_e set (certainty
  ##                    equivalence: TT/RR/ZZ/DD/d unchanged).
  ##   * num_names  -> solution_derivatives() blocks (dTT/dRR/dZZ/dDD/dd),
  ##                    dSigma_e = NULL (one shared Sylvester/QR factorization
  ##                    across all of num_names).
  ## Per-parameter robustness: any parameter with ok == FALSE or a non-finite
  ## tangent-filter gradient entry falls back to .fd_loglik_grad1().
  ##
  ## Non-stationary init / missing data: the dense tangent and adjoint kernels
  ## assume a stationary P0 = solve_lyapunov(TT, QQ) and complete data. Such
  ## draws are detected per-draw (mirrors kalman_filter(..., lik_init = "auto"))
  ## and take the exact FD-"hybrid" fallback. (Previously only unit roots were
  ## caught; missing data reached the dense adjoint, which errors on NA.)
  warned_special <- FALSE

  ## F5: lightweight fallback-usage counters, so a chain silently degrading
  ## from the analytic Kalman-adjoint/tangent kernel to the FD-hybrid path
  ## leaves a trace (fallback draws do NOT register as NUTS divergences and,
  ## after the first warning(), the `warned_special` latch above goes silent
  ## for every later occurrence in the chain). Attached to the returned
  ## closure as attr(grad_fn, "kernel_stats") -- read with as.list(). A
  ## healthy chain has every n_fallback_*/n_score_nonfinite count at 0 and
  ## n_calls > 0. Uses env-slot assignment (kernel_stats$x <- ...), NOT
  ## `<<-` -- test-function-uniqueness.R's allow-list caps `<<-` uses in
  ## this file at 3 (the three `warned_special <<- TRUE` latches).
  kernel_stats <- new.env(parent = emptyenv())
  kernel_stats$n_calls                 <- 0L
  kernel_stats$n_fallback_special_init <- 0L  # non-stationary init/missing-data combo no kernel covers
  kernel_stats$n_fallback_diffuse      <- 0L  # exact-diffuse adjoint kernel unavailable for the draw
  kernel_stats$n_fallback_nondiffuse   <- 0L  # tangent/adjoint (non-diffuse) kernel threw/failed
  kernel_stats$n_score_nonfinite       <- 0L  # per-parameter non-finite analytic score -> FD fallback

  .impl_adj_grad_fn <- function(theta) {
    kernel_stats$n_calls <- kernel_stats$n_calls + 1L
    theta <- .theta_by_name(theta, par_names, "make_posterior_grad")
    g <- .dlog_prior(theta, prior_spec)        # analytic prior score (all params)

    sol <- .solve_dr(theta)
    if (is.null(sol)) return(g)                # infeasible: prior-only (rare on-path)
    dr <- sol$dr; params <- sol$params
    si <- dr$state_idx; oi <- obs_in_endo(dr)
    TT <- dr$ghx[si, , drop = FALSE]; RR <- dr$ghu[si, , drop = FALSE]
    ZZ <- dr$ghx[oi, , drop = FALSE]; DD <- dr$ghu[oi, , drop = FALSE]
    dvec <- dr$ys[obs_vars]
    Sigma_e <- .get_shock_cov(model, exo, params)
    Y <- if (is.null(dim(data))) matrix(data, nrow = length(obs_vars)) else data
    if (nrow(Y) != length(obs_vars)) Y <- t(Y)

    ## Per-draw guard (mirrors kalman_filter's lik_init = "auto" eigenvalue test):
    ## the dense tangent/adjoint kernels assume a stationary Lyapunov P0 and
    ## complete data. Unit-root models take the exact FD-"hybrid" fallback
    ## (which differentiates the forward likelihood directly, so it is always
    ## consistent). Missing data on a STATIONARY model is handled by the
    ## missing-data univariate adjoint kernel (.kf_loglik_adjoint_uni): it does
    ## per-period observed-row subsetting on the stationary Lyapunov P0 and is
    ## forward-consistent with the production KF to machine precision (Tier 14;
    ## the orphaning rationale -- a ~0.5 nat loglik gap -- was a production
    ## steady-state-lock-vs-fully-missing bug, now fixed). That kernel supports
    ## neither me_extra nor shock_scale, so those (and unit-root draws) still
    ## take the exact FD-hybrid path.
    has_missing <- anyNA(Y)
    ## Mirror kalman_filter(lik_init = "auto") EXACTLY: a near-/at-unit root only
    ## forces the diffuse path when the stationary Lyapunov P0 is non-finite or
    ## non-PSD (a true unit/explosive root). For roots in [1-1e-6, 1) the filter
    ## now uses the stationary init, so the gradient MUST use the stationary
    ## adjoint too -- otherwise the gradient would be d(diffuse loglik)/dtheta
    ## while the loglik is the stationary one, breaking newrat's line search.
    ## Explicit lik_init (W68): "diffuse" is the diffuse init on every root
    ## within 1e-6 of the unit circle (not only where the Lyapunov P0 fails,
    ## as under "auto"); "stationary" on a unit root is a -Inf posterior, so
    ## the prior score is all there is. ("kappa" never reaches this closure:
    ## refused at build time.)
    init_now <- .grad_init_in_force(TT, RR, Sigma_e, lik_init)
    if (identical(init_now, "reject")) return(g)
    has_unit_root <- identical(init_now, "diffuse")
    ## Specialized kernels cover two of the special-init cases exactly:
    ##  * stationary + missing data            -> .kf_loglik_adjoint_uni
    ##    (supports shock_scale since 0.9.0.0008 -- per-period Se_t in the R
    ##    reference, sandwiched Sigma_e adjoint; validated vs numDeriv in
    ##    test-gradient-shock-scale-wired.R. me_extra remains excluded.)
    ##  * unit-root (diffuse) + complete data  -> .kf_loglik_adjoint_diffuse
    ##    (no me_extra/shock_scale; needs complete data through the diffuse phase)
    ## Any remaining special case (unit-root + missing; special-init + me_extra;
    ## diffuse + shock_scale) takes the exact FD-hybrid path, which
    ## differentiates the forward likelihood directly.
    use_uni_adjoint <- has_missing && !has_unit_root &&
      is.null(me_extra)
    use_diffuse_adjoint <- has_unit_root && !has_missing &&
      is.null(me_extra) && is.null(shock_scale)
    ## == (has_unit_root || has_missing) && !use_uni_adjoint &&
    ##    !use_diffuse_adjoint; shared with the "auto" resolution.
    if (.grad_no_kernel_covers(has_unit_root, has_missing, me_extra, shock_scale)) {
      kernel_stats$n_fallback_special_init <- kernel_stats$n_fallback_special_init + 1L
      if (!warned_special) {
        .dynhr_warn("make_posterior_grad: ", grad_method, " gradient with a ",
                "non-stationary init and/or (me_extra/shock_scale +) missing ",
                "data falls back to the exact FD-hybrid path for this case.",
                call. = FALSE)
        warned_special <<- TRUE
      }
      return(hybrid_grad_fn(theta))
    }

    ## Construct dsge_ss at the KF-gradient boundary (timing = "lagged").
    ## Add TT/RR/ZZ/DD aliases after construction so the gradient KF functions
    ## (.kf_loglik_tangent / .kf_loglik_adjoint) can read ss$TT etc. without
    ## touching their field-access code or the C++ backends.
    ss <- new_dsge_ss(
      T_mat   = TT,
      R_mat   = RR,
      Z_mat   = ZZ,
      D_mat   = DD,
      Sigma_e = Sigma_e,
      d       = dvec,
      timing  = "lagged"
    )
    ss$TT <- ss$T_mat; ss$RR <- ss$R_mat
    ss$ZZ <- ss$Z_mat; ss$DD <- ss$D_mat

    ## Full reverse mode ("adjoint_solution"): the dense adjoint kernel exports
    ## its bar matrices and .solution_adjoint() replaces the per-parameter
    ## forward Sylvester solves. The uni/diffuse special-case kernels do not
    ## export bars, so those draws use the "adjoint" construction instead.
    use_sol_adjoint <- grad_method == "adjoint_solution" &&
      !use_uni_adjoint && !use_diffuse_adjoint

    ## Build d_ss_list: sigma params get dSigma_e only; all other params get
    ## solution_derivatives() blocks (one shared call/factorization). Under
    ## use_sol_adjoint the structural blocks are NOT built -- num_names carry
    ## only their dSigma_e channel (stderr-expression coupling; exact-zero for
    ## parameters absent from the shocks block) and the structural gradient
    ## comes from .solution_adjoint() below.
    d_ss_list <- vector("list", np)
    names(d_ss_list) <- par_names

    for (nm in sig_names)
      d_ss_list[[nm]] <- list(dSigma_e = .dSigma_e_d(theta, params, nm))

    sd_res <- NULL
    if (use_sol_adjoint) {
      for (nm in num_names)
        d_ss_list[[nm]] <- list(dSigma_e = .dSigma_e_d(theta, params, nm))
    } else if (length(num_names) > 0) {
      sd_res <- tryCatch(
        solution_derivatives(model, compiled, dr, params,
                              param_names = num_names, obs_vars = obs_vars),
        error = function(e) .dynhr_reraise_bug(e, NULL)
      )
      for (nm in num_names) {
        d <- if (!is.null(sd_res)) sd_res$derivs[[nm]] else NULL
        if (!is.null(d) && isTRUE(d$ok)) {
          ## dSigma_e is NOT necessarily zero for solution-moving parameters:
          ## since shocks-block stderr/variance expressions re-evaluate
          ## against the current params, a parameter like sig_i in
          ## `i = ... + sig_i*eps_i;` PLUS `stderr sig_i;` moves both the
          ## decision rule and Sigma_e. The FD is exact-zero for parameters
          ## absent from the shocks block, and cheap either way.
          d_ss_list[[nm]] <- list(dTT = d$dTT, dRR = d$dRR, dZZ = d$dZZ,
                                   dDD = d$dDD, dd = d$dd,
                                   dSigma_e = .dSigma_e_d(theta, params, nm))
        } else {
          ## Keep the slot (x[[nm]] <- NULL would DELETE it and misalign the
          ## positional grad vector for every later parameter).
          d_ss_list[nm] <- list(NULL)   # zero block; FD fallback applied below
        }
      }
    }

    if (use_diffuse_adjoint) {
      ## Unit-root draw: exact-diffuse adjoint (Durbin-Koopman ch.5). Returns the
      ## exact within-regime gradient for unit-root-preserving TT directions
      ## (estimated params move TT within the smooth regime). If the diffuse
      ## init/recursion gradient is unavailable for this draw (stage flags or a
      ## non-finite loglik), degrade to the FD-hybrid path rather than trust a
      ## partial gradient.
      tang <- tryCatch(
        .kf_loglik_adjoint_diffuse(Y, ss, d_ss_list, me_variance = me_variance),
        error = function(e) .dynhr_reraise_bug(e, NULL))
      if (is.null(tang) || !isTRUE(tang$stage1_ok) || !isTRUE(tang$stage2_ok) ||
          !is.finite(tang$loglik)) {
        kernel_stats$n_fallback_diffuse <- kernel_stats$n_fallback_diffuse + 1L
        if (!warned_special) {
          .dynhr_warn("make_posterior_grad: exact-diffuse adjoint gradient ",
                  "unavailable for this draw; using the FD-hybrid fallback.",
                  call. = FALSE)
          warned_special <<- TRUE
        }
        return(hybrid_grad_fn(theta))
      }
    } else {
      ## Non-diffuse adjoint/tangent kernels: the dense KF kernels already
      ## degrade gracefully (ok = FALSE -> loglik = -Inf, grad = NA) for an
      ## infeasible draw. But a sufficiently extreme theta (e.g. a NUTS
      ## step-size search trial) can still throw from deeper in the C++
      ## kernel (Sylvester/solution-derivative solves in .solution_adjoint /
      ## solution_derivatives below, or an unanticipated numerical failure);
      ## wrap the whole dispatch in tryCatch, mirroring use_diffuse_adjoint
      ## above, so every grad_method degrades to the FD-hybrid fallback
      ## instead of crashing the sampler. Programming errors (subscript /
      ## argument / object errors) are re-raised, not degraded: the 0.9.3.85
      ## adjoint-uni "[[<- NULL" bug was masked by this handler.
      tang <- tryCatch({
        if (use_uni_adjoint) {
          ## Missing-data multivariate adjoint with per-period observed-row
          ## subsetting; stationary Lyapunov P0 (P0 = NULL). me_variance is scalar.
          ## shock_scale (when present) routes to the kernel's R reference,
          ## which substitutes Se_t per period (see gradient-adjoint-uni.R).
          .kf_loglik_adjoint_uni(Y, ss, d_ss_list, me_variance = me_variance,
                                 P0 = NULL, shock_scale = shock_scale)
        } else if (use_sol_adjoint) {
          ## Dense stationary adjoint WITH bar export (the _ss long-T variant does
          ## not export bars; the O(T) n^2 storage is accepted here).
          .kf_loglik_adjoint(Y, ss, d_ss_list, me_variance = me_variance,
                             me_extra = me_extra, shock_scale = shock_scale,
                             return_bars = TRUE)
        } else if (grad_method == "adjoint" || grad_method == "adjoint_solution") {
          ## Stationary dense adjoint. For LONG samples the steady-state-aware
          ## variant (.kf_loglik_adjoint_ss) cuts n^2-matrix storage from O(T) to
          ## O(t_conv) once the Riccati recursion converges; it is identical to the
          ## dense kernel to <1e-9 (test-gradient-adjoint-ss.R) but does not support
          ## me_extra/shock_scale. Use it only when it actually saves memory (T large)
          ## and those features are off; the dense path stays the default otherwise.
          use_ss_adjoint <- ncol(Y) >= 1000L &&
            is.null(me_extra) && is.null(shock_scale)
          if (use_ss_adjoint) {
            .kf_loglik_adjoint_ss(Y, ss, d_ss_list, me_variance = me_variance)
          } else {
            .kf_loglik_adjoint(Y, ss, d_ss_list, me_variance = me_variance,
                               me_extra = me_extra, shock_scale = shock_scale)
          }
        } else {
          .kf_loglik_tangent(Y, ss, d_ss_list, me_variance = me_variance,
                             me_extra = me_extra, shock_scale = shock_scale)
        }
      }, error = function(e) .dynhr_reraise_bug(e, NULL))

      if (is.null(tang)) {
        kernel_stats$n_fallback_nondiffuse <- kernel_stats$n_fallback_nondiffuse + 1L
        if (!warned_special) {
          .dynhr_warn("make_posterior_grad: ", grad_method, " gradient kernel ",
                  "failed for this draw (numerical KF failure); using the ",
                  "FD-hybrid fallback.", call. = FALSE)
          warned_special <<- TRUE
        }
        return(hybrid_grad_fn(theta))
      }
    }
    base_ll <- tang$loglik
    ## The analytic kernel's loglik IS the forward filter's (every kernel that
    ## reaches here is validated forward-consistent; the sampler-side fusion
    ## re-checks it against the objective at its start point). A non-finite
    ## one is replaced by the objective's own below, or noted by
    ## .fused_of()'s fallback.
    if (is.finite(base_ll)) .fuse_note(base_ll)

    fd_names <- character(0)
    if (use_sol_adjoint && length(num_names) > 0) {
      ## Reverse mode through the solve: contract the exported bars against
      ## the analytic primitive derivatives. tang$grad[nm] holds ONLY the
      ## Sigma_e channel for these parameters (their d_ss_list entries carry
      ## just dSigma_e); the structural piece comes from .solution_adjoint.
      sol_adj <- if (!is.null(tang$bars)) tryCatch(
        .solution_adjoint(model, compiled, dr, params,
                          param_names = num_names, obs_vars = obs_vars,
                          bars = tang$bars),
        error = function(e) .dynhr_reraise_bug(e, NULL)
      ) else NULL
      for (nm in num_names) {
        gj_sig <- tang$grad[[match(nm, par_names)]]
        gj_str <- if (!is.null(sol_adj) && isTRUE(sol_adj$ok[[nm]]))
          sol_adj$grad[[nm]] else NA_real_
        if (is.finite(gj_sig) && is.finite(gj_str)) {
          g[nm] <- g[nm] + gj_sig + gj_str
        } else {
          fd_names <- c(fd_names, nm)
        }
      }
    } else if (length(num_names) > 0) {
      for (nm in num_names) {
        d <- if (!is.null(sd_res)) sd_res$derivs[[nm]] else NULL
        ok <- !is.null(d) && isTRUE(d$ok)
        gj <- if (ok) tang$grad[[match(nm, par_names)]] else NA_real_
        if (ok && is.finite(gj)) {
          g[nm] <- g[nm] + gj
        } else {
          fd_names <- c(fd_names, nm)
        }
      }
    }

    if (length(sig_names) > 0) {
      for (nm in sig_names) {
        gj <- tang$grad[[match(nm, par_names)]]
        if (is.finite(gj)) {
          g[nm] <- g[nm] + gj
        } else {
          fd_names <- c(fd_names, nm)
        }
      }
    }

    ## Per-parameter FD fallback (only for parameters that need it).
    if (length(fd_names) > 0) {
      kernel_stats$n_score_nonfinite <- kernel_stats$n_score_nonfinite + length(fd_names)
      if (!is.finite(base_ll)) {
        base_ll <- lp_fn(theta)$loglik
        .fuse_note(base_ll)
      }
      for (nm in fd_names) {
        d1 <- .fd_loglik_grad1(theta, nm, base_ll)
        if (is.finite(d1)) g[nm] <- g[nm] + d1
      }
    }

    g
  }

  attr(.impl_adj_grad_fn, "kernel_stats") <- kernel_stats
  .attach_fused(.finish_grad(.impl_adj_grad_fn))
}


## TRUE when the stationary Lyapunov P0 of (TT, RR Sigma_e RR') does not exist
## or is not PSD -- a true unit/explosive root. Mirrors kalman_filter(lik_init =
## "auto") exactly: a root in [1 - 1e-6, 1) with a valid stationary P0 keeps
## the stationary init, so the gradient kernels must too. (Moved verbatim from
## the tangent/adjoint closure of make_posterior_grad(); the "auto" resolution
## applies the same test at theta_ref.)
.grad_needs_diffuse_init <- function(TT, RR, Sigma_e) {
  if (nrow(TT) == 0 ||
      !any(Mod(eigen(TT, symmetric = FALSE, only.values = TRUE)$values) > 1 - 1e-6))
    return(FALSE)
  QQ_g <- RR %*% Sigma_e %*% t(RR)
  P0_g <- tryCatch(solve_lyapunov(TT, QQ_g), error = function(e) NULL)
  ## the forward filter's (relative) validity rule, shared (W77)
  !.kf_stationary_P0_ok(P0_g)
}


## The P0 initialization make_log_posterior(lik_init = ...) runs at a draw
## with state transition TT, for the gradient kernels to follow:
##   "stationary" -- the Lyapunov P0 (stationary tangent/adjoint kernels);
##   "diffuse"    -- the exact diffuse (P_inf, P_star) init (diffuse adjoint);
##   "kappa"      -- the big-kappa P0 (no analytic kernel; FD only);
##   "reject"     -- lik_init = "stationary" on a unit root: the posterior is
##                   -Inf there (make_log_posterior's stationarity guard).
## "auto" is kalman_filter()'s rule (diffuse only where the Lyapunov P0 is not
## a valid covariance, .grad_needs_diffuse_init); an explicit "diffuse" is
## .kf_diffuse_P0()'s unit-root block, every root with ||lambda| - 1| < 1e-6,
## which includes the still-stationary roots in (1 - 1e-6, 1) that "auto"
## keeps stationary.
.grad_init_in_force <- function(TT, RR, Sigma_e, lik_init = "auto") {
  if (identical(lik_init, "kappa")) return("kappa")
  if (identical(lik_init, "auto"))
    return(if (.grad_needs_diffuse_init(TT, RR, Sigma_e)) "diffuse"
           else "stationary")
  if (nrow(TT) == 0) return("stationary")
  mods <- Mod(eigen(TT, symmetric = FALSE, only.values = TRUE)$values)
  if (identical(lik_init, "stationary"))
    return(if (max(mods) >= 1) "reject" else "stationary")
  if (any(abs(mods - 1) < 1e-6)) "diffuse" else "stationary"
}


## TRUE when no analytic (tangent / adjoint) Kalman kernel covers a draw, so
## the implicit/adjoint/adjoint_solution closures take the exact FD-"hybrid"
## path for it: the missing-data kernel needs a stationary init and no
## me_extra; the exact-diffuse kernel needs complete data and neither me_extra
## nor shock_scale.
.grad_no_kernel_covers <- function(has_unit_root, has_missing, me_extra,
                                   shock_scale) {
  use_uni     <- has_missing && !has_unit_root && is.null(me_extra)
  use_diffuse <- has_unit_root && !has_missing &&
    is.null(me_extra) && is.null(shock_scale)
  (has_unit_root || has_missing) && !use_uni && !use_diffuse
}


## Does the Gaussian forward drop observation components at the base solve
## (kalman_filter's singular-F univariate fallback skips every component whose
## conditional variance is <= the ABSOLUTE kalman_tol; $diagnostics$n_dropped)
## while the dense analytic kernel evaluates a FINITE, different likelihood
## there? Returns list(n_dropped, ll_forward, ll_kernel) in that case, NULL
## otherwise. Two cases are deliberately NOT flagged:
##  * an EXACT singularity (an observable that is an exact combination of
##    others, no measurement error): the dense kernel fails (non-finite) on
##    it, so every draw already takes the exact FD fallback of the forward,
##    whose dropped set is stable (test-gradient-shock-scale-wired.R W3);
##  * a draw outside the dense stationary kernel's domain (diffuse init,
##    missing data): those draws take the FD-hybrid path anyway.
## Flagged: a singularity by SCALE that Dynare's rule also routes to the
## univariate filter (rcond(F) < kalman_tol with a component variance below
## kalman_tol, e.g. one observable of variance ~1e-12 beside O(1) ones) --
## chol succeeds, the kernel returns the undropped likelihood, and the
## forward drops the tiny component. (A well-conditioned F of small overall
## scale -- art_zlb_mcp, observable variances ~1e-9 -- is no longer routed
## there since W74, so it is not flagged.) One filter pass (plus one kernel pass
## when components are dropped) at build time; the forward's own fallback
## warning is superseded by the builder's classed one.
.grad_singular_F_mismatch <- function(model, data, obs_vars, base,
                                      me_variance = 0, lik_init = "auto",
                                      me_extra = NULL, shock_scale = NULL) {
  if (is.null(base)) return(NULL)
  dr <- base$dr; si <- dr$state_idx; oi <- match(obs_vars, dr$endo_names)
  TT <- dr$ghx[si, , drop = FALSE]; RR <- dr$ghu[si, , drop = FALSE]
  Sigma_e <- .get_shock_cov(model, model$varexo_names, base$params)
  init <- .grad_init_in_force(TT, RR, Sigma_e, lik_init)
  if (!identical(init, "stationary")) return(NULL)
  Y <- if (is.null(dim(data))) matrix(data, nrow = length(obs_vars)) else data
  if (nrow(Y) != length(obs_vars)) Y <- t(Y)
  if (anyNA(Y)) return(NULL)
  kf <- suppressWarnings(
    kalman_filter(data, dr, model, base$params, obs_vars,
                  return_filtered = FALSE, me_variance = me_variance,
                  lik_init = lik_init, me_extra = me_extra,
                  shock_scale = shock_scale))
  nd <- as.integer(kf$diagnostics$n_dropped %||% 0L)
  if (nd == 0L || !is.finite(kf$loglik)) return(NULL)
  ss <- list(TT = TT, RR = RR, ZZ = dr$ghx[oi, , drop = FALSE],
             DD = dr$ghu[oi, , drop = FALSE], d = dr$ys[obs_vars],
             Sigma_e = Sigma_e)
  zero <- list(dSigma_e = matrix(0, ncol(RR), ncol(RR)))
  kern <- .kf_loglik_adjoint(Y, ss, list(zero), me_variance = me_variance,
                             me_extra = me_extra, shock_scale = shock_scale)
  if (!is.finite(kern$loglik) ||
      abs(kern$loglik - kf$loglik) <= 1e-6 * max(1, abs(kf$loglik)))
    return(NULL)
  list(n_dropped = nd, ll_forward = kf$loglik, ll_kernel = kern$loglik)
}


## First-order solve at theta for the gradient builders: list(dr, params), or
## NULL when the steady state or the Blanchard-Kahn condition fails.
.grad_solve_dr <- function(model, compiled, sys_cache, theta) {
  params <- .apply_theta_to_params(model, theta)
  ss <- solve_steady_state(model, compiled, params, verbose = FALSE)
  if (is.null(ss) || !isTRUE(ss$converged)) return(NULL)
  ## Re-derive any steady_state_model-computed parameter so the base decision
  ## rule (and the returned params) are consistent with the re-solved steady
  ## state (no-op for non-SSM-parameter models; Tier 13 #1).
  params <- ss$params %||% params
  sys <- extract_system_matrices_fast(sys_cache, ss$ss, params)
  dr  <- .solve_from_system(sys, model, compiled, ss$ss, params, FALSE)
  if (is.null(dr) || !isTRUE(dr$bk_satisfied)) return(NULL)
  list(dr = dr, params = params)
}


## The method make_posterior_grad() would build for these inputs, without
## building it: `grad_method` itself unless it is "auto", else the
## resolution at theta_ref (default: the prior means, as make_posterior_grad
## uses) -- one first-order solve. Used where the closure is built elsewhere
## (the mirai NUTS daemons) but the method must still be reported/recorded.
.grad_method_resolved <- function(model, data, prior_spec, obs_vars, compiled,
                                  grad_method = "auto",
                                  likelihood = "gaussian", me_extra = NULL,
                                  shock_scale = NULL, theta_ref = NULL,
                                  cumulant_orders = 1:4,
                                  lik_init = "auto", me_variance = 0) {
  if (!identical(grad_method, "auto")) return(grad_method)
  if (is.null(theta_ref))
    theta_ref <- setNames(prior_spec$mean, prior_spec$name)
  base <- .grad_solve_dr(model, compiled, cache_system_structure(compiled),
                         theta_ref)
  singular_F <- identical(likelihood, "gaussian") &&
    !is.null(.grad_singular_F_mismatch(model, data, obs_vars, base,
                                       me_variance, lik_init, me_extra,
                                       shock_scale))
  .resolve_auto_grad_method(likelihood, base = base, model = model,
                            data = data, me_extra = me_extra,
                            shock_scale = shock_scale,
                            cumulant_orders = cumulant_orders,
                            lik_init = lik_init, singular_F = singular_F)
}


## "adjoint_solution", or "auto -> adjoint_solution" when it was resolved.
.grad_method_label <- function(resolved, requested) {
  if (is.null(resolved) || identical(resolved, requested)) return(requested)
  paste0(requested, " -> ", resolved)
}


## Resolve make_posterior_grad(grad_method = "auto") to a concrete method.
##
## Rule (W62 benchmark, 2026-09-26; median of >= 20 interleaved calls per
## method, accuracy against a five-point central difference):
##   * gaussian -> "adjoint_solution": fastest (or tied within timing noise)
##     of the exact methods on EVERY model measured -- AR(1) (2 params,
##     1 state), rbc2shock (7, 3), nk_demo (9, 4), fs2000 (9, 4), ireland_2004
##     (12, 5), ls2003 (17, 6), SW2007 (36, 20: 16 ms vs adjoint 32, implicit
##     91, hybrid 97) and NZSIM (68, 37: 45 ms vs 120 / 391 / 316), and at
##     T = 2000 (where "adjoint" switches to its slower steady-state kernel).
##     No size threshold is warranted: at the smallest size all four methods
##     cost the same R-level overhead (~1.4 ms).
##     Exception -> "hybrid": when at theta_ref no analytic kernel covers the
##     data / init combination (.grad_no_kernel_covers()), every draw would
##     take the FD-hybrid fallback anyway (with a warning); hybrid is then the
##     same gradient without the wasted kernel attempt.
##   * cumulant -> "adjoint_solution" when the third cumulant is matched and
##     the fourth is not (3 in cumulant_orders, 4 not): the reverse-mode path
##     covers exactly orders 1-3 on the order-2 rule the forward solves there.
##     W65 benchmark (median of 5-7 interleaved calls, T = 150-200, accuracy
##     against a five-point central difference of the log-posterior):
##     rbc2shock 4 structural params, orders 1:3: 22 ms vs implicit 1921 ms
##     (rel err 2.9e-8 vs 1.9e-7); rbc2shock with 2 estimated shock stds,
##     orders 2:3: 20 vs 67 ms (1.6e-8 vs 4.2e-6); nk_demo (9), orders 1:3:
##     28 vs 121 ms (9.9e-11 vs 8.4e-7); ls2003 (17), orders 1:3: 97 vs
##     178 ms (3.2e-9 vs 1.2e-7). Otherwise "implicit": with order 4 the
##     reverse path declines wholesale and runs the implicit path anyway
##     (rbc2shock 1:4: 2127 vs 2188 ms, identical gradient), and without
##     order 3 the forward solves a first-order rule the reverse path does not
##     take (same code path, nk_demo 1:2: 13 vs 12 ms). Since W68 that shared
##     path for orders within 1:2 is EXACT (.cumulant_loglik_grad_order1:
##     solution_derivatives + derivative Lyapunov), not FD of the forward.
##   * whittle, pruned -> "implicit". Whittle runs one analytic path for every
##     non-hybrid method; pruned does not consult grad_method.
##   * gaussian with lik_init = "kappa" -> "hybrid" (W68): the big-kappa P0 has
##     no analytic kernel, and the other methods refuse it. An explicit
##     "diffuse" / "stationary" moves the unit-root test to that init's rule
##     (.grad_init_in_force()).
##   * gaussian where the forward drops observation components at theta_ref
##     (singular_F: the singular-F univariate fallback) -> "hybrid" (W68):
##     the objective is discontinuous there and no analytic kernel follows
##     the dropping; make_posterior_grad() warns and uses FD for every
##     parameter.
## `base` is make_posterior_grad()'s solve at theta_ref (NULL if it failed:
## the stationary default then applies -- per-draw fallbacks stay exact).
.resolve_auto_grad_method <- function(likelihood, base, model, data,
                                      me_extra = NULL, shock_scale = NULL,
                                      cumulant_orders = 1:4,
                                      lik_init = "auto",
                                      singular_F = FALSE) {
  if (identical(likelihood, "cumulant"))
    return(if (3L %in% cumulant_orders && !(4L %in% cumulant_orders))
      "adjoint_solution" else "implicit")
  if (!identical(likelihood, "gaussian")) return("implicit")
  if (identical(lik_init, "kappa") || isTRUE(singular_F)) return("hybrid")
  has_missing <- anyNA(data)
  has_unit_root <- FALSE
  if (!is.null(base)) {
    dr <- base$dr; si <- dr$state_idx
    has_unit_root <- identical(.grad_init_in_force(
      dr$ghx[si, , drop = FALSE], dr$ghu[si, , drop = FALSE],
      .get_shock_cov(model, model$varexo_names, base$params),
      lik_init), "diffuse")
  }
  if (.grad_no_kernel_covers(has_unit_root, has_missing, me_extra, shock_scale))
    return("hybrid")
  "adjoint_solution"
}
