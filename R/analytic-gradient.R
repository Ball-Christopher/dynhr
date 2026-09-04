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
.dlog_prior_density1 <- function(x, dist, p1, p2) {
  switch(.normalize_dist(dist),
    "inv_gamma" =, "inv_gamma1" = {
      if (x <= 0) return(0)
      alpha <- .ig1_alpha(p1, p2)
      theta <- (alpha - 1) * (p2^2 + p1^2)
      # log f = log2 + a*log(theta) - lgamma(a) - (2a+1) log x - theta/x^2
      -(2 * alpha + 1) / x + 2 * theta / x^3
    },
    "inv_gamma2" = {
      if (x <= 0) return(0)
      shape <- (p1 / p2)^2 + 2
      scale <- p1 * (shape - 1)
      if (shape <= 2) return(-1 / x)
      # log f = dgamma(1/x; shape, rate=scale, log) - 2 log x
      #       = (shape-1) log(1/x) - scale/x + const - 2 log x
      #       = -(shape+1) log x - scale/x + const
      -(shape + 1) / x + scale / x^2
    },
    "beta" = {
      if (x <= 0 || x >= 1) return(0)
      v <- p2^2
      a <- p1 * (p1 * (1 - p1) / v - 1)
      b <- (1 - p1) * (p1 * (1 - p1) / v - 1)
      if (a <= 0 || b <= 0) return(0)
      (a - 1) / x - (b - 1) / (1 - x)
    },
    "gamma" = {
      if (x <= 0) return(0)
      shape <- (p1 / p2)^2
      rate  <- p1 / p2^2
      (shape - 1) / x - rate
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
  for (i in seq_len(nrow(prior_spec))) {
    nm <- prior_spec$name[i]
    if (!(nm %in% names(theta))) next
    g[nm] <- .dlog_prior_density1(theta[[nm]], prior_spec$distribution[i],
                                  prior_spec$p1[i], prior_spec$p2[i])
  }
  g
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
    out <- kf_score_sigma_cpp(Yd, TT, RR, ZZ, DD, Sigma_e, HH, SS, me_diag, P,
                              dSigma_list, dH, dS, dPk)
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
#' @param grad_method \code{"hybrid"} (default), \code{"implicit"},
#'   \code{"adjoint"}, or \code{"adjoint_solution"}.
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
#'   that draw; whittle/cumulant/pruned likelihoods treat it as
#'   \code{"implicit"}. Benchmark before preferring it as a default: on
#'   small/medium models the forward layer's shared factorization is already
#'   cheap and the savings may not clear the R-level overhead (see the
#'   E-wave lesson in the NZSIM records).
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
#' @return \code{function(theta)} returning the gradient vector, suitable for the
#'   \code{grad_fn} argument of \code{\link{nuts}}.
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
#' @references Childers, D., Fernandez-Villaverde, J., Perla, J., Rackauckas,
#'   C., & Wu, P. (2022). \emph{Differentiable State-Space Models and
#'   Hamiltonian Monte Carlo Estimation}. NBER Working Paper No. 30573.
#' @seealso \code{\link{make_posterior}}, \code{\link{nuts}},
#'   \code{solution_derivatives}
#' @export
make_posterior_grad <- function(model, data, prior_spec, obs_vars, compiled,
                                me_variance = 0, theta_ref = NULL,
                                verbose = FALSE,
                                grad_method = c("hybrid", "implicit",
                                                "adjoint",
                                                "adjoint_solution"),
                                me_extra = NULL, shock_scale = NULL,
                                likelihood = "gaussian",
                                freq_band = c(0, pi),
                                cumulant_orders = 1:4,
                                cumulant_weight = "identity",
                                debias = TRUE,
                                pruned_order = 2L) {
  grad_method <- match.arg(grad_method)
  .validate_prior_spec_dist(prior_spec)  # F3: fail loud at BUILD time, not mid-chain
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
  if (!likelihood %in% c("gaussian", "whittle", "cumulant", "pruned"))
    stop("make_posterior_grad(): no gradient path for likelihood '", likelihood,
         "'. Supported: gaussian, whittle, cumulant, pruned. ",
         "(pskf/student_t/tpf/ppf/copf have no gradient path.)")
  lp_fn <- if (use_whittle) {
    make_log_posterior_whittle(model, data, prior_spec, obs_vars, compiled,
                               me_variance = me_variance,
                               freq_band = freq_band,
                               debias = debias)
  } else if (use_cumulant) {
    make_log_posterior_cumulant(model, data, prior_spec, obs_vars, compiled,
                                me_variance = me_variance,
                                cumulant_orders = cumulant_orders,
                                cumulant_weight = cumulant_weight)
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
                                 me_variance = me_variance)
    } else {
      make_log_posterior_pruned(model, data, prior_spec, obs_vars, compiled,
                                me_variance = me_variance)
    }
  } else {
    make_log_posterior(model, data, prior_spec, obs_vars, compiled,
                       me_variance = me_variance,
                       me_extra = me_extra, shock_scale = shock_scale)
  }
  ## Time-varying inputs present? .kf_loglik_score_sigma (the analytic sigma
  ## score used by the "hybrid" closure) is NOT tv-aware: its loglik would be
  ## of a different likelihood than lp_fn, and mixing the two inside the FD
  ## step produces garbage gradients. Under tv, hybrid treats sigma params by
  ## consistent FD instead (tangent/adjoint paths are fully tv-aware).
  has_tv <- !is.null(me_extra) ||
            (!is.null(shock_scale) && !all(shock_scale == 1))

  ## Solve the model at a parameter vector, returning the decision rule (or NULL).
  .solve_dr <- function(theta) {
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

  ## --- One-time classification: which params leave the decision rule fixed? ---
  if (is.null(theta_ref)) theta_ref <- setNames(prior_spec$mean, par_names)
  base <- .solve_dr(theta_ref)
  is_sigma <- rep(FALSE, np); names(is_sigma) <- par_names
  if (!is.null(base)) {
    g0 <- list(ghx = base$dr$ghx, ghu = base$dr$ghu)
    for (k in seq_len(np)) {
      th <- theta_ref; h <- 1e-5 * max(abs(th[k]), 1e-3)
      th[k] <- th[k] + h
      d2 <- .solve_dr(th)
      if (!is.null(d2)) {
        dchg <- max(abs(d2$dr$ghx - g0$ghx), abs(d2$dr$ghu - g0$ghu))
        is_sigma[k] <- dchg < 1e-10        # decision rule invariant -> sigma-like
      }
    }
  }
  sig_names <- par_names[is_sigma]
  num_names <- par_names[!is_sigma]
  if (verbose)
    cat(sprintf("  Analytic gradient: %d analytic (Kalman score: %s), %d numerical (%s)\n",
                length(sig_names), paste(sig_names, collapse = ","),
                length(num_names), paste(num_names, collapse = ",")))

  ## State-space matrices from a decision rule (mirrors kalman_filter()).
  obs_in_endo <- function(dr) match(obs_vars, dr$endo_names)
  me_diag <- me_variance * diag(length(obs_vars))

  ## dSigma_e/dtheta_nm by central FD of .get_shock_cov, for sigma-like params
  ## (and reused as the certainty-equivalence dSigma_e block under "implicit").
  .dSigma_e_fd <- function(theta, params, nm) {
    h <- 1e-6 * max(abs(theta[[nm]]), 1e-3)
    tp <- theta; tm <- theta; tp[nm] <- tp[nm] + h; tm[nm] <- tm[nm] - h
    pp <- params; pm <- params
    pp[nm] <- tp[[nm]]; pm[nm] <- tm[[nm]]
    (.get_shock_cov(model, exo, pp) - .get_shock_cov(model, exo, pm)) / (2 * h)
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
      names(theta) <- par_names
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
          error = function(e) NULL)
        if (is.null(dr_use)) dr_use <- d1   # order-2 failed: orders 3-4 drop out
      } else {
        dr_use <- d1                        # first-order: matches the forward
      }
      ## grad_method == "adjoint_solution": reverse-mode (O(1) in P) path for
      ## orders 1-3 via .cumulant_loglik_grad_adjoint; any NA entry (e.g. order
      ## 4 requested, or a not-ok adjoint block) FD-fallbacks below, so the
      ## returned gradient stays exactly consistent with lp_fn. All other
      ## grad_method values keep the (unchanged) implicit path.
      gll <- if (grad_method == "adjoint_solution") {
        g_adj <- tryCatch(
          .cumulant_loglik_grad_adjoint(model, compiled, dr_use, params,
                                        par_names, obs_vars, data,
                                        orders = cumulant_orders,
                                        me_variance = me_variance),
          error = function(e) setNames(rep(NA_real_, np), par_names))
        ## If the reverse path declined wholesale (all NA — e.g. order 4 or a
        ## non-DecisionRules2), fall back to the implicit path rather than a
        ## full per-param FD, matching the accuracy of the other methods.
        if (all(is.na(g_adj))) {
          tryCatch(
            cumulant_loglik_grad(model, compiled, dr_use, params, par_names,
                                 obs_vars, data, orders = cumulant_orders,
                                 me_variance = me_variance, deriv = "implicit"),
            error = function(e) setNames(rep(NA_real_, np), par_names))
        } else g_adj
      } else {
        tryCatch(
          cumulant_loglik_grad(model, compiled, dr_use, params, par_names,
                               obs_vars, data, orders = cumulant_orders,
                               me_variance = me_variance, deriv = "implicit"),
          error = function(e) setNames(rep(NA_real_, np), par_names))
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
    return(cumulant_grad_fn)
  }

  ## --- "pruned" gradient closure ---------------------------------------------
  ## Order 2: analytic adjoint-chain gradient (phase a: R/pruned-kf-adjoint.R,
  ## phase b: R/pruned-grad-chain.R). One order-2 solve +
  ## solution_derivatives_order2() call covers ALL structural parameters via a
  ## shared factorization (Childers et al. efficiency point, same as the
  ## Gaussian "implicit" path); each parameter also needs a dSigma_e (central
  ## FD of .get_shock_cov, same convention as the Gaussian hybrid path's
  ## .dSigma_e_fd). Any parameter whose solution_derivatives_order2() block
  ## comes back not-ok, or whose chained gradient is non-finite, falls back to
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
      names(theta) <- par_names
      g <- .dlog_prior(theta, prior_spec)        # analytic prior score (all params)

      if (is_pruned3) {
        Y3 <- if (is.null(dim(data))) matrix(data, nrow = length(obs_vars)) else data
        if (nrow(Y3) != length(obs_vars)) Y3 <- t(Y3)

        chain_res3 <- tryCatch(
          .pgo3_grad_chain(theta, model, compiled, par_names, Y3, obs_vars,
                           me_variance = me_variance),
          error = function(e) NULL)

        base_ll <- if (!is.null(chain_res3)) chain_res3$loglik else
          tryCatch(lp_fn(theta)$loglik, error = function(e) -Inf)
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
        error = function(e) NULL)
      if (is.null(ss_result) || !isTRUE(ss_result$converged)) return(g)
      params <- ss_result$params %||% params

      sys <- extract_system_matrices_fast(sys_cache, ss_result$ss, params)
      dr1 <- tryCatch(.solve_from_system(sys, model, compiled, ss_result$ss,
                                         params, FALSE), error = function(e) NULL)
      if (is.null(dr1) || !isTRUE(dr1$bk_satisfied)) return(g)

      dr2 <- tryCatch(
        solve_perturbation_order2(model, compiled, ss_result$ss, params, dr1,
                                  verbose = FALSE),
        error = function(e) NULL)
      if (is.null(dr2)) return(g)

      pss <- tryCatch(pruned_state_space(dr2, model, params), error = function(e) NULL)
      if (is.null(pss)) return(g)

      Y <- if (is.null(dim(data))) matrix(data, nrow = length(obs_vars)) else data
      if (nrow(Y) != length(obs_vars)) Y <- t(Y)

      sd2 <- tryCatch(
        solution_derivatives_order2(model, compiled, dr2, params, par_names),
        error = function(e) NULL)
      if (is.null(sd2)) {
        base_ll <- tryCatch(pruned_ss_loglik(pss, Y, obs_vars, me_variance = me_variance),
                            error = function(e) -Inf)
        if (!is.finite(base_ll)) return(g)
        for (nm in par_names) {
          d1 <- .fd_loglik_grad1(theta, nm, base_ll)
          if (is.finite(d1)) g[nm] <- g[nm] + d1
        }
        return(g)
      }

      dSigma_e_list <- setNames(vector("list", length(par_names)), par_names)
      for (nm in par_names) dSigma_e_list[[nm]] <- .dSigma_e_fd(theta, params, nm)

      chain_res <- tryCatch(
        .pruned_ss_loglik_grad_chain(pss, Y, obs_vars, sd2, dSigma_e_list,
                                     me_variance = me_variance,
                                     model = model, compiled = compiled,
                                     ss = ss_result$ss, dr1 = dr1, params = params),
        error = function(e) NULL)

      base_ll <- if (!is.null(chain_res)) chain_res$loglik else
        tryCatch(pruned_ss_loglik(pss, Y, obs_vars, me_variance = me_variance),
                 error = function(e) -Inf)
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
    return(pruned_grad_fn)
  }

  ## --- "hybrid" gradient closure (existing behaviour, bit-identical) -------
  hybrid_grad_fn <- function(theta) {
    names(theta) <- par_names
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

    ## Analytic Kalman score for the sigma-like params. The recursion also
    ## returns the loglik at theta, reused below as the central value for the
    ## numerical structural-parameter differences (saves one solve+filter).
    base_ll <- NULL
    if (length(sig_names) > 0 && !has_tv) {
      dS_list <- lapply(sig_names, function(nm) .dSigma_e_fd(theta, params, nm))
      names(dS_list) <- sig_names
      sc <- .kf_loglik_score_sigma(Y, TT, RR, ZZ, DD, Sigma_e, dvec, me_diag, dS_list)
      finite_sc <- is.finite(sc$score)
      g[sig_names[finite_sc]] <- g[sig_names[finite_sc]] + sc$score[finite_sc]
      if (is.finite(sc$loglik)) base_ll <- sc$loglik
      ## Sigma params whose analytic score came back non-finite (e.g. a
      ## (near-)unit-root TT makes the stationary Lyapunov P0 undefined):
      ## fall back to relative-step FD for just those parameters.
      if (!all(finite_sc)) {
        if (is.null(base_ll)) base_ll <- lp_fn(theta)$loglik
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
    fd_names <- if (has_tv) c(sig_names, num_names) else num_names
    if (length(fd_names) > 0) {
      if (is.null(base_ll)) base_ll <- lp_fn(theta)$loglik
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
      names(theta) <- par_names
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
        d_ss_list[[nm]] <- list(dSigma_e = .dSigma_e_fd(theta, params, nm))

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
          error = function(e) NULL
        )
        for (nm in num_names) {
          d <- if (!is.null(sd_res)) sd_res$derivs[[nm]] else NULL
          if (!is.null(d) && isTRUE(d$ok)) {
            d_ss_list[[nm]] <- list(
              dTT = d$dTT, dRR = d$dRR, dZZ = d$dZZ, dDD = d$dDD,
              dSigma_e = .dSigma_e_fd(theta, params, nm)
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
              dEI_per_param[[k_nm]] <- NULL  ## FD fallback
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
        }, error = function(e) FALSE)
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
    return(whittle_grad_fn)
  }

  if (grad_method == "hybrid") return(hybrid_grad_fn)

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
    names(theta) <- par_names
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
    has_unit_root <- FALSE
    if (nrow(TT) > 0 &&
        any(Mod(eigen(TT, symmetric = FALSE, only.values = TRUE)$values) > 1 - 1e-6)) {
      QQ_g <- RR %*% Sigma_e %*% t(RR)
      P0_g <- tryCatch(solve_lyapunov(TT, QQ_g), error = function(e) NULL)
      stat_ok <- !is.null(P0_g) && all(is.finite(P0_g)) &&
        min(Re(eigen((P0_g + t(P0_g)) / 2, symmetric = TRUE,
                     only.values = TRUE)$values)) > -1e-8
      has_unit_root <- !stat_ok
    }
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
    if ((has_unit_root || has_missing) &&
        !use_uni_adjoint && !use_diffuse_adjoint) {
      kernel_stats$n_fallback_special_init <- kernel_stats$n_fallback_special_init + 1L
      if (!warned_special) {
        warning("make_posterior_grad: ", grad_method, " gradient with a ",
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
      d_ss_list[[nm]] <- list(dSigma_e = .dSigma_e_fd(theta, params, nm))

    sd_res <- NULL
    if (use_sol_adjoint) {
      for (nm in num_names)
        d_ss_list[[nm]] <- list(dSigma_e = .dSigma_e_fd(theta, params, nm))
    } else if (length(num_names) > 0) {
      sd_res <- tryCatch(
        solution_derivatives(model, compiled, dr, params,
                              param_names = num_names, obs_vars = obs_vars),
        error = function(e) NULL
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
                                   dSigma_e = .dSigma_e_fd(theta, params, nm))
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
        error = function(e) NULL)
      if (is.null(tang) || !isTRUE(tang$stage1_ok) || !isTRUE(tang$stage2_ok) ||
          !is.finite(tang$loglik)) {
        kernel_stats$n_fallback_diffuse <- kernel_stats$n_fallback_diffuse + 1L
        if (!warned_special) {
          warning("make_posterior_grad: exact-diffuse adjoint gradient ",
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
      ## instead of crashing the sampler.
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
      }, error = function(e) NULL)

      if (is.null(tang)) {
        kernel_stats$n_fallback_nondiffuse <- kernel_stats$n_fallback_nondiffuse + 1L
        if (!warned_special) {
          warning("make_posterior_grad: ", grad_method, " gradient kernel ",
                  "failed for this draw (numerical KF failure); using the ",
                  "FD-hybrid fallback.", call. = FALSE)
          warned_special <<- TRUE
        }
        return(hybrid_grad_fn(theta))
      }
    }
    base_ll <- tang$loglik

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
        error = function(e) NULL
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
      if (!is.finite(base_ll)) base_ll <- lp_fn(theta)$loglik
      for (nm in fd_names) {
        d1 <- .fd_loglik_grad1(theta, nm, base_ll)
        if (is.finite(d1)) g[nm] <- g[nm] + d1
      }
    }

    g
  }

  attr(.impl_adj_grad_fn, "kernel_stats") <- kernel_stats
  .impl_adj_grad_fn
}
