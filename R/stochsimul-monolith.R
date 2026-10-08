################################################################################
## dynhr_stochsimul.R v0.2 (refactored 2026-05-08)
## Phase 4 -- First-order perturbation solution (stoch_simul)
## Depends on: dynhr_parser.R, dynhr_jacobian.R, dynhr_steady.R
##
## Changelog v0.2 (refactor):
##   - Removed extract_jacobian_partitions() -- superseded by
##     extract_system_matrices() in dynhr_perturbation.R
##   - Removed .compute_decision_rules() -- superseded by Villemot solver
##     in dynhr_perturbation.R
##   - Removed conditional source() guards -- use dynhr_load_all.R
##   - Fixed compute_irfs(): params argument was silently overwritten
##   - Fixed compute_moments(): params argument was silently overwritten
################################################################################

## Helper: extract shock standard deviations from model$shocks

## Package-private memo for parse()d expressions. Shock stderr/variance specs
## are constant strings across MCMC draws; only the parameter values they are
## evaluated against change. Caching the parse tree (keyed on the literal
## string) removes a per-draw parse() that showed up hot in profiling. eval()
## still runs every call against the current parameter environment.
##
## A-SEC: this cache is a security choke point (see the note above
## `.ssm_expr_cache` in steady-monolith.R).  An expression is AST-checked
## against the .mod allowlist when it is first inserted, so a model built
## programmatically -- which never went through parse_shocks_block()'s sandbox
## -- still cannot smuggle `file.create()` into a variance expression.  The
## check aborts with `dynhr_error_unsafe_mod_expression`; callers that map
## ordinary evaluation errors to NA must re-raise that class.  `envir` should
## be a child of `.dynhr_mod_sandbox_parent()` (see `.dynhr_param_eval_env`).
.dynhr_expr_cache <- new.env(parent = emptyenv())

.eval_cached_expr <- function(expr_str, envir) {
  ex <- .dynhr_expr_cache[[expr_str]]
  if (is.null(ex)) {
    ex <- parse(text = expr_str)
    .dynhr_check_mod_expr(ex, expr_str, "the shocks-block expression")
    .dynhr_expr_cache[[expr_str]] <- ex
  }
  eval(ex, envir = envir)
}

## Parameter environment for re-evaluating shocks-block *_expr text: the
## params bound in a child of the session allowlist sandbox (never baseenv()).
.dynhr_param_eval_env <- function(params) {
  list2env(as.list(params), parent = .dynhr_mod_sandbox_parent())
}

.get_shock_stderr <- function(model, exo_names, params = NULL) {
  stderr <- setNames(rep(0, length(exo_names)), exo_names)
  if (is.null(params)) params <- model$param_values

  ## Eval environment from parameters, built lazily: models whose shocks are
  ## all resolved by Priority 0 (estimated shock-named params) or have no
  ## *_expr never need it, so we skip the per-call list2env for them.
  penv <- NULL
  get_penv <- function() {
    if (is.null(penv)) penv <<- .dynhr_param_eval_env(params)
    penv
  }

  sv <- if (!is.null(model$shocks) && is.data.frame(model$shocks$variances))
    model$shocks$variances else NULL

  for (i in seq_along(exo_names)) {
    nm <- exo_names[i]
    se <- NA_real_

    ## Priority 0: a parameter named EXACTLY the shock (Dynare `stderr <shock>`
    ## estimated_params convention). When a shock std is estimated this way the
    ## draw is injected into `params` under the shock name (see
    ## .apply_theta_to_params); it must override the frozen parse-time stderr in
    ## the shocks block (Priority 1 would otherwise short-circuit on the
    ## calibration constant, e.g. 0.01, giving the estimated std ZERO effect on
    ## the likelihood -- the P0 bug). Inert for every model that does not inject
    ## a shock-named param (shock names are never structural params).
    if (!is.null(params) && nm %in% names(params) && is.finite(params[[nm]])) {
      se <- params[[nm]]
    }

    ## [BODY PRESERVED FROM ORIGINAL -- shock stderr lookup logic]
    ## Search in shocks$variances for matching shock name, evaluate
    ## stderr expression in parameter environment, handle symbolic
    ## references (e.g. sig_a), fall back to param_values lookup.
    ## ----------------------------------------------------------
    ## Priority 1: shocks block variance/stderr. The *_expr columns hold the
    ## raw expression text from the .mod file; re-evaluating them against the
    ## CURRENT params (not the parse-time snapshot in the numeric columns) is
    ## what lets an estimated shock std like `stderr sig_a;` track theta
    ## during MCMC/SBC. Numeric columns remain as calibration fallbacks.
    ## Skipped entirely when Priority 0 already resolved an estimated shock std.
    if ((is.na(se) || !is.finite(se)) && !is.null(sv)) {
      idx <- which(sv$name == nm)
      if (length(idx) > 0) {
        row <- sv[idx[1], ]
        safe_eval <- function(expr_str) {
          tryCatch(.eval_cached_expr(expr_str, get_penv()),
                   error = function(e) .dynhr_reraise_unsafe(e, NA_real_))
        }
        ## Priority: raw expression columns first (re-evaluate against current
        ## params), then fall back to parse-time numeric snapshots.  This
        ## ordering matters because the var/variance syntax stores the
        ## expression only in variance_expr (stderr_expr is NA), while the
        ## var+stderr two-statement syntax stores it only in stderr_expr.  If
        ## we checked the numeric stderr column before variance_expr, a
        ## `var eps_a = sig_a^2;` shock would return the frozen calibration
        ## value (0.007) instead of sqrt(sig_a^2) at the current theta.
        ##
        ## Correct priority: stderr_expr → variance_expr → numeric stderr → numeric variance.

        ## 1. Re-evaluate the stderr expression against current params
        if ("stderr_expr" %in% names(row) && !is.na(row$stderr_expr)) {
          se <- safe_eval(as.character(row$stderr_expr))
        }
        ## 2. Re-evaluate the variance expression against current params
        if ((is.na(se) || !is.finite(se)) && "variance_expr" %in% names(row) && !is.na(row$variance_expr)) {
          v <- safe_eval(as.character(row$variance_expr))
          if (!is.na(v) && is.finite(v) && v >= 0) se <- sqrt(v)
        }
        ## 3. Fall back to the parse-time stderr snapshot
        if ((is.na(se) || !is.finite(se)) && "stderr" %in% names(row) && !is.na(row$stderr)) {
          se <- safe_eval(as.character(row$stderr))
        }
        ## 4. Fall back to the parse-time variance snapshot
        if ((is.na(se) || !is.finite(se)) && "variance" %in% names(row) && !is.na(row$variance)) {
          v <- safe_eval(as.character(row$variance))
          if (!is.na(v) && is.finite(v) && v >= 0) se <- sqrt(v)
        }
      }
    }

    ## There is deliberately NO name-matching fallback here.  Until
    ## 0.9.4 a shock missing from the shocks block took the value of any
    ## parameter named sig_/stderr_/sigma_<core> -- so `sigma_c = 2` (risk
    ## aversion) silently gave eps_c a variance of 4.  Dynare semantics: a
    ## shock the shocks block does not declare has ZERO variance.

    ## Default: 0 (shock has no variance)
    if (is.na(se) || !is.finite(se)) se <- 0
    stderr[nm] <- se
  }
  stderr
}

## Helper: the shock covariance an IRF routine should scale its impulses by.
##
## SHARED by compute_irfs() (order 1) and compute_irfs_order2() (order 2) so the
## two cannot drift apart again: before 0.9.4 the order-1 path preferred
## `dr$Sigma_e` (with correlations) while the order-2 path always rebuilt a
## DIAGONAL covariance from `params`, so the same model gave different IRFs at
## the two orders.
##
## Precedence follows the 0.9.3.7 Kalman-path rule -- `params` WINS over
## `dr$Sigma_e`:
##   1. caller-supplied `params`  -> .get_shock_cov(model, exo, params)
##   2. else `dr$Sigma_e`, when it carries the right dimension (set by
##      solve_perturbation(Sigma_e = ))
##   3. else .get_shock_cov(model, exo, model$param_values)
## `params = NULL` therefore keeps the old dr$Sigma_e behaviour.
## @noRd
.irf_shock_scale <- function(dr, model, params = NULL) {
  exo   <- dr$exo_names
  n_exo <- length(exo)
  if (!is.null(params)) return(.get_shock_cov(model, exo, params))
  if (!is.null(dr$Sigma_e) && all(dim(dr$Sigma_e) == c(n_exo, n_exo)))
    return(dr$Sigma_e)
  .get_shock_cov(model, exo, model$param_values)
}

## Helper: lower-triangular Cholesky factor L of a shock covariance,
## Sigma_e = L %*% t(L), in the DECLARED shock order (Dynare's orthogonalisation
## convention for IRFs, variance decompositions and simulation draws).
##
## - Diagonal Sigma_e  -> diag(sqrt(diag(Sigma_e))), i.e. byte-identical to the
##   old per-shock `sigma_k` scaling, so nothing changes for uncorrelated shocks.
## - Positive-definite   -> t(chol(.)).
## - PSD but singular (a shock with zero declared variance, a perfectly
##   correlated pair) -> an explicit outer-product Cholesky that emits a zero
##   column where the pivot vanishes.  chol() would error there; this keeps the
##   factorisation exact (L L' == Sigma_e) instead of silently dropping the
##   off-diagonals.  No tryCatch: the branch is decided by the pivot/eigenvalue.
## @noRd
.sigma_e_chol_lower <- function(Sigma_e) {
  n <- nrow(Sigma_e)
  if (is.null(n) || n == 0L) return(matrix(0, 0, 0))
  S <- (Sigma_e + t(Sigma_e)) / 2
  off <- S
  diag(off) <- 0
  if (all(off == 0)) return(diag(sqrt(pmax(diag(S), 0)), nrow = n))

  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) > 1e-12 * max(1, max(ev))) return(t(chol(S)))

  ## Outer-product Cholesky with zero pivots tolerated (PSD / singular case).
  L <- matrix(0, n, n)
  for (j in seq_len(n)) {
    prev <- seq_len(j - 1L)
    d <- S[j, j] - sum(L[j, prev]^2)
    if (d <= 0) next                       # zero pivot -> zero column
    L[j, j] <- sqrt(d)
    if (j < n) {
      rows <- seq.int(j + 1L, n)
      L[rows, j] <- (S[rows, j] -
                     as.numeric(L[rows, prev, drop = FALSE] %*% L[j, prev])) /
                    L[j, j]
    }
  }
  L
}

## Helper: extract shock skewness shape parameters (alpha) from model$shocks
##
## The skew column in model$shocks$variances holds a parse-time numeric
## snapshot; skew_expr holds the raw expression text so the alpha can be
## re-evaluated against the current parameter vector theta during MCMC (mirrors
## the stderr_expr / .get_shock_stderr pattern above).  alpha = 0 (default)
## gives a symmetric Gaussian shock; alpha != 0 gives a skew-normal with
## E[eps_i] = sigma_i * delta_i * sqrt(2/pi), delta_i = alpha_i/sqrt(1+alpha_i^2).
## alpha can be NEGATIVE.
#' Draw zero-mean shock vectors from the joint closed skew-normal
#'
#' The law is the one the package's own skew likelihood evaluates
#' (\code{.get_csn_shock_params()} / \code{.csn_state_noise_lift()} in
#' \code{R/pskf-likelihood.R}):
#' \deqn{e \sim CSN(\mu_e, \Sigma_e, \Gamma_e, 0, I_q)}
#' with \eqn{\Gamma_e} the rows \eqn{\alpha_i/\sigma_i} of the skewed shocks and
#' the FULL \eqn{\Sigma_e} (off-diagonals included) as the seed covariance. All
#' cross-shock coupling is carried by \eqn{\Sigma_e}: the truncation latents
#' inherit its correlation, which is what makes the cross-shock co-skewness
#' nonzero and sign-matched to \eqn{\rho}.
#'
#' Sampling uses the stochastic (conditioning) representation of
#' Gonzalez-Farias, Dominguez-Molina & Gupta (2004) / Arellano-Valle & Azzalini
#' (2006): draw \eqn{[e; U]} jointly Gaussian with
#' \deqn{Cov = [[\Sigma_e, \Sigma_e \Gamma_e'], [\Gamma_e \Sigma_e,
#'              I + \Gamma_e \Sigma_e \Gamma_e']]}
#' and keep the draws with \eqn{U \ge 0} componentwise. Acceptance is
#' \eqn{\Phi_q(0; 0, I + \Gamma_e \Sigma_e \Gamma_e')}, at least \eqn{2^{-q}},
#' so the vectorised over-draw below needs no per-period loop. Only shocks with
#' \eqn{\alpha_i \ne 0} contribute a latent, keeping \eqn{q} as small as
#' possible.
#'
#' The returned draws are shifted by the SAME mean correction the filter
#' applies, the exact mean of the joint CSN law (\code{.csn_shock_mean()} in
#' \code{R/pskf-likelihood.R}), so the draws are mean-zero and the simulated
#' data sit where the filter places the steady state. With no skewed shock
#' correlated to another shock this is the per-shock
#' \eqn{E[e_i] = \sigma_i \delta_i \sqrt{2/\pi}},
#' \eqn{\delta_i = \alpha_i/\sqrt{1+\alpha_i^2}}. (Before 0.9.4 that
#' per-shock shift was applied at \eqn{\rho \ne 0} too, in both the filter
#' and here, so the draws were NOT mean-zero and the circular DGP/likelihood
#' agreement hid it.)
#' Likewise the realised sample CORRELATION of the draws is NOT \eqn{\rho}:
#' \eqn{\rho} parameterises the SEED covariance of the pre-truncation Gaussian,
#' and the conditioning step changes it.
#'
#' @param n_draws  Number of shock vectors (rows) to return.
#' @param Sigma_e  n_exo x n_exo seed shock covariance (from .get_shock_cov).
#' @param alpha    length-n_exo skewness shape vector (from .get_shock_skewness).
#' @param sigma_e  length-n_exo shock standard deviations (sqrt(diag(Sigma_e))).
#' @return n_draws x n_exo numeric matrix.
#' @noRd
.draw_csn_shocks <- function(n_draws, Sigma_e, alpha, sigma_e) {
  n_exo <- length(alpha)
  alpha <- as.numeric(alpha); sigma_e <- as.numeric(sigma_e)

  ## Only genuinely skewed shocks get a truncation latent. A zero-stderr shock
  ## cannot be skewed (alpha/sigma is undefined), so it is excluded too.
  sk <- which(alpha != 0 & is.finite(alpha) & sigma_e > 0)
  q  <- length(sk)
  if (q == 0L)
    return(matrix(rnorm(n_draws * n_exo), n_draws, n_exo) %*%
             t(.sigma_e_chol_lower(Sigma_e)))

  Gam <- matrix(0, q, n_exo)
  Gam[cbind(seq_len(q), sk)] <- alpha[sk] / sigma_e[sk]

  SG  <- Sigma_e %*% t(Gam)                       # n_exo x q
  M_U <- diag(q) + Gam %*% Sigma_e %*% t(Gam)     # q x q
  J   <- rbind(cbind(Sigma_e, SG), cbind(t(SG), M_U))
  J   <- 0.5 * (J + t(J))
  Lj  <- chol(J + diag(1e-12, n_exo + q))         # upper factor: z %*% Lj

  acc <- matrix(NA_real_, n_draws, n_exo)
  got <- 0L
  ## Acceptance >= 2^-q, so over-drawing by 2^q * 1.25 clears the remainder in
  ## one or two passes for any practical n_exo.
  repeat {
    need <- n_draws - got
    if (need <= 0L) break
    m  <- max(128L, as.integer(ceiling(need * (2^q) * 1.25)))
    Zc <- matrix(rnorm(m * (n_exo + q)), m, n_exo + q) %*% Lj
    U  <- Zc[, n_exo + seq_len(q), drop = FALSE]
    keep <- .rowSums(U >= 0, m, q) == q
    if (any(keep)) {
      ok   <- Zc[keep, seq_len(n_exo), drop = FALSE]
      take <- min(nrow(ok), need)
      acc[got + seq_len(take), ] <- ok[seq_len(take), , drop = FALSE]
      got <- got + take
    }
  }

  ## Filter-matching mean correction (see roxygen above): the exact joint-CSN
  ## mean. The seed covariance's own diagonal is used for sigma inside, which
  ## is what sigma_e is at every call site.
  sweep(acc, 2L, .csn_shock_mean(Sigma_e, alpha), `-`)
}


.get_shock_skewness <- function(model, exo_names, params = NULL) {
  alpha <- setNames(rep(0, length(exo_names)), exo_names)
  if (is.null(params)) params <- model$param_values

  penv <- NULL
  get_penv <- function() {
    if (is.null(penv)) penv <<- .dynhr_param_eval_env(params)
    penv
  }

  sv <- if (!is.null(model$shocks) && is.data.frame(model$shocks$variances))
    model$shocks$variances else NULL

  ## Priority 0 (mirrors .get_shock_stderr's Priority 0 and .get_shock_cov's
  ## correlation Priority 0): an ESTIMATED `skew <shock>` arrives in `params`
  ## under the canonical "skew <shock>" key (injected by
  ## apply_theta_to_params) and must override the shocks block -- INCLUDING a
  ## `skew_expr` row, which would otherwise freeze alpha at its calibrated
  ## value and leave the PSKF likelihood flat in the estimated skewness.
  ## Applied first, and the block lookup below skips any shock it resolved.
  skidx <- .skew_params_index(params, exo_names)
  resolved <- logical(length(exo_names))
  for (k in seq_len(nrow(skidx))) {
    v <- params[[skidx$k[k]]]
    if (!is.finite(v)) next
    alpha[skidx$i[k]] <- v
    resolved[skidx$i[k]] <- TRUE
  }

  if (is.null(sv)) return(alpha)
  if (!("skew" %in% names(sv))) return(alpha)

  for (i in seq_along(exo_names)) {
    if (resolved[i]) next
    nm  <- exo_names[i]
    idx <- which(sv$name == nm)
    if (length(idx) == 0L) next
    row <- sv[idx[1], ]

    al <- NA_real_
    if ("skew_expr" %in% names(row) && !is.na(row$skew_expr) &&
        nzchar(row$skew_expr)) {
      al <- tryCatch(.eval_cached_expr(as.character(row$skew_expr), get_penv()),
                     error = function(e) .dynhr_reraise_unsafe(e, NA_real_))
    }
    if (is.na(al) || !is.finite(al)) {
      if (!is.na(row$skew) && is.finite(row$skew)) al <- row$skew
    }
    if (!is.na(al) && is.finite(al)) alpha[nm] <- al
  }
  alpha
}

## ============================================================================
## 1. IMPULSE RESPONSE FUNCTIONS
## ============================================================================

#' Compute impulse response functions for each shock
#'
#' @param dr DecisionRules object from solve_perturbation
#' @param model dynhr_mod object (for shock variances)
#' @param n_periods Number of IRF periods (default 40)
#' @param shock_size Size of shock in std dev units (default 1)
#' @param params Named numeric parameter vector (default: model$param_values)
#' @return An \code{IRFCollection} object: a named list with one entry per
#'   shock.  Each entry is a \code{T x n_endo} numeric matrix where rows are
#'   horizons 1…T and columns are endogenous variables in model order.
#'   Access individual shocks and variables as
#'   \code{irfs[["shock_name"]][t, "var_name"]}.
#'   This is \strong{not} a 3-D array; it is a plain R list of matrices.
#'
#' @references
#'   Koop, G., Pesaran, M. H., & Potter, S. M. (1996). Impulse response
#'     analysis in nonlinear multivariate models. \emph{Journal of Econometrics},
#'     74(1), 119-147.
#'   Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez, J. F. (2018).
#'     The pruned state-space system for non-linear DSGE models.
#'     \emph{Review of Economic Studies}, 85(1), 1-49.
#' @export
compute_irfs <- function(dr, model, n_periods = 40L, shock_size = 1, params = NULL) {
  ghx <- dr$ghx
  ghu <- dr$ghu
  endo <- dr$endo_names
  exo  <- dr$exo_names
  state_idx <- dr$state_idx
  n_endo <- length(endo)
  n_exo  <- length(exo)
  n_state <- length(state_idx)

  ## Build full shock covariance and its lower Cholesky factor.
  ## Dynare convention: IRF for shock k uses the k-th column of chol(Sigma_e)
  ## (lower triangular), so correlated shocks propagate to all impact responses.
  ## For a diagonal Sigma_e this reduces to ghu[,k]*sigma_k (backward-compat).
  ## The precedence between caller `params` and a covariance
  ## carried on `dr` now lives in the shared .irf_shock_scale(), which
  ## compute_irfs_order2() uses too -- `params` wins when supplied.
  Sigma_e_irf <- .irf_shock_scale(dr, model, params)
  if (is.null(params)) params <- model$param_values
  L_chol <- .sigma_e_chol_lower(Sigma_e_irf)

  irfs <- list()
  for (k in seq_along(exo)) {
    shock_name <- exo[k]
    irf_mat <- matrix(0, nrow = n_periods, ncol = n_endo)
    colnames(irf_mat) <- endo
    rownames(irf_mat) <- paste0("t", seq_len(n_periods))

    ## Shock vector: k-th column of lower Cholesky of Sigma_e, scaled by
    ## shock_size.  Reduces to shock_stderr[k]*shock_size when Sigma_e diagonal.
    eps <- L_chol[, k] * shock_size

    ## A shock with no declared variance (all-zero Cholesky column) would
    ## otherwise give an identically-zero IRF even though the policy function
    ## ghu[,k] is correct. Drive a UNIT impulse of magnitude `shock_size`
    ## instead, so e.g. monetary-policy IRFs work when the shock variance is set
    ## to 0 in the .mod (common in IRF-only / IRF-matching exercises).
    if (all(eps == 0) && shock_size != 0) {
      eps <- numeric(n_exo)
      eps[k] <- shock_size
    }

    ## Period 1: impact
    y <- as.numeric(ghu %*% eps)
    irf_mat[1, ] <- y

    ## Periods 2..n_periods: propagation through state transition.
    ## seq_len(n_periods - 1L) + 1L avoids the 2:1 = c(2,1) pitfall when
    ## n_periods == 1.
    if (n_periods > 1L) {
      for (t in seq.int(2L, n_periods)) {
        y_state <- y[state_idx]
        y <- as.numeric(ghx %*% y_state)
        irf_mat[t, ] <- y
      }
    }
    irfs[[shock_name]] <- irf_mat
  }

  class(irfs) <- "IRFCollection"
  attr(irfs, "n_periods")  <- n_periods
  attr(irfs, "endo_names") <- endo
  attr(irfs, "exo_names")  <- exo
  irfs
}

## ============================================================================
## 2. SIMULATION
## ============================================================================

#' Simulate the model forward from steady state
#'
#' @param dr DecisionRules object
#' @param n_periods Number of simulation periods
#' @param shocks Matrix of shocks with \code{n_periods + burn_in} rows and
#'   \code{n_exo} columns (the burn-in rows are simulated and then discarded;
#'   pass \code{burn_in = 0} to supply exactly \code{n_periods} rows). If the
#'   matrix has column names they must be the model's exogenous names and the
#'   columns are matched by name. If NULL, shocks are
#'   drawn from the model's full shock covariance \eqn{\Sigma_e} (from
#'   \code{\link{shock_cov}}), \strong{including} the \code{corr} /
#'   \code{var a, b} entries of the shocks block, via its lower Cholesky
#'   factor. Skewed shocks (\code{alpha != 0}) are still drawn independently
#'   from their skew-normal marginals and warn when correlations are declared.
#' @param n_replications Number of replications for stochastic simulation.
#'   When \code{1} (default) a plain \code{n_periods x n_endo} matrix is
#'   returned (backward-compatible).  When \code{> 1} a 3-D array of dimension
#'   \code{(n_periods, n_endo, n_replications)} is returned; each slice
#'   \code{[,,r]} is one independent replication drawn with a different random
#'   seed.  The pre-existing \code{attr(., "levels")} attribute is attached to
#'   the first slice only (2-D case) or omitted (3-D case).
#' @param model dynhr_mod (for shock variances)
#' @param burn_in Burn-in periods to discard
#' @param init_state Optional named numeric vector of initial state deviations
#'   (over endogenous names) loaded into the period-1 state; pair with
#'   \code{burn_in = 0}. \code{NULL} starts at the steady state.
#' @param linear Logical (default \code{FALSE}). For a second- or third-order
#'   rule (\code{DecisionRules2} / \code{DecisionRules3}) the simulation uses
#'   the rule's native order (pruned, via \code{simulate_model_order2} /
#'   \code{simulate_model_order3}); \code{linear = TRUE} forces the
#'   first-order (\code{ghx}/\code{ghu}) simulation. No effect on a
#'   first-order rule.
#' @return Matrix (n_periods x n_endo) when \code{n_replications = 1};
#'   3-D array (n_periods x n_endo x n_replications) when \code{> 1}.
#'
#'   The values are \strong{deviations from the steady state}; the same paths
#'   in LEVELS are attached as \code{attr(., "levels")}. Every filtering and
#'   smoothing entry point -- \code{\link{kalman_filter}},
#'   \code{\link{kalman_smoother}}, \code{\link{forecast_backtest}},
#'   \code{\link{realtime_decomposition}} -- takes LEVELS and subtracts the
#'   model's steady state itself, so feed them
#'   \code{attr(sim, "levels")[, obs_vars]}, not the raw return value, unless
#'   the observables' steady states are zero.
#' @export
simulate_model <- function(dr, n_periods = 200L, shocks = NULL,
                           n_replications = 1L, model = NULL,
                           burn_in = 100L, init_state = NULL, linear = FALSE) {
  if (length(n_periods) != 1L || !is.numeric(n_periods) || !is.finite(n_periods) ||
      n_periods < 1 || n_periods != round(n_periods))
    stop("simulate_model: n_periods must be a positive whole number.", call. = FALSE)
  if (length(burn_in) != 1L || !is.numeric(burn_in) || !is.finite(burn_in) ||
      burn_in < 0 || burn_in != round(burn_in))
    stop("simulate_model: burn_in must be a non-negative whole number.", call. = FALSE)
  n_periods <- as.integer(n_periods); burn_in <- as.integer(burn_in)
  ## n_replications > 1: run the single-path body in a loop and stack results
  ## into a 3-D array.  Each replication draws independently (shocks = NULL).
  if (!is.null(n_replications) && n_replications > 1L) {
    n_replications <- as.integer(n_replications)
    endo <- dr$endo_names
    n_endo <- length(endo)
    out <- array(NA_real_, dim = c(n_periods, n_endo, n_replications),
                 dimnames = list(
                   paste0("t", seq_len(n_periods)),
                   endo,
                   paste0("rep", seq_len(n_replications))
                 ))
    for (r in seq_len(n_replications)) {
      out[, , r] <- simulate_model(dr, n_periods = n_periods,
                                   shocks = shocks,
                                   n_replications = 1L,
                                   model = model,
                                   burn_in = burn_in,
                                   init_state = init_state,
                                   linear = linear)
    }
    return(out)
  }

  ## Supplied shocks: validated (and matched by name) BEFORE the order
  ## dispatch so the order-2 / order-3 simulators get the same contract.
  exo <- dr$exo_names
  if (!is.null(shocks)) {
    if (!is.matrix(shocks) || !is.numeric(shocks))
      stop("simulate_model: shocks must be a numeric matrix with n_periods + ",
           "burn_in rows and one column per exogenous shock.", call. = FALSE)
    if (ncol(shocks) != length(exo))
      stop(sprintf("simulate_model: supplied shocks has %d columns; the model has %d shocks (%s).",
                   ncol(shocks), length(exo), paste(exo, collapse = ", ")), call. = FALSE)
    if (!is.null(colnames(shocks))) {
      if (anyDuplicated(colnames(shocks)) || !setequal(colnames(shocks), exo))
        stop("simulate_model: shocks column names (", paste(colnames(shocks), collapse = ", "),
             ") do not match the model's exogenous names (", paste(exo, collapse = ", "),
             ").", call. = FALSE)
      shocks <- shocks[, exo, drop = FALSE]
    }
    if (nrow(shocks) != n_periods + burn_in)
      stop(sprintf(paste0("simulate_model: supplied shocks has %d rows; needs ",
                          "n_periods + burn_in = %d rows, or pass burn_in = 0."),
                   nrow(shocks), n_periods + burn_in), call. = FALSE)
  }

  ## Native-order dispatch (DR3 inherits DR2, so test DR3 first).
  if (!isTRUE(linear)) {
    if (inherits(dr, "DecisionRules3"))
      return(simulate_model_order3(dr, n_periods = n_periods, model = model,
                                   burn_in = burn_in, shocks = shocks,
                                   init_state = init_state))
    if (inherits(dr, "DecisionRules2"))
      return(simulate_model_order2(dr, n_periods = n_periods, model = model,
                                   burn_in = burn_in, shocks = shocks,
                                   init_state = init_state))
  }

  ghx <- dr$ghx
  ghu <- dr$ghu
  endo <- dr$endo_names
  exo  <- dr$exo_names
  state_idx <- dr$state_idx
  n_endo <- length(endo)
  n_exo  <- length(exo)
  params <- if (!is.null(model)) model$param_values else NULL

  total_periods <- n_periods + burn_in

  ## Get shock standard deviations
  shock_stderr <- .get_shock_stderr(model, exo, params)

  if (is.null(shocks)) {
    alpha   <- .get_shock_skewness(model, exo, params)
    sigma_e <- shock_stderr[exo]   # named vector, length n_exo

    ## The drawn shocks must carry the model's FULL covariance --
    ## `corr a, b` / `var a, b` entries in the shocks block were previously
    ## ignored here and every path was drawn with independent innovations,
    ## which silently contradicted compute_moments()/compute_irfs() and the
    ## Kalman filter (all of which use .get_shock_cov()).
    Sigma_e_sim <- .get_shock_cov(model, exo, params)
    off_sim <- Sigma_e_sim
    diag(off_sim) <- 0
    has_corr <- any(abs(off_sim) > 0)

    if (any(alpha != 0) && has_corr) {
      ## L2: CORRELATED + SKEWED shocks.  This branch used to warn that
      ## the declared `corr` / `var a,b` entries were IGNORED and draw each
      ## shock independently from its skew-normal marginal -- a simulation whose
      ## DGP silently differed from the model's declared covariance AND from the
      ## law the package's own skew likelihood evaluates.  It now draws from the
      ## joint closed skew-normal that .get_csn_shock_params() (R/pskf-likelihood.R)
      ## builds, so simulate_model() and the PSKF likelihood share one DGP.
      shocks <- .draw_csn_shocks(total_periods, Sigma_e_sim, alpha, sigma_e)
    } else if (any(alpha != 0)) {
      ## Diagonal Sigma_e: the per-shock construction below is the SAME law as
      ## .draw_csn_shocks() (with Sigma_e diagonal the truncation latents are
      ## independent, so the joint CSN factorises into its marginals), and it is
      ## kept because it consumes a fixed RNG stream that seeded regression
      ## tests depend on, and needs no rejection step.
      ## CSN draw: e_i = sigma_i * (delta_i*|z1| + sqrt(1-delta_i^2)*z2) - mu_i
      ## where mu_i = sigma_i * delta_i * sqrt(2/pi)  (subtract mean
      ## so shocks are mean-zero; sign error here shifts the steady state).
      delta <- alpha / sqrt(1 + alpha^2)            # length n_exo
      mu_e  <- sigma_e * delta * sqrt(2 / pi)       # E[e_i] before correction

      z1 <- matrix(abs(rnorm(total_periods * n_exo)), total_periods, n_exo)
      z2 <- matrix(    rnorm(total_periods * n_exo),  total_periods, n_exo)

      ## delta is per-SHOCK (length n_exo): apply along columns via sweep --
      ## `delta * z1` would recycle delta down the ROWS (column-major).
      shocks_std <- sweep(z1, 2L, delta, `*`) +
                    sweep(z2, 2L, sqrt(1 - delta^2), `*`)
      ## Column k scaled by sigma_e[k], then subtract mean correction mu_e[k]
      shocks <- sweep(sweep(shocks_std, 2L, sigma_e, `*`), 2L, mu_e, `-`)
    } else if (!has_corr) {
      ## Gaussian, diagonal Sigma_e: unchanged (regression-safe, and keeps the
      ## exact RNG stream of every existing seeded test).
      shocks <- matrix(rnorm(total_periods * n_exo), ncol = n_exo)
      for (k in seq_along(exo)) shocks[, k] <- shocks[, k] * sigma_e[k]
    } else {
      ## Gaussian, correlated Sigma_e: e_t = L z_t with Sigma_e = L L'
      ## (lower Cholesky, declared shock order -- the same factor
      ## compute_irfs() / compute_moments() orthogonalise with).  MASS-free.
      L_sim  <- .sigma_e_chol_lower(Sigma_e_sim)
      z      <- matrix(rnorm(total_periods * n_exo), ncol = n_exo)
      shocks <- z %*% t(L_sim)
    }
  }

  sim <- matrix(0, nrow = total_periods, ncol = n_endo)
  colnames(sim) <- endo

  ## State entering period 1.  Defaults to the steady state (deviation 0).
  ## A supplied `init_state` (named deviations over endo names) conditions the
  ## path on a custom starting state -- only the state_idx entries matter, the
  ## rest of y is overwritten in period 1.  Pair with burn_in = 0 to keep the
  ## initial state from being washed out.
  y <- rep(0, n_endo)
  if (!is.null(init_state)) {
    idx <- match(names(init_state), endo)
    ok  <- !is.na(idx)
    if (any(ok)) y[idx[ok]] <- as.numeric(init_state)[ok]
  }
  for (t in seq_len(total_periods)) {
    e <- shocks[t, ]
    y_state <- y[state_idx]
    y <- as.numeric(ghx %*% y_state + ghu %*% e)
    sim[t, ] <- y
  }

  ## Discard burn-in
  sim <- sim[(burn_in + 1):total_periods, , drop = FALSE]

  ## Add steady state for levels
  sim_levels <- sim
  for (j in seq_along(endo)) {
    sim_levels[, j] <- sim[, j] + dr$ys[endo[j]]
  }

  attr(sim, "levels") <- sim_levels
  sim
}


## Simulate a decision rule at its native perturbation order.
##
## simulate_model() uses ONLY ghx/ghu (first-order); calling it on a
## DecisionRules2/3 silently drops the second-order risk/volatility correction
## (ghss/ghxx/ghuu) -- a recurring bug in the welfare/Ramsey paths where the
## code claimed to "work for any order". This dispatcher routes an order-3 rule
## to simulate_model_order3(), an order-2 rule to simulate_model_order2(), and
## a first-order rule to simulate_model(). DecisionRules3 MUST be caught before
## DecisionRules2 because inherits(dr3, "DecisionRules2") is TRUE (class
## hierarchy: DR3 > DR2 > DR). All three return an object carrying
## attr(., "levels"), so callers are unchanged.
## @noRd
.simulate_dr_any_order <- function(dr, n_periods = 200L, model = NULL,
                                   burn_in = 100L, shocks = NULL,
                                   init_state = NULL) {
  if (inherits(dr, "DecisionRules3")) {
    simulate_model_order3(dr, n_periods = as.integer(n_periods),
                          model = model, burn_in = as.integer(burn_in),
                          shocks = shocks, init_state = init_state)
  } else if (inherits(dr, "DecisionRules2")) {
    simulate_model_order2(dr, n_periods = as.integer(n_periods),
                          model = model, burn_in = as.integer(burn_in),
                          shocks = shocks, init_state = init_state)
  } else {
    simulate_model(dr, n_periods = as.integer(n_periods),
                   model = model, burn_in = as.integer(burn_in),
                   shocks = shocks, init_state = init_state)
  }
}

## ============================================================================
## 3. THEORETICAL MOMENTS (Lyapunov equation)
## ============================================================================

## NOTE: solve_lyapunov() used to live here.  It is now the
## package's SINGLE discrete-Lyapunov solver and lives in R/solve-helpers.R
## (together with the `.solve_lyapunov()` alias that R/backend-monolith.R used
## to define separately).  Do not re-add a local copy: the near-unit-root
## relative-tolerance / stability-gate / NaN contract only exists in one place
## on purpose.

## ---------------------------------------------------------------------------
## Stationary-subspace projection via the modal (eigen) decomposition.
##
## For a state transition A with unit/explosive eigenvalues, the unconditional
## state covariance Sigma_s = sum_{k>=0} A^k Q (A')^k DIVERGES.  Dynare reports
## NaN only for the nonstationary variables and FINITE moments for the rest.
## We reproduce that by working in the modal basis A = V Lambda V^{-1}:
##
##   z = V^{-1} x,  z_t = Lambda z_{t-1} + V^{-1} (shock loading) e_t
##   modal innovation cov   M = V^{-1} Q V^{-H}
##   modal stationary cov   S[i,j] = M[i,j] / (1 - lam_i conj(lam_j))
##
## The stationary-SUBSPACE projection keeps only mode pairs where BOTH modes
## are stable (|lam| < 1); pairs touching a nonstationary mode are dropped
## (their variance is infinite).  Sigma_s_finite = Re(V S V^H) is then the
## covariance of the projection of the state onto the stable invariant
## subspace -- finite, and exact for stationary directions.
##
## A variable is nonstationary iff it loads on a nonstationary mode; the test
## is basis-correct (it uses the modal loadings ghx %*% V), unlike a test on
## the raw state-coordinate columns of ghx, which is wrong when the transition
## matrix is dense (the usual case -- ghx_state is generally NOT triangular,
## so abs(diag(ghx_state)) does NOT equal the eigenvalues).
## ---------------------------------------------------------------------------

## Eigendecomposition of the state transition + mode classification.
## Returns ok = FALSE when V is too ill-conditioned to invert reliably
## (near-defective at a repeated unit root), so the caller can fall back.
.modal_decomp <- function(A, tol = 1e-6) {
  eg  <- eigen(A)
  lam <- eg$values
  V   <- eg$vectors
  rc  <- tryCatch(rcond(V), error = function(e) 0)
  if (!is.finite(rc) || rc < 1e-12) return(list(ok = FALSE))
  W <- tryCatch(solve(V), error = function(e) NULL)
  if (is.null(W)) return(list(ok = FALSE))
  list(ok = TRUE, lam = lam, V = V, W = W,
       stat    = which(Mod(lam) <  1 - tol),
       nonstat = which(Mod(lam) >= 1 - tol))
}

## Stationary-subspace covariance for innovation cov Q, given a modal decomp.
## Sigma = Re( V[,stat] %*% (M_ss / (1 - lam_i conj(lam_j))) %*% V[,stat]^H ).
.modal_project_cov <- function(dec, Q) {
  n   <- nrow(Q)
  stat <- dec$stat
  S <- matrix(0 + 0i, n, n)
  if (length(stat) > 0L) {
    M  <- dec$W %*% Q %*% Conj(t(dec$W))          # modal innovation cov
    li <- dec$lam[stat]
    denom <- 1 - outer(li, Conj(li))              # |stat| x |stat|
    S[stat, stat] <- M[stat, stat, drop = FALSE] / denom
  }
  Sigma <- Re(dec$V %*% S %*% Conj(t(dec$V)))
  (Sigma + t(Sigma)) / 2                          # symmetrize FP noise
}

#' Compute theoretical moments of the model
#'
#' Solves the Lyapunov equation for the unconditional variance-covariance
#' matrix and returns correlations, autocorrelations, and variance
#' decompositions.
#'
#' \strong{Dynare name mapping:} the field \code{var_cov} corresponds to
#' Dynare's \code{oo_.var}; there is no \code{variance_covariance} field.
#' The full list of returned fields is:
#' \describe{
#'   \item{var_cov}{n_endo x n_endo unconditional variance-covariance matrix
#'     (Dynare: \code{oo_.var})}
#'   \item{std_dev}{named vector of unconditional standard deviations}
#'   \item{correlation}{contemporaneous correlation matrix}
#'   \item{autocorr}{n_endo x n_endo x n_ar array of autocorrelation matrices}
#'   \item{var_decomp}{n_endo x n_exo variance-decomposition (shares, not
#'     percentages)}
#'   \item{var_decomp_pct}{same, scaled to 100}
#' }
#'
#' \strong{Correlated shocks (ordering dependence).} When \eqn{\Sigma_e} has
#' non-zero off-diagonal entries the shocks are orthogonalised with its lower
#' Cholesky factor, \eqn{\Sigma_e = L L'}, taken in the shock order the model
#' DECLARES them (\code{dr$exo_names}) -- Dynare's \code{stoch_simul}
#' convention. Shock \eqn{k} is credited the variance generated by column
#' \eqn{k} of \eqn{g_u L}, so the per-shock contributions sum exactly to each
#' variable's total variance. Like every Cholesky decomposition this is
#' \emph{order dependent}: reordering \code{varexo} reallocates the shared
#' (covariance) variance between the correlated shocks. With a diagonal
#' \eqn{\Sigma_e} the decomposition is order-invariant and identical to the
#' naive per-shock contribution \eqn{g_{u,k} \sigma_k^2 g_{u,k}'}.
#'
#' @param dr DecisionRules object
#' @param model dynhr_mod object (for shock covariance)
#' @param n_ar Number of autocorrelation lags to compute
#' @param params Named numeric parameter vector (default: model$param_values)
#' @return Named list; see Details.
#' @export
compute_moments <- function(dr, model, n_ar = 5L, params = NULL) {
  ghx <- dr$ghx
  ghu <- dr$ghu
  endo <- dr$endo_names
  exo  <- dr$exo_names
  state_idx <- dr$state_idx
  n_endo <- length(endo)
  n_exo  <- length(exo)
  n_state <- length(state_idx)

  ## FIX: respect caller-supplied params (was unconditionally overwritten)
  if (is.null(params)) params <- model$param_values

  ## Build shock covariance matrix Sigma_e (handles cross-shock correlations).
  ## Honour a covariance carried on the decision rules (set when the caller
  ## passed solve_perturbation(Sigma_e=)) -- including off-diagonal correlations --
  ## in preference to re-deriving from model$shocks. Falls back to the parsed
  ## shocks block when dr$Sigma_e is absent (common case, byte-identical to before).
  Sigma_e <- if (!is.null(dr$Sigma_e) &&
                 all(dim(dr$Sigma_e) == c(n_exo, n_exo))) {
    dr$Sigma_e
  } else {
    .get_shock_cov(model, exo, params)
  }

  ## State transition matrix (for state vars only)
  ghx_state <- ghx[state_idx, , drop = FALSE]   # n_state x n_state

  ## Shock impact on states
  ghu_state <- ghu[state_idx, , drop = FALSE]    # n_state x n_exo

  ## ---- stationary-subspace projection ------------------------------------
  ## Detect unit/explosive roots via the EIGENVALUES of ghx_state (exact and
  ## basis-independent).  When all roots are stable, solve the Lyapunov
  ## equation as usual.  When a unit/explosive root is present, project onto
  ## the stable invariant subspace (modal decomposition) so stationary
  ## variables keep finite moments and only nonstationary ones become NaN --
  ## matching Dynare.  See .modal_decomp / .modal_project_cov above.
  nonstat_tol <- getOption("dynhr.unit_root_tol", 1e-6)
  ## A model with NO state variables (purely static in the shocks) has a 0x0
  ## transition; eigen() / the Lyapunov solver fail on it, and its state
  ## covariance is simply empty.
  ev <- if (n_state > 0L) eigen(ghx_state, only.values = TRUE)$values else complex(0)
  has_unit_root <- any(Mod(ev) >= 1 - nonstat_tol)

  Q_state <- ghu_state %*% Sigma_e %*% t(ghu_state)   # innovation cov of states
  dec <- NULL
  nonstat_var <- integer(0)

  if (n_state == 0L) {
    Sigma_state <- matrix(0, 0L, 0L)
  } else if (!has_unit_root) {
    Sigma_state <- solve_lyapunov(ghx_state, Q_state)
  } else {
    if (getOption("dynhr.warn_unit_root", TRUE)) {
      .dynhr_warn(sprintf(
        paste0("compute_moments(): model has %d nonstationary root(s) ",
               "(|eigenvalue| >= %.6f). Finite moments are returned for the ",
               "stationary subspace; variables loading on a nonstationary mode ",
               "have NaN variance. Suppress with ",
               "options(dynhr.warn_unit_root = FALSE)."),
        sum(Mod(ev) >= 1 - nonstat_tol), 1 - nonstat_tol
      ), call. = FALSE)
    }
    dec <- .modal_decomp(ghx_state, nonstat_tol)
    if (!isTRUE(dec$ok)) {
      ## Near-defective transition (e.g. repeated unit root): cannot project
      ## reliably -> report all-NaN rather than a silently-wrong finite number.
      Sigma_state <- matrix(NaN, n_state, n_state)
      nonstat_var <- seq_len(n_endo)
    } else {
      Sigma_state <- .modal_project_cov(dec, Q_state)
      ## A variable is nonstationary iff it loads on a nonstationary mode.
      ## Modal loadings of endo vars: ghx %*% V (basis-correct).
      GV <- ghx %*% dec$V
      load_ns <- rowSums(Mod(GV[, dec$nonstat, drop = FALSE]))
      nonstat_var <- which(load_ns > 1e-8)
    }
  }

  ## Full variance-covariance of all endogenous variables.  Sigma_state is
  ## finite (stationary-subspace projection), so the product is finite
  ## everywhere; we then NaN-out only the nonstationary variables.
  Sigma_y <- ghx %*% Sigma_state %*% t(ghx) + ghu %*% Sigma_e %*% t(ghu)
  if (length(nonstat_var) > 0L) {
    Sigma_y[nonstat_var, ] <- NaN
    Sigma_y[, nonstat_var] <- NaN
  }

  rownames(Sigma_y) <- endo
  colnames(Sigma_y) <- endo

  ## Standard deviations
  variances <- pmax(diag(Sigma_y), 0)
  std_dev <- sqrt(variances)
  names(std_dev) <- endo

  ## Correlation matrix
  sd_outer <- outer(std_dev, std_dev)
  sd_outer[sd_outer == 0] <- Inf
  corr_mat <- Sigma_y / sd_outer
  diag(corr_mat) <- 1
  rownames(corr_mat) <- endo
  colnames(corr_mat) <- endo

  ## Autocorrelations
  autocorr <- array(NaN, dim = c(n_endo, n_endo, n_ar))
  dimnames(autocorr) <- list(endo, endo, paste0("lag", 1:n_ar))

  S_sel <- matrix(0, nrow = n_state, ncol = n_endo)
  for (i in seq_along(state_idx)) S_sel[i, state_idx[i]] <- 1

  ## SGU-2b fix: restrict the lag-covariance recursion to the STATIONARY
  ## subspace so that NaN from unit-root rows/cols of Sigma_y does not bleed
  ## into the stationary-variable autocorrelations.
  ##
  ## Strategy:
  ##   - stat_var : indices into 1:n_endo that are stationary (complement of
  ##                nonstat_var, which was determined above by the SAME
  ##                modal-eigenvalue criterion used for var_cov).
  ##   - The transition operator on the full space is G = ghx %*% S_sel
  ##     (n_endo x n_endo).  Restricted to the stationary block it is
  ##     G_ss = G[stat_var, stat_var] (stat x stat).
  ##   - Seed the recursion from the stationary block of Sigma_y, which is
  ##     guaranteed finite by the M2 fix.
  ##   - Results are placed into autocorr[stat_var, stat_var, lag]; all entries
  ##     touching a nonstationary variable remain NaN (the initial value).
  stat_var <- if (length(nonstat_var) > 0L) {
    setdiff(seq_len(n_endo), nonstat_var)
  } else {
    seq_len(n_endo)
  }

  if (length(stat_var) > 0L) {
    G       <- ghx %*% S_sel                        # n_endo x n_endo
    G_ss    <- G[stat_var, stat_var, drop = FALSE]   # n_stat x n_stat
    sd_ss   <- sd_outer[stat_var, stat_var, drop = FALSE]
    sd_ss[sd_ss == 0] <- Inf
    ## Seed: stationary block of Sigma_y (finite by construction)
    Gamma_ss <- Sigma_y[stat_var, stat_var, drop = FALSE]
    for (lag in seq_len(n_ar)) {
      Gamma_ss    <- G_ss %*% Gamma_ss
      autocorr[stat_var, stat_var, lag] <- Gamma_ss / sd_ss
    }
  }

  ## Variance decomposition: contribution of each shock to each variable.
  ## Mirrors the stationary-subspace logic above for per-shock Lyapunov solves.
  ##
  ## Correlated shocks (Dynare convention): the innovations are orthogonalised
  ## with the lower Cholesky factor of Sigma_e in the DECLARED shock order,
  ## Sigma_e = L L'.  Shock k is credited the variance generated by column k of
  ## ghu %*% L (unit-variance orthogonal innovation), so the per-shock
  ## contributions sum EXACTLY to the total variance and no covariance term is
  ## left unattributed.  For a diagonal Sigma_e, L = diag(sd) and this reduces
  ## to the old ghu[, k] * sigma_k contribution byte-for-byte.
  L_e        <- .sigma_e_chol_lower(Sigma_e)
  ghu_L      <- ghu %*% L_e                       # n_endo x n_exo
  ghu_L_st   <- ghu_L[state_idx, , drop = FALSE]  # n_state x n_exo

  var_decomp <- matrix(0, nrow = n_endo, ncol = n_exo)
  rownames(var_decomp) <- endo
  colnames(var_decomp) <- exo

  for (k in seq_along(exo)) {
    ghu_state_k <- ghu_L_st[, k, drop = FALSE]
    ghu_k       <- ghu_L[, k, drop = FALSE]
    Q_state_k   <- ghu_state_k %*% t(ghu_state_k)

    if (!has_unit_root) {
      Sigma_state_k <- solve_lyapunov(ghx_state, Q_state_k)
    } else if (isTRUE(dec$ok)) {
      ## Reuse the eigendecomposition; project shock k's innovation cov onto
      ## the stable subspace (same machinery as the full covariance above).
      Sigma_state_k <- .modal_project_cov(dec, Q_state_k)
    } else {
      Sigma_state_k <- matrix(0, n_state, n_state)   # cannot decompose
    }

    ## Full contribution of shock k (finite); nonstationary vars are zeroed so
    ## the per-variable shares still sum sensibly for the stationary block.
    Sigma_y_k <- ghx %*% Sigma_state_k %*% t(ghx) + tcrossprod(ghu_k)
    dk <- diag(Sigma_y_k)
    if (length(nonstat_var) > 0L) dk[nonstat_var] <- 0
    var_decomp[, k] <- pmax(replace(dk, is.nan(dk), 0), 0)
  }

  ## Normalise to percentages
  total_var <- rowSums(var_decomp)
  total_var[total_var == 0] <- 1
  var_decomp_pct <- var_decomp / total_var * 100

  list(
    var_cov = Sigma_y,
    std_dev = std_dev,
    correlation = corr_mat,
    autocorr = autocorr,
    var_decomp = var_decomp,
    var_decomp_pct = var_decomp_pct,
    Sigma_e = Sigma_e,
    Sigma_state = Sigma_state
  )
}

## ============================================================================
## 3a-ii. CONDITIONAL VARIANCE DECOMPOSITION (FEVD at finite and infinite horizons)
## ============================================================================

#' Forecast-error variance decomposition at finite and infinite horizons
#'
#' Computes the h-step forecast-error variance decomposition (FEVD) for each
#' endogenous variable and shock at the requested horizons. The decomposition
#' follows the Dynare convention: horizon \eqn{h} gives the variance of the
#' h-step-ahead forecast error of \eqn{y_{t+h}} given information at \eqn{t}.
#'
#' @details
#' The state-space representation (lagged-state convention) is:
#' \deqn{s_t = A s_{t-1} + B \varepsilon_t, \quad
#'       y_t = C s_{t-1} + D \varepsilon_t}
#' where \eqn{A = \code{ghx[state\_idx, ]}}, \eqn{B = \code{ghu[state\_idx, ]}},
#' \eqn{C = \code{ghx}}, \eqn{D = \code{ghu}}.
#'
#' The h-step forecast-error covariance is:
#' \deqn{V(h) = \sum_{j=0}^{h-1} \Psi_j \Sigma_e \Psi_j'}
#' with \eqn{\Psi_0 = D} and \eqn{\Psi_j = C A^{j-1} B} for \eqn{j \ge 1}.
#' Writing \eqn{\Sigma_e = L L'} for the lower Cholesky factor in the DECLARED
#' shock order (Dynare's convention), the contribution of shock \eqn{k} to
#' variable \eqn{i} at horizon \eqn{h} is
#' \eqn{\sum_{j=0}^{h-1} (\Psi_j L)_{ik}^2}, so the shares sum exactly to
#' \eqn{V(h)_{ii}} even with correlated shocks. With a diagonal \eqn{\Sigma_e}
#' this is \eqn{\sum_j (\Psi_j)_{ik}^2 \Sigma_e[k,k]}, as before. As with any
#' Cholesky orthogonalisation the split of the shared variance between
#' correlated shocks depends on the \code{varexo} declaration order.
#'
#' For \code{horizon = Inf}, the per-shock Lyapunov solution is used (same
#' machinery as \code{\link{compute_moments}}), which equals the limit of the
#' finite-horizon recursion.
#'
#' @param dr     DecisionRules object (from \code{\link{solve_perturbation}}).
#' @param model  dynhr_mod object (for shock covariance).
#' @param horizons Integer or numeric vector of forecast horizons.
#'   \code{Inf} is allowed and triggers the Lyapunov-based unconditional
#'   decomposition. Default: \code{c(1L, 4L, 8L, 16L, 40L, Inf)}.
#' @param params Named numeric parameter vector (default: \code{model$param_values}).
#'
#' @return A named list with:
#'   \describe{
#'     \item{fevd}{n_endo x n_exo x n_horizons array of raw forecast-error
#'       variance contributions (absolute, not percentages).}
#'     \item{fevd_pct}{Same array normalised so each variable's shares sum to
#'       100 over shocks (i.e. row-normalised percentages).}
#'     \item{horizons}{The horizon vector used (matching the third dimension of
#'       \code{fevd}).}
#'   }
#'   \code{dimnames(fevd)} = \code{list(endo_names, exo_names, as.character(horizons))}.
#'
#' @seealso \code{\link{compute_moments}} for unconditional moments.
#' @export
conditional_variance_decomposition <- function(dr, model,
                                               horizons = c(1L, 4L, 8L, 16L, 40L, Inf),
                                               params = NULL) {
  ghx       <- dr$ghx
  ghu       <- dr$ghu
  endo      <- dr$endo_names
  exo       <- dr$exo_names
  state_idx <- dr$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)
  n_state   <- length(state_idx)
  n_h       <- length(horizons)

  if (is.null(params)) params <- model$param_values

  ## Shock covariance
  Sigma_e   <- .get_shock_cov(model, exo, params)

  ## State matrices
  ghx_state <- ghx[state_idx, , drop = FALSE]  # A: n_state x n_state
  ghu_state <- ghu[state_idx, , drop = FALSE]  # B: n_state x n_exo

  ## Unit-root detection (reuse same logic as compute_moments)
  nonstat_tol <- getOption("dynhr.unit_root_tol", 1e-6)
  ev <- eigen(ghx_state, only.values = TRUE)$values
  has_unit_root <- any(Mod(ev) >= 1 - nonstat_tol)
  dec <- if (has_unit_root) .modal_decomp(ghx_state, nonstat_tol) else NULL

  ## Orthogonalise the innovations with the lower Cholesky factor of Sigma_e in
  ## the DECLARED shock order (Dynare convention; see compute_moments()).  The
  ## impulse matrices are then post-multiplied by L and each shock's
  ## contribution is the squared column, so the shares sum exactly to the
  ## h-step forecast-error variance even when Sigma_e is not diagonal.  For a
  ## diagonal Sigma_e, L = diag(sd) and this is the old Psi[, k]^2 * sigma_k^2.
  L_e <- .sigma_e_chol_lower(Sigma_e)

  ## Output array: n_endo x n_exo x n_h (raw FEV contributions)
  fevd <- array(0, dim = c(n_endo, n_exo, n_h),
                dimnames = list(endo, exo, as.character(horizons)))

  ## Sort finite horizons; handle Inf separately
  h_vals   <- horizons[is.finite(horizons)]
  h_inf    <- which(!is.finite(horizons))
  h_finite <- which(is.finite(horizons))

  ## ---- Finite horizons via iterative Psi recursion --------------------------
  if (length(h_finite) > 0L) {
    max_h <- max(h_vals)

    ## Psi_0 = D = ghu;  for j>=1: Psi_j = C * (A^{j-1} B) = ghx * phi_{j-1}
    ## phi_0 = B = ghu_state; phi_{j} = A * phi_{j-1} = ghx_state * phi_{j-1}
    Psi_cur <- ghu %*% L_e         # n_endo x n_exo  (Psi_0 L)
    phi     <- ghu_state %*% L_e   # n_state x n_exo (phi_0 = B L)

    ## Running FEV accumulator
    FEV <- Psi_cur^2

    for (hidx in h_finite) {
      if (horizons[hidx] == 1L || horizons[hidx] == 1) {
        fevd[,, hidx] <- FEV
      }
    }

    for (j in seq_len(max_h - 1L)) {
      ## Advance: phi_j = ghx_state * phi_{j-1}; Psi_j = ghx * phi_{j-1}
      ## (phi currently holds phi_{j-1} at the START of this iteration,
      ##  which gives Psi_j = ghx * phi_{j-1})
      Psi_cur <- ghx %*% phi
      phi     <- ghx_state %*% phi

      ## Accumulate FEV
      FEV <- FEV + Psi_cur^2

      ## Store whenever we've hit a requested horizon (j+1 steps accumulated)
      h_done <- j + 1L
      for (hidx in h_finite) {
        if (horizons[hidx] == h_done) {
          fevd[,, hidx] <- FEV
        }
      }
    }
  }

  ## ---- Infinite horizon via per-shock Lyapunov solve ------------------------
  if (length(h_inf) > 0L) {
    ghu_L    <- ghu %*% L_e
    ghu_L_st <- ghu_state %*% L_e
    for (k in seq_along(exo)) {
      ghu_state_k <- ghu_L_st[, k, drop = FALSE]
      ghu_k       <- ghu_L[, k, drop = FALSE]
      Q_state_k   <- tcrossprod(ghu_state_k)

      if (!has_unit_root) {
        Sigma_state_k <- solve_lyapunov(ghx_state, Q_state_k)
      } else if (isTRUE(dec$ok)) {
        Sigma_state_k <- .modal_project_cov(dec, Q_state_k)
      } else {
        Sigma_state_k <- matrix(0, n_state, n_state)
      }

      Sigma_y_k <- ghx %*% Sigma_state_k %*% t(ghx) + tcrossprod(ghu_k)
      dk <- diag(Sigma_y_k)
      dk <- pmax(replace(dk, is.nan(dk), 0), 0)

      for (hidx in h_inf) {
        fevd[, k, hidx] <- dk
      }
    }
  }

  ## ---- Normalise to percentages (row-sum to 100 per variable per horizon) ---
  fevd_pct <- array(0, dim = dim(fevd), dimnames = dimnames(fevd))
  for (hidx in seq_len(n_h)) {
    row_totals <- rowSums(fevd[,, hidx, drop = FALSE])
    row_totals[row_totals == 0] <- 1  # avoid 0/0
    fevd_pct[,, hidx] <- fevd[,, hidx] / row_totals * 100
  }

  list(
    fevd     = fevd,
    fevd_pct = fevd_pct,
    horizons = horizons
  )
}

## ============================================================================
## 3a-iii. HP-FILTERED THEORETICAL MOMENTS (spectral integration)
## ============================================================================

#' HP-filtered theoretical second moments via spectral integration
#'
#' Computes the theoretical variance-covariance matrix, standard deviations,
#' and correlations of all endogenous variables after applying the
#' Hodrick-Prescott filter with smoothing parameter \eqn{\lambda}.
#'
#' @details
#' The HP filter has squared gain function:
#' \deqn{G_{\rm hp}(\omega) =
#'   \left(\frac{4\lambda(1-\cos\omega)^2}{1 + 4\lambda(1-\cos\omega)^2}\right)^2}
#' The HP-filtered variance of variable \eqn{i} is:
#' \deqn{\mathrm{Var}_{\rm hp}(y_i) =
#'   \frac{1}{\pi} \int_0^\pi G_{\rm hp}(\omega)\, S_{ii}(\omega)\, d\omega}
#' where \eqn{S_{yy}(\omega)} is the one-sided power spectral density matrix
#' (returned by \code{.spectral_density_core} without a normalisation factor,
#' so \eqn{(1/\pi)\int_0^\pi S_{ii}(\omega)\,d\omega = \mathrm{Var}(y_i)}).
#'
#' Numerical integration uses a uniform \code{n_freq}-point grid on
#' \eqn{(0, \pi]} (Riemann mid-point rule). The grid avoids \eqn{\omega=0}
#' where \eqn{S(\omega)} may diverge for near-unit-root models.
#'
#' @param dr      DecisionRules object (from \code{\link{solve_perturbation}}).
#' @param model   dynhr_mod object (for shock covariance).
#' @param lambda  HP smoothing parameter. Default 1600 (quarterly data).
#' @param n_freq  Number of quadrature points on \eqn{(0,\pi]}. Default 512.
#' @param params  Named numeric parameter vector (default: \code{model$param_values}).
#'
#' @return A named list:
#'   \describe{
#'     \item{var_cov}{n_endo x n_endo HP-filtered variance-covariance matrix.}
#'     \item{std_dev}{Named vector of HP-filtered standard deviations.}
#'     \item{correlation}{HP-filtered correlation matrix.}
#'     \item{lambda}{The \eqn{\lambda} value used.}
#'   }
#'
#' @seealso \code{\link{compute_moments}}, \code{\link{spectral_density}}
#' @export
hp_filtered_moments <- function(dr, model, lambda = 1600, n_freq = 512L,
                                 params = NULL) {
  if (is.null(params)) params <- model$param_values
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda < 0)
    stop("hp_filtered_moments: 'lambda' must be a single non-negative number.")

  endo      <- dr$endo_names
  exo       <- dr$exo_names
  state_idx <- dr$state_idx
  n_endo    <- length(endo)

  ## Build state-space matrices directly from dr (avoids needing m$lead_lag_incidence)
  ghx_state <- dr$ghx[state_idx, , drop = FALSE]  # TT: n_state x n_state
  ghu_state <- dr$ghu[state_idx, , drop = FALSE]  # RR: n_state x n_exo
  ZZ        <- dr$ghx                              # n_endo x n_state (observe all)
  DD        <- dr$ghu                              # n_endo x n_exo

  Sigma_e   <- .get_shock_cov(model, exo, params)

  ## Uniform quadrature grid on (0, pi]  (avoid omega=0 for near-unit-root safety)
  omegas <- seq(pi / n_freq, pi, length.out = n_freq)
  dw     <- pi / n_freq   # width of each bin (Riemann rule)

  ## HP squared gain at each frequency
  G_hp <- function(w) {
    x <- 4 * lambda * (1 - cos(w))^2
    (x / (1 + x))^2
  }
  gains <- G_hp(omegas)   # length n_freq

  ## Accumulate (1/pi) * sum_j G_hp(w_j) * Re(S(w_j)) * dw
  ## = (dw/pi) * sum_j G_hp(w_j) * Re(S(w_j))
  Sigma_hp <- matrix(0, n_endo, n_endo)
  for (idx in seq_along(omegas)) {
    S_w <- .spectral_density_core(omegas[idx],
                                   TT      = ghx_state,
                                   RR      = ghu_state,
                                   ZZ      = ZZ,
                                   DD      = DD,
                                   Sigma_e = Sigma_e)
    Sigma_hp <- Sigma_hp + gains[idx] * Re(S_w)
  }
  Sigma_hp <- Sigma_hp * (dw / pi)

  ## Symmetrize (floating-point rounding can break exact symmetry)
  Sigma_hp <- (Sigma_hp + t(Sigma_hp)) / 2
  rownames(Sigma_hp) <- endo
  colnames(Sigma_hp) <- endo

  ## Standard deviations and correlations
  variances <- pmax(diag(Sigma_hp), 0)
  std_dev   <- sqrt(variances)
  names(std_dev) <- endo

  sd_outer <- outer(std_dev, std_dev)
  sd_outer[sd_outer == 0] <- Inf
  corr_mat <- Sigma_hp / sd_outer
  diag(corr_mat) <- 1
  rownames(corr_mat) <- endo
  colnames(corr_mat) <- endo

  list(
    var_cov     = Sigma_hp,
    std_dev     = std_dev,
    correlation = corr_mat,
    lambda      = lambda
  )
}

## ============================================================================
## 3b. ORDER-2 THEORETICAL MOMENTS (pruned second-order decision rules)
## ============================================================================

#' Compute theoretical unconditional moments of a pruned second-order DSGE
#'
#' Returns the analytic unconditional mean, variance-covariance, standard
#' deviations, correlations, autocorrelations, and variance decomposition for
#' the full endogenous variable vector under a pruned second-order perturbation
#' solution. Gaussian shocks only.
#'
#' @details
#' The pruned second-order output equation is:
#' \deqn{y_t = y_s + g_{hx} x_{t-1}^{(1)} + g_{hu} \varepsilon_t
#'           + g_{hx} x_{t-1}^{(2)}
#'           + \tfrac{1}{2} g_{hxx} (x_{t-1}^{(1)} \otimes x_{t-1}^{(1)})
#'           + g_{hxu} (\varepsilon_t \otimes x_{t-1}^{(1)})
#'           + \tfrac{1}{2} g_{huu} (\varepsilon_t \otimes \varepsilon_t)
#'           + \tfrac{1}{2} g_{hss}}
#'
#' The unconditional mean (in level space) is:
#' \deqn{E[y] = y_s + g_{hx} E[x^{(2)}]
#'            + \tfrac{1}{2} g_{hxx} \mathrm{vec}(\Sigma_x)
#'            + \tfrac{1}{2} g_{huu} \mathrm{vec}(\Sigma_e)
#'            + \tfrac{1}{2} g_{hss}}
#'
#' where \eqn{E[x^{(2)}] = \tfrac{1}{2}(I - h_x)^{-1}
#'   (h_{xx}\,\mathrm{vec}(\Sigma_x) + h_{uu}\,\mathrm{vec}(\Sigma_e) + h_{ss})}.
#'
#' The variance uses centered fourth moments (Gaussian: \eqn{M4 - \mathrm{vec}(\Sigma)\mathrm{vec}(\Sigma)'})
#' and the third cross-cumulant \eqn{C_3^{211}} to account for the covariance
#' between the second-order state \eqn{x^{(2)}} and the quadratic term
#' \eqn{x^{(1)} \otimes x^{(1)}}.
#'
#' **Scope:** Gaussian shocks only. For skew-normal shocks the fourth-moment
#' formulas differ and are not yet implemented; a warning is issued.
#'
#' \strong{Correlated shocks (ordering dependence).} As in
#' \code{\link{compute_moments}}, \code{var_decomp} orthogonalises the shocks
#' with the lower Cholesky factor of \eqn{\Sigma_e} in the shock order the model
#' DECLARES them, \eqn{\Sigma_e = L L'} (Dynare's convention). Shock \eqn{k} is
#' credited the order-2 variance generated by \eqn{L_{\cdot k} L_{\cdot k}'}
#' alone, so no covariance is left unattributed, and the split is \emph{order
#' dependent}: reordering \code{varexo} reallocates the shared variance between
#' correlated shocks. With a diagonal \eqn{\Sigma_e} it is identical to the
#' per-shock \eqn{\Sigma_e[k,k]} contribution. Note the order-2 map
#' \eqn{\Sigma_e \mapsto \mathrm{Var}(y)} is quadratic, so the shares sum to the
#' total only up to genuine cross-shock interaction terms (exactly, for a linear
#' model).
#'
#' @param dr DecisionRules2 object (from \code{solve_perturbation(order=2)}).
#' @param model dynhr_mod object (for shock covariance).
#' @param n_ar Number of autocorrelation lags to compute (default 5).
#' @param params Named numeric parameter vector (default: \code{model$param_values}).
#' @return A list with slots:
#'   \describe{
#'     \item{mean}{Named n_endo vector: full unconditional mean in level space.}
#'     \item{var_cov}{n_endo x n_endo full order-2 variance-covariance.}
#'     \item{std_dev}{Named n_endo vector of standard deviations.}
#'     \item{correlation}{n_endo x n_endo correlation matrix.}
#'     \item{autocorr}{n_endo x n_endo x n_ar autocorrelation array.}
#'     \item{var_decomp}{n_endo x n_exo raw variance contributions by shock.}
#'     \item{var_decomp_pct}{n_endo x n_exo percent variance decomposition.}
#'     \item{Sigma_e}{n_exo x n_exo shock covariance (passed through).}
#'     \item{Sigma_x}{n_s x n_s first-order state covariance.}
#'     \item{Var_x2}{n_s x n_s second-order state variance.}
#'     \item{mean_x2}{n_s vector: \eqn{E[x_t^{(2)}]} (for diagnostics).}
#'   }
#' @export
compute_moments_order2 <- function(dr, model, n_ar = 5L, params = NULL) {
  if (!inherits(dr, "DecisionRules2")) {
    stop("compute_moments_order2: dr must be a DecisionRules2 object.")
  }

  if (is.null(params)) params <- model$param_values

  ## --- Gaussian-scope guard --------------------------------------------------
  exo <- dr$exo_names
  alpha <- .get_shock_skewness(model, exo, params)
  if (!all(alpha == 0)) {
    .dynhr_warn(
      "compute_moments_order2: one or more shocks have non-zero skewness. ",
      "The formulas assume Gaussian shocks; results will be approximate."
    )
  }

  ## --- Extract matrices ------------------------------------------------------
  ghx  <- dr$ghx;  ghu  <- dr$ghu
  ghxx <- dr$ghxx; ghxu <- dr$ghxu
  ghuu <- dr$ghuu; ghss <- dr$ghss
  ys   <- dr$ys

  endo      <- dr$endo_names
  state_idx <- dr$state_idx
  n_endo    <- length(endo)
  n_exo     <- length(exo)
  n_s       <- length(state_idx)

  ## State sub-matrices (rows = states only)
  hx  <- ghx[state_idx, , drop = FALSE]   # n_s x n_s
  hu  <- ghu[state_idx, , drop = FALSE]   # n_s x n_exo
  hxx <- ghxx[state_idx, , drop = FALSE]  # n_s x n_s²
  hxu <- ghxu[state_idx, , drop = FALSE]  # n_s x (n_s·n_exo)
  huu <- ghuu[state_idx, , drop = FALSE]  # n_s x n_exo²
  hss <- ghss[state_idx]                  # n_s

  ## Shock covariance
  Sigma_e <- .get_shock_cov(model, exo, params)

  ## --- Exact order-2 moments via the augmented pruned state-space -----------
  ## xi_t = [x1_t; x2_t; x1_t (x) x1_t] is linear in the contemporaneous shock,
  ## so mean and Var(y) follow from a deterministic Lyapunov solve (Gaussian
  ## shocks).  This captures the FULL x2 <-> (x1 (x) x1) cross-covariance,
  ## including the hxu/huu shock contributions the earlier per-block formula
  ## dropped.  Validated against cross-sectional simulation to MC precision.
  sys <- .order2_aug_system(dr, Sigma_e)
  st  <- .order2_stationary_moments(sys)
  Sigma_x <- st$Sigma_x
  Var_x2  <- st$Var_x2
  mean_x2 <- st$mean_x2
  mn      <- st$mean
  names(mn) <- endo
  Sigma_y <- st$var_cov
  rownames(Sigma_y) <- endo
  colnames(Sigma_y) <- endo

  ## --- Standard deviations and correlations ----------------------------------
  variances <- pmax(diag(Sigma_y), 0)
  std_dev   <- sqrt(variances)
  names(std_dev) <- endo

  sd_outer <- outer(std_dev, std_dev)
  sd_outer[sd_outer == 0] <- Inf
  corr_mat <- Sigma_y / sd_outer
  diag(corr_mat) <- 1
  rownames(corr_mat) <- endo
  colnames(corr_mat) <- endo

  ## --- Autocovariances (lag τ ≥ 1) ------------------------------------------
  ## Exact augmented-system formula; see .order2_autocov() for the
  ## derivation and for what the old ghx/S_sel recursion got wrong.
  autocov  <- .order2_autocov(sys, st, n_ar)
  autocorr <- array(0, dim = c(n_endo, n_endo, n_ar))
  dimnames(autocorr) <- list(endo, endo, paste0("lag", seq_len(n_ar)))
  for (lag in seq_len(n_ar)) autocorr[, , lag] <- autocov[, , lag] / sd_outer

  ## --- Variance decomposition (exact per-shock via the augmented system) ----
  ## Shock k's contribution = order-2 Var(y) driven by orthogonal innovation k
  ## ALONE, i.e. with Sigma_e replaced by the rank-1 block L[, k] L[, k]', where
  ## Sigma_e = L L' is the lower Cholesky factor in the DECLARED shock order
  ## (Dynare's convention; identical to what compute_moments() and
  ## conditional_variance_decomposition() now do at order 1).
  ##
  ## 0.9.4 : this used to zero every Sigma_e entry
  ## except the diagonal (k,k), which simply DISCARDS the off-diagonal
  ## covariance -- with correlated shocks the per-shock pieces then failed to
  ## add up, and the order-1 and order-2 decompositions of the SAME linear model
  ## disagreed.  Because sum_k L[, k] L[, k]' = L L' = Sigma_e exactly, the
  ## rank-1 split attributes all of the covariance.  For a diagonal Sigma_e,
  ## L = diag(sd) and L[, k] L[, k]' is the old (k,k)-only matrix, so nothing
  ## changes for uncorrelated shocks.
  ##
  ## Like every Cholesky orthogonalisation this is ORDER DEPENDENT: reordering
  ## `varexo` reallocates the shared variance between correlated shocks.  And
  ## note the order-2 map Sigma_e -> Var(y) is quadratic, so the pieces sum to
  ## the total only up to genuine cross-shock interaction terms (they sum
  ## EXACTLY when the model is linear).
  L_e <- .sigma_e_chol_lower(Sigma_e)
  var_decomp <- matrix(0, nrow = n_endo, ncol = n_exo)
  rownames(var_decomp) <- endo
  colnames(var_decomp) <- exo
  for (k in seq_len(n_exo)) {
    Sigma_ek <- tcrossprod(L_e[, k])
    sys_k <- .order2_aug_system(dr, Sigma_ek)
    var_decomp[, k] <- pmax(diag(.order2_stationary_moments(sys_k)$var_cov), 0)
  }

  ## Normalise to percentages
  total_var <- rowSums(var_decomp)
  total_var[total_var == 0] <- 1
  var_decomp_pct <- var_decomp / total_var * 100

  ## --- Return ----------------------------------------------------------------
  list(
    mean           = mn,
    var_cov        = Sigma_y,
    std_dev        = std_dev,
    correlation    = corr_mat,
    autocorr       = autocorr,
    var_decomp     = var_decomp,
    var_decomp_pct = var_decomp_pct,
    Sigma_e        = Sigma_e,
    Sigma_x        = Sigma_x,
    Var_x2         = Var_x2,
    mean_x2        = mean_x2
  )
}

## ============================================================================
## Augmented pruned-state moment machinery (order 2) -- shared by
## compute_moments_order2() (stationary, unconditional) and conditional_welfare()
## (transient, s0-conditional).
##
## The pruned 2nd-order system is LINEAR in the augmented state
##   xi_t = [ x1_t ; x2_t ; x1_t (x) x1_t ]   (dim d = 2*n_s + n_s^2),
## driven only by the contemporaneous shock eps_t.  Its mean and covariance
## therefore propagate by EXACT deterministic recursions:
##   mu_t   = T mu_{t-1} + c + c_u
##   Sxi_t  = T Sxi_{t-1} T' + G Cov(r) G',   Cov(xi_{t-1}, u_t) = 0
## where r = [eps ; eps(x)x1 ; x1(x)eps ; eps(x)eps].  Cov(xi_{t-1}, u_t) = 0
## because every innovation term is odd in eps or eps-independent.  This gives
## the unconditional moments (Lyapunov fixed point) and the s0-conditional
## moments (transient from a known state) with NO Monte Carlo, exact for
## Gaussian shocks.  Validated against cross-sectional simulation to MC
## precision.  (Supersedes the earlier per-block formula, which built the
## x2/x1(x)x1 cross-cumulant from the hxx term only -- omitting the hxu/huu
## shock contributions -- biasing Var(y) by up to ~2% for jump variables.)
## ============================================================================

## Cov(r), r = [eps ; eps(x)x1 ; x1(x)eps ; eps(x)eps] (centered), given
## x1 ~ N(a, P) and eps ~ N(0, Sigma_e) independent.  M = P + a a'.
.order2_cov_r <- function(a, P, Sigma_e) {
  n_u <- nrow(Sigma_e); n_s <- length(a)
  M <- P + outer(a, a)
  vecSe <- as.numeric(Sigma_e)
  C_ee_ee <- .fourth_moment_gaussian(Sigma_e) - outer(vecSe, vecSe)
  ## cross block (eps(x)x1, x1(x)eps): entry (i,j),(k,l) = Sigma_e[i,l]*M[j,k]
  C_ex1_x1e <- matrix(0, n_u * n_s, n_s * n_u)
  for (i in seq_len(n_u)) for (j in seq_len(n_s))
    for (k in seq_len(n_s)) for (l in seq_len(n_u))
      C_ex1_x1e[(i - 1L) * n_s + j, (k - 1L) * n_u + l] <- Sigma_e[i, l] * M[j, k]
  d1 <- n_u; d2 <- n_u * n_s; d3 <- n_s * n_u; d4 <- n_u * n_u
  D  <- d1 + d2 + d3 + d4
  Cr <- matrix(0, D, D)
  ## seq_len offsets, not `a:b`: with no states (n_s = 0) the x1 blocks are
  ## EMPTY, and `(d1 + 1L):(d1 + 0L)` is the DESCENDING pair c(d1 + 1, d1).
  i1 <- seq_len(d1); i2 <- d1 + seq_len(d2)
  i3 <- d1 + d2 + seq_len(d3); i4 <- d1 + d2 + d3 + seq_len(d4)
  Cr[i1, i1] <- Sigma_e
  Cr[i1, i2] <- kronecker(Sigma_e, t(a))
  Cr[i1, i3] <- kronecker(t(a), Sigma_e)
  Cr[i2, i2] <- kronecker(Sigma_e, M)
  Cr[i2, i3] <- C_ex1_x1e
  Cr[i3, i3] <- kronecker(M, Sigma_e)
  Cr[i4, i4] <- C_ee_ee
  Cr[lower.tri(Cr)] <- t(Cr)[lower.tri(Cr)]
  Cr
}

## Constant augmented-system matrices for a DecisionRules2 object + Sigma_e.
.order2_aug_system <- function(dr, Sigma_e) {
  ghx <- dr$ghx; ghu <- dr$ghu; ghxx <- dr$ghxx; ghxu <- dr$ghxu
  ghuu <- dr$ghuu; ghss <- dr$ghss; ys <- dr$ys
  sidx <- dr$state_idx; endo <- dr$endo_names
  n_endo <- length(endo); n_u <- ncol(ghu); n_s <- length(sidx)
  hx <- ghx[sidx, , drop = FALSE]; hu <- ghu[sidx, , drop = FALSE]
  hxx <- ghxx[sidx, , drop = FALSE]; hxu <- ghxu[sidx, , drop = FALSE]
  huu <- ghuu[sidx, , drop = FALSE]; hss <- ghss[sidx]
  vecSe <- as.numeric(Sigma_e)
  d <- 2L * n_s + n_s * n_s
  ## seq_len offsets (see .order2_cov_r): a stateless model (n_s = 0, d = 0)
  ## must give empty blocks, not the descending `(n_s + 1L):(2L * n_s)` = 1:0.
  ix1 <- seq_len(n_s); ix2 <- n_s + seq_len(n_s); ik <- 2L * n_s + seq_len(n_s * n_s)
  Tlin <- matrix(0, d, d)
  Tlin[ix1, ix1] <- hx; Tlin[ix2, ix2] <- hx
  Tlin[ix2, ik] <- 0.5 * hxx; Tlin[ik, ik] <- kronecker(hx, hx)
  cc <- numeric(d); cc[ix2] <- 0.5 * hss
  c_u <- numeric(d)
  c_u[ix2] <- 0.5 * as.numeric(huu %*% vecSe)
  c_u[ik]  <- as.numeric(kronecker(hu, hu) %*% vecSe)
  D1 <- n_u; D2 <- n_u * n_s; D3 <- n_s * n_u; D4 <- n_u * n_u
  Dr <- D1 + D2 + D3 + D4
  j1 <- seq_len(D1); j2 <- D1 + seq_len(D2)
  j3 <- D1 + D2 + seq_len(D3); j4 <- D1 + D2 + D3 + seq_len(D4)
  G <- matrix(0, d, Dr)
  G[ix1, j1] <- hu; G[ix2, j2] <- hxu; G[ix2, j4] <- 0.5 * huu
  G[ik, j2] <- kronecker(hu, hx); G[ik, j3] <- kronecker(hx, hu)
  G[ik, j4] <- kronecker(hu, hu)
  Dxi <- cbind(ghx, ghx, 0.5 * ghxx)
  Gv  <- matrix(0, n_endo, Dr)
  Gv[, j1] <- ghu; Gv[, j2] <- ghxu; Gv[, j4] <- 0.5 * ghuu
  r_mean <- numeric(Dr); r_mean[j4] <- vecSe
  c_v <- as.numeric(Gv %*% r_mean)
  list(Tlin = Tlin, cc = cc, c_u = c_u, G = G, Dxi = Dxi, Gv = Gv,
       c_v = c_v, ghss = ghss, ys = ys, hx = hx, hu = hu, Sigma_e = Sigma_e,
       n_s = n_s, n_u = n_u, n_endo = n_endo, d = d,
       ix1 = ix1, ix2 = ix2, ik = ik, endo = endo)
}

## Stationary (unconditional) order-2 output moments.
.order2_stationary_moments <- function(sys) {
  Sigma_x <- solve_lyapunov(sys$hx, sys$hu %*% sys$Sigma_e %*% t(sys$hu))
  Cr0 <- .order2_cov_r(numeric(sys$n_s), Sigma_x, sys$Sigma_e)
  Sxi <- solve_lyapunov(sys$Tlin, sys$G %*% Cr0 %*% t(sys$G))
  var_cov <- sys$Dxi %*% Sxi %*% t(sys$Dxi) + sys$Gv %*% Cr0 %*% t(sys$Gv)
  mu_xi <- if (sys$d == 0L) numeric(0) else
    as.numeric(solve(diag(sys$d) - sys$Tlin, sys$cc + sys$c_u))
  mean_dev <- as.numeric(sys$Dxi %*% mu_xi) + 0.5 * sys$ghss + sys$c_v
  list(var_cov = var_cov, mean = sys$ys + mean_dev, Sigma_x = Sigma_x,
       ## The FULL augmented-state covariance and the innovation
       ## covariance are needed for the exact lag-tau autocovariance; they were
       ## computed here already but thrown away, and the callers reconstructed a
       ## wrong recursion from Sigma_x / Var_x2 alone.  See .order2_autocov().
       Sxi = Sxi, Cr0 = Cr0,
       Var_x2 = Sxi[sys$ix2, sys$ix2, drop = FALSE], mean_x2 = mu_xi[sys$ix2],
       ## Cov(x2_t, x1_t (x) x1_t): exact order-2 cross block (n_s x n_s^2),
       ## x1(x)x1 column (d,e) with e fastest.  Needed by the order-3 pruned-SS
       ## Cr0 to build the j5(=eps(x)x2) cross-category blocks correctly (the
       ## connected non-Gaussian moment E[x2c x1 x1]); see .order3_cov_r.
       Cov_x2_x11 = Sxi[sys$ix2, sys$ik, drop = FALSE])
}

## EXACT lag-tau autocovariances of the AFVRR order-2 pruned state space.
##
## A13a. Both compute_moments_order2() and pruned_ss_moments() used
##   Gamma(tau) = ghx S_sel Gamma(tau-1),  seeded at
##   Gamma(0)_hat = ghx (Sigma_x + Var_x2) ghx'
## which is wrong twice over:
##   (i) the seed omits the contemporaneous shock term (ghu Sigma_e ghu' and the
##       order-2 Gv Cr0 Gv' blocks), so even a plain AR(1) state came out with
##       autocorrelation rho^(tau+1) instead of rho^tau -- e.g. rho = 0.85 was
##       reported as 0.614 at lag 1 against a 4e5-period MC value of 0.8486
##       (MCSE ~ 0.0016, i.e. ~140 MCSE out);
##   (ii) y_t is NOT a function of the TOTAL state x1+x2 alone -- its quadratic
##       term 0.5 ghxx (x1 (x) x1) sees x1 only -- so no S_sel recursion on
##       Gamma can be right for the pruned system.
##
## The pruned system IS linear in the augmented state
##   xi_t = [x1; x2; x1 (x) x1],   xi_{t+1} = Tlin xi_t + c + G r_t
##   y_t  = Dxi xi_t + const + Gv r_t
## and r_t is a martingale difference given xi_t (E[r_t | xi_t] = (0,0,0,vecSe)
## is constant because eps_t is independent of x1_t).  Hence for tau >= 1
##   Cov(xi_t, r_{t-tau}) = Tlin^{tau-1} G Cr0 ,  Cov(r_t, xi_{t-tau}) = 0
## and therefore
##   Gamma(tau) = Dxi Tlin^{tau-1} [ Tlin Sxi Dxi' + G Cr0 Gv' ].
## At tau = 0 the same algebra gives Dxi Sxi Dxi' + Gv Cr0 Gv' = var_cov, which
## is exactly what .order2_stationary_moments() already returns -- a free
## internal consistency check on the derivation.
##
## Returns an n_endo x n_endo x n_ar array of AUTOCOVARIANCES (not scaled).
## @noRd
.order2_autocov <- function(sys, st, n_ar) {
  n_endo <- sys$n_endo
  out <- array(0, dim = c(n_endo, n_endo, n_ar))
  if (n_ar < 1L) return(out)
  ## A = Tlin^{tau-1} [ Tlin Sxi Dxi' + G Cr0 Gv' ]   (d x n_endo)
  A <- sys$Tlin %*% st$Sxi %*% t(sys$Dxi) + sys$G %*% st$Cr0 %*% t(sys$Gv)
  for (lag in seq_len(n_ar)) {
    out[, , lag] <- sys$Dxi %*% A
    A <- sys$Tlin %*% A
  }
  out
}

## Transient s0-conditional order-2 output moments (levels), t = 1..n_periods.
## Returns mean (n_periods x n_endo) and cov (length-n_periods list).
.order2_conditional_moments <- function(sys, s0_state, n_periods) {
  d <- sys$d; ix1 <- sys$ix1; ik <- sys$ik
  mu <- numeric(d); mu[ix1] <- s0_state
  mu[ik] <- as.numeric(kronecker(s0_state, s0_state))
  Sx <- matrix(0, d, d)
  mean_mat <- matrix(0, n_periods, sys$n_endo, dimnames = list(NULL, sys$endo))
  cov_list <- vector("list", n_periods)
  for (t in seq_len(n_periods)) {
    a <- mu[ix1]; P <- Sx[ix1, ix1, drop = FALSE]
    Cr <- .order2_cov_r(a, P, sys$Sigma_e)
    mean_mat[t, ] <- as.numeric(sys$Dxi %*% mu) + 0.5 * sys$ghss + sys$c_v + sys$ys
    cov_list[[t]] <- sys$Dxi %*% Sx %*% t(sys$Dxi) + sys$Gv %*% Cr %*% t(sys$Gv)
    mu <- as.numeric(sys$Tlin %*% mu + sys$cc + sys$c_u)
    Sx <- sys$Tlin %*% Sx %*% t(sys$Tlin) + sys$G %*% Cr %*% t(sys$G)
  }
  list(mean = mean_mat, cov = cov_list)
}


## ============================================================================
## 4. MAIN stoch_simul INTERFACE
## ============================================================================

#' Run stoch_simul: the main entry point for perturbation solution
#'
#' @param model dynhr_mod or filename
#' @param compiled dynhr_compiled (optional)
#' @param ss Steady state (optional)
#' @param params Parameter values (optional)
#' @param order Perturbation order (only 1 supported)
#' @param irf Number of IRF periods (0 = no IRFs).  Also accepted as
#'   \code{n_periods} for consistency with \code{compute_irfs()}.  When both
#'   are supplied, \code{n_periods} takes precedence with a message.
#' @param n_periods Alias for \code{irf} (preferred name; mirrors
#'   \code{compute_irfs(n_periods=)}).
#' @param periods Simulation periods (0 = no simulation)
#' @param verbose Print Dynare-style output
#' @return List with dr, irfs, moments, simulation results
#' @export
stoch_simul <- function(model, compiled = NULL, ss = NULL, params = NULL,
                        order = 1L, irf = 40L, n_periods = NULL, periods = 0L,
                        verbose = TRUE) {
  ## n_periods= alias for irf= (L4 fix: stoch_simul should accept n_periods
  ## for consistency with compute_irfs(); when supplied it overrides irf=).
  if (!is.null(n_periods)) {
    if (!missing(irf) && irf != 40L) {
      .dynhr_inform("stoch_simul: both 'irf' and 'n_periods' supplied; using 'n_periods'.")
    }
    irf <- n_periods
  }
  ## Parse if filename
  if (is.character(model) && length(model) == 1 && file.exists(model)) {
    model <- parse_mod(model)
  }

  if (is.null(params)) params <- model$param_values
  if (is.null(compiled)) {
    if (verbose) .dynhr_cat("Compiling model...\n")
    compiled <- compile_model(model, verbose = verbose)
  }

  ## Steady state
  if (is.null(ss)) {
    if (verbose) .dynhr_cat("Computing steady state...\n")
    ss_result <- solve_steady_state(model, compiled, params, verbose = verbose)
    ss <- ss_result$ss
  } else {
    ## L31: catch the positional-arg trap stoch_simul(model, compiled, params, .)
    ## where `params` lands in the `ss` slot — a valid steady state is named over
    ## endogenous variables; if `ss` instead matches parameter names the result
    ## would be a wrong Jacobian + spurious BK violation with no other signal.
    ss_nm <- names(ss)
    if (!is.null(ss_nm) &&
        length(intersect(ss_nm, model$var_names)) == 0L &&
        length(intersect(ss_nm, model$param_names)) > 0L) {
      stop("stoch_simul: `ss` names no endogenous variable but matches parameter ",
           "names -- you likely called stoch_simul(model, compiled, params, ...) ",
           "positionally. Pass ss = ss$values and params = p by name.",
           call. = FALSE)
    }
  }

  ## First-order perturbation
  if (order != 1L) stop("Only order=1 perturbation is supported")
  if (verbose) .dynhr_cat("Solving first-order perturbation...\n")
  dr <- solve_perturbation(model, compiled, ss, params, verbose = verbose)

  ## Print eigenvalues
  if (verbose) {
    .dynhr_cat("\nEigenvalues:\n")
    eig <- dr$eigenvalues
    eig_mod <- Mod(eig)
    for (i in seq_along(eig)) {
      flag <- if (eig_mod[i] < 1) "stable" else "UNSTABLE"
      .dynhr_cat(sprintf("  %3d: %8.4f + %8.4fi (|lambda| = %7.4f) %s\n",
                  i, Re(eig[i]), Im(eig[i]), eig_mod[i], flag))
    }
    .dynhr_cat(sprintf("\n%d stable, %d unstable, %d forward-looking\n",
                dr$n_stable, dr$n_unstable, length(dr$fwd_vars)))
    if (dr$bk_satisfied) {
      .dynhr_cat("Blanchard-Kahn conditions are satisfied.\n\n")
    } else {
      .dynhr_cat("WARNING: Blanchard-Kahn conditions NOT satisfied!\n\n")
    }
  }

  ## IRFs
  irfs <- NULL
  if (irf > 0) {
    if (verbose) .dynhr_cat("Computing IRFs (", irf, " periods)...\n")
    irfs <- compute_irfs(dr, model, n_periods = irf, params = params)
  }

  ## Moments
  if (verbose) .dynhr_cat("Computing theoretical moments...\n")
  moments <- compute_moments(dr, model, params = params)

  ## Progress output inside a run, not a user-invoked render: emit it on the
  ## levelled message stream so `dynhr_set_verbosity()` can silence it. One
  ## call, so the tables cannot split across messages.
  if (verbose) .dynhr_cat(format_moments(moments, model), sep = "\n")

  ## Simulation
  sim <- NULL
  if (periods > 0) {
    if (verbose) .dynhr_cat("Simulating ", periods, " periods...\n")
    sim <- simulate_model(dr, n_periods = periods, model = model)
  }

  result <- list(
    model = model,
    compiled = compiled,
    ss = ss,
    params = params,
    dr = dr,
    irfs = irfs,
    moments = moments,
    sim = sim
  )
  class(result) <- "StochSimulResult"
  invisible(result)
}

## ============================================================================
## 5. DISPLAY FUNCTIONS
## ============================================================================

#' Format theoretical moments as Dynare-style report lines
#'
#' Pure: builds and returns the lines, writes nothing. See R/format-report.R
#' for why the reports are split this way -- in short, a report that emits ONCE
#' works on any stream, while one that emits per table cell is bound to
#' \code{cat()} forever.
#'
#' @param moments A moments list from [compute_moments()].
#' @param model Optional model; unused, kept for signature compatibility.
#' @return Character vector of lines.
#' @noRd
format_moments <- function(moments, model = NULL) {
  endo <- names(moments$std_dev)
  rule <- .fmt_rule(70L)

  tbl1 <- .fmt_cols_lines(endo,
                          Mean = rep(0, length(endo)),
                          `Std. Dev.` = as.numeric(moments$std_dev[endo]),
                          label_width = 20L, col_width = 12L, digits = 6L)
  out <- c("THEORETICAL MOMENTS", rule, tbl1[1L], rule, tbl1[-1L], rule, "")

  out <- c(out, "CORRELATION MATRIX", rule,
           .fmt_matrix_lines(moments$correlation,
                             row_labels = endo, col_labels = endo,
                             row_width = 10L, col_width = 9L, digits = 4L,
                             max_label = 8L),
           "")

  n_ar <- dim(moments$autocorr)[3L]
  ac   <- vapply(seq_len(n_ar),
                 function(lag) moments$autocorr[cbind(seq_along(endo),
                                                      seq_along(endo), lag)],
                 numeric(length(endo)))
  ac <- matrix(ac, nrow = length(endo), ncol = n_ar)
  out <- c(out, "AUTOCORRELATION (diagonal)", rule,
           .fmt_matrix_lines(ac, row_labels = endo,
                             col_labels = paste0("lag", seq_len(n_ar)),
                             row_width = 15L, col_width = 10L, digits = 4L,
                             corner = "Variable"),
           "")

  exo <- colnames(moments$var_decomp_pct)
  c(out, "VARIANCE DECOMPOSITION (percent)", rule,
    .fmt_matrix_lines(moments$var_decomp_pct, row_labels = endo,
                      col_labels = exo, row_width = 15L, col_width = 10L,
                      digits = 2L, max_label = 9L, corner = "Variable"),
    rule, "")
}

#' Print theoretical moments in Dynare-style format
#'
#' One emit call, so the whole report is a single write on whichever stream it
#' goes to.
#'
#' @inheritParams format_moments
#' @return The formatted lines, invisibly.
#' @noRd
print_moments <- function(moments, model = NULL) {
  lines <- format_moments(moments, model)
  cat(paste(lines, collapse = "\n"), "\n", sep = "")
  invisible(lines)
}


#' Print first-order decision rules
#'
#' @param x A \code{DecisionRules} object from [solve_perturbation()].
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.DecisionRules <- function(x, ...) {
  cat("=== Decision Rules (first order) ===\n")
  cat("State variables:", paste(x$state_vars, collapse = ", "), "\n")
  cat("Eigenvalues:", x$n_stable, "stable,", x$n_unstable, "unstable\n")
  cat("Blanchard-Kahn:", if (x$bk_satisfied) "satisfied" else "VIOLATED", "\n\n")

  cat("ghx (state transition): ", nrow(x$ghx), "x", ncol(x$ghx), "\n")
  if (ncol(x$ghx) <= 10 && nrow(x$ghx) <= 20) {
    print(round(x$ghx, 6))
  }
  cat("\nghu (shock impact): ", nrow(x$ghu), "x", ncol(x$ghu), "\n")
  if (ncol(x$ghu) <= 10 && nrow(x$ghu) <= 20) {
    print(round(x$ghu, 6))
  }
  cat("\n")
  invisible(x)
}


#' Print a stoch_simul result
#'
#' @param x A \code{StochSimulResult} object from [stoch_simul()], [compute_irfs()] or [compute_moments()].
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.StochSimulResult <- function(x, ...) {
  cat("=== stoch_simul Results ===\n")
  cat("Endogenous vars:", length(x$dr$endo_names), "\n")
  cat("Exogenous vars: ", length(x$dr$exo_names), "\n")
  cat("State variables:", x$dr$n_state, "\n")
  cat("Eigenvalues: ", x$dr$n_stable, "stable,", x$dr$n_unstable, "unstable\n")
  cat("BK conditions: ", if (x$dr$bk_satisfied) "satisfied" else "VIOLATED", "\n")
  if (!is.null(x$irfs)) {
    cat("IRFs computed:  yes (", attr(x$irfs, "n_periods"), " periods)\n")
  }
  if (!is.null(x$sim)) {
    cat("Simulation:     yes (", nrow(x$sim), " periods)\n")
  }
  cat("===========================\n")
  invisible(x)
}
