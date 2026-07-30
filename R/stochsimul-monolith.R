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
.dynhr_expr_cache <- new.env(parent = emptyenv())

.eval_cached_expr <- function(expr_str, envir) {
  ex <- .dynhr_expr_cache[[expr_str]]
  if (is.null(ex)) {
    ex <- parse(text = expr_str)
    .dynhr_expr_cache[[expr_str]] <- ex
  }
  eval(ex, envir = envir)
}

## Memoized shock-name -> core stripping (eps_/e_ prefix and trailing _).
## Depends only on the shock name, never on parameter values, so it is safe to
## cache across MCMC draws; replaces three sub() regex calls per shock per call.
.dynhr_shock_core_cache <- new.env(parent = emptyenv())

.shock_core <- function(nm) {
  core <- .dynhr_shock_core_cache[[nm]]
  if (is.null(core)) {
    core <- sub("^eps_", "", nm)
    core <- sub("^e_",   "", core)
    core <- sub("_$",    "", core)   # strip trailing "_" (e.g. eps_pref_ -> "pref")
    .dynhr_shock_core_cache[[nm]] <- core
  }
  core
}

.get_shock_stderr <- function(model, exo_names, params = NULL) {
  stderr <- setNames(rep(0, length(exo_names)), exo_names)
  if (is.null(params)) params <- model$param_values

  ## Eval environment from parameters, built lazily: models that resolve every
  ## shock via the param-name path (Priority 2) never need it, so we skip the
  ## per-call list2env entirely for them.
  penv <- NULL
  get_penv <- function() {
    if (is.null(penv)) penv <<- list2env(as.list(params), parent = baseenv())
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
                   error = function(e) NA_real_)
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

    ## Priority 2: parameter named sig_<shock> or stderr_<shock>
    if (is.na(se) || !is.finite(se)) {
      core <- .shock_core(nm)             # memoized regex strip (nm-only, value-free)
      for (prefix in c("sig_", "stderr_", "sigma_")) {
        pnm <- paste0(prefix, core)
        if (pnm %in% names(params) && is.finite(params[[pnm]])) {
          se <- params[[pnm]]
          break
        }
      }
    }

    ## Default: 0 (shock has no variance)
    if (is.na(se) || !is.finite(se)) se <- 0
    stderr[nm] <- se
  }
  stderr
}

## Helper: extract shock skewness shape parameters (alpha) from model$shocks
##
## The skew column in model$shocks$variances holds a parse-time numeric
## snapshot; skew_expr holds the raw expression text so the alpha can be
## re-evaluated against the current parameter vector theta during MCMC (mirrors
## the stderr_expr / .get_shock_stderr pattern above).  alpha = 0 (default)
## gives a symmetric Gaussian shock; alpha != 0 gives a skew-normal with
## E[eps_i] = sigma_i * delta_i * sqrt(2/pi), delta_i = alpha_i/sqrt(1+alpha_i^2).
## alpha can be NEGATIVE (brief Landmine 7).
.get_shock_skewness <- function(model, exo_names, params = NULL) {
  alpha <- setNames(rep(0, length(exo_names)), exo_names)
  if (is.null(params)) params <- model$param_values

  penv <- NULL
  get_penv <- function() {
    if (is.null(penv)) penv <<- list2env(as.list(params), parent = baseenv())
    penv
  }

  sv <- if (!is.null(model$shocks) && is.data.frame(model$shocks$variances))
    model$shocks$variances else NULL

  if (is.null(sv)) return(alpha)
  if (!("skew" %in% names(sv))) return(alpha)

  for (i in seq_along(exo_names)) {
    nm  <- exo_names[i]
    idx <- which(sv$name == nm)
    if (length(idx) == 0L) next
    row <- sv[idx[1], ]

    al <- NA_real_
    if ("skew_expr" %in% names(row) && !is.na(row$skew_expr) &&
        nzchar(row$skew_expr)) {
      al <- tryCatch(.eval_cached_expr(as.character(row$skew_expr), get_penv()),
                     error = function(e) NA_real_)
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

  ## FIX: respect caller-supplied params (was unconditionally overwritten)
  if (is.null(params)) params <- model$param_values

  ## Build full shock covariance and its lower Cholesky factor.
  ## Dynare convention: IRF for shock k uses the k-th column of chol(Sigma_e)
  ## (lower triangular), so correlated shocks propagate to all impact responses.
  ## For a diagonal Sigma_e this reduces to ghu[,k]*sigma_k (backward-compat).
  ## M28: honour a shock covariance carried on the decision rules (set when the
  ## caller passed solve_perturbation(Sigma_e=)) -- including off-diagonal
  ## correlations -- in preference to re-deriving a covariance from model$shocks.
  ## Falls back to the parsed shocks block when dr$Sigma_e is absent (the common
  ## case, byte-identical to before).
  Sigma_e_irf <- if (!is.null(dr$Sigma_e) &&
                     all(dim(dr$Sigma_e) == c(n_exo, n_exo))) {
    dr$Sigma_e
  } else {
    .get_shock_cov(model, exo, params)
  }
  L_chol <- tryCatch(
    t(chol(Sigma_e_irf)),            # lower triangular: L %*% t(L) = Sigma_e
    error = function(e) {
      ## Fallback: diagonal (square-root of diagonal); preserves old behaviour
      ## for degenerate / near-singular Sigma_e (e.g. shocks with zero stderr)
      diag(sqrt(pmax(diag(Sigma_e_irf), 0)), nrow = n_exo)
    }
  )

  irfs <- list()
  for (k in seq_along(exo)) {
    shock_name <- exo[k]
    irf_mat <- matrix(0, nrow = n_periods, ncol = n_endo)
    colnames(irf_mat) <- endo
    rownames(irf_mat) <- paste0("t", seq_len(n_periods))

    ## Shock vector: k-th column of lower Cholesky of Sigma_e, scaled by
    ## shock_size.  Reduces to shock_stderr[k]*shock_size when Sigma_e diagonal.
    eps <- L_chol[, k] * shock_size

    ## M26: a shock with no declared variance (all-zero Cholesky column) would
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
#' @param shocks Matrix of shocks (n_periods x n_exo). If NULL, draws random.
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
#' @return Matrix (n_periods x n_endo) when \code{n_replications = 1};
#'   3-D array (n_periods x n_endo x n_replications) when \code{> 1}.
#' @export
simulate_model <- function(dr, n_periods = 200L, shocks = NULL,
                           n_replications = 1L, model = NULL,
                           burn_in = 100L, init_state = NULL) {
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
                                   init_state = init_state)
    }
    return(out)
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

    if (any(alpha != 0)) {
      ## CSN draw: e_i = sigma_i * (delta_i*|z1| + sqrt(1-delta_i^2)*z2) - mu_i
      ## where mu_i = sigma_i * delta_i * sqrt(2/pi)  (Landmine 6: subtract mean
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
    } else {
      ## Gaussian path (unchanged; regression-safe at alpha = 0)
      shocks <- matrix(rnorm(total_periods * n_exo), ncol = n_exo)
      for (k in seq_along(exo)) shocks[, k] <- shocks[, k] * sigma_e[k]
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

#' Solve the discrete Lyapunov equation X = A X A' + B
#'
#' Uses the doubling algorithm for efficiency.
#'
#' @param A Square matrix
#' @param B Symmetric positive semi-definite matrix
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' @return Solution matrix X
#' @export
solve_lyapunov <- function(A, B, max_iter = 500L, tol = 1e-14) {
  n <- nrow(A)

  ## Fast path: if B is all zeros, solution is zero
  if (all(B == 0)) return(matrix(0, n, n))

  ## Doubling algorithm
  X <- B
  A_pow <- A
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    X_new <- X + A_pow %*% X %*% t(A_pow)
    if (any(!is.finite(X_new))) break
    diff <- max(abs(X_new - X))
    if (!is.finite(diff)) break
    ## RELATIVE convergence: a near-unit root gives a huge stationary covariance
    ## (entries ~ 1/(1-rho^2)), so the per-step increment can never fall below an
    ## ABSOLUTE 1e-14 (it plateaus at ~max|X| * machine-eps). An absolute test
    ## therefore never converges for near-unit-root systems -> the loop runs all
    ## max_iter steps and falls through to the O(n^6) kronecker solve (~0.5 s for
    ## n = 37). Scaling by max|X| makes it converge in the proper ~log2(mixing)
    ## steps for ANY stable A. (Harmless for well-damped A where max|X| ~ O(1).)
    if (diff < tol * max(1, max(abs(X_new)))) { converged <- TRUE; break }
    A_pow <- A_pow %*% A_pow
    if (any(!is.finite(A_pow))) break
    X <- X_new
  }
  if (converged) return(X)

  ## Stability gate before the O(n^6) vec/kronecker fallback.
  ##
  ## The discrete Lyapunov equation X = A X A' + B has a finite (PSD) solution
  ## only when A is stable (spectral radius < 1). The doubling loop above only
  ## fails to converge when A has a unit/explosive root: A^k does not decay, so
  ## the stationary covariance diverges and NO valid X exists. In that case the
  ## kronecker fallback is doubly bad -- it spends O(n^6) solving an (n^2 x n^2)
  ## system (e.g. ~0.5 s for n = 37) AND, because M = I - A (x) A is only
  ## *near*-singular for an explosive root (rcond just above machine-eps, so the
  ## guard below misses it), it returns a garbage non-PSD X instead of NaN.
  ##
  ## A cheap eigenvalue check (~0.2 ms for n = 37) short-circuits this: signal
  ## non-stationarity with NaN (the same contract as the singular-M guard
  ## below), letting the caller fall back to the exact-diffuse Kalman init /
  ## simulation-based moments. This is the dominant hot path when an estimation
  ## sampler probes BK-boundary draws whose decision rule is near-explosive
  ## (verified on NZSIM: ~473 ms/eval of wasted kronecker solves).
  if (max(Mod(eigen(A, only.values = TRUE)$values)) >= 1)
    return(matrix(NaN, n, n))

  ## Fallback: vec method  vec(X) = (I - A (x) A)^{-1} vec(B)
  I_n2 <- diag(n^2)
  AkA <- kronecker(A, A)
  M <- I_n2 - AkA
  # Check for singular system (unit roots from NN1 placeholders, etc.)
  if (rcond(M) < .Machine$double.eps) {
    # Return a matrix with NaN to signal non-stationarity; caller can
    # fall back to simulation-based welfare.
    # Unit root detected (e.g. NN1 placeholder equations). Not an error;
    # the caller can fall back to simulation-based computations.
    if (getOption("dynhr.warn_lyapunov", FALSE)) {
      warning("solve_lyapunov: system is singular (unit root detected). Returning NaN.")
    }
    return(matrix(NaN, n, n))
  }
  x_vec <- solve(M, as.vector(B))
  matrix(x_vec, nrow = n, ncol = n)
}

## ---------------------------------------------------------------------------
## M2: stationary-subspace projection via the modal (eigen) decomposition.
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
  ## M28: honour a covariance carried on the decision rules (set when the caller
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

  ## ---- M2: stationary-subspace projection ------------------------------------
  ## Detect unit/explosive roots via the EIGENVALUES of ghx_state (exact and
  ## basis-independent).  When all roots are stable, solve the Lyapunov
  ## equation as usual.  When a unit/explosive root is present, project onto
  ## the stable invariant subspace (modal decomposition) so stationary
  ## variables keep finite moments and only nonstationary ones become NaN --
  ## matching Dynare.  See .modal_decomp / .modal_project_cov above.
  nonstat_tol <- getOption("dynhr.unit_root_tol", 1e-6)
  ev <- eigen(ghx_state, only.values = TRUE)$values
  has_unit_root <- any(Mod(ev) >= 1 - nonstat_tol)

  Q_state <- ghu_state %*% Sigma_e %*% t(ghu_state)   # innovation cov of states
  dec <- NULL
  nonstat_var <- integer(0)

  if (!has_unit_root) {
    Sigma_state <- solve_lyapunov(ghx_state, Q_state)
  } else {
    if (getOption("dynhr.warn_unit_root", TRUE)) {
      warning(sprintf(
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
  ##                modal-eigenvalue criterion used for var_cov / M2).
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
    ## Seed: stationary block of Sigma_y (finite by M2)
    Gamma_ss <- Sigma_y[stat_var, stat_var, drop = FALSE]
    for (lag in seq_len(n_ar)) {
      Gamma_ss    <- G_ss %*% Gamma_ss
      autocorr[stat_var, stat_var, lag] <- Gamma_ss / sd_ss
    }
  }

  ## Variance decomposition: contribution of each shock to each variable.
  ## Mirrors the stationary-subspace logic above for per-shock Lyapunov solves.
  var_decomp <- matrix(0, nrow = n_endo, ncol = n_exo)
  rownames(var_decomp) <- endo
  colnames(var_decomp) <- exo

  for (k in seq_along(exo)) {
    ghu_state_k <- ghu_state[, k, drop = FALSE]
    ghu_k       <- ghu[, k, drop = FALSE]
    Q_state_k   <- ghu_state_k %*% (Sigma_e[k, k]) %*% t(ghu_state_k)

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
    Sigma_y_k <- ghx %*% Sigma_state_k %*% t(ghx) +
                 ghu_k %*% (Sigma_e[k, k]) %*% t(ghu_k)
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
#' The contribution of shock \eqn{k} to variable \eqn{i} at horizon \eqn{h} is
#' \eqn{\sum_{j=0}^{h-1} (\Psi_j)_{ik}^2 \Sigma_e[k,k]}.
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

  ## Pre-extract diagonal of Sigma_e (per-shock variances)
  sigma_k2 <- diag(Sigma_e)   # length n_exo

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
    Psi_cur <- ghu             # n_endo x n_exo  (Psi_0)
    phi     <- ghu_state       # n_state x n_exo (phi_0 = B)

    ## Running FEV accumulator
    FEV <- matrix(0, n_endo, n_exo)
    for (k in seq_len(n_exo)) {
      FEV[, k] <- Psi_cur[, k]^2 * sigma_k2[k]
    }

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
      for (k in seq_len(n_exo)) {
        FEV[, k] <- FEV[, k] + Psi_cur[, k]^2 * sigma_k2[k]
      }

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
    for (k in seq_along(exo)) {
      ghu_state_k <- ghu_state[, k, drop = FALSE]
      ghu_k       <- ghu[, k, drop = FALSE]
      Q_state_k   <- ghu_state_k %*% (Sigma_e[k, k]) %*% t(ghu_state_k)

      if (!has_unit_root) {
        Sigma_state_k <- solve_lyapunov(ghx_state, Q_state_k)
      } else if (isTRUE(dec$ok)) {
        Sigma_state_k <- .modal_project_cov(dec, Q_state_k)
      } else {
        Sigma_state_k <- matrix(0, n_state, n_state)
      }

      Sigma_y_k <- ghx %*% Sigma_state_k %*% t(ghx) +
                   ghu_k %*% (Sigma_e[k, k]) %*% t(ghu_k)
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
    warning(
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
  ## For Gaussian pruned system the lag-τ autocovariance is:
  ##   Γ(τ) = ghx · hx^(τ-1) · (Σ_x + Var(x²)) · ghx'    for τ ≥ 1
  ## (quadratic innovations are independent of lagged y under Gaussianity at lag≥1)
  ## τ = 0 uses the full Sigma_y computed above.
  autocorr <- array(0, dim = c(n_endo, n_endo, n_ar))
  dimnames(autocorr) <- list(endo, endo, paste0("lag", seq_len(n_ar)))

  ## Selection matrix S_sel maps states back to full endo space
  S_sel <- matrix(0, nrow = n_s, ncol = n_endo)
  for (i in seq_along(state_idx)) S_sel[i, state_idx[i]] <- 1

  ## Propagation recurrence: Γ(τ) = ghx · S_sel · Γ(τ-1)
  ## initialised from Γ(0)_x1 = ghx*Sigma_x*ghx', Γ(0)_x2 = ghx*Var_x2*ghx'.
  Gamma_prev_x1 <- ghx %*% Sigma_x %*% t(ghx)
  Gamma_prev_x2 <- ghx %*% Var_x2  %*% t(ghx)

  for (lag in seq_len(n_ar)) {
    Gamma_lag_x1 <- ghx %*% S_sel %*% Gamma_prev_x1
    Gamma_lag_x2 <- ghx %*% S_sel %*% Gamma_prev_x2
    Gamma_lag <- Gamma_lag_x1 + Gamma_lag_x2
    autocorr[, , lag] <- Gamma_lag / sd_outer
    Gamma_prev_x1 <- Gamma_lag_x1
    Gamma_prev_x2 <- Gamma_lag_x2
  }

  ## --- Variance decomposition (exact per-shock via the augmented system) ----
  ## Shock k's contribution = order-2 Var(y) with Sigma_e zeroed except entry
  ## (k,k).  Same machinery as the total, so per-shock attribution is internally
  ## consistent (sums to the total up to genuine cross-shock interaction terms).
  var_decomp <- matrix(0, nrow = n_endo, ncol = n_exo)
  rownames(var_decomp) <- endo
  colnames(var_decomp) <- exo
  for (k in seq_len(n_exo)) {
    Sigma_ek <- matrix(0, n_exo, n_exo)
    Sigma_ek[k, k] <- Sigma_e[k, k]
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
## (transient, s0-conditional, Tier 15 B4).
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
  i1 <- seq_len(d1); i2 <- (d1 + 1L):(d1 + d2)
  i3 <- (d1 + d2 + 1L):(d1 + d2 + d3); i4 <- (d1 + d2 + d3 + 1L):D
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
  ix1 <- seq_len(n_s); ix2 <- (n_s + 1L):(2L * n_s); ik <- (2L * n_s + 1L):d
  Tlin <- matrix(0, d, d)
  Tlin[ix1, ix1] <- hx; Tlin[ix2, ix2] <- hx
  Tlin[ix2, ik] <- 0.5 * hxx; Tlin[ik, ik] <- kronecker(hx, hx)
  cc <- numeric(d); cc[ix2] <- 0.5 * hss
  c_u <- numeric(d)
  c_u[ix2] <- 0.5 * as.numeric(huu %*% vecSe)
  c_u[ik]  <- as.numeric(kronecker(hu, hu) %*% vecSe)
  D1 <- n_u; D2 <- n_u * n_s; D3 <- n_s * n_u; D4 <- n_u * n_u
  Dr <- D1 + D2 + D3 + D4
  j1 <- seq_len(D1); j2 <- (D1 + 1L):(D1 + D2)
  j3 <- (D1 + D2 + 1L):(D1 + D2 + D3); j4 <- (D1 + D2 + D3 + 1L):Dr
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
  mu_xi <- as.numeric(solve(diag(sys$d) - sys$Tlin, sys$cc + sys$c_u))
  mean_dev <- as.numeric(sys$Dxi %*% mu_xi) + 0.5 * sys$ghss + sys$c_v
  list(var_cov = var_cov, mean = sys$ys + mean_dev, Sigma_x = Sigma_x,
       Var_x2 = Sxi[sys$ix2, sys$ix2, drop = FALSE], mean_x2 = mu_xi[sys$ix2],
       ## Cov(x2_t, x1_t (x) x1_t): exact order-2 cross block (n_s x n_s^2),
       ## x1(x)x1 column (d,e) with e fastest.  Needed by the order-3 pruned-SS
       ## Cr0 to build the j5(=eps(x)x2) cross-category blocks correctly (the
       ## connected non-Gaussian moment E[x2c x1 x1]); see .order3_cov_r.
       Cov_x2_x11 = Sxi[sys$ix2, sys$ik, drop = FALSE])
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
      message("stoch_simul: both 'irf' and 'n_periods' supplied; using 'n_periods'.")
    }
    irf <- n_periods
  }
  ## Parse if filename
  if (is.character(model) && length(model) == 1 && file.exists(model)) {
    model <- parse_mod(model)
  }

  if (is.null(params)) params <- model$param_values
  if (is.null(compiled)) {
    if (verbose) cat("Compiling model...\n")
    compiled <- compile_model(model, verbose = verbose)
  }

  ## Steady state
  if (is.null(ss)) {
    if (verbose) cat("Computing steady state...\n")
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
  if (verbose) cat("Solving first-order perturbation...\n")
  dr <- solve_perturbation(model, compiled, ss, params, verbose = verbose)

  ## Print eigenvalues
  if (verbose) {
    cat("\nEigenvalues:\n")
    eig <- dr$eigenvalues
    eig_mod <- Mod(eig)
    for (i in seq_along(eig)) {
      flag <- if (eig_mod[i] < 1) "stable" else "UNSTABLE"
      cat(sprintf("  %3d: %8.4f + %8.4fi (|lambda| = %7.4f) %s\n",
                  i, Re(eig[i]), Im(eig[i]), eig_mod[i], flag))
    }
    cat(sprintf("\n%d stable, %d unstable, %d forward-looking\n",
                dr$n_stable, dr$n_unstable, length(dr$fwd_vars)))
    if (dr$bk_satisfied) {
      cat("Blanchard-Kahn conditions are satisfied.\n\n")
    } else {
      cat("WARNING: Blanchard-Kahn conditions NOT satisfied!\n\n")
    }
  }

  ## IRFs
  irfs <- NULL
  if (irf > 0) {
    if (verbose) cat("Computing IRFs (", irf, " periods)...\n")
    irfs <- compute_irfs(dr, model, n_periods = irf, params = params)
  }

  ## Moments
  if (verbose) cat("Computing theoretical moments...\n")
  moments <- compute_moments(dr, model, params = params)

  if (verbose) {
    print_moments(moments, model)
  }

  ## Simulation
  sim <- NULL
  if (periods > 0) {
    if (verbose) cat("Simulating ", periods, " periods...\n")
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

#' Print theoretical moments in Dynare-style format
#' @noRd
print_moments <- function(moments, model = NULL) {
  endo <- names(moments$std_dev)

  cat("THEORETICAL MOMENTS\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("%-20s %12s %12s\n", "Variable", "Mean", "Std. Dev."))
  cat(strrep("-", 70), "\n")
  for (nm in endo) {
    cat(sprintf("%-20s %12.6f %12.6f\n", nm, 0, moments$std_dev[nm]))
  }
  cat(strrep("-", 70), "\n\n")

  ## Correlation matrix
  cat("CORRELATION MATRIX\n")
  cat(strrep("-", 70), "\n")
  n <- length(endo)
  short <- substr(endo, 1, 8)
  cat(sprintf("%-10s", ""))
  for (nm in short) cat(sprintf("%9s", nm))
  cat("\n")
  for (i in seq_along(endo)) {
    cat(sprintf("%-10s", short[i]))
    for (j in seq_along(endo)) {
      cat(sprintf("%9.4f", moments$correlation[i, j]))
    }
    cat("\n")
  }
  cat("\n")

  ## Autocorrelation
  n_ar <- dim(moments$autocorr)[3]
  cat("AUTOCORRELATION (diagonal)\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("%-15s", "Variable"))
  for (lag in seq_len(n_ar)) cat(sprintf("%10s", paste0("lag", lag)))
  cat("\n")
  for (i in seq_along(endo)) {
    cat(sprintf("%-15s", endo[i]))
    for (lag in seq_len(n_ar)) {
      cat(sprintf("%10.4f", moments$autocorr[i, i, lag]))
    }
    cat("\n")
  }
  cat("\n")

  ## Variance decomposition
  exo <- colnames(moments$var_decomp_pct)
  cat("VARIANCE DECOMPOSITION (percent)\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("%-15s", "Variable"))
  for (nm in exo) cat(sprintf("%10s", substr(nm, 1, 9)))
  cat("\n")
  for (i in seq_along(endo)) {
    cat(sprintf("%-15s", endo[i]))
    for (j in seq_along(exo)) {
      cat(sprintf("%10.2f", moments$var_decomp_pct[i, j]))
    }
    cat("\n")
  }
  cat(strrep("-", 70), "\n\n")
}

#' Print decision rules
#' @noRd
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

#' Plot impulse response functions
#'
#' @param irfs IRFCollection from compute_irfs
#' @param vars Character vector of variables to plot (NULL = all)
#' @param shocks Character vector of shocks to plot (NULL = all)
#' @param ncol Number of columns in plot grid
#' @noRd
plot_irfs <- function(irfs, vars = NULL, shocks = NULL, ncol = 3L) {
  endo <- attr(irfs, "endo_names")
  exo  <- attr(irfs, "exo_names")
  n_periods <- attr(irfs, "n_periods")

  if (is.null(vars)) vars <- endo
  if (is.null(shocks)) shocks <- exo

  for (shock_name in shocks) {
    irf_mat <- irfs[[shock_name]]
    if (is.null(irf_mat)) next

    plot_vars <- intersect(vars, colnames(irf_mat))
    if (length(plot_vars) == 0) next

    nrow_plot <- ceiling(length(plot_vars) / ncol)
    par(mfrow = c(nrow_plot, ncol), mar = c(3, 3, 2, 1))
    for (v in plot_vars) {
      plot(seq_len(n_periods), irf_mat[, v], type = "l",
           main = paste(v, "<-", shock_name),
           xlab = "", ylab = "", col = "steelblue", lwd = 2)
      abline(h = 0, lty = 2, col = "grey60")
    }
  }
  par(mfrow = c(1, 1))
  invisible(irfs)
}

#' Print summary for StochSimulResult
#' @noRd
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
