## R/global-solution-methods.R
## --------------------------------------------------------------------------
## S3 methods for GlobalSolution objects returned by solve_global().
##
## Methods:
##   print.GlobalSolution
##   predict.GlobalSolution   -- evaluate policy at arbitrary state points
##   simulate.GlobalSolution  -- stochastic simulation using global policy
##   euler_errors             -- generic and default method
##   euler_errors.GlobalSolution
## --------------------------------------------------------------------------

#' Print a GlobalSolution object
#' @param x   A \code{GlobalSolution} object.
#' @param ... Ignored.
#' @export
print.GlobalSolution <- function(x, ...) {
  cat("GlobalSolution (Chebyshev projection)\n")
  cat(sprintf("  States:   %s\n", paste(x$state_names, collapse = ", ")))
  cat(sprintf("  Endo:     %s\n", paste(x$all_endo_names, collapse = ", ")))
  cat(sprintf("  Shocks:   %s\n", paste(x$shock_names, collapse = ", ")))
  cat(sprintf("  Degree:   %d  (n_basis = %d)\n",
              x$poly_degree, ncol(x$coefs)))
  cat(sprintf("  Converged: %s (in %d iterations, delta = %.2e)\n",
              x$converged, x$n_iter, x$last_delta))
  ## The domain-clipping diagnostic is printed only when it is non-zero:
  ## a clip-free solution is the normal case and does not need a line.
  if (isTRUE(is.finite(x$domain_clip_frac) && x$domain_clip_frac > 0))
    cat(sprintf(
      "  Domain:   %.2f%% of the Euler-quadrature t+1 state coordinates are CLIPPED over the ergodic +-3 sd box\n",
      100 * x$domain_clip_frac))
  invisible(x)
}

#' Evaluate the global policy function at given state (lag) points
#'
#' @param object  A \code{GlobalSolution} object.
#' @param newdata A named matrix or data.frame with columns for each state
#'   variable (lag values), or a named numeric vector for a single point.
#' @param ...     Ignored.
#' @return A matrix with one row per input point and one column per
#'   endogenous variable (in declaration order).
#' @export
predict.GlobalSolution <- function(object, newdata, ...) {
  state_names <- object$state_names
  n_state     <- length(state_names)

  ## Coerce newdata to a matrix
  if (is.vector(newdata) && !is.list(newdata)) {
    newdata <- matrix(newdata, nrow = 1L,
                      dimnames = list(NULL, names(newdata)))
  }
  newdata <- as.matrix(newdata)

  if (!all(state_names %in% colnames(newdata)))
    stop(sprintf(
      "predict.GlobalSolution: newdata must contain columns for states: %s",
      paste(state_names, collapse = ", ")))

  ## Normalize each state column
  n_pts <- nrow(newdata)
  x_norm <- matrix(0.0, nrow = n_pts, ncol = n_state)
  colnames(x_norm) <- state_names
  for (j in seq_len(n_state)) {
    nm  <- state_names[j]
    dom <- object$state_domain[[nm]]
    x_norm[, j] <- pmax(-1.0, pmin(1.0,
      cheb_normalize(newdata[, nm], dom[1], dom[2])))
  }

  ## Evaluate basis and multiply by coefficients
  Phi    <- cheb_basis(x_norm, object$poly_degree)
  result <- tcrossprod(Phi, object$coefs)   # n_pts x n_endo
  colnames(result) <- object$all_endo_names
  result
}

#' Stochastic simulation from the global policy function
#'
#' Period \eqn{t} draws \eqn{eps_t}, feeds it into the state lag through the
#' solver's own AR(1) injection (\code{compute_next_lag()}), and evaluates
#' the policy:
#' \deqn{y_t = policy(s_{t-1} + \psi \epsilon_t / \rho), \quad
#'       s_t = y_t[state\_names].}
#'
#' \strong{Timing (changed in F1-C).}  This is the timing
#' \code{make_log_posterior_global_pf()}'s particle filter uses, and the two
#' now agree period-for-period.  Previously the shock drawn in period
#' \eqn{t} was applied to the state lag entering period \eqn{t+1}, which
#' made the FIRST returned row deterministic (the shock-free policy value at
#' \code{init_state}) and dropped the last drawn shock entirely: the old
#' output's row \eqn{t+1} equals the new output's row \eqn{t}.  Any use that
#' relied on a zero-shock first period must now drop it explicitly.
#'
#' \code{init_state} defaults to the steady state -- a COLD start.  For a
#' draw from the model's stationary distribution, either discard a burn-in
#' or pass an \code{init_state} drawn from the stationary state covariance
#' (what the particle filter initialises its cloud with).
#'
#' @param object  A \code{GlobalSolution} object.
#' @param nsim    Number of periods to simulate (default 100).
#' @param seed    Random seed (default NULL = don't set).
#' @param init_state Named numeric; initial state lag values.
#'   Defaults to steady state.
#' @param ...     Ignored.
#' @return A matrix of size \code{nsim x n_endo} with simulated paths.
#' @export
simulate.GlobalSolution <- function(object, nsim = 100L, seed = NULL,
                                    init_state = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)

  state_names <- object$state_names
  endo        <- object$all_endo_names
  exo         <- object$shock_names
  n_state     <- length(state_names)
  n_endo      <- length(endo)
  n_exo       <- length(exo)
  shock_sds   <- object$shock_sds

  ## Initial state lag
  if (is.null(init_state)) {
    state_lag <- object$ss_vals[state_names]
  } else {
    state_lag <- init_state[state_names]
  }

  out <- matrix(NA_real_, nrow = nsim, ncol = n_endo)
  colnames(out) <- endo

  for (t in seq_len(nsim)) {
    ## Draw the CURRENT period's shock and inject it into the lag the policy
    ## is evaluated at -- the solver's own AR(1) feed. `compute_next_lag()`
    ## indexes its first argument by `state_names`, so a state-lag vector is
    ## exactly what it wants; the result is the feed point
    ## s_{t-1} + psi*eps_t/rho, whose policy value has z_t = rho*z_{t-1} +
    ## psi*eps_t. This is the recursion .global_pf_loglik() filters (F1-C
    ## timing audit); the previous ordering applied eps_t to period t+1 and
    ## so returned a shock-free first row.
    eps_t <- rnorm(n_exo, mean = 0, sd = shock_sds)
    names(eps_t) <- exo
    feed <- object$compute_next_lag(state_lag, eps_t)

    new_state_norm <- matrix(
      vapply(seq_len(n_state), function(j) {
        nm  <- state_names[j]
        dom <- object$state_domain[[nm]]
        max(-1.0, min(1.0, cheb_normalize(feed[nm], dom[1], dom[2])))
      }, numeric(1L)),
      nrow = 1L)

    Phi_t <- cheb_basis(new_state_norm, object$poly_degree)
    y_t   <- as.vector(Phi_t %*% t(object$coefs))
    names(y_t) <- endo
    out[t, ] <- y_t

    ## The time-t state (the lag entering t+1) is the policy output itself.
    state_lag <- y_t[state_names]
  }

  out
}

#' Euler equation errors generic
#'
#' @param x    An object with a global solution (e.g., \code{GlobalSolution}).
#' @param ...  Additional arguments.
#' @export
euler_errors <- function(x, ...) UseMethod("euler_errors")

## =========================================================================
## Model-agnostic accuracy battery (E1-B)
## -------------------------------------------------------------------------
## `euler_errors()` used to hard-code the stochastic-growth-with-full-
## depreciation Euler equation, keyed on the PARAMETER NAMES `alpha`/`beta`:
##
##   ee = E_t[beta * alpha * exp(z_{t+1}) * k_t^(alpha-1) / c_{t+1}] * c_t - 1
##
## On any model that happens to declare parameters called `alpha` and `beta`
## -- an NK model with a Calvo `alpha`, say -- that formula silently computed
## something that is not an accuracy measure at all, and on every other model
## it returned all-NA.  It is now replaced by residuals of the COMPILED
## model's own equations:
##
##   R_i(s) = sum_k w_k * F_i(y_lag = s, y_cur = g(s, 0),
##                            y_lead = g(g(s,0)[states], eps_k), eps_cur = 0)
##
## with (w_k, eps_k) Gauss-Hermite over the shocks (`gauss_hermite()`, the
## same rule solve_global() integrates with).  `g` is the policy: a
## `GlobalSolution`'s Chebyshev interpolant or a perturbation `DecisionRules`
## / `DecisionRules2` / `DecisionRules3` Taylor rule (evaluated UNPRUNED, the
## Aruoba-Fernandez-Villaverde-Rubio-Ramirez 2006 accuracy convention).
##
## NORMALISATION.  A raw residual carries the units of its equation, so the
## numbers are not comparable across equations or models.  Where the model
## has a natural numeraire variable `nu` (consumption by default) the
## residual is converted to CONSUMPTION-EQUIVALENT units -- the fractional
## change in nu_t that would zero the residual:
##
##   ee_i = R_i / (dE[F_i]/d nu_t * nu_t)
##
## For the RBC Euler equation dF_1/dc_t = -1/c_t^2, so ee_1 = -R_1 * c_t,
## which is EXACTLY (up to a sign that |.| discards, and up to floating-point
## association) the Maliar-Maliar unit-free error the old hard-wired formula
## produced -- see test-global-accuracy.R, which pins the pre-E1-B numbers.
## Equations in which nu_t does not appear (an AR(1) shock process, an
## accounting identity) have no such normalisation and are reported RAW; the
## `$normalized` flag says which is which.
## =========================================================================


## Column map + dy assembler for a compiled model.  ONLY +-1 timing is
## supported: the residual is evaluated from (y_{t-1}, y_t, y_{t+1}, eps_t),
## and a model carrying y(-2) or y(+2) must be run through the auxiliary-
## variable expansion first (which is what compile_model() does anyway).
#' @noRd
.acc_dy_map <- function(compiled) {
  dyn <- compiled$dynamic
  cmp <- dyn$dyn_col_map
  if (any(abs(cmp$lead_lag) > 1L))
    stop("euler_errors: the compiled model carries |lead_lag| > 1 columns; ",
         "the accuracy battery evaluates residuals from (t-1, t, t+1) only.",
         call. = FALSE)

  tot  <- dyn$total_cols
  endo <- dyn$endo_names
  exo  <- dyn$exo_names

  keys <- character(tot)
  for (k in seq_len(nrow(cmp))) {
    ll  <- cmp$lead_lag[k]
    sfx <- if (ll == 0L) "__0" else if (ll > 0L) paste0("__p", ll) else
             paste0("__m", abs(ll))
    keys[cmp$col[k]] <- paste0(cmp$name[k], sfx)
  }

  is_exo <- cmp$name %in% exo
  sel <- function(ll, exo_side) {
    take <- cmp$lead_lag == ll & (is_exo == exo_side)
    list(col = cmp$col[take], name = cmp$name[take])
  }
  m1 <- sel(-1L, FALSE); c0 <- sel(0L, FALSE); p1 <- sel(1L, FALSE)
  ux <- sel(0L, TRUE)

  list(
    keys     = keys,
    template = setNames(numeric(tot), keys),
    endo     = endo,
    exo      = exo,
    n_eq     = dyn$n_eq,
    m1 = m1, c0 = c0, p1 = p1, ux = ux,
    ## dy column of each endogenous variable at t (for the numeraire scale)
    col0_of  = setNames(c0$col, c0$name),
    ## equations that carry at least one t+1 term are the expectational ones
    lead_cols = p1$col
  )
}


## Assemble the dy vector the compiled residual/Jacobian functions expect.
#' @noRd
.acc_assemble_dy <- function(map, y_lag, y_cur, y_lead, eps_cur) {
  dy <- map$template
  if (length(map$m1$col)) dy[map$m1$col] <- y_lag[map$m1$name]
  if (length(map$c0$col)) dy[map$c0$col] <- y_cur[map$c0$name]
  if (length(map$p1$col)) dy[map$p1$col] <- y_lead[map$p1$name]
  if (length(map$ux$col)) dy[map$ux$col] <- eps_cur[map$ux$name]
  dy
}


## Tensor-product Gauss-Hermite rule over `n_exo` independent shocks with
## standard deviations `sds`.  Identical construction (and therefore identical
## node/weight ORDER) to solve_global()'s `expected_residual()`.
#' @noRd
.acc_quad <- function(n_quad, sds) {
  gh    <- gauss_hermite(n_quad)
  n_exo <- length(sds)
  if (n_exo == 1L)
    return(list(nodes = matrix(gh$nodes * sds[1L], ncol = 1L),
                weights = gh$weights))
  qn <- lapply(seq_len(n_exo), function(k) gh$nodes * sds[k])
  nodes <- as.matrix(expand.grid(qn))
  wm    <- as.matrix(expand.grid(replicate(n_exo, gh$weights,
                                           simplify = FALSE)))
  w     <- apply(wm, 1L, prod)
  list(nodes = nodes, weights = w / sum(w))
}


## ---- Policy adapters -----------------------------------------------------
## Both adapters expose the SAME two primitives, so the accuracy engine below
## never branches on the solution type:
##
##   y_at(state_lag, eps)  ->  named length-n_endo vector of LEVELS y_t
##   next_lag(y_t)         ->  named vector of the state values that are the
##                             lag at t+1  (= y_t[state_names] for both)
##
## The GlobalSolution's policy is a function of the state lag ALONE (its
## time-iteration convention sets eps_t = 0 and injects the shock through the
## AR(1) "feed" lag, see solve_global()'s compute_next_lag()); the adapter
## therefore routes `eps` through compute_next_lag() before interpolating.
#' @noRd
.acc_policy_global <- function(g) {
  sn <- g$state_names
  list(
    state_names  = sn,
    endo         = g$all_endo_names,
    exo          = g$shock_names,
    ss_vals      = g$ss_vals,
    shock_sds    = g$shock_sds[g$shock_names],
    state_domain = g$state_domain,
    y_at = function(state_lag, eps) {
      lag <- if (any(eps != 0)) g$compute_next_lag(state_lag, eps) else
               state_lag[sn]
      out <- predict(g, matrix(lag[sn], nrow = 1L,
                               dimnames = list(NULL, sn)))[1L, ]
      setNames(out, g$all_endo_names)
    },
    next_lag = function(y_t) y_t[sn]
  )
}


## Unpruned Taylor evaluation of a perturbation decision rule at (x, u), where
## x is the state deviation from the steady state and u the current shock.
## Kronecker conventions follow simulate_model_order2/3 EXACTLY (ghxu columns
## are state-FAST/exo-SLOW, hence `u %x% x`).
#' @noRd
.acc_dr_taylor <- function(dr, x, u, order) {
  y <- as.numeric(dr$ghx %*% x + dr$ghu %*% u)
  if (order >= 2L) {
    y <- y +
      0.5 * as.numeric(dr$ghxx %*% (x %x% x)) +
      as.numeric(dr$ghxu %*% (u %x% x)) +
      0.5 * as.numeric(dr$ghuu %*% (u %x% u)) +
      0.5 * dr$ghss
  }
  if (order >= 3L) {
    y <- y +
      (1 / 6) * as.numeric(dr$ghxxx %*% (x %x% x %x% x)) +
      0.5     * as.numeric(dr$ghxxu %*% (u %x% x %x% x)) +
      0.5     * as.numeric(dr$ghxuu %*% (u %x% u %x% x)) +
      (1 / 6) * as.numeric(dr$ghuuu %*% (u %x% u %x% u))
    if (!is.null(dr$ghxss)) y <- y + 0.5 * as.numeric(dr$ghxss %*% x)
    if (!is.null(dr$ghuss)) y <- y + 0.5 * as.numeric(dr$ghuss %*% u)
    if (!is.null(dr$ghs3))  y <- y + (1 / 6) * dr$ghs3
  }
  y + dr$ys
}


#' @noRd
.acc_policy_dr <- function(dr, model, params, shock_sds = NULL) {
  order <- if (inherits(dr, "DecisionRules3")) 3L else
           if (inherits(dr, "DecisionRules2")) 2L else 1L
  sn    <- dr$state_vars
  endo  <- dr$endo_names
  exo   <- dr$exo_names
  ss    <- dr$ys
  ss_s  <- ss[sn]
  if (is.null(shock_sds)) shock_sds <- .get_shock_sds(model, params)[exo]

  list(
    state_names  = sn,
    endo         = endo,
    exo          = exo,
    ss_vals      = ss,
    shock_sds    = shock_sds,
    state_domain = NULL,
    y_at = function(state_lag, eps) {
      x <- as.numeric(state_lag[sn] - ss_s)
      setNames(.acc_dr_taylor(dr, x, as.numeric(eps[exo]), order), endo)
    },
    next_lag = function(y_t) y_t[sn]
  )
}


## Default state grid for a policy that carries no explicit domain: the
## steady state +- `n_sd` unconditional standard deviations of each state,
## taken from the FIRST-ORDER stationary covariance (solve_lyapunov via
## kf_stationary_init).  A policy with an explicit `state_domain` (a
## GlobalSolution) uses that instead, so the perturbation and projection
## batteries can be pointed at the SAME box.
#' @noRd
.acc_state_domain_dr <- function(dr, model, params, n_sd = 2) {
  sn <- dr$state_vars
  s  <- .dr_state_sd(dr, model, params)
  sd_s <- if (is.null(s)) rep(0.1, length(sn)) else unname(s[sn])
  dom <- vector("list", length(sn)); names(dom) <- sn
  for (j in seq_along(sn)) {
    h <- n_sd * sd_s[j]
    if (!is.finite(h) || h <= 0) h <- max(1e-3 * abs(dr$ys[sn[j]]), 1e-3)
    dom[[sn[j]]] <- c(dr$ys[sn[j]] - h, dr$ys[sn[j]] + h)
  }
  dom
}


## Resolve the numeraire variable used for the consumption-equivalent
## normalisation.  "auto" picks the first of c/C/cons/consumption present
## among the endogenous variables; NA/NULL/FALSE disables normalisation.
#' @noRd
.acc_numeraire <- function(normalize_by, endo) {
  if (is.null(normalize_by) || isFALSE(normalize_by) ||
      (length(normalize_by) == 1L && is.na(normalize_by)))
    return(NA_character_)
  if (identical(normalize_by, "auto")) {
    hit <- intersect(c("c", "C", "cons", "consumption"), endo)
    return(if (length(hit)) hit[1L] else NA_character_)
  }
  if (!normalize_by %in% endo)
    stop("euler_errors: normalize_by = '", normalize_by,
         "' is not an endogenous variable of the model.", call. = FALSE)
  normalize_by
}


## The engine: expected residuals + expected Jacobians on a grid of state
## lags.  Returns the per-point n_pts x n_eq matrices the public methods
## summarise.
#' @noRd
.acc_residual_grid <- function(pol, compiled, params, grid_nat, n_quad,
                               normalize_by = "auto") {
  map  <- .acc_dy_map(compiled)
  dyn  <- compiled$dynamic
  n_eq <- map$n_eq
  exo  <- map$exo
  sn   <- pol$state_names
  ssv  <- pol$ss_vals

  if (!all(map$m1$name %in% sn))
    stop("euler_errors: the model has lagged variables (",
         paste(setdiff(map$m1$name, sn), collapse = ", "),
         ") that are not states of the supplied solution.", call. = FALSE)

  quad    <- .acc_quad(n_quad, pol$shock_sds[exo])
  n_combo <- nrow(quad$nodes)
  zero_e  <- setNames(numeric(length(exo)), exo)

  nu      <- .acc_numeraire(normalize_by, map$endo)
  nu_col  <- if (is.na(nu)) NA_integer_ else unname(map$col0_of[nu])

  n_pts <- nrow(grid_nat)
  R     <- matrix(NA_real_, n_pts, n_eq)
  S     <- matrix(NA_real_, n_pts, n_eq)   # normalisation scale
  Mg    <- matrix(NA_real_, n_pts, n_eq)   # equation term magnitude

  for (j in seq_len(n_pts)) {
    state_lag <- setNames(grid_nat[j, sn], sn)
    y_t <- tryCatch(pol$y_at(state_lag, zero_e), error = function(e) NULL)
    if (is.null(y_t) || anyNA(y_t) || !all(is.finite(y_t))) next
    nl <- pol$next_lag(y_t)

    acc_r <- numeric(n_eq)
    acc_j <- numeric(n_eq)
    acc_m <- numeric(n_eq)
    ok    <- TRUE
    for (ki in seq_len(n_combo)) {
      eps_next <- setNames(quad$nodes[ki, ], exo)
      y_lead <- tryCatch(pol$y_at(nl, eps_next), error = function(e) NULL)
      if (is.null(y_lead) || !all(is.finite(y_lead))) { ok <- FALSE; break }
      dy   <- .acc_assemble_dy(map, state_lag, y_t, y_lead, zero_e)
      resk <- tryCatch(dyn$residuals_fn(dy, params, ssv),
                       error = function(e) NULL)
      if (is.null(resk) || !all(is.finite(resk))) { ok <- FALSE; break }
      acc_r <- acc_r + quad$weights[ki] * resk
      if (!is.na(nu_col)) {
        Jk <- tryCatch(dyn$jacobian_fn(dy, params, ssv),
                       error = function(e) NULL)
        if (is.null(Jk) || !all(is.finite(Jk))) { ok <- FALSE; break }
        acc_j <- acc_j + quad$weights[ki] * Jk[, nu_col]
        acc_m <- acc_m + quad$weights[ki] *
          apply(abs(Jk) * rep(abs(dy), each = n_eq), 1L, max)
      }
    }
    if (!ok) next
    R[j, ] <- acc_r
    if (!is.na(nu_col)) {
      S[j, ] <- acc_j * y_t[[nu]]
      Mg[j, ] <- acc_m
    }
  }

  ## An equation is reported in consumption-equivalent units only when the
  ## numeraire channel is a non-negligible part of it at EVERY grid point --
  ## otherwise (an AR(1) shock process, in which c_t does not appear at all)
  ## dividing by ~0 would manufacture an enormous fake "error".
  normalized <- rep(FALSE, n_eq)
  if (!is.na(nu_col)) {
    tolS <- sqrt(.Machine$double.eps)
    fin  <- stats::complete.cases(S)
    if (any(fin))
      normalized <- apply(abs(S[fin, , drop = FALSE]) >
                            tolS * pmax(Mg[fin, , drop = FALSE], 1e-300),
                          2L, all)
  }

  E <- R
  for (i in seq_len(n_eq)) if (normalized[i]) E[, i] <- R[, i] / S[, i]

  list(residuals = R, errors = E, scale = S,
       normalized = normalized, numeraire = nu,
       has_lead = .acc_lead_equations(compiled, params, ssv),
       n_eq = n_eq, map = map, grid = grid_nat)
}


## Assemble the public return value from an .acc_residual_grid() result.
#' @noRd
.acc_summarise <- function(res, eq_labels = NULL) {
  E   <- res$errors
  L10 <- suppressWarnings(log10(abs(E)))
  colnm <- eq_labels %||% paste0("eq", seq_len(res$n_eq))
  colnames(E) <- colnm; colnames(L10) <- colnm
  colnames(res$residuals) <- colnm
  names(res$normalized)   <- colnm

  safe <- function(v) if (!any(is.finite(v))) NA_real_ else max(v[is.finite(v)])
  mn   <- function(v) if (!any(is.finite(v))) NA_real_ else mean(v[is.finite(v)])

  max_by_eq  <- apply(L10, 2L, safe)
  mean_by_eq <- apply(L10, 2L, mn)

  ## `max_log10_error` is the EULER error: the max over the expectational
  ## equations (those carrying a t+1 term).  A static identity -- the resource
  ## constraint, an AR(1) shock process -- is imposed exactly at the
  ## collocation nodes and its off-node residual measures interpolation, not
  ## the accuracy of the conditional expectation; folding it into the headline
  ## number would report a projection solution as three orders of magnitude
  ## worse than it is.  `max_log10_all` keeps the across-everything figure.
  hl <- res$has_lead
  if (!any(hl)) hl <- rep(TRUE, res$n_eq)
  names(hl) <- colnm
  Lh <- L10[, hl, drop = FALSE]

  list(
    errors           = E,
    log10_errors     = L10,
    residuals        = res$residuals,
    max_log10_by_eq  = max_by_eq,
    mean_log10_by_eq = mean_by_eq,
    max_log10_error  = safe(as.numeric(Lh)),
    mean_log10_error = mn(as.numeric(Lh)),
    max_log10_all    = safe(as.numeric(L10)),
    mean_log10_all   = mn(as.numeric(L10)),
    has_lead         = hl,
    normalized       = res$normalized,
    numeraire        = res$numeraire,
    grid             = res$grid
  )
}


## Equations that carry at least one t+1 column (the expectational block).
#' @noRd
.acc_lead_equations <- function(compiled, params, ss_vals) {
  map <- .acc_dy_map(compiled)
  if (!length(map$p1$col)) return(logical(map$n_eq))
  dyn <- compiled$dynamic
  dy  <- map$template
  ## Probe the sparsity pattern at a generic point (the steady state jittered
  ## so that a term whose derivative vanishes AT the SS is still detected).
  for (k in seq_along(map$keys)) {
    nm <- sub("__(0|p1|m1)$", "", map$keys[k])
    v  <- if (nm %in% names(ss_vals)) ss_vals[[nm]] else 0
    dy[k] <- v + 1e-3 * (1 + abs(v))
  }
  J <- tryCatch(dyn$jacobian_fn(dy, params, ss_vals), error = function(e) NULL)
  if (is.null(J)) return(rep(TRUE, map$n_eq))
  rowSums(abs(J[, map$p1$col, drop = FALSE]) > 0) > 0
}


## Build the regular (non-Chebyshev) evaluation grid from a state domain.
#' @noRd
.acc_grid <- function(state_domain, state_names, n_grid) {
  g1 <- lapply(state_names, function(nm)
    seq(state_domain[[nm]][1], state_domain[[nm]][2], length.out = n_grid))
  names(g1) <- state_names
  args <- rev(g1); names(args) <- rev(state_names)
  out <- as.matrix(expand.grid(args))
  out[, state_names, drop = FALSE]
}


#' Model-agnostic Euler-equation accuracy for a GlobalSolution
#'
#' Evaluates the compiled model's OWN dynamic residuals under the projection
#' policy on a regular grid of state lags, taking the \eqn{t+1} expectation by
#' Gauss-Hermite quadrature over the shocks.  Residuals of equations in which
#' the numeraire variable appears are reported in consumption-equivalent
#' units (the fractional change in \eqn{\nu_t} that would zero the residual);
#' the rest are reported raw.  See the file header for the definition.
#'
#' @param x       A \code{GlobalSolution} object.
#' @param n_grid  Number of test points per state dimension (default 10).
#' @param n_quad  Gauss-Hermite nodes per shock dimension (default 20).
#' @param state_domain Named list of \code{c(lo, hi)} per state giving the
#'   evaluation box.  \code{NULL} (default) uses the solution's own
#'   \code{state_domain}.  Supply one to compare a projection and a
#'   perturbation policy on the SAME box, or to stay strictly inside the
#'   collocation nodes: Chebyshev-Gauss nodes stop at
#'   \eqn{\cos(\pi/(2n))} of the half-width, so the default box's corners are
#'   a (short) EXTRAPOLATION and their error need not fall with
#'   \code{poly_degree}.
#' @param normalize_by Numeraire variable for the unit-free normalisation.
#'   \code{"auto"} (default) picks the first of \code{c}, \code{C},
#'   \code{cons}, \code{consumption} present in the model; a variable name
#'   forces one; \code{NA} reports every equation raw.
#' @param ...     Ignored.
#' @return A list with \code{errors} and \code{log10_errors} (both
#'   \code{n_points x n_eq} matrices), the raw \code{residuals},
#'   \code{max_log10_by_eq} / \code{mean_log10_by_eq}, the overall
#'   \code{max_log10_error} / \code{mean_log10_error}, the per-equation
#'   \code{normalized} flag, the \code{numeraire} used and the evaluation
#'   \code{grid}.
#' @seealso \code{\link{den_haan_marcet}}
#' @export
euler_errors.GlobalSolution <- function(x, n_grid = 10L, n_quad = 20L,
                                        state_domain = NULL,
                                        normalize_by = "auto", ...) {
  pol  <- .acc_policy_global(x)
  grid <- .acc_grid(state_domain %||% x$state_domain, pol$state_names, n_grid)
  res  <- .acc_residual_grid(pol, x$compiled, x$params, grid, n_quad,
                             normalize_by)
  .acc_summarise(res)
}


#' Model-agnostic Euler-equation accuracy for a perturbation decision rule
#'
#' The perturbation counterpart of \code{\link{euler_errors.GlobalSolution}}:
#' the Taylor policy is evaluated UNPRUNED (the standard accuracy convention)
#' at each grid point and the model's own residuals are integrated over the
#' \eqn{t+1} shocks by Gauss-Hermite quadrature.  Dispatches for
#' \code{DecisionRules}, \code{DecisionRules2} and \code{DecisionRules3}, so
#' the order-1 / order-2 / order-3 accuracy ordering is a single call each.
#'
#' @param x        A \code{DecisionRules}, \code{DecisionRules2} or
#'   \code{DecisionRules3} object.
#' @param compiled The \code{dynhr_compiled} the rule was solved from.
#' @param params   Named numeric parameter vector (defaults to
#'   \code{compiled$model$param_values}).
#' @param model    The \code{dynhr_mod} (defaults to \code{compiled$model});
#'   used only for the shock standard deviations.
#' @param n_grid   Test points per state dimension (default 10).
#' @param n_quad   Gauss-Hermite nodes per shock dimension (default 20).
#' @param state_domain Named list of \code{c(lo, hi)} per state.  \code{NULL}
#'   (default) uses steady state \eqn{\pm} \code{n_sd} unconditional standard
#'   deviations from the first-order stationary covariance.
#' @param n_sd     Half-width of the default box in state standard deviations
#'   (default 2).
#' @param normalize_by Numeraire variable; see
#'   \code{\link{euler_errors.GlobalSolution}}.
#' @param ...      Ignored.
#' @return The same list \code{\link{euler_errors.GlobalSolution}} returns.
#' @export
euler_errors.DecisionRules <- function(x, compiled, params = NULL,
                                       model = NULL,
                                       n_grid = 10L, n_quad = 20L,
                                       state_domain = NULL, n_sd = 2,
                                       normalize_by = "auto", ...) {
  if (missing(compiled) || is.null(compiled))
    stop("euler_errors: a `compiled` model is required for a perturbation ",
         "decision rule (the residuals are the model's own).", call. = FALSE)
  model  <- model  %||% compiled$model
  params <- params %||% model$param_values
  pol    <- .acc_policy_dr(x, model, params)
  dom    <- state_domain %||% .acc_state_domain_dr(x, model, params, n_sd)
  grid   <- .acc_grid(dom, pol$state_names, n_grid)
  res    <- .acc_residual_grid(pol, compiled, params, grid, n_quad,
                               normalize_by)
  .acc_summarise(res)
}


## ---- Den Haan - Marcet statistic ----------------------------------------

## Moore-Penrose inverse of a small symmetric PSD matrix (base R only; MASS
## is not a dynhr dependency).  The DHM weighting matrix is a sample second-
## moment matrix and can be numerically singular when an instrument is
## collinear with another, so a plain solve() is not enough.
#' @noRd
.acc_pinv <- function(A, rtol = sqrt(.Machine$double.eps)) {
  sv   <- svd(A)
  keep <- sv$d > rtol * max(sv$d, 0)
  if (!any(keep)) return(matrix(0, nrow(A), ncol(A)))
  sv$v[, keep, drop = FALSE] %*%
    (t(sv$u[, keep, drop = FALSE]) / sv$d[keep])
}



#' Den Haan-Marcet accuracy statistic
#'
#' @param x   A solved policy: a \code{GlobalSolution} or a perturbation
#'   \code{DecisionRules*} object.
#' @param ... Method arguments; see
#'   \code{\link{den_haan_marcet.GlobalSolution}}.
#' @export
den_haan_marcet <- function(x, ...) UseMethod("den_haan_marcet")


## The shared statistic.  `pol` is a policy adapter, `compiled`/`params` the
## model.  Simulates the policy, forms the REALISED one-step residuals
## u_t = F(y_{t-1}, y_t, y_{t+1}, eps_t) (normalised exactly as in
## euler_errors()), and tests E[u_t (x) h_t] = 0 with the Den Haan-Marcet
## (1994) chi-square statistic
##
##   J_T = T * M' A^{-1} M,   M = mean_t (u_t (x) h_t),
##                            A = mean_t (u_t (x) h_t)(u_t (x) h_t)'
##
## which is chi-square with dim(u) * dim(h) degrees of freedom under the null
## that the policy solves the model exactly.
#' @noRd
.acc_dhm <- function(pol, compiled, params, n_sim, burn_in, instruments,
                     normalize_by, eq_idx, init_state) {
  map  <- .acc_dy_map(compiled)
  dyn  <- compiled$dynamic
  sn   <- pol$state_names
  exo  <- map$exo
  ssv  <- pol$ss_vals
  sds  <- pol$shock_sds[exo]

  if (!all(map$m1$name %in% sn))
    stop("den_haan_marcet: the model has lagged variables that are not ",
         "states of the supplied solution.", call. = FALSE)
  if (!length(eq_idx))
    stop("den_haan_marcet: no equations to test -- the model carries no ",
         "expectational (t+1) term. Pass `equations =` explicitly.",
         call. = FALSE)

  total <- n_sim + burn_in + 1L
  lag_m <- matrix(NA_real_, total, length(sn), dimnames = list(NULL, sn))
  y_m   <- matrix(NA_real_, total, length(pol$endo),
                  dimnames = list(NULL, pol$endo))
  e_m   <- matrix(NA_real_, total, length(exo), dimnames = list(NULL, exo))

  lag <- if (is.null(init_state)) setNames(ssv[sn], sn) else
           setNames(init_state[sn], sn)
  for (t in seq_len(total)) {
    eps <- setNames(stats::rnorm(length(exo), 0, sds), exo)
    y_t <- pol$y_at(lag, eps)
    if (!all(is.finite(y_t)))
      stop("den_haan_marcet: the simulated policy left the finite domain at ",
           "period ", t, ".", call. = FALSE)
    lag_m[t, ] <- lag; y_m[t, ] <- y_t; e_m[t, ] <- eps
    lag <- pol$next_lag(y_t)
  }

  keep   <- (burn_in + 1L):(total - 1L)
  nu     <- .acc_numeraire(normalize_by, map$endo)
  nu_col <- if (is.na(nu)) NA_integer_ else unname(map$col0_of[nu])

  n_eq <- map$n_eq
  U    <- matrix(NA_real_, length(keep), n_eq)
  for (ii in seq_along(keep)) {
    t  <- keep[ii]
    dy <- .acc_assemble_dy(map, lag_m[t, ], y_m[t, ], y_m[t + 1L, ], e_m[t, ])
    r  <- dyn$residuals_fn(dy, params, ssv)
    if (!is.na(nu_col)) {
      J <- dyn$jacobian_fn(dy, params, ssv)
      s <- J[, nu_col] * y_m[t, nu]
      r <- ifelse(abs(s) > 0, r / s, r)
    }
    U[ii, ] <- r
  }

  U <- U[, eq_idx, drop = FALSE]
  H <- if (is.null(instruments))
         cbind(1, sweep(lag_m[keep, , drop = FALSE], 2L, ssv[sn], "-"))
       else instruments
  H <- as.matrix(H)
  if (nrow(H) != nrow(U))
    stop("den_haan_marcet: `instruments` must have one row per retained ",
         "simulation period (", nrow(U), ").", call. = FALSE)

  T_eff <- nrow(U)
  nq    <- ncol(U) * ncol(H)
  G     <- matrix(NA_real_, T_eff, nq)
  for (t in seq_len(T_eff)) G[t, ] <- as.numeric(U[t, ] %x% H[t, ])

  M <- colMeans(G)
  A <- crossprod(G) / T_eff
  Ai <- .acc_pinv(A)
  J_stat <- T_eff * drop(crossprod(M, Ai %*% M))

  list(J_stat   = J_stat,
       df       = nq,
       p_value  = stats::pchisq(J_stat, nq, lower.tail = FALSE),
       crit_95  = stats::qchisq(0.95, nq),
       reject   = J_stat > stats::qchisq(0.95, nq),
       n_used   = T_eff,
       equations = eq_idx,
       numeraire = nu,
       residual_sd = apply(U, 2L, stats::sd))
}


#' Den Haan-Marcet statistic for a GlobalSolution
#'
#' Simulates the projection policy, forms the realised one-step model
#' residuals in the same (consumption-equivalent where available) units as
#' \code{\link{euler_errors}}, and tests their orthogonality to a set of
#' time-\eqn{t} instruments.  Under the null that the policy solves the model
#' exactly the statistic is \eqn{\chi^2} with \code{df} degrees of freedom, so
#' \code{J_stat > crit_95} rejects the policy.
#'
#' @param x        A \code{GlobalSolution}.
#' @param n_sim    Simulation length AFTER burn-in (default 10000).
#' @param burn_in  Discarded initial periods (default 200).
#' @param instruments Optional \code{n_sim x q} instrument matrix.
#'   \code{NULL} (default) uses \code{(1, state lags in deviation)}.
#' @param normalize_by Numeraire variable; see \code{\link{euler_errors}}.
#' @param equations Integer positions of the equations to test.  \code{NULL}
#'   (default) uses the EXPECTATIONAL equations (those carrying a \eqn{t+1}
#'   term); a static identity has an identically-zero residual and would make
#'   the weighting matrix singular.
#' @param seed     Optional integer seed.  Applied through
#'   \code{.with_local_seed()}, so the caller's RNG stream is untouched.
#' @param init_state Optional named initial state lag (default: steady state).
#' @param ...      Ignored.
#' @return A list with \code{J_stat}, \code{df}, \code{p_value},
#'   \code{crit_95}, \code{reject}, \code{n_used}, \code{equations},
#'   \code{numeraire} and the per-equation \code{residual_sd}.
#' @export
den_haan_marcet.GlobalSolution <- function(x, n_sim = 10000L, burn_in = 200L,
                                           instruments = NULL,
                                           normalize_by = "auto",
                                           equations = NULL, seed = NULL,
                                           init_state = NULL, ...) {
  pol <- .acc_policy_global(x)
  eqi <- equations %||% which(.acc_lead_equations(x$compiled, x$params,
                                                  x$ss_vals))
  .with_local_seed(seed,
    .acc_dhm(pol, x$compiled, x$params, n_sim, burn_in, instruments,
             normalize_by, eqi, init_state))
}


#' Den Haan-Marcet statistic for a perturbation decision rule
#'
#' @inheritParams euler_errors.DecisionRules
#' @param n_sim,burn_in,instruments,normalize_by,equations,seed,init_state
#'   As in \code{\link{den_haan_marcet.GlobalSolution}}.
#' @return As \code{\link{den_haan_marcet.GlobalSolution}}.
#' @export
den_haan_marcet.DecisionRules <- function(x, compiled, params = NULL,
                                          model = NULL,
                                          n_sim = 10000L, burn_in = 200L,
                                          instruments = NULL,
                                          normalize_by = "auto",
                                          equations = NULL, seed = NULL,
                                          init_state = NULL, ...) {
  if (missing(compiled) || is.null(compiled))
    stop("den_haan_marcet: a `compiled` model is required for a perturbation ",
         "decision rule.", call. = FALSE)
  model  <- model  %||% compiled$model
  params <- params %||% model$param_values
  pol    <- .acc_policy_dr(x, model, params)
  eqi    <- equations %||% which(.acc_lead_equations(compiled, params, x$ys))
  .with_local_seed(seed,
    .acc_dhm(pol, compiled, params, n_sim, burn_in, instruments,
             normalize_by, eqi, init_state))
}
