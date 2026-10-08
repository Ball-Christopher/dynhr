## R/balanced-growth.R
## --------------------------------------------------------------------------
## Balanced-growth steady state of a linear (affine) model with unit roots
## and drift -- the IRIS `sstate(m, 'growth=', true)` statement.
##
## A `model(linear)` such as
##     a - a(-1) = rho*(a(-1) - a(-2)) + (1 - rho)*g + e;
## has no STATIC steady state when g != 0: the static system A y = -c is
## inconsistent. It does have a balanced-growth path y_t = y0 + g t on which
## every dynamic equation holds exactly with the shocks at zero. Writing the
## period-t residual as B y_{t-1} + A y_t + C y_{t+1} + c (aux lag/lead
## variables are ordinary endogenous variables, so one lag and one lead are
## all there is), substituting the path and matching the terms in t^0 and t^1
## gives
##     S g = 0,   S y0 + (C - B) g = -c,   S = A + B + C.
## The stacked system is solved in the MINIMUM-NORM sense: a unit-root level
## (and, for an I(2) block, a growth direction) is indeterminate, and the
## minimum-norm member is as good a representative as any -- the deviations
## the filter estimates absorb it under a diffuse initialisation. A candidate
## is ACCEPTED only when the dynamic residual of y0 + g t is ~0 at two
## different t, so an inconsistent stacked system is still reported as
## converged = FALSE.
##
## Downstream the deviations x_t = y_t - (y0 + g t) obey the homogeneous
## system, so the linear decision rule (ghx, ghu) is unchanged; the growth
## becomes a deterministic observation trend (see .obs_trend_slopes() in
## R/state-space.R) and a reported state growth path.
## --------------------------------------------------------------------------


#' Solve the balanced-growth path of a linear model
#'
#' @param compiled Compiled model (its `dynamic` block is used).
#' @param params Named parameter vector.
#' @param endo_names Endogenous names (declaration order, incl. aux).
#' @param tol Absolute residual tolerance.
#' @return list(ok, y0, growth, residuals, max_residual)
#' @noRd
.linear_balanced_growth <- function(compiled, params, endo_names, tol) {
  dyn <- compiled$dynamic
  dcm <- dyn$dyn_col_map
  n   <- length(endo_names)
  zero_ss <- stats::setNames(numeric(n), endo_names)
  key_of <- function(nm, ll)
    paste0(nm, if (ll == 0L) "__0" else if (ll > 0L) paste0("__p", ll)
               else paste0("__m", abs(ll)))
  keys <- vapply(seq_len(nrow(dcm)),
                 function(k) key_of(dcm$name[k], dcm$lead_lag[k]), "")
  ## dy at a deterministic path value(name, lead_lag); exogenous at zero.
  dy_at <- function(val) {
    v <- numeric(nrow(dcm))
    for (k in seq_len(nrow(dcm)))
      if (dcm$name[k] %in% endo_names) v[k] <- val(dcm$name[k], dcm$lead_lag[k])
    stats::setNames(v, keys)
  }
  dy0  <- dy_at(function(nm, ll) 0)
  cvec <- as.numeric(dyn$residuals_fn(dy0, params, zero_ss))
  J    <- as.matrix(dyn$jacobian_fn(dy0, params, zero_ss))
  n_eq <- length(cvec)
  if (!all(is.finite(cvec)) || !all(is.finite(J)) || ncol(J) != nrow(dcm))
    return(list(ok = FALSE, max_residual = Inf))

  ## Collapse the dynamic Jacobian onto the endogenous variables:
  ## S = sum_ll J_ll, Dm = sum_ll ll * J_ll (= C - B for one lag/lead).
  S  <- matrix(0, n_eq, n)
  Dm <- matrix(0, n_eq, n)
  for (k in seq_len(nrow(dcm))) {
    j <- match(dcm$name[k], endo_names)
    if (is.na(j)) next
    ## dyn_col_map rows are ordered by Jacobian column (see
    ## .build_dy_ss_o2), so column dcm$col[k] belongs to row k.
    col <- dcm$col[k]
    S[, j]  <- S[, j]  + J[, col]
    Dm[, j] <- Dm[, j] + dcm$lead_lag[k] * J[, col]
  }
  K   <- rbind(cbind(S, Dm), cbind(matrix(0, n_eq, n), S))
  rhs <- c(-cvec, numeric(n_eq))
  sol <- as.numeric(.safe_inv(K) %*% rhs)
  y0  <- stats::setNames(sol[seq_len(n)], endo_names)
  g   <- stats::setNames(sol[n + seq_len(n)], endo_names)
  ## Exact zeros for components that are zero to round-off: a stationary
  ## variable's growth is 0 by construction, and a 1e-17 "growth" would
  ## otherwise switch on the trend machinery for it.
  g_scale <- max(1, abs(g))
  g[abs(g) <= 64 * .Machine$double.eps * g_scale] <- 0

  ## Accept only if y0 + g t solves the DYNAMIC system at two different t.
  path_res <- function(tt)
    as.numeric(dyn$residuals_fn(
      dy_at(function(nm, ll) y0[[nm]] + g[[nm]] * (tt + ll)), params, y0))
  r0 <- path_res(0)
  r1 <- path_res(1)
  max_r <- max(abs(c(r0, r1)))
  scale <- max(abs(cvec), max(abs(J)) * max(abs(c(y0, g)), 1))
  ok <- all(is.finite(c(r0, r1))) && max_r <= max(tol, 1e-8 * scale)
  list(ok = ok, y0 = y0, growth = g, residuals = stats::setNames(r0, NULL),
       max_residual = max_r)
}


## ---- Growth on the model / decision-rule objects ---------------------------
##
## The growth vector is a function of the parameters, so it lives on the
## decision rule (`dr$growth`, solved at the same parameters as dr$ys). The
## MODEL carries only the request (`model$balanced_growth = TRUE`, set by
## solve_model() when the solved growth is non-zero): every steady-state
## re-solve on that model -- a posterior evaluation at theta, a diagnostic --
## then solves the balanced-growth path again at its own parameters, and every
## site that refuses deterministic observation trends sees the flag on the
## model object it already receives.

#' Did the caller (or the model) request the balanced-growth solve?
#' @noRd
.ss_growth_requested <- function(growth, model) {
  if (!is.logical(growth) || length(growth) != 1L || is.na(growth))
    .dynhr_abort("solve_steady_state: `growth` must be TRUE or FALSE.",
                 class = "dynhr_error_invalid_argument")
  isTRUE(growth) || isTRUE(model$balanced_growth)
}

#' Non-zero growth carried by a decision rule, or NULL
#' @noRd
.dr_growth <- function(dr) {
  g <- if (is.list(dr)) dr$growth else NULL
  if (is.null(g) || !any(g != 0)) NULL else g
}

#' The model and the decision rule must agree about balanced growth.
#'
#' A growth `dr` with an unflagged model would be honoured by the filter but
#' NOT refused by the sites that see only the model (they would score
#' untrended data); a flagged model with a `dr` that carries no growth would
#' filter with the trend switched off. Both are refused.
#' @noRd
.check_growth_pair <- function(model, dr) {
  if (is.null(dr) || !is.list(dr)) return(invisible(TRUE))
  flagged <- isTRUE(model$balanced_growth)
  if (flagged && is.null(dr$growth))
    .dynhr_abort(
      "The model is a balanced-growth model (solved with ",
      "steady_options = list(growth = TRUE)) but the decision rule carries ",
      "no `growth`: it was not solved on the balanced-growth path. Use the ",
      "`dr` that solve_model() returned (or set dr$growth to the steady ",
      "state's `$growth` at the same parameters).",
      class = "dynhr_error_balanced_growth_mismatch")
  if (!flagged && !is.null(.dr_growth(dr)))
    .dynhr_abort(
      "The decision rule carries a non-zero balanced growth (`dr$growth`) ",
      "but the model object does not: pass the `model` that solve_model() ",
      "returned alongside that `dr` (sol$model, which carries ",
      "balanced_growth = TRUE), so that every routine can see the trend.",
      class = "dynhr_error_balanced_growth_mismatch")
  invisible(TRUE)
}

#' Deterministic balanced-growth path of named variables over sample periods
#'
#' Column t is y0 + g * (first_obs + t - 1): the same period index the
#' observation-trend machinery uses. Rows not in `dr$growth` (e.g. a
#' mixed-frequency augmentation state) are NA.
#' @return n_var x n_T matrix, or NULL when `dr` carries no growth.
#' @noRd
.growth_state_path <- function(dr, model, var_names, n_T) {
  g <- .dr_growth(dr)
  if (is.null(g)) return(NULL)
  first_obs <- model$observation_trends$first_obs %||% 1L
  tt <- first_obs - 1 + seq_len(n_T)
  ys <- dr$ys
  out <- matrix(NA_real_, length(var_names), n_T,
                dimnames = list(var_names, NULL))
  hit <- var_names %in% names(g) & var_names %in% names(ys)
  if (any(hit)) {
    v <- var_names[hit]
    out[hit, ] <- as.numeric(ys[v]) + outer(as.numeric(g[v]), tt)
  }
  out
}
