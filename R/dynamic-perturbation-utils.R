## R/dynamic-perturbation-utils.R
## --------------------------------------------------------------------------
## Helper utilities for dynamic_path_perturbation():
##   .dp_build_dy()       -- dy-format adapter: raw [xp; y; xt] -> named dy
##   .dp_split_jacobian() -- IFT block split: jacobian -> (Eaa, Ebb)
##
## The dy-format adapter is the critical engineering piece (scope C7 §8):
## dynhr's jacobian_fn / residuals_fn take a named vector `dy` with keys
## like "c__m1", "pi__0", "v__p1", whereas the standalone DP code works
## with raw vectors z = [xp; y] and xcur = xt.  This module converts
## between the two representations using the model's dyn_col_map.
##
## Column ordering convention (from extract_system_matrices):
##   f_minus  = df/d(state at t-1) = df/d xt   [<-> Ebb in DP notation]
##   f_zero   = df/d(all vars at t)             [<-> Eaa restricted to
##                                                state__0 + jump__0]
##   f_plus   = df/d(jump at t+1)  = df/d y'   [handled by ghat policy]
##
## IFT split (matching standalone grad_dynamic):
##   Eaa = df/d[xp, y] = f_zero[, c(state_cols, jump_cols)]
##   Ebb = df/d xt     = f_minus[, state_cols]
##
## The IFT update is then: tmp = -solve(Eaa, Ebb)  (same as standalone).
## --------------------------------------------------------------------------

#' Build a named dy vector for jacobian_fn / residuals_fn from raw DP variables.
#'
#' Maps the standalone DP representation [xp; y; xt] (where xp = next-period
#' states, y = current jumps, xt = current states) to dynhr's named `dy`
#' vector using the model's `dyn_col_map`.
#'
#' The mapping follows dynhr's lead/lag convention:
#'   xt  (current state x_t)   -> state variables at lead_lag = -1  (__m1)
#'   xp  (next-period state)   -> state variables at lead_lag =  0  (__0)
#'   y   (current jump)        -> jump variables  at lead_lag =  0  (__0)
#'   yp  (next-period jump from ghat policy) -> jump variables at lead_lag = +1 (__p1)
#'
#' @param xp       Numeric vector: next-period states (length nx).
#' @param y        Numeric vector: current jumps (length ny).
#' @param xt       Numeric vector: current states (length nx).
#' @param yp       Numeric vector: next-period jumps from ghat policy (length ny).
#' @param state_vars Character vector: names of state (predetermined) variables.
#' @param jump_vars  Character vector: names of jump (non-predetermined) variables.
#' @param dyn_col_map Data frame with columns name, lead_lag, col.
#' @param ss       Numeric vector: steady state (for SS-oracle arg).
#' @return Named numeric vector `dy` suitable for jacobian_fn / residuals_fn.
#' @noRd
.dp_build_dy <- function(xp, y, xt, yp, state_vars, jump_vars, dyn_col_map, ss) {
  n_cols <- nrow(dyn_col_map)
  dy <- numeric(n_cols)
  names(dy) <- character(n_cols)  # will be named below

  for (k in seq_len(n_cols)) {
    nm  <- dyn_col_map$name[k]
    ll  <- dyn_col_map$lead_lag[k]
    sfx <- if (ll == 0L) "__0"
           else if (ll > 0L) paste0("__p", ll)
           else paste0("__m", abs(ll))
    key <- paste0(nm, sfx)
    names(dy)[k] <- key

    if (nm %in% state_vars) {
      if (ll == -1L) {
        # state at t-1 = current state xt
        idx <- match(nm, state_vars)
        dy[k] <- xt[idx]
      } else if (ll == 0L) {
        # state at t = next-period state xp (the choice made this period)
        idx <- match(nm, state_vars)
        dy[k] <- xp[idx]
      } else {
        # state at t+n (n>=1): approximate via ss (should not occur in standard models)
        dy[k] <- if (nm %in% names(ss)) ss[[nm]] else 0
      }
    } else if (nm %in% jump_vars) {
      if (ll == 0L) {
        # jump at t = current jump y
        idx <- match(nm, jump_vars)
        dy[k] <- y[idx]
      } else if (ll == 1L) {
        # jump at t+1 = ghat prediction yp
        idx <- match(nm, jump_vars)
        dy[k] <- yp[idx]
      } else {
        # exo shocks or other leads: use ss
        dy[k] <- if (nm %in% names(ss)) ss[[nm]] else 0
      }
    } else {
      # Exogenous (varexo) at t: zero shock (deterministic path)
      dy[k] <- 0
    }
  }
  dy
}

#' Split the full Jacobian from jacobian_fn into Eaa and Ebb blocks.
#'
#' After solving the nonlinear system at [xp; y] given xt, the IFT update
#' requires:
#'   Eaa = df/d[xp, y]  (the Jacobian block w.r.t. what we solved for)
#'   Ebb = df/d xt      (the Jacobian block w.r.t. the current state)
#'
#' Using dynhr's column ordering via lead_lag_incidence or dyn_col_map:
#'   Eaa <-> f_zero columns for state + jump variables (all at t=0)
#'   Ebb <-> f_minus columns for state variables (at t=-1)
#'
#' The full Jacobian J = jacobian_fn(dy, params, ss) has columns in the order
#' given by dyn_col_map.  This function reads columns directly from J using
#' the column indices from dyn_col_map.
#'
#' @param J          Full Jacobian matrix from jacobian_fn (n_eq x n_dy_cols).
#' @param state_vars Character: names of state variables.
#' @param jump_vars  Character: names of jump variables.
#' @param dyn_col_map Data frame with name, lead_lag, col (1-indexed into J).
#' @param gg_df Numeric matrix (ny x nx): slope of the current local jump
#'   policy, d yp / d xp (the backward sweep's \code{gg$df}).  REQUIRED.  The
#'   next-period jumps are not free -- yp = ghat(xp) -- so the total derivative
#'   of f w.r.t. xp picks up a chain-rule term J[, jump@+1] \%*\% gg_df that is
#'   added to the xp block of Eaa.  Omitting it zeroes the forward-looking
#'   jumps' response and gives the wrong IFT local rule.
#' @return List with Eaa (n_eq x (nx+ny)) and Ebb (n_eq x nx).
#' @noRd
.dp_split_jacobian <- function(J, state_vars, jump_vars, dyn_col_map, gg_df) {
  nx <- length(state_vars)
  ny <- length(jump_vars)

  ## Find column indices in J for each block
  ## Eaa: columns for state variables at lead_lag=0 (xp) + jump variables at lead_lag=0 (y)
  ## Ebb: columns for state variables at lead_lag=-1 (xt)

  # Build lookup: (name, lead_lag) -> position k in dyn_col_map
  dcm <- dyn_col_map

  # Eaa: [xp block (states at 0), y block (jumps at 0)]
  Eaa <- matrix(0, nrow(J), nx + ny)
  for (i in seq_along(state_vars)) {
    nm <- state_vars[i]
    k  <- which(dcm$name == nm & dcm$lead_lag == 0L)
    if (length(k) == 1L) Eaa[, i] <- J[, k]
    # else: state doesn't appear at t=0 in Jacobian => column stays 0
  }
  for (j in seq_along(jump_vars)) {
    nm <- jump_vars[j]
    k  <- which(dcm$name == nm & dcm$lead_lag == 0L)
    if (length(k) == 1L) Eaa[, nx + j] <- J[, k]
  }

  ## Chain-rule correction for the xp block.  The next-period jumps are not
  ## free: yp = ghat(xp) = gg$f0 + gg_df %*% (xp - gg$x0), so dyp/dxp = gg_df.
  ## jacobian_fn reports df/dxp (state@0) and df/dyp (jump@+1) as if
  ## independent; the TOTAL derivative w.r.t. xp is
  ##   df/dxp|_total = J[, state@0] + (df/dyp) %*% (dyp/dxp)
  ##                 = J[, state@0] + J[, jump@+1] %*% gg_df.
  ## Without this term the forward-looking jumps' response is dropped and the
  ## IFT yields the wrong local rule (paper ZLB-NK interior error 1e-4 -> 2e-8
  ## once included).
  if (ny > 0L && nx > 0L) {
    J_yp <- matrix(0, nrow(J), ny)
    for (j in seq_along(jump_vars)) {
      nm <- jump_vars[j]
      k  <- which(dcm$name == nm & dcm$lead_lag == 1L)
      if (length(k) == 1L) J_yp[, j] <- J[, k]
    }
    Eaa[, seq_len(nx)] <- Eaa[, seq_len(nx)] + J_yp %*% gg_df
  }

  # Ebb: [xt block (states at -1)]
  Ebb <- matrix(0, nrow(J), nx)
  for (i in seq_along(state_vars)) {
    nm <- state_vars[i]
    k  <- which(dcm$name == nm & dcm$lead_lag == -1L)
    if (length(k) == 1L) Ebb[, i] <- J[, k]
    # else: state has no lag in equations => column stays 0
  }

  list(Eaa = Eaa, Ebb = Ebb)
}
