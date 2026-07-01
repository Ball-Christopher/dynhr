## R/solve-perturbation-order3-sigma-helpers.R
## --------------------------------------------------------------------------
## Compound-derivative builders for Phase 7c.2 / 7c.3 (ghuss, ghxss).
##
## Mutschler (2022) perturbation_solver_nonsymmetric_order3.m operates on
## a stacked compound vector  z = [y_{t-1}; y_t; y_{t+1}; u_t]  in DR order
## and writes the sigma-correction RHS in terms of partial derivatives of
## z with respect to (x, u, u', s):
##
##   zx     = dz/dx     (already = T_x in dynhr)
##   zu     = dz/du     (already = T_u in dynhr)
##   zup    = dz/du'                                  -- THIS FILE
##   zss    = dz/d(s^2)/2  (the ghs2-driven sigma piece) -- THIS FILE
##   zxup   = d^2 z / dx du'                          -- THIS FILE
##   zuup   = d^2 z / du du'                          -- THIS FILE
##   zupup  = d^2 z / du' du'                         -- THIS FILE
##
## In dynhr's column convention z is indexed by `total_cols` (one column
## per (variable name, lead_lag) pair in dyn_col_map).  These builders
## return matrices of shape (total_cols x n_*), so they can be left- or
## right-multiplied into the dense second- and third-order Hessians sliced
## at compound-column indices -- exactly the pattern used by ghs3.
##
## Layout for kron columns: (a_outer slow, k_shock fast), matching MATLAB
## kron(A, I_u) ordering used by Mutschler.
## --------------------------------------------------------------------------


#' Build T_up = dz/du' (total_cols x n_u).
#'
#' Non-zero only at lead-block jumper compound-cols, where row equals the
#' corresponding jumper row of `ghu`.
#'
#' @noRd
.build_T_up <- function(dyn, ghu, endo_names, exo_names, has_lead) {
  total_cols <- dyn$total_cols
  n_u        <- length(exo_names)
  T_up       <- matrix(0, total_cols, n_u)

  jm  <- .identify_jumpers(dyn, endo_names, has_lead)
  for (i in seq_along(jm$jumper_idx)) {
    j <- jm$jumper_idx[i]
    c <- jm$jumper_compound_c[i]
    T_up[c, ] <- ghu[j, ]
  }
  T_up
}


#' Build T_ss = dz/d(s^2)/2 (total_cols x 1).
#'
#' Mutschler line 128-131:
#'   lag block       : 0
#'   current block   : ghs2[j]
#'   lead jumper blk : ghs2[j] + (ghx[j, ]) %*% ghs2[state_idx]
#'   shock block     : 0
#'
#' @noRd
.build_T_ss <- function(dyn, ghx, ghs2, state_idx, endo_names, exo_names,
                        has_lead) {
  total_cols <- dyn$total_cols
  T_ss       <- matrix(0, total_cols, 1L)
  dcm        <- dyn$dyn_col_map

  ghs2_state <- ghs2[state_idx]

  for (k in seq_len(nrow(dcm))) {
    c  <- dcm$col[k]
    nm <- dcm$name[k]
    ll <- dcm$lead_lag[k]
    if (nm %in% exo_names) next
    j <- which(endo_names == nm)
    if (length(j) != 1L) next

    if (ll == 0L) {
      T_ss[c, 1] <- ghs2[j]
    } else if (ll == 1L && isTRUE(has_lead[j])) {
      T_ss[c, 1] <- ghs2[j] + sum(ghx[j, ] * ghs2_state)
    }
  }
  T_ss
}


#' Build W_xup = d^2 z / dx du' (total_cols x (n_s * n_u)).
#'
#' Lead jumper row j: W_xup[c, (a'-1)*n_u + k'] = sum_a ghxu[j,(a-1)*n_u+k'] * hx[a,a'].
#' Column layout: (s_outer=a' slow, k_shock=k' fast).
#'
#' `ghxu` is expected in STANDARD Kron convention (state SLOW, shock FAST):
#' col = (state-1)*n_u + shock  (as stored in dr3$ghxu).
#' Reshaping ghxu[j,] as (n_u x n_s) gives a [shock, state] matrix; transposing
#' gives ghxu_mat[state=a, shock=k].  Then t(ghxu_mat) %*% hx is (n_u x n_s)
#' with element [k', a'] = sum_a ghxu_mat[a,k']*hx[a,a'].
#' as.vector of that matrix: k' FAST (rows), a' SLOW (cols) → correct layout.
#'
#' @noRd
.build_W_xup <- function(dyn, ghxu, hx, endo_names, exo_names, has_lead) {
  total_cols <- dyn$total_cols
  n_s        <- ncol(hx)
  n_u        <- length(exo_names)
  W_xup      <- matrix(0, total_cols, n_s * n_u)

  if (n_s == 0L || n_u == 0L) return(W_xup)

  jm  <- .identify_jumpers(dyn, endo_names, has_lead)
  for (i in seq_along(jm$jumper_idx)) {
    j <- jm$jumper_idx[i]
    c <- jm$jumper_compound_c[i]
    # dr3$ghxu uses STANDARD convention: col = (state-1)*n_u + shock.
    # matrix(., n_u, n_s) → [shock, state]; t(.) → ghxu_mat[state=a, shock=k].
    ghxu_mat <- t(matrix(ghxu[j, , drop = FALSE], n_u, n_s))  # [a, k]
    W_xup[c, ] <- as.vector(t(ghxu_mat) %*% hx)               # (n_u x n_s): k' fast, a' slow
  }
  W_xup
}


#' Build W_uup = d^2 z / du du' (total_cols x (n_u * n_u)).
#'
#' Lead jumper row j: W_uup[c, (l-1)*n_u + k'] = sum_a ghxu[j,(a-1)*n_u+k'] * hu[a,l].
#' Column layout: (l=current-shock slow, k'=lead-shock fast).
#'
#' `ghxu` is in STANDARD Kron convention (state SLOW, shock FAST) as stored in
#' dr3$ghxu.  Same reshape trick as W_xup: t(matrix(., n_u, n_s)) gives
#' ghxu_mat[a, k], then t(ghxu_mat) %*% hu is (n_u x n_u) with
#' element [k', l] = sum_a ghxu_mat[a,k']*hu[a,l]. as.vector: k' FAST, l SLOW.
#'
#' @noRd
.build_W_uup <- function(dyn, ghxu, hu, endo_names, exo_names, has_lead) {
  total_cols <- dyn$total_cols
  n_s        <- nrow(hu)
  n_u        <- length(exo_names)
  W_uup      <- matrix(0, total_cols, n_u * n_u)

  if (n_u == 0L) return(W_uup)

  jm  <- .identify_jumpers(dyn, endo_names, has_lead)
  for (i in seq_along(jm$jumper_idx)) {
    j <- jm$jumper_idx[i]
    c <- jm$jumper_compound_c[i]
    # dr3$ghxu uses STANDARD convention: col = (state-1)*n_u + shock.
    # matrix(., n_u, n_s) → [shock, state]; t(.) → ghxu_mat[state=a, shock=k].
    ghxu_mat <- t(matrix(ghxu[j, , drop = FALSE], n_u, n_s))  # [a, k]
    W_uup[c, ] <- as.vector(t(ghxu_mat) %*% hu)               # (n_u x n_u): k' fast, l slow
  }
  W_uup
}


#' Build W_upup = d^2 z / du' du' (total_cols x n_u^2).
#'
#' Lead jumper rows: ghuu[j, ] directly (no chain through h: both u's are
#' future shocks).
#'
#' @noRd
.build_W_upup <- function(dyn, ghuu, endo_names, exo_names, has_lead) {
  total_cols <- dyn$total_cols
  n_u        <- length(exo_names)
  W_upup     <- matrix(0, total_cols, n_u * n_u)

  if (n_u == 0L) return(W_upup)

  jm  <- .identify_jumpers(dyn, endo_names, has_lead)
  for (i in seq_along(jm$jumper_idx)) {
    j <- jm$jumper_idx[i]
    c <- jm$jumper_compound_c[i]
    W_upup[c, ] <- ghuu[j, ]
  }
  W_upup
}


# =====================================================================
# FD verification helpers (used by tests, not by the solvers)
# =====================================================================

#' Evaluate the policy y_t = g(x_t, u_t, sigma) via the dr2 second-order
#' decision rule (including the ghss * sigma^2 / 2 term).
#'
#' Used for FD verification of T_up, T_ss, W_xup, W_uup, W_upup: by
#' constructing the compound vector z(x_t, u_t, sigma, u_{t+1}) =
#' (y_{t-1}=SS_state, y_t = g, y_{t+1} = g(h(x_t,u_t,sigma), u_{t+1}, sigma),
#'  u_t), we can numerically differentiate z and compare against the
#' analytic builders.
#'
#' @noRd
.policy_eval_o2 <- function(dr2, x_dev, u, sigma) {
  ss   <- numeric(length(dr2$endo_names))
  ss[] <- 0   # we work in deviations
  n_s  <- length(dr2$state_idx)
  n_u  <- length(dr2$exo_names)
  # ghxu cols are (state FAST, exo SLOW); the matching Kron is (u %x% x_dev).
  y <- as.numeric(dr2$ghx %*% x_dev) +
       as.numeric(dr2$ghu %*% u) +
       0.5 * as.numeric(dr2$ghxx %*% (x_dev %x% x_dev)) +
       as.numeric(dr2$ghxu %*% (u %x% x_dev)) +
       0.5 * as.numeric(dr2$ghuu %*% (u %x% u)) +
       0.5 * sigma^2 * dr2$ghss
  y
}


#' Evaluate the compound vector z in dynhr's total_cols layout for a
#' given (x_dev, u, sigma, up).  Returns a numeric vector of length
#' total_cols.
#'
#' The state-lag block (ll=-1) takes the SS values of state vars
#' (deviations = 0).  The current block (ll=0) takes y_t = g(x_dev, u, sigma).
#' The lead block (ll=+1) takes y_{t+1} = g(h(x_dev, u, sigma), up, sigma)
#' where h = state-row slice of g.  The shock block takes u.
#'
#' @noRd
.compound_z <- function(dyn, dr2, x_dev, u, sigma, up) {
  total_cols <- dyn$total_cols
  dcm        <- dyn$dyn_col_map
  state_idx  <- dr2$state_idx
  endo_names <- dr2$endo_names
  exo_names  <- dr2$exo_names

  y_t   <- .policy_eval_o2(dr2, x_dev, u,  sigma)
  x_new <- y_t[state_idx]            # h(x,u,sigma) deviations
  y_tp1 <- .policy_eval_o2(dr2, x_new, up, sigma)

  z <- numeric(total_cols)
  for (k in seq_len(nrow(dcm))) {
    c  <- dcm$col[k]
    nm <- dcm$name[k]
    ll <- dcm$lead_lag[k]
    if (nm %in% exo_names) {
      k_u <- which(exo_names == nm)
      if (length(k_u) == 1L) z[c] <- u[k_u]
      next
    }
    j <- which(endo_names == nm)
    if (length(j) != 1L) next
    if (ll == -1L) {
      # State-lag block: equals x_dev[s] where s is the position of j
      # within state_idx (deviations of predetermined state at t).
      s <- which(state_idx == j)
      if (length(s) == 1L) z[c] <- x_dev[s]
    } else if (ll == 0L) {
      z[c] <- y_t[j]
    } else if (ll == 1L) {
      z[c] <- y_tp1[j]
    }
  }
  z
}
