## R/hank-diagnostics.R
## --------------------------------------------------------------------------
## Diagnostics and plots for heterogeneous-agent (HANK) blocks: marginal
## propensity to consume (MPC) by wealth, sequence-space Jacobian heatmaps, and
## impulse-response plots.
## --------------------------------------------------------------------------


#' Marginal propensity to consume from a het-block steady state
#'
#' Computes the one-period MPC out of a small transitory cash-on-hand transfer
#' at each \code{(e, a)} gridpoint, the wealth-averaged MPC profile, and the
#' distribution-weighted aggregate MPC.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param dcash Numeric > 0: size of the transitory cash transfer used in the
#'   finite-difference (in cash-on-hand units).
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{mpc}}{\code{n_e x n_a} matrix of MPCs.}
#'     \item{\code{aggregate}}{Distribution-weighted aggregate MPC.}
#'     \item{\code{by_wealth}}{\code{data.frame(a, mpc)} — MPC averaged over
#'       income states with the conditional wealth distribution.}
#'   }
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' ag  <- hank_asset_grid(100, 60, 0)
#' blk <- hank_het_block(ag, inc$Pi, inc$e, beta = 0.96, eis = 1,
#'                       r = 0.02, w = 1)
#' hank_mpc(blk)$aggregate
#' @export
hank_mpc <- function(block, dcash = 1e-4) {
  a_grid <- block$a_grid; r <- block$r
  n_e <- block$n_e; n_a <- block$n_a
  mpc <- matrix(0, n_e, n_a)
  ## A cash transfer dcash raises coh by dcash <=> raises assets by dcash/(1+r).
  da <- dcash / (1 + r)
  for (e in seq_len(n_e)) {
    c_new <- .hank_interp1(a_grid, block$c[e, ], a_grid + da)
    mpc[e, ] <- (c_new - block$c[e, ]) / dcash
  }
  agg <- hank_aggregate(block$D, mpc)

  ## MPC by wealth: average over income states using the conditional
  ## distribution at each asset gridpoint.
  Dmat <- .hank_vec_to_mat(block$D, n_e, n_a)       # n_e x n_a
  col_mass <- colSums(Dmat)
  col_mass[col_mass <= 0] <- NA_real_
  mpc_by_a <- colSums(Dmat * mpc) / col_mass
  by_wealth <- data.frame(a = a_grid, mpc = mpc_by_a)

  ## MPC by income state: average over wealth using the conditional wealth
  ## distribution within each income state (per-income-state cut; satisfies the
  ## multi-value vectorization convention across the >= 2 income states).
  row_mass <- rowSums(Dmat)
  row_mass[row_mass <= 0] <- NA_real_
  mpc_by_e <- rowSums(Dmat * mpc) / row_mass
  by_income <- data.frame(e = block$e, mpc = mpc_by_e)

  list(mpc = mpc, aggregate = agg, by_wealth = by_wealth, by_income = by_income)
}


#' Gini coefficient of a discrete distribution
#'
#' @param values Numeric vector of outcomes (e.g. assets), any order.
#' @param mass Numeric vector of probability masses (same length; need not sum
#'   to exactly 1 — it is renormalized).
#' @return The Gini coefficient in \[0, 1) (0 = perfect equality).
#' @examples
#' hank_gini(c(0, 1), c(0.5, 0.5))   # 0.5
#' @export
hank_gini <- function(values, mass) {
  keep <- mass > 0
  v <- values[keep]; p <- mass[keep] / sum(mass[keep])
  o <- order(v); v <- v[o]; p <- p[o]
  mu <- sum(p * v)
  if (mu <= 0) return(0)                     # no positive wealth -> define 0
  ## Lorenz-curve (trapezoidal) Gini: G = 1 - sum p_k (L_{k-1} + L_k), where L
  ## is the cumulative wealth share.
  cw <- cumsum(p * v) / mu                    # cumulative wealth share L_k
  L_prev <- c(0, cw[-length(cw)])
  1 - sum(p * (L_prev + cw))
}


#' Wealth/consumption distribution statistics for a het block
#'
#' Summarizes the stationary cross-sectional distribution: Gini, top wealth
#' shares, the hand-to-mouth fraction (mass at the borrowing constraint), and
#' wealth percentiles.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param top Numeric vector of top-share fractions (default \code{c(0.1, 0.01)}
#'   for the top-10\% and top-1\% wealth shares).
#' @param probs Numeric vector of percentile levels for the wealth quantiles.
#' @param constraint_tol Agents with assets within this of \code{amin} are
#'   counted as hand-to-mouth.
#'
#' @return A list with \code{gini} (wealth), \code{top_shares} (named by
#'   \code{top}), \code{hand_to_mouth} (mass fraction at the constraint),
#'   \code{mean_wealth}, \code{quantiles}, and \code{gini_consumption}.
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' ag  <- hank_asset_grid(100, 80, 0)
#' blk <- hank_het_block(ag, inc$Pi, inc$e, beta = 0.96, eis = 1, r = 0.02, w = 1)
#' hank_distribution_stats(blk)$gini
#' @export
hank_distribution_stats <- function(block, top = c(0.1, 0.01),
                                    probs = c(0.5, 0.9, 0.99),
                                    constraint_tol = 1e-9) {
  a_grid <- block$a_grid; n_e <- block$n_e; n_a <- block$n_a
  Dmat <- .hank_vec_to_mat(block$D, n_e, n_a)      # n_e x n_a
  a_mass <- colSums(Dmat)                          # marginal wealth distribution
  a_mass <- a_mass / sum(a_mass)

  gini <- hank_gini(a_grid, a_mass)

  ## Top wealth shares: order by assets descending, take the top fraction of
  ## population, report its share of total wealth.
  o <- order(a_grid, decreasing = TRUE)
  p_desc <- a_mass[o]; w_desc <- a_grid[o] * a_mass[o]
  cum_pop <- cumsum(p_desc)
  total_w <- sum(a_grid * a_mass)
  top_shares <- vapply(top, function(f) {
    ## interpolate the wealth share at the population cutoff f
    idx <- which(cum_pop >= f)[1]
    if (is.na(idx)) return(NA_real_)
    ## share of wealth held by the richest f of the population
    w_below <- if (idx > 1) sum(w_desc[seq_len(idx - 1)]) else 0
    frac_in <- (f - (if (idx > 1) cum_pop[idx - 1] else 0)) / p_desc[idx]
    (w_below + frac_in * w_desc[idx]) / total_w
  }, numeric(1))
  names(top_shares) <- paste0("top_", format(top * 100, trim = TRUE), "pct")

  ## Hand-to-mouth: mass with savings policy at the borrowing constraint.
  htm <- sum(block$D[.hank_mat_to_vec(block$a) <= a_grid[1L] + constraint_tol])

  ## Wealth quantiles from the marginal CDF.
  cdf <- cumsum(a_mass)
  quantiles <- vapply(probs, function(q) a_grid[which(cdf >= q)[1]], numeric(1))
  names(quantiles) <- paste0("p", format(probs * 100, trim = TRUE))

  ## Consumption Gini (consumption is nonnegative).
  c_vec <- .hank_mat_to_vec(block$c)
  gini_c <- hank_gini(c_vec, block$D)

  list(gini = gini, top_shares = top_shares, hand_to_mouth = htm,
       mean_wealth = total_w, quantiles = quantiles,
       gini_consumption = gini_c)
}


#' Sequence-space determinacy / well-posedness of a HANK model
#'
#' Reports the conditioning of the GE Jacobian \code{H_U}.  A well-posed
#' (locally determinate) linear rational-expectations equilibrium requires
#' \code{H_U} to be invertible; a near-singular \code{H_U} means
#' \code{\link{hank_model_irf}} would return an unreliable solution.
#'
#' @param model A \code{\link{hank_model}}.
#' @param tol Reciprocal-condition threshold below which \code{H_U} is flagged
#'   ill-conditioned (default 1e-10).
#'
#' @return A list with \code{rcond} (reciprocal condition number), \code{cond},
#'   \code{min_sv}/\code{max_sv} (extreme singular values), \code{determinate}
#'   (logical: \code{rcond >= tol}), and a \code{message}.
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' ag  <- hank_asset_grid(60, 40, 0)
#' nk  <- hank_nk_hank(ag, inc$Pi, inc$e, beta = 0.96, eis = 1,
#'                     r_ss = 0.005, phi = 1.5, kappa = 0.1, T_h = 80)
#' hank_determinacy(nk$model)$determinate
#' @export
hank_determinacy <- function(model, tol = 1e-10) {
  H <- model$H_U
  sv <- svd(H, nu = 0, nv = 0)$d
  min_sv <- min(sv); max_sv <- max(sv)
  rc <- if (max_sv > 0) min_sv / max_sv else 0
  determinate <- rc >= tol
  msg <- if (determinate)
    sprintf("H_U well-conditioned (rcond = %.2e): locally determinate.", rc)
  else
    sprintf(paste0("H_U ill-conditioned (rcond = %.2e < %.0e): the GE solution ",
                   "is unreliable (near-singular / indeterminate)."), rc, tol)
  list(rcond = rc, cond = if (min_sv > 0) max_sv / min_sv else Inf,
       min_sv = min_sv, max_sv = max_sv,
       determinate = determinate, message = msg)
}


#' Intertemporal MPCs (iMPC matrix) of a het block
#'
#' Computes the intertemporal marginal propensities to consume
#' \eqn{M_{t,s} = dC_t / d\tau_s} — the response of aggregate consumption at date
#' \eqn{t} to a one-time UNIFORM lump-sum income transfer at date \eqn{s} — via
#' the fake-news algorithm (the Auclert-Rognlie-Straub "intertemporal Keynesian
#' cross").  The impact column is the notch response to a transitory transfer;
#' \code{M[1, 1]} is the contemporaneous aggregate MPC.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param T_h Integer horizon.
#' @param delta_tr,delta_va,delta_d Finite-difference steps for the transfer
#'   shock, the value-function propagation, and the distributional response.
#'
#' @return A \code{T x T} matrix \code{M} with \code{[t, s] = dC_t/dtau_s}.
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' ag  <- hank_asset_grid(60, 60, 0)
#' blk <- hank_het_block(ag, inc$Pi, inc$e, beta = 0.96, eis = 1, r = 0.02, w = 1)
#' M <- hank_impc(blk, 40)
#' M[1, 1]   # contemporaneous aggregate MPC
#' @export
hank_impc <- function(block, T_h, delta_tr = 1e-5, delta_va = 1e-6,
                      delta_d = 1e-6) {
  a_grid <- block$a_grid; Pi <- block$Pi; D_ss <- block$D
  Va_ss <- block$Va; a_ss <- block$a; Lam <- block$Lambda
  y_ss <- block$w * block$e                      # steady-state income by state
  step <- function(Va, dtr) .hank_egm_step(Va, a_grid, y = y_ss + dtr,
                                           r = block$r, beta = block$beta,
                                           eis = block$eis, Pi = Pi)

  ## Expectation vectors for the consumption outcome: E_s = Lambda^s c_ss.
  E <- vector("list", T_h)
  E[[1L]] <- .hank_mat_to_vec(block$c)
  for (s in 2L:T_h) E[[s]] <- as.numeric(Lam %*% E[[s - 1L]])

  curlyD_from_dA <- function(dA) {
    Lp <- hank_forward_operator(a_ss + delta_d * dA, a_grid, Pi)
    Lm <- hank_forward_operator(a_ss - delta_d * dA, a_grid, Pi)
    (as.numeric(Matrix::t(Lp) %*% D_ss) -
       as.numeric(Matrix::t(Lm) %*% D_ss)) / (2 * delta_d)
  }

  curlyY <- numeric(T_h)
  curlyD <- matrix(0, block$n_e * block$n_a, T_h)

  ## s = 1: a lump-sum transfer NOW (adds delta_tr to every agent's income).
  sp <- step(Va_ss, delta_tr); sm <- step(Va_ss, -delta_tr)
  dA  <- (sp$a  - sm$a)  / (2 * delta_tr)
  dC  <- (sp$c  - sm$c)  / (2 * delta_tr)
  dVa <- (sp$Va - sm$Va) / (2 * delta_tr)
  curlyY[1L]  <- hank_aggregate(D_ss, dC)
  curlyD[, 1L] <- curlyD_from_dA(dA)

  ## s >= 2: anticipation of a future transfer, propagated via the value fn.
  dVa_prev <- dVa
  for (s in 2L:T_h) {
    h  <- delta_va / max(1, max(abs(dVa_prev)))
    sp <- .hank_egm_step(Va_ss + h * dVa_prev, a_grid, y = y_ss, r = block$r,
                         beta = block$beta, eis = block$eis, Pi = Pi)
    sm <- .hank_egm_step(Va_ss - h * dVa_prev, a_grid, y = y_ss, r = block$r,
                         beta = block$beta, eis = block$eis, Pi = Pi)
    dA  <- (sp$a  - sm$a)  / (2 * h)
    dC  <- (sp$c  - sm$c)  / (2 * h)
    dVa <- (sp$Va - sm$Va) / (2 * h)
    curlyY[s]  <- hank_aggregate(D_ss, dC)
    curlyD[, s] <- curlyD_from_dA(dA)
    dVa_prev <- dVa
  }

  ## Fake-news matrix F then iMPC matrix M (diagonal cumulative sum).
  Fm <- matrix(0, T_h, T_h); Fm[1L, ] <- curlyY
  for (tt in 2L:T_h) Fm[tt, ] <- as.numeric(crossprod(curlyD, E[[tt - 1L]]))
  M <- matrix(0, T_h, T_h); M[1L, ] <- Fm[1L, ]
  for (tt in 2L:T_h) {
    M[tt, 1L] <- Fm[tt, 1L]
    M[tt, 2L:T_h] <- M[tt - 1L, 1L:(T_h - 1L)] + Fm[tt, 2L:T_h]
  }
  M
}


#' Heatmap of a sequence-space Jacobian block
#'
#' @param Jblock A \code{T x T} Jacobian matrix (e.g. \code{J[["C"]][["r"]]}).
#' @param main Plot title.
#' @param ... Passed to \code{graphics::image}.
#' @return Invisibly, \code{Jblock}.
#' @export
hank_plot_jacobian <- function(Jblock, main = "Sequence-space Jacobian", ...) {
  Th <- nrow(Jblock)
  ## image() plots columns as x; transpose+flip so [t,s] reads t down, s across.
  graphics::image(x = seq_len(Th), y = seq_len(Th),
                  z = t(Jblock[Th:1L, , drop = FALSE]),
                  xlab = "shock date s", ylab = "response date t",
                  main = main, ...)
  invisible(Jblock)
}


#' Plot impulse-response paths from a linear GE IRF
#'
#' @param irf A list of named deviation paths (e.g. from
#'   \code{\link{hank_ks_linear_irf}}).
#' @param vars Character: which paths to plot (defaults to all \code{d*} paths).
#' @param main Plot title.
#' @return Invisibly, a matrix of the plotted paths (columns = variables).
#' @export
hank_plot_irf <- function(irf, vars = NULL,
                          main = "HANK GE impulse responses") {
  if (is.null(vars)) vars <- grep("^d", names(irf), value = TRUE)
  M <- sapply(vars, function(v) irf[[v]])
  graphics::matplot(M, type = "l", lty = 1, xlab = "period",
                    ylab = "deviation", main = main)
  graphics::legend("topright", legend = vars, col = seq_along(vars), lty = 1,
                   bty = "n")
  graphics::abline(h = 0, col = "grey70", lty = 3)
  invisible(M)
}
