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
  .hank_reject_het2(block, "hank_mpc", use = NULL)
  .hank_reject_wedge(block, "hank_mpc")
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
  ## Especially important here: without this the (e, b, a) distribution reaches
  ## .hank_vec_to_mat() and matrix() RECYCLES it with a warning rather than an
  ## error, so the failure is one arithmetic coincidence away from returning a
  ## plausible, wrong Gini.
  .hank_reject_het2(block, "hank_distribution_stats", use = NULL)
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


#' Sequence-space determinacy of a HANK model (winding-number criterion)
#'
#' Decides local determinacy of the linearized general equilibrium by the
#' winding-number criterion of Auclert, Rognlie and Straub ("Determinacy and
#' Existence in the Sequence Space", Prop. 1; ARS below) as extended to
#' the stacked (matrix-valued) GE Jacobian by Rion (Dynare WP 88, Prop. 7):
#' the untruncated stacked Jacobian \eqn{\tilde J_X} is block quasi-Toeplitz,
#' \eqn{\tilde J_X = T(j_X) + C}, and when \eqn{\det j_X(z) \neq 0} on the unit
#' circle it is Fredholm with
#' \deqn{\dim\ker \tilde J_X - \mathrm{codim\,ran}\, \tilde J_X =
#'   -\,\mathrm{wind}(\det j_X).}
#' The conditioning of the TRUNCATED \code{H_U} alone cannot decide this: a
#' \code{T x T} corner may be perfectly well conditioned while the infinite
#' system is indeterminate (e.g. a Taylor rule with \code{phi < 1}).
#'
#' \strong{Symbol.} Conventions follow both papers: \eqn{j_k = \lim_s
#' J_{s+k,s}} (date-\eqn{t} response to a date-\eqn{s} movement,
#' \eqn{k = t - s}) and \eqn{j(z) = \sum_k j_k z^k}, so the lag operator has
#' symbol \eqn{z} (winding \eqn{+1}) and the lead operator \eqn{z^{-1}}
#' (winding \eqn{-1}). For every (target, unknown) block of \code{H_U} the
#' coefficients are read at an interior date \eqn{\tau} (0-based) as in ARS
#' (2023 version, Sec. 5.1): \eqn{j_k = J_{\tau+k,\tau}} for \eqn{k \le 0} (column
#' \eqn{\tau}, anticipation) and \eqn{j_k = J_{\tau,\tau-k}} for \eqn{k \ge 0}
#' (row \eqn{\tau}, propagation). Every entry used lies in the leading
#' \eqn{(\tau+1)\times(\tau+1)} corner, far from the truncation artefacts
#' that composing truncated Jacobians leaves near date \code{T_h}. The default
#' is \eqn{\tau = \lfloor (T_h-1)/2 \rfloor}; a second estimate from
#' \eqn{\tau_2 = \lfloor 2\tau/3 \rfloor} is used as a convergence
#' diagnostic (\code{symbol_discrepancy}).
#'
#' \strong{Winding.} \eqn{\det j_X} is evaluated (FFT) on \code{n_grid}
#' counter-clockwise roots of unity; the winding number is the sum of the
#' principal-value phase increments over \eqn{2\pi} (unwrapped arg). The grid
#' is doubled (up to \eqn{2^{16}}) until no phase step exceeds \eqn{\pi/4}.
#'
#' \strong{Classification} (\code{status}):
#' \describe{
#'   \item{\code{"determinate"}}{\eqn{w = 0} and the truncated \code{H_U} is
#'     numerically nonsingular: the index is zero and the kernel is trivial,
#'     so a unique bounded response exists. \eqn{w = 0} is the generic
#'     certificate (invertible on an open dense set; Rion Prop. 7 item 3).}
#'   \item{\code{"indeterminate"}}{\eqn{w < 0}: at least \eqn{-w} independent
#'     square-summable solutions of the homogeneous system -- the bounded
#'     response is not unique.}
#'   \item{\code{"nonexistence"}}{\eqn{w > 0}: the range has codimension at
#'     least \eqn{w}; for some forcings no bounded response exists (whether
#'     THIS model's shocks hit an unreachable direction is a separate check,
#'     Rion Prop. 8).}
#'   \item{\code{"singular_index_zero"}}{\eqn{w = 0} but \code{H_U} is
#'     numerically singular (\code{rcond < tol}): the non-generic index-zero
#'     case with a nontrivial kernel, which couples indeterminacy with
#'     non-existence (Rion, discussion after Prop. 8).}
#'   \item{\code{"boundary"}}{\eqn{\det j_X} vanishes on the unit circle
#'     (\code{min_abs_det <= boundary_tol * max_abs_det}): the operator is not
#'     Fredholm and \code{winding} is \code{NA}.}
#' }
#' \code{determinate} is \code{TRUE} only for \code{"determinate"}.
#'
#' \strong{Convergence.} Rion (Prop. 10) shows the verdict read off an
#' approximate symbol is exact once the approximation error stays below the
#' determinant floor \eqn{\min_{|z|=1} |\det j_X|} (the curves cannot be
#' pulled across the origin). \code{converged} applies that argument between
#' the two readings: it is \code{TRUE} when
#' \eqn{|\det \hat j^{(\tau_2)}(z) - \det \hat j^{(\tau)}(z)| <
#' |\det \hat j^{(\tau)}(z)|} at every grid point, which by Rouche's theorem
#' forces both readings to wind identically. \code{FALSE} means the symbol has
#' not visibly converged relative to the determinant floor: increase
#' \code{T_h}. It does not change \code{determinate}.
#'
#' @param model A \code{\link{hank_model}} (anything with a square stacked
#'   \code{H_U}; \code{T_h} and \code{unknowns} are read when present,
#'   otherwise one unknown with \code{T_h = nrow(H_U)} is assumed).
#' @param tol Reciprocal-condition threshold below which the truncated
#'   \code{H_U} is flagged numerically singular (default 1e-10). It no longer
#'   decides \code{determinate} on its own: it only separates
#'   \code{"determinate"} from \code{"singular_index_zero"} when \eqn{w = 0}.
#' @param tau Optional 0-based interior date at which the symbol is read
#'   (default \code{floor((T_h - 1) / 2)}).
#' @param n_grid Integer: initial number of unit-circle points (default
#'   4096; raised to exceed \eqn{4\tau} and doubled adaptively).
#' @param boundary_tol Relative floor for \eqn{\min |\det j_X|} below which the
#'   symbol is treated as vanishing on the circle (default 1e-8).
#'
#' @return A list with the headline \code{determinate} (logical, from the
#'   winding number -- see above), \code{status}, \code{winding} (integer, or
#'   \code{NA} on the boundary), \code{index} (\eqn{= -w}),
#'   \code{min_abs_det}/\code{max_abs_det} over the circle, \code{converged},
#'   \code{symbol_discrepancy} (\eqn{\sup_z \|\hat j^{(\tau)} -
#'   \hat j^{(\tau_2)}\|_F} relative to \eqn{\sup_z \|\hat j^{(\tau)}\|_F}), \code{symbol_tail} (largest
#'   \eqn{|j_k|} in the outer tenth of the band, relative), \code{tau},
#'   \code{n_grid} (final), \code{symbol} (\code{n x n x (2 tau + 1)} array of
#'   the blocks \eqn{j_k}, \eqn{k = -\tau..\tau}); the numerical-conditioning
#'   fields of the truncated \code{H_U}: \code{rcond}, \code{cond},
#'   \code{min_sv}/\code{max_sv}, \code{well_conditioned}
#'   (\code{rcond >= tol}); and a \code{message}.
#' @references
#' Auclert, A., Rognlie, M., Straub, L. (2023; cited by Rion as 2025).
#' Determinacy and Existence in the Sequence Space. Working paper.
#'
#' Rion, N. (2026). Rich Heterogeneity in Dynare 7: A Practical Description.
#' Dynare Working Paper 88.
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' ag  <- hank_asset_grid(60, 40, 0)
#' nk  <- hank_nk_hank(ag, inc$Pi, inc$e, beta = 0.96, eis = 1,
#'                     r_ss = 0.005, phi = 1.5, kappa = 0.1, T_h = 80)
#' d <- hank_determinacy(nk$model)
#' d$determinate
#' d$winding
#' @export
hank_determinacy <- function(model, tol = 1e-10, tau = NULL, n_grid = 4096L,
                             boundary_tol = 1e-8) {
  H <- model$H_U
  if (!is.matrix(H) || nrow(H) != ncol(H))
    .dynhr_abort("hank_determinacy(): `model$H_U` must be a square matrix.",
                 class = "dynhr_error_input")
  n_u <- if (is.null(model$unknowns)) 1L else length(model$unknowns)
  T_h <- if (is.null(model$T_h)) nrow(H) %/% n_u else as.integer(model$T_h)
  if (n_u * T_h != nrow(H))
    .dynhr_abort(sprintf(paste0("hank_determinacy(): H_U is %d x %d, not ",
                                "(n_unknowns * T_h) = %d square."),
                         nrow(H), ncol(H), n_u * T_h),
                 class = "dynhr_error_input")
  if (T_h < 7L)
    .dynhr_abort("hank_determinacy(): T_h must be at least 7 to read an ",
                 "interior symbol.", class = "dynhr_error_input")

  ## ---- secondary: conditioning of the truncated corner --------------------
  sv <- svd(H, nu = 0, nv = 0)$d
  min_sv <- min(sv); max_sv <- max(sv)
  rc <- if (max_sv > 0) min_sv / max_sv else 0
  well_conditioned <- is.finite(rc) && rc >= tol

  ## ---- headline: winding number of det j_X --------------------------------
  if (is.null(tau)) tau <- (T_h - 1L) %/% 2L
  tau <- as.integer(tau)
  if (tau < 3L || tau > T_h - 1L)
    .dynhr_abort(sprintf("hank_determinacy(): `tau` must lie in [3, %d].",
                         T_h - 1L), class = "dynhr_error_input")
  tau2 <- max(2L, (tau * 2L) %/% 3L)
  j1 <- .hank_symbol_blocks(H, n_u, T_h, tau)
  j2 <- .hank_symbol_blocks(H, n_u, T_h, tau2)

  N <- as.integer(max(n_grid, 2^ceiling(log2(4 * tau + 2))))
  repeat {
    Jz <- .hank_symbol_eval(j1, tau, N)
    d  <- .hank_cdet_grid(Jz)
    abs_d <- Mod(d)
    min_abs <- min(abs_d); max_abs <- max(abs_d)
    on_boundary <- !is.finite(min_abs) || max_abs <= 0 ||
      min_abs <= boundary_tol * max_abs
    if (on_boundary) break
    ## phase increments between consecutive counter-clockwise nodes, closing
    ## the loop; principal values are exact while |step| < pi
    dphi <- Arg(d[c(seq_len(N)[-1L], 1L)] / d)
    if (max(abs(dphi)) <= pi / 4 || N >= 65536L) break
    N <- 2L * N
  }
  w <- if (on_boundary) NA_integer_ else as.integer(round(sum(dphi) / (2 * pi)))

  ## ---- convergence diagnostic: a second reading date ----------------------
  ## The verdict is trusted when the symbol read at tau2 winds identically,
  ## certified pointwise by Rouche: |det j2 - det j1| < |det j1| on the whole
  ## circle forces wind(det j2) = wind(det j1) (Rion Prop. 10's argument,
  ## with the second reading standing in for the unknown exact symbol).
  Jz2 <- .hank_symbol_eval(j2, tau2, N)
  d2  <- .hank_cdet_grid(Jz2)
  disc_abs <- max(sqrt(apply(Mod(Jz - Jz2)^2, 1L, sum)))
  scale_j <- max(sqrt(apply(Mod(Jz)^2, 1L, sum)))
  sym_disc <- if (scale_j > 0) disc_abs / scale_j else Inf
  converged <- !on_boundary && all(Mod(d2 - d) < abs_d)
  kk <- seq(-tau, tau)
  outer_band <- abs(kk) >= ceiling(0.9 * tau)
  jmax <- max(abs(j1))
  sym_tail <- if (jmax > 0) max(abs(j1[, , outer_band])) / jmax else 0

  status <- if (on_boundary) "boundary"
    else if (w < 0L) "indeterminate"
    else if (w > 0L) "nonexistence"
    else if (!well_conditioned) "singular_index_zero"
    else "determinate"
  determinate <- identical(status, "determinate")

  msg <- switch(status,
    boundary = sprintf(paste0(
      "det j_X(z) vanishes on the unit circle (min |det| = %.2e, relative ",
      "%.2e): the stacked GE Jacobian is not Fredholm -- boundary between ",
      "regimes; no determinacy verdict."), min_abs,
      if (max_abs > 0) min_abs / max_abs else NA_real_),
    indeterminate = sprintf(paste0(
      "winding(det j_X) = %d < 0: at least %d independent bounded solution(s) ",
      "of the homogeneous system -- INDETERMINATE (the truncated H_U has ",
      "rcond = %.2e; the corner hides this)."), w, -w, rc),
    nonexistence = sprintf(paste0(
      "winding(det j_X) = %d > 0: the range has codimension >= %d, so for ",
      "some forcings no bounded response exists (non-existence risk; whether ",
      "this model's shocks hit an unreachable direction is not checked)."),
      w, w),
    singular_index_zero = sprintf(paste0(
      "winding(det j_X) = 0 but the truncated H_U is numerically singular ",
      "(rcond = %.2e < %.0e): the non-generic index-zero case with a ",
      "nontrivial kernel (indeterminacy together with non-existence)."),
      rc, tol),
    determinate = sprintf(paste0(
      "winding(det j_X) = 0 and H_U nonsingular (rcond = %.2e): locally ",
      "determinate (index zero, trivial kernel)."), rc))
  if (!converged && !on_boundary)
    msg <- paste0(msg, sprintf(paste0(
      " Not converged: the symbol read at date %d does not provably wind ",
      "like the one read at date %d (|det j_X| floor %.2e vs discrepancy); ",
      "increase T_h."), tau2, tau, min_abs))

  list(determinate = determinate, status = status, winding = w,
       index = if (is.na(w)) NA_integer_ else -w,
       min_abs_det = min_abs, max_abs_det = max_abs,
       converged = converged, symbol_discrepancy = sym_disc,
       symbol_tail = sym_tail, tau = tau, n_grid = N, symbol = j1,
       rcond = rc, cond = if (min_sv > 0) max_sv / min_sv else Inf,
       min_sv = min_sv, max_sv = max_sv, well_conditioned = well_conditioned,
       message = msg)
}


#' Asymptotic Toeplitz symbol blocks of a stacked sequence-space Jacobian
#'
#' Reads \eqn{j_k} for every (row-block, column-block) pair of the stacked
#' \code{(n T) x (n T)} matrix \code{H} at the 0-based interior date
#' \code{tau}: \eqn{j_k = H_{\tau+k,\tau}} for \eqn{k \le 0},
#' \eqn{j_k = H_{\tau,\tau-k}} for \eqn{k \ge 0} (ARS 2023 version, Sec. 5.1).
#'
#' @param H Stacked square matrix.
#' @param n Number of blocks per side.
#' @param T_h Block size (horizon).
#' @param tau 0-based interior date.
#' @return Real array \code{n x n x (2 tau + 1)}; slice \code{k + tau + 1}
#'   holds \eqn{j_k}.
#' @keywords internal
.hank_symbol_blocks <- function(H, n, T_h, tau) {
  out <- array(0, c(n, n, 2L * tau + 1L))
  k_neg <- seq(-tau, 0L)                  # anticipation: column tau
  k_pos <- seq_len(tau)                   # propagation: row tau
  for (a in seq_len(n)) for (b in seq_len(n)) {
    r0 <- (a - 1L) * T_h; c0 <- (b - 1L) * T_h
    out[a, b, k_neg + tau + 1L] <- H[r0 + tau + k_neg + 1L, c0 + tau + 1L]
    out[a, b, k_pos + tau + 1L] <- H[r0 + tau + 1L, c0 + tau - k_pos + 1L]
  }
  out
}


#' Evaluate a matrix symbol on the unit circle by FFT
#'
#' @param jk Array \code{n x n x (2 tau + 1)} from
#'   \code{\link{.hank_symbol_blocks}}.
#' @param tau Band half-width.
#' @param N Number of nodes \eqn{z_m = e^{2\pi i m/N}}, \eqn{m = 0..N-1}
#'   (counter-clockwise); needs \code{N > 2 tau}.
#' @return Complex array \code{N x n x n} with \code{[m + 1, , ] = j(z_m)}.
#' @keywords internal
.hank_symbol_eval <- function(jk, tau, N) {
  n <- dim(jk)[1L]
  idx <- (seq(-tau, tau) %% N) + 1L       # j_k at position k mod N
  out <- array(0 + 0i, c(N, n, n))
  for (a in seq_len(n)) for (b in seq_len(n)) {
    x <- numeric(N)
    x[idx] <- jk[a, b, ]
    ## unnormalized inverse DFT: sum_k j_k exp(+2 pi i k m / N) = j(z_m)
    out[, a, b] <- stats::fft(x, inverse = TRUE)
  }
  out
}


#' Determinant of a complex matrix at every grid point
#'
#' Closed form for \code{n <= 2}; Gaussian elimination with partial pivoting
#' otherwise (base \code{det()} does not accept complex input).
#'
#' @param Jz Complex array \code{N x n x n}.
#' @return Complex vector of length \code{N}.
#' @keywords internal
.hank_cdet_grid <- function(Jz) {
  n <- dim(Jz)[2L]
  if (n == 1L) return(Jz[, 1L, 1L])
  if (n == 2L) return(Jz[, 1L, 1L] * Jz[, 2L, 2L] - Jz[, 1L, 2L] * Jz[, 2L, 1L])
  vapply(seq_len(dim(Jz)[1L]), function(m) {
    A <- Jz[m, , ]
    d <- 1 + 0i
    for (col in seq_len(n)) {
      p <- which.max(Mod(A[col:n, col])) + col - 1L
      if (Mod(A[p, col]) == 0) return(0 + 0i)
      if (p != col) { A[c(col, p), ] <- A[c(p, col), ]; d <- -d }
      d <- d * A[col, col]
      if (col < n) {
        rows <- (col + 1L):n
        f <- A[rows, col] / A[col, col]
        A[rows, ] <- A[rows, , drop = FALSE] - outer(f, A[col, ])
      }
    }
    d
  }, complex(1))
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
  .hank_reject_het2(block, "hank_impc", use = NULL)
  .hank_reject_wedge(block, "hank_impc")
  a_grid <- block$a_grid; Pi <- block$Pi; D_ss <- block$D
  Va_ss <- block$Va; a_ss <- block$a; Lam <- block$Lambda
  ## Steady-state income by state, INCLUDING any lump-sum transfer the block
  ## was built with -- otherwise the iMPC of a Tr > 0 block would be
  ## differenced around the wrong baseline.
  ##
  ## The notch is scaled by the block's incidence weight, so what this returns
  ## is the iMPC out of a transfer distributed by THAT rule. That keeps it the
  ## correct independent cross-oracle for the "Tr" column of
  ## hank_het_jacobian() under any incidence, and it reduces to the uniform
  ## iMPC (omega == 1) byte-for-byte.
  omega <- .hank_block_omega(block)
  y_ss <- block$w * block$e + .hank_block_tr(block) * omega
  step <- function(Va, dtr) .hank_egm_step(Va, a_grid, y = y_ss + dtr * omega,
                                           r = block$r, beta = block$beta,
                                           eis = block$eis, Pi = Pi)

  ## Expectation vectors for the consumption outcome: E_s = Lambda^s c_ss.
  E <- vector("list", T_h)
  E[[1L]] <- .hank_mat_to_vec(block$c)
  ## s = 2 .. T_h, empty if T_h < 2
  for (s in seq_len(T_h - 1L) + 1L) E[[s]] <- as.numeric(Lam %*% E[[s - 1L]])

  ## Matrix-free, and it must STAY matrix-free in lockstep with
  ## .hank_curly_sweep: test-hank-transfer.R asserts J[C][Tr] equals this
  ## iMPC matrix to 1e-12, which holds only because the two routes difference
  ## the SAME steps around the SAME baseline. Round-off differences here are
  ## amplified by the 1/(2*delta_d) division to ~1e-10, so a mismatched pair
  ## of implementations breaks that cross-oracle (it did, exactly once).
  curlyD_from_dA <- function(dA)
    (.hank_forward_push(a_ss + delta_d * dA, a_grid, Pi, D_ss) -
       .hank_forward_push(a_ss - delta_d * dA, a_grid, Pi, D_ss)) / (2 * delta_d)

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
  for (s in seq_len(T_h - 1L) + 1L) {              # s = 2 .. T_h, empty if T_h < 2
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
  ## tt = 2 .. T_h, empty if T_h < 2
  for (tt in seq_len(T_h - 1L) + 1L)
    Fm[tt, ] <- as.numeric(crossprod(curlyD, E[[tt - 1L]]))
  M <- matrix(0, T_h, T_h); M[1L, ] <- Fm[1L, ]
  for (tt in seq_len(T_h - 1L) + 1L) {            # tt = 2 .. T_h, empty if T_h < 2
    M[tt, 1L] <- Fm[tt, 1L]
    ## Body only runs when T_h >= 2, so the 2L:T_h slices are in range here.
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
#' The two IRF producers in this package name their paths differently:
#' \code{\link{hank_ks_linear_irf}} returns \code{d}-prefixed names
#' (\code{dK}, \code{dC}, ...), while the general \code{\link{hank_model_irf}}
#' returns the model's own bare variable names (\code{K}, \code{C}, ...). The
#' default therefore takes the \code{d*} paths when there are any and every
#' path otherwise, and a \code{vars} entry that is not in \code{irf} is an
#' error naming what IS available -- previously it silently produced a list
#' column and died inside \code{matplot} with "'list' object cannot be coerced
#' to type 'double'".
#'
#' @param irf A list of named deviation paths (e.g. from
#'   \code{\link{hank_ks_linear_irf}} or \code{\link{hank_model_irf}}).
#' @param vars Character: which paths to plot. Defaults to all \code{d*} paths
#'   if the IRF has any, otherwise to every path of the right length.
#' @param main Plot title.
#' @return Invisibly, a matrix of the plotted paths (columns = variables).
#' @export
hank_plot_irf <- function(irf, vars = NULL,
                          main = "HANK GE impulse responses") {
  paths <- names(irf)[vapply(irf, function(v)
    is.numeric(v) && length(v) > 1L, logical(1))]
  if (is.null(vars)) {
    vars <- grep("^d", paths, value = TRUE)
    if (!length(vars)) vars <- paths
  }
  if (!length(vars))
    stop("hank_plot_irf(): the IRF has no numeric deviation paths to plot.")
  bad <- setdiff(vars, paths)
  if (length(bad))
    stop("hank_plot_irf(): no such path(s) in this IRF: ",
         paste0("'", bad, "'", collapse = ", "), ". Available: ",
         paste0("'", paths, "'", collapse = ", "), ".")
  M <- sapply(vars, function(v) irf[[v]])
  graphics::matplot(M, type = "l", lty = 1, xlab = "period",
                    ylab = "deviation", main = main)
  graphics::legend("topright", legend = vars, col = seq_along(vars), lty = 1,
                   bty = "n")
  graphics::abline(h = 0, col = "grey70", lty = 3)
  invisible(M)
}


#' Preflight validation of a solved HANK household block
#'
#' One cheap, interpretable gate to run BEFORE spending a costly likelihood,
#' Jacobian, or posterior evaluation on a \code{\link{hank_het_block}}
#' (adversarial review 2026-07-13, extension #3). Reports -- rather than
#' silently assumes -- the invariants every downstream consumer relies on:
#' grid monotonicity, the Markov contract on \code{Pi}, forward-operator
#' row-stochasticity, distribution mass/positivity/stationarity, policy
#' feasibility against the block's own borrowing limit \code{amin},
#' consumption positivity, the stored-aggregate identities
#' \code{A == hank_aggregate(D, a)} / \code{C == hank_aggregate(D, c)}
#' (catching hand-spliced blocks that edit policies without re-aggregating),
#' and stationary-distribution convergence. Optionally re-solves the EGM on
#' BOTH backends and reports their policy gap (\code{check_parity}).
#'
#' A LIST of blocks (e.g. \code{mks$blocks} from a mixture steady state) is
#' validated per type, with a \code{type} column prepended.
#'
#' @param block A \code{\link{hank_het_block}}, or a list of them.
#' @param check_parity Logical: additionally re-solve the household EGM under
#'   both \code{backend = "R"} and \code{"cpp"} and report the max absolute
#'   savings-policy gap (costs two EGM solves; default \code{FALSE}).
#' @param tol Numerical tolerance for the residual checks (default
#'   \code{1e-8}); the stationarity residual uses \code{max(tol, 1e-10)}.
#'
#' @return A data frame with one row per check: \code{check}, \code{value}
#'   (the measured residual, or \code{NA} for boolean checks), \code{threshold},
#'   \code{pass}. The attribute \code{"ok"} is \code{TRUE} iff every row
#'   passes. All checks are evaluated even when earlier ones fail, so a
#'   corrupted block yields the full damage report.
#' @seealso \code{\link{hank_het_block}}, \code{\link{hank_theta_boundary_check}},
#'   \code{\link{hank_mixture_ks_assemble}}
#' @export
validate_hank_block <- function(block, check_parity = FALSE, tol = 1e-8) {
  if (!inherits(block, "hank_het_block") && is.list(block) &&
      length(block) >= 1L &&
      all(vapply(block, inherits, logical(1), "hank_het_block"))) {
    per <- lapply(seq_along(block), function(k) {
      r <- validate_hank_block(block[[k]], check_parity = check_parity,
                               tol = tol)
      cbind(data.frame(type = k), r)
    })
    out <- do.call(rbind, per)
    rownames(out) <- NULL
    attr(out, "ok") <- all(out$pass)
    return(out)
  }
  if (inherits(block, "hank_het2_block"))
    stop("validate_hank_block(): this is a two-asset block ",
         "(hank_het2_block), whose policies are n_e x n_b x n_a arrays over a ",
         "joint (e, b, a) cell space -- the checks here assume the one-asset ",
         "n_e x n_a shape and its single borrowing limit, so they would be ",
         "meaningless rather than merely wrong. No two-asset preflight exists ",
         "yet; validate the block's own invariants directly (Lambda ",
         "row-stochastic, D stationary, B/A/C/CHI vs hank_aggregate2).")
  if (!inherits(block, "hank_het_block"))
    stop("validate_hank_block(): 'block' must be a hank_het_block ",
         "or a list of them.")

  ag <- block$a_grid
  rows <- list()
  row <- function(check, value, threshold, pass)
    data.frame(check = check, value = value, threshold = threshold,
               pass = pass)

  ## grid: finite, strictly increasing
  grid_ok <- is.numeric(ag) && all(is.finite(ag)) &&
    (length(ag) < 2L || all(diff(ag) > 0))
  rows$grid <- row("grid_strictly_increasing", NA_real_, NA_real_, grid_ok)

  ## Markov contract on Pi (same validator as the public boundaries)
  pi_ok <- tryCatch({
    .hank_check_markov(block$Pi, block$n_e, caller = "validate_hank_block",
                       tol = tol)
    TRUE
  }, error = function(e) FALSE)
  pi_res <- if (is.matrix(block$Pi) && is.numeric(block$Pi))
    max(abs(rowSums(block$Pi) - 1), -min(block$Pi, 0)) else NA_real_
  rows$markov <- row("Pi_markov_contract", pi_res, tol, pi_ok)

  ## forward operator row-stochastic
  lam_res <- max(abs(Matrix::rowSums(block$Lambda) - 1))
  rows$lambda <- row("Lambda_row_stochastic", lam_res, tol, lam_res <= tol)

  ## distribution: shape, positivity, unit mass, stationarity
  n_cell <- block$n_e * block$n_a
  d_shape <- length(block$D) == n_cell && all(is.finite(block$D))
  rows$dshape <- row("D_length_and_finite", NA_real_, NA_real_, d_shape)
  d_neg <- if (d_shape) max(0, -min(block$D)) else NA_real_
  rows$dneg <- row("D_nonnegative", d_neg, tol,
                   isTRUE(d_neg <= tol))
  d_mass <- if (d_shape) abs(sum(block$D) - 1) else NA_real_
  rows$dmass <- row("D_unit_mass", d_mass, tol, isTRUE(d_mass <= tol))
  stat_tol <- max(tol, 1e-10)
  d_stat <- if (d_shape)
    max(abs(as.numeric(Matrix::t(block$Lambda) %*% block$D) - block$D))
  else NA_real_
  rows$dstat <- row("D_stationary_residual", d_stat, stat_tol,
                    isTRUE(d_stat <= stat_tol))

  ## policies: finite, feasible against the block's own amin, c > 0
  amin <- if (!is.null(block$amin)) block$amin else ag[1L]
  pol_fin <- all(is.finite(block$a)) && all(is.finite(block$c))
  rows$pfin <- row("policies_finite", NA_real_, NA_real_, pol_fin)
  a_feas <- if (pol_fin) max(0, amin - min(block$a)) else NA_real_
  rows$afeas <- row("a_policy_respects_amin", a_feas, tol,
                    isTRUE(a_feas <= tol))
  c_min <- if (pol_fin) min(block$c) else NA_real_
  rows$cpos <- row("c_policy_positive", c_min, 0,
                   isTRUE(c_min > 0))

  ## stored aggregates match distribution-weighted policies
  agg_A <- tryCatch(abs(block$A - hank_aggregate(block$D, block$a)),
                    error = function(e) NA_real_)
  agg_C <- tryCatch(abs(block$C - hank_aggregate(block$D, block$c)),
                    error = function(e) NA_real_)
  rows$aggA <- row("A_matches_aggregated_policy", agg_A, tol,
                   isTRUE(agg_A <= tol))
  rows$aggC <- row("C_matches_aggregated_policy", agg_C, tol,
                   isTRUE(agg_C <= tol))

  ## stationary-distribution solver convergence flag
  conv <- isTRUE(block$dist_converged)
  rows$conv <- row("stationary_dist_converged", NA_real_, NA_real_, conv)

  ## optional R/C++ backend parity on a fresh EGM re-solve
  if (isTRUE(check_parity)) {
    gap <- tryCatch({
      hhR <- hank_egm_solve(ag, y = block$w * block$e, r = block$r,
                            beta = block$beta, eis = block$eis,
                            Pi = block$Pi, amin = amin, backend = "R")
      hhC <- hank_egm_solve(ag, y = block$w * block$e, r = block$r,
                            beta = block$beta, eis = block$eis,
                            Pi = block$Pi, amin = amin, backend = "cpp")
      max(abs(hhR$a - hhC$a))
    }, error = function(e) NA_real_)
    rows$parity <- row("egm_backend_parity", gap, 1e-9,
                       isTRUE(gap <= 1e-9))
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "ok") <- all(out$pass)
  out
}
