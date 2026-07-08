## R/hank-partial-id.R
## --------------------------------------------------------------------------
## Partial-identification frontier for the het-beta HANK mixture's discount
## SPREAD: a REDO of an earlier robustness-efficiency study whose
## contamination model was wrong (a MULTIPLICATIVE distortion, which
## perturbs the net cross-section reweighting too, so it never showcased
## the response channel's robustness niche and misleadingly found the LEVEL
## channel dominates everywhere).
##
## The CORRECT contamination model here is a FIXED, mean-zero, ADDITIVE
## cross-sectional distortion \code{f} on the asset cells (a date-invariant,
## parameter-invariant wealth-survey misreporting pattern), applied
## IDENTICALLY to the pre-shock and post-shock cross-sections:
##   D0_obs = D0 + eps*f,   D1_obs = D1 + eps*f,   sum(f) == 0.
## Because the RESPONSE channel's observable is the DIFFERENCE of the two
## cross-sections, the additive f cancels EXACTLY:
##   dD_obs = D1_obs - D0_obs = (D1 + eps f) - (D0 + eps f) = D1 - D0 = dD,
## i.e. bias_response(eps) === 0 for ANY f and ANY eps -- the response is
## EXACTLY robust to this class of contamination. The LEVEL channel's
## observable is the single (uncancelled) cross-section D0_obs, so its
## implied spread estimate IS biased by eps*f, with the bias magnitude
## depending on how much f aligns with the spread's identifying direction
## in D0 (see hank_partial_id_contamination_f()'s "adversarial" choice).
##
## Two channels, reusing R/hank-reweighting.R's grid-invariant functional
## metrics (hank_reweight_level_metric / hank_reweight_functional_metric)
## and the SAME cell-space convention hank_mixture_joint_logpost()
## (R/hank-mixture-estimation.R) uses for its m_lev_hat / dD_hat channels:
##   - LEVEL:    the stationary wealth cross-section D0 (sharp but fragile
##               under this contamination).
##   - RESPONSE: the net cross-section reweighting dD = D1 - D0 under a
##               price shock (weak but EXACTLY robust to this contamination).
##
## hank_partial_id_level_response() bundles: (1) the two channels' Fisher
## information for the spread (g' Sigma^{-1} g by central FD, both
## verified step-robust); (2) the bias each channel suffers under a supplied
## additive f across an eps grid; (3) the resulting honest RMSE frontier
## RMSE(eps, n) = sqrt(bias(eps)^2 + 1/(I*n)) and its crossover -- the eps
## (or n) at which the robust-but-weak response overtakes the sharp-but-
## fragile level.
## --------------------------------------------------------------------------


#' One shared GE solve: stationary distribution + snapshot reweighting
#'
#' Internal helper shared by \code{\link{hank_partial_id_level_response}}'s
#' Fisher-information and bias computations: solves the mixture GE steady
#' state at \code{(centre, spread, omega)} and returns both the stationary
#' distribution \code{D0} (the LEVEL channel's object) and the date-
#' \code{t_star} snapshot distribution-Jacobian response to \code{rw_shock}
#' (the RESPONSE channel's object, BEFORE the \code{rw_scale} multiplier is
#' applied -- mirrors \code{\link{hank_mixture_joint_logpost}}'s internal
#' \code{dD_mod <- config$rw_scale * hank_dist_response_snapshot(...)}
#' construction exactly).
#'
#' @param a_grid,Pi,e Household asset grid, income transition, income levels.
#' @param centre,spread Mixture centre/spread (\code{betas = centre + spread
#'   * beta_shape}).
#' @param omega Length-K mixture weights.
#' @param beta_shape Length-K heterogeneity pattern (see
#'   \code{\link{hank_partial_id_level_response}}).
#' @param eis,alpha,delta,Z GE calibration (see
#'   \code{\link{hank_mixture_ks_steady}}).
#' @param T_h,rw_shock,t_star Distribution-Jacobian horizon, shock path, and
#'   snapshot date (see \code{\link{hank_dist_response_snapshot}}).
#'
#' @return A list with \code{D0} (length-\code{n_cell} stationary
#'   distribution) and \code{dD_snap} (length-\code{n_cell} UNSCALED
#'   snapshot reweighting; multiply by \code{rw_scale} before use).
#' @keywords internal
.hank_partial_id_ge <- function(a_grid, Pi, e, centre, spread, omega,
                                 beta_shape, eis, alpha, delta, Z, T_h,
                                 rw_shock, t_star) {
  betas <- centre + spread * beta_shape
  mks <- hank_mixture_ks_steady(a_grid, Pi, e, betas, omega, eis = eis,
                                alpha = alpha, delta = delta, Z = Z)
  b <- hank_mixture_blocks(a_grid, Pi, e, betas = betas, eis = eis,
                           r = mks$r, w = mks$w)
  D0 <- hank_mixture_dist(b, omega)$D
  JD <- hank_mixture_dist_jacobian(b, omega, T_h, inputs = "r")
  dD_snap <- hank_dist_response_snapshot(JD, rw_shock, t_star)
  list(D0 = D0, dD_snap = dD_snap)
}


#' Central-FD Fisher information of a channel's spread moment
#'
#' \eqn{I = g' \Sigma^{-1} g}, \code{g = d(moment)/d(spread)} by CENTRAL
#' finite difference, holding \code{centre} and \code{omega} fixed at the
#' evaluation point. Shared by the LEVEL (\code{moment_fn} projects the
#' stationary distribution) and RESPONSE (\code{moment_fn} projects the
#' scaled snapshot reweighting) channels in
#' \code{\link{hank_partial_id_level_response}}.
#'
#' @param moment_fn Function of a single argument \code{spread} returning
#'   the length-\code{K} projected moment (\code{t(B) \%*\% <cell object>}).
#' @param Sigma_inv The channel metric's \code{K x K} GLS precision (e.g.
#'   \code{hank_reweight_level_metric(...)$Sigma_inv}).
#' @param spread0 Evaluation point.
#' @param h Central-FD step size.
#'
#' @return A list with \code{g} (the length-\code{K} gradient) and \code{I}
#'   (the scalar Fisher information).
#' @keywords internal
.hank_partial_id_fisher <- function(moment_fn, Sigma_inv, spread0, h) {
  g <- (moment_fn(spread0 + h) - moment_fn(spread0 - h)) / (2 * h)
  list(g = g, I = as.numeric(crossprod(g, Sigma_inv %*% g)))
}


#' Canonical fixed additive cross-sectional contamination directions
#'
#' Builds one of three PINNED choices of a fixed, mean-zero, additive
#' cross-sectional distortion \code{f} on the distribution cells (see the
#' file header comment for the contamination model this supports), each
#' normalized so \code{0.5 * sum(abs(f)) <= 1} (an L1/mass normalization
#' under which \code{eps * f} moves AT MOST a fraction \code{eps} of total
#' unit mass -- \code{eps = 0.10} reads as "at most 10% of mass distorted"),
#' and then uniformly shrunk (same direction/shape, just rescaled) so that
#' \strong{both} \code{D0 + eps_max*f >= 0} and \code{D1 + eps_max*f >= 0}
#' cell-by-cell, guaranteeing every observed cross-section in the eps range
#' \code{[0, eps_max]} stays a valid (non-negative) distribution. The
#' realized \code{0.5 * sum(abs(f))} after this safety shrink is reported
#' (strictly less than 1 whenever nonnegativity binds on a thin cell --
#' e.g. a household grid's thin extreme-wealth cells).
#'
#' \describe{
#'   \item{\code{"random"}}{A fixed mean-zero smooth-ish random cell
#'     distortion: iid noise, mildly smoothed across cells sorted by asset
#'     level, scaled proportional to each cell's own \code{D0} mass (the
#'     natural "misreporting" shape), then mean-zeroed. Seeded by
#'     \code{seed}.}
#'   \item{\code{"adversarial"}}{The WORST-CASE direction: proportional to
#'     \code{d(D0)/d(spread)} (by central FD) -- the cell-space SIGNATURE of
#'     a spread change, i.e. the fixed distortion that is maximally
#'     confusable with a genuine spread shift under the LEVEL channel's own
#'     loss. Winsorized at the 99th percentile of \code{|d(D0)/d(spread)|}
#'     before normalizing, so a single extreme cell (typically the
#'     borrowing-constraint cell, where patient/impatient households pile up
#'     most differently) does not force an excessive safety shrink.}
#'   \item{\code{"bracket"}}{A plausible survey misreport: moves half of one
#'     asset-wealth DECILE's own per-cell mass into the adjacent decile
#'     (mass-shaped within the receiving decile), a stylized "some
#'     respondents round their decile up" distortion.}
#' }
#'
#' @param a_cell Length-\code{n_cell} asset level of each distribution cell
#'   (asset-fast cell order; see \code{\link{hank_mixture_dist}}).
#' @param D0,D1 Length-\code{n_cell} pre-/post-shock TRUTH distributions
#'   (non-negative, summing to 1) the safety rescale is calibrated against.
#' @param which One of \code{"random"}, \code{"adversarial"},
#'   \code{"bracket"}.
#' @param eps_max The largest eps the safety rescale must keep non-negative
#'   for (default \code{0.10}).
#' @param d_D0_d_spread Only used when \code{which == "adversarial"}: the
#'   length-\code{n_cell} vector \code{d(D0)/d(spread)} by central FD (the
#'   caller supplies this since it requires a GE re-solve at
#'   \code{spread0 +- h}, which \code{\link{hank_partial_id_level_response}}
#'   already computes for its Fisher-information step).
#' @param seed RNG seed, only used when \code{which == "random"}.
#' @param margin Safety-rescale margin in \code{(0, 1]} (default \code{0.98}):
#'   leaves headroom below the exact zero-crossing rather than pinning the
#'   single thinnest cell to exactly 0.
#'
#' @return A length-\code{n_cell} numeric vector \code{f}, with
#'   \code{sum(f) == 0} (to numerical precision) and
#'   \code{D0 + eps_max*f >= 0}, \code{D1 + eps_max*f >= 0} cell-by-cell.
#' @seealso \code{\link{hank_partial_id_level_response}}
#' @export
hank_partial_id_contamination_f <- function(a_cell, D0, D1,
                                             which = c("random", "adversarial", "bracket"),
                                             eps_max = 0.10, d_D0_d_spread = NULL,
                                             seed = 42L, margin = 0.98) {
  which <- match.arg(which)
  n_cell <- length(D0)
  if (length(a_cell) != n_cell || length(D1) != n_cell)
    stop("hank_partial_id_contamination_f(): 'a_cell', 'D0', 'D1' must have the same length.")

  safe_nonneg_rescale <- function(f) {
    Dmin <- pmin(D0, D1)
    neg_idx <- which(f < 0)
    if (length(neg_idx) == 0L) return(f)
    ratios <- Dmin[neg_idx] / (eps_max * (-f[neg_idx]))
    f * min(1, margin * min(ratios))
  }

  f <- switch(which,
    random = {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
      set.seed(seed)
      raw <- stats::rnorm(n_cell)
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      ord <- order(a_cell)
      sm <- function(x, k = 2) {
        n <- length(x); y <- x
        for (i in seq_len(n)) { lo <- max(1L, i - k); hi <- min(n, i + k); y[i] <- mean(x[lo:hi]) }
        y
      }
      shape <- numeric(n_cell); shape[ord] <- sm(raw[ord])
      shape * D0
    },
    adversarial = {
      if (is.null(d_D0_d_spread))
        stop("hank_partial_id_contamination_f(): which = 'adversarial' requires 'd_D0_d_spread'.")
      wins_cap <- stats::quantile(abs(d_D0_d_spread), 0.99)
      pmax(pmin(d_D0_d_spread, wins_cap), -wins_cap)
    },
    bracket = {
      deciles <- cut(rank(a_cell, ties.method = "first"),
                     breaks = stats::quantile(rank(a_cell, ties.method = "first"), probs = seq(0, 1, 0.1)),
                     include.lowest = TRUE, labels = FALSE)
      out <- numeric(n_cell)
      dec3 <- which(deciles == 3L); dec4 <- which(deciles == 4L)
      mass_moved <- 0.5 * D0[dec3]
      out[dec3] <- -mass_moved
      total_moved <- sum(mass_moved)
      if (sum(D0[dec4]) > 0) out[dec4] <- total_moved * D0[dec4] / sum(D0[dec4])
      out
    })

  f <- f - mean(f)                          # exact mean-zero (sums to 0 over cells)
  l1 <- 0.5 * sum(abs(f))
  if (l1 > 0) f <- f / l1                   # L1/mass normalization (see Details)
  safe_nonneg_rescale(f)
}


#' Partial-identification frontier: LEVEL vs. RESPONSE spread estimation
#' under a fixed additive cross-sectional contamination
#'
#' Bundles the full partial-ID comparison of the het-beta HANK mixture's two
#' spread-identifying channels (see the file header comment):
#' \enumerate{
#'   \item \strong{Efficiency}: \code{I_level}, \code{I_response} -- the
#'     spread's Fisher information under each channel's grid-invariant
#'     functional metric (\code{\link{hank_reweight_level_metric}},
#'     \code{\link{hank_reweight_functional_metric}}), via CENTRAL finite
#'     difference at TWO step sizes (\code{fd_steps}), so step-robustness is
#'     verified rather than assumed.
#'   \item \strong{Robustness}: for a supplied additive contamination
#'     \code{f} (\code{sum(f) == 0}; see
#'     \code{\link{hank_partial_id_contamination_f}} for three canonical
#'     choices) and an \code{eps_grid}, the argmin-over-spread bias each
#'     channel's estimate suffers from observing \code{D0 + eps*f} (LEVEL)
#'     or \code{(D1 + eps*f) - (D0 + eps*f)} (RESPONSE) instead of the clean
#'     \code{D0}/\code{dD}. By construction \code{bias_response(eps) == 0}
#'     for every \code{eps} (the additive \code{f} cancels EXACTLY in the
#'     response's own difference, before the argmin search even runs) --
#'     this is verified, not assumed, via \code{dD_cancel_err}.
#'   \item \strong{Honest frontier}: \code{RMSE_channel(eps, n) =
#'     sqrt(bias_channel(eps)^2 + 1/(I_channel * n))} for
#'     \code{n in n_surveys_grid}, and the crossover \code{eps} at which
#'     \code{RMSE_response(n)} (flat in \code{eps}, since
#'     \code{bias_response == 0}) first drops below \code{RMSE_level(eps,
#'     n)} -- reported both on the literal \code{eps_grid} (\code{NA} if no
#'     crossover occurs there) and via a linear extrapolation of
#'     \code{bias_level(eps)} fit through the \code{eps_grid} (whether or
#'     not that extrapolated crossover falls within a physically plausible
#'     \code{eps <= 1} range is reported explicitly, not asserted).
#' }
#'
#' @param a_grid,Pi,e Household asset grid, income transition, income
#'   levels, shared by every mixture type (see
#'   \code{\link{hank_mixture_ks_steady}}).
#' @param centre,spread TRUTH mixture centre/spread: \code{betas = centre +
#'   spread * beta_shape}. At the default \code{beta_shape} with \code{K = 2}
#'   types this is exactly \code{\link{hank_mixture_joint_logpost}}'s
#'   \code{theta} convention, \code{betas = c(centre - spread, centre +
#'   spread)}. \code{spread} stays the SCALAR estimand at any \code{K}: it
#'   scales the whole heterogeneity pattern, and both channels' Fisher
#'   information is with respect to it.
#' @param omega Length-K TRUTH mixture weights (K >= 2, non-negative, summing
#'   to 1; the mixture GE machinery is K-generic).
#' @param beta_shape Length-K heterogeneity pattern \code{u} in \code{betas =
#'   centre + spread * u}, or \code{NULL} (default) for the equispaced
#'   \code{seq(-1, 1, length.out = K)} (which reduces to the classic
#'   two-type \code{c(-1, 1)} at \code{K = 2}). NOT internally renormalized:
#'   the pattern's scale defines \code{spread}'s units (the default has
#'   half-range 1, so \code{spread} keeps its two-type meaning of "half the
#'   patient-impatient beta gap"). Must contain at least two distinct values
#'   (otherwise \code{spread} moves every type identically with \code{centre}
#'   and is not separately identified).
#' @param eis,alpha,delta,Z GE calibration (see
#'   \code{\link{hank_mixture_ks_steady}}).
#' @param T_h Integer distribution-Jacobian horizon.
#' @param rw_shock Named list, one length-\code{T_h} shock path per
#'   reweighting input (see \code{\link{hank_dist_response_snapshot}}).
#' @param t_star Integer output date in \code{1..T_h} for the RESPONSE
#'   channel's snapshot (see \code{\link{hank_dist_response_snapshot}}).
#' @param rw_scale Scalar multiplier applied to the snapshot response before
#'   it is treated as the RESPONSE channel's \code{dD} (mirrors
#'   \code{\link{hank_mixture_joint_logpost}}'s \code{config$rw_scale}).
#' @param N Reference survey sample size the functional metrics
#'   (\code{\link{hank_reweight_level_metric}},
#'   \code{\link{hank_reweight_functional_metric}}) are built at.
#' @param degree,tail_trim Functional-basis degree / tail trim (see
#'   \code{\link{hank_reweight_functional_basis}}).
#' @param f Length-\code{n_cell} additive contamination vector (must sum to
#'   0; see \code{\link{hank_partial_id_contamination_f}}), or \code{NULL}
#'   (default) to skip the robustness/frontier sections and return only
#'   \code{I_level}/\code{I_response}.
#' @param eps_grid Numeric vector of contamination fractions to evaluate the
#'   bias/frontier at (default \code{c(0, 0.02, 0.05, 0.10)}).
#' @param spread_search Numeric vector of candidate spreads the argmin
#'   search scans over (default a 61-point grid spanning \code{spread +-
#'   0.006}); every point requires one GE solve, precomputed ONCE and reused
#'   across every \code{(f, eps)} bias evaluation (there is only one
#'   \code{f}/vector of \code{eps} per call, so this is a single sweep).
#' @param n_surveys_grid Numeric vector of survey-count multipliers for the
#'   RMSE frontier (default \code{c(1, 100, 1000)}; \code{n} multiplies each
#'   channel's Fisher information linearly, \code{I(n) = n * I(1)}).
#' @param fd_steps Length-2 central-FD step sizes for the Fisher-information
#'   step-robustness check (default \code{c(1e-4, 5e-5)}; the SMALLER step's
#'   information is what is returned as \code{I_level}/\code{I_response}).
#'
#' @return A list with:
#'   \item{I_level, I_response}{Scalars, the Fisher information at the
#'     smaller \code{fd_steps} entry.}
#'   \item{fisher_by_step}{A \code{data.frame} with one row per
#'     \code{fd_steps} entry (\code{h}, \code{I_level}, \code{I_response}),
#'     for inspecting step-robustness directly.}
#'   \item{bias}{\code{NULL} if \code{f} was \code{NULL}, else a
#'     \code{data.frame} with one row per \code{eps_grid} entry:
#'     \code{eps}, \code{bias_level}, \code{bias_response} (identically 0),
#'     \code{dD_cancel_err} (the numeric verification that the additive
#'     \code{f} cancels in the response's difference).}
#'   \item{frontier}{\code{NULL} if \code{f} was \code{NULL}, else a
#'     \code{data.frame} with one row per \code{(eps, n)} pair in
#'     \code{eps_grid x n_surveys_grid}: \code{eps}, \code{n},
#'     \code{RMSE_level}, \code{RMSE_response}.}
#'   \item{crossover}{\code{NULL} if \code{f} was \code{NULL}, else a
#'     \code{data.frame} with one row per \code{n_surveys_grid} entry:
#'     \code{n}, \code{eps_star_grid} (first \code{eps_grid} entry at which
#'     \code{RMSE_level > RMSE_response}, or \code{NA} if none), and
#'     \code{eps_star_extrap} (the linear-extrapolation crossover solving
#'     \code{bias_level(eps)^2 == (1/I_response - 1/I_level)/n} using a
#'     through-origin OLS slope fit to \code{bias_level(eps_grid[-1])}; may
#'     exceed 1 -- i.e. fall outside any physically achievable contamination
#'     -- reported as-is, not clamped).}
#' @seealso \code{\link{hank_reweight_level_metric}},
#'   \code{\link{hank_reweight_functional_metric}},
#'   \code{\link{hank_reweight_functional_loss}},
#'   \code{\link{hank_partial_id_contamination_f}},
#'   \code{\link{hank_mixture_joint_logpost}}
#' @export
hank_partial_id_level_response <- function(a_grid, Pi, e, centre, spread, omega,
                                            beta_shape = NULL,
                                            eis = 1, alpha, delta, Z = 1, T_h,
                                            rw_shock, t_star, rw_scale, N,
                                            degree = 8L, tail_trim = 0,
                                            f = NULL, eps_grid = c(0, 0.02, 0.05, 0.10),
                                            spread_search = seq(spread - 0.006, spread + 0.006, length.out = 61L),
                                            n_surveys_grid = c(1, 100, 1000),
                                            fd_steps = c(1e-4, 5e-5)) {
  K <- length(omega)
  if (K < 2L)
    stop("hank_partial_id_level_response(): 'omega' must have length >= 2 (one weight per mixture type).")
  if (any(!is.finite(omega)) || any(omega < 0))
    stop("hank_partial_id_level_response(): 'omega' must be finite and non-negative.")
  if (abs(sum(omega) - 1) > 1e-8)
    stop("hank_partial_id_level_response(): 'omega' must sum to 1.")
  if (is.null(beta_shape)) beta_shape <- seq(-1, 1, length.out = K)
  if (length(beta_shape) != K)
    stop(sprintf(
      "hank_partial_id_level_response(): length(beta_shape) (%d) must equal length(omega) (%d).",
      length(beta_shape), K))
  if (any(!is.finite(beta_shape)) || diff(range(beta_shape)) <= 0)
    stop("hank_partial_id_level_response(): 'beta_shape' must be finite with at least two distinct values (otherwise 'spread' is not separately identified from 'centre').")
  n_cell <- length(a_grid) * length(e)
  a_cell <- rep(a_grid, times = length(e))

  ge_at <- function(cc, ss) .hank_partial_id_ge(a_grid, Pi, e, cc, ss, omega,
                                                beta_shape, eis,
                                                alpha, delta, Z, T_h, rw_shock, t_star)

  ge0 <- ge_at(centre, spread)
  D0  <- ge0$D0
  D1  <- D0 + rw_scale * ge0$dD_snap
  if (!isTRUE(all.equal(sum(D0), 1, tolerance = 1e-8)))
    stop("hank_partial_id_level_response(): the mixture stationary distribution D0 does not sum to 1.")

  level_metric <- hank_reweight_level_metric(a_cell, D0, N, degree = degree, tail_trim = tail_trim)
  resp_metric  <- hank_reweight_functional_metric(a_cell, D0, D1, N, degree = degree, tail_trim = tail_trim)

  ## ---- 1. Efficiency: Fisher information at each fd_steps entry ----------
  moment_level    <- function(ss) as.numeric(crossprod(level_metric$B, ge_at(centre, ss)$D0))
  moment_response <- function(ss) as.numeric(crossprod(resp_metric$B, rw_scale * ge_at(centre, ss)$dD_snap))

  fisher_by_step <- do.call(rbind, lapply(fd_steps, function(h) {
    fl <- .hank_partial_id_fisher(moment_level, level_metric$Sigma_inv, spread, h)
    fr <- .hank_partial_id_fisher(moment_response, resp_metric$Sigma_inv, spread, h)
    data.frame(h = h, I_level = fl$I, I_response = fr$I)
  }))
  I_level    <- fisher_by_step$I_level[which.min(fd_steps)]
  I_response <- fisher_by_step$I_response[which.min(fd_steps)]

  out <- list(I_level = I_level, I_response = I_response, fisher_by_step = fisher_by_step,
              bias = NULL, frontier = NULL, crossover = NULL)
  if (is.null(f)) return(out)

  if (length(f) != n_cell)
    stop(sprintf("hank_partial_id_level_response(): length(f) (%d) must equal n_cell (%d).",
                 length(f), n_cell))
  if (abs(sum(f)) > 1e-6 * max(1, sum(abs(f))))
    stop("hank_partial_id_level_response(): 'f' must sum to (approximately) 0.")

  ## ---- 2. Bias: precompute the argmin search grid's GE solves ONCE -------
  grid_ge     <- lapply(spread_search, function(ss) ge_at(centre, ss))
  grid_dD_mod <- lapply(grid_ge, function(g) rw_scale * g$dD_snap)

  parabolic_refine <- function(sg, losses, idx) {
    if (idx <= 1L || idx >= length(sg)) return(sg[idx])
    s1 <- sg[idx - 1L]; s2 <- sg[idx]; s3 <- sg[idx + 1L]
    y1 <- losses[idx - 1L]; y2 <- losses[idx]; y3 <- losses[idx + 1L]
    denom <- y1 - 2 * y2 + y3
    if (abs(denom) < 1e-300) return(s2)
    s2 + 0.5 * (s2 - s1) * (y1 - y3) / denom
  }
  argmin_over_grid <- function(loss_of_idx, refine = TRUE) {
    losses <- vapply(seq_along(spread_search), loss_of_idx, numeric(1))
    idx <- which.min(losses)
    if (refine) parabolic_refine(spread_search, losses, idx) else spread_search[idx]
  }

  ## RESPONSE bias is IDENTICAL for every eps (the additive f cancels
  ## EXACTLY in dD_obs = D1_obs - D0_obs -- computed once via the RAW grid
  ## argmin, since the noiseless response loss surface is a near-perfect
  ## quadratic bowl whose true minimum sits so sharply at the truth that a
  ## 3-point parabolic refit amplifies floating-point asymmetry noise rather
  ## than recovering genuine signal; see the analogous note in the scratch
  ## study this function promotes).
  dD_clean <- D1 - D0
  spread_hat_resp <- argmin_over_grid(function(i)
    hank_reweight_functional_loss(dD_clean, grid_dD_mod[[i]], resp_metric)$rho, refine = FALSE)
  bias_response_val <- spread_hat_resp - spread

  bias_rows <- lapply(eps_grid, function(eps) {
    D0_obs <- D0 + eps * f
    D1_obs <- D1 + eps * f
    dD_obs <- D1_obs - D0_obs
    dD_cancel_err <- max(abs(dD_obs - dD_clean))

    spread_hat_lev <- argmin_over_grid(function(i)
      hank_reweight_functional_loss(D0_obs, grid_ge[[i]]$D0, level_metric)$rho, refine = TRUE)
    data.frame(eps = eps, bias_level = spread_hat_lev - spread,
               bias_response = bias_response_val, dD_cancel_err = dD_cancel_err)
  })
  bias_df <- do.call(rbind, bias_rows)

  ## ---- 3. Honest frontier + crossover --------------------------------
  rmse_level_fn    <- function(bias, n) sqrt(bias^2 + 1 / (I_level * n))
  rmse_response_fn <- function(n) sqrt(1 / (I_response * n))

  frontier_rows <- lapply(seq_len(nrow(bias_df)), function(i) {
    do.call(rbind, lapply(n_surveys_grid, function(n) {
      data.frame(eps = bias_df$eps[i], n = n,
                 RMSE_level = rmse_level_fn(bias_df$bias_level[i], n),
                 RMSE_response = rmse_response_fn(n))
    }))
  })
  frontier_df <- do.call(rbind, frontier_rows)

  ## extrapolated bias-vs-eps slope (through-origin OLS on eps_grid > 0
  ## points; bias_level(0) is a pure argmin-search noise floor, not a
  ## contamination effect, so it is excluded from the slope fit)
  pos_idx <- which(bias_df$eps > 0)
  slope <- if (length(pos_idx) >= 1L)
    sum(bias_df$eps[pos_idx] * bias_df$bias_level[pos_idx]) / sum(bias_df$eps[pos_idx]^2)
  else NA_real_

  crossover_rows <- lapply(n_surveys_grid, function(n) {
    sub <- frontier_df[frontier_df$n == n, ]
    exceeds <- sub$RMSE_level > sub$RMSE_response
    eps_star_grid <- if (any(exceeds)) sub$eps[which(exceeds)[1L]] else NA_real_
    target <- sqrt(max(0, 1 / I_response - 1 / I_level) / n)
    eps_star_extrap <- if (!is.na(slope) && abs(slope) > 0) target / abs(slope) else NA_real_
    data.frame(n = n, eps_star_grid = eps_star_grid, eps_star_extrap = eps_star_extrap)
  })
  crossover_df <- do.call(rbind, crossover_rows)

  out$bias <- bias_df
  out$frontier <- frontier_df
  out$crossover <- crossover_df
  out
}
