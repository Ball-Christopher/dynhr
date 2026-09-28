## R/diag-deep-d36-calibration.R
## --------------------------------------------------------------------------
## D36. Calibration deepness / validation -- "is this *fixed* parameter
## actually deep, or is the calibration hiding misspecification?"
##
## Most DSGE deep parameters are calibrated, not estimated, so the rest of the
## deep-parameter suite (D1/D6/D33/D34/D35 -- which work off estimated draws or
## the estimated mode) reports them as "?" (unassessed). But a calibrated
## parameter still leaves a footprint in the likelihood: even though we fix it,
## we can ask the data what value IT would choose. D36 sweeps the likelihood
## along each calibrated deep parameter -- either a one-dimensional SLICE
## (the rest of the model held fixed; a conditional likelihood) or the
## CONCENTRATED/PROFILE likelihood (the `estimated` block re-optimised at every
## grid point) -- and reports three things:
##
##   1. IDENTIFICATION -- is the curve curved at all? A flat curve means the
##      data cannot speak to the parameter; the calibration is then a pure
##      assumption ("untestable") and its deepness cannot be verified or refuted.
##   2. TENSION -- where does the likelihood peak relative to the calibrated
##      value? A likelihood-ratio test of H0: c = c_cal,
##      LR = 2 (l_max - l(c_cal)) vs chi-square(1) (3.84 at 95%, i.e. a drop of
##      1.92 nats), plus the matching LR interval {c : l(c) >= l_max - 1.92}.
##   3. CONSTANCY (optional) -- profiling on different sub-samples, does the
##      likelihood-implied value drift? A "deep" constant should not. This is the
##      Lucas-critique constancy test extended to *calibrated* parameters.
##
## Verdict per calibrated deep parameter: untestable / consistent / tension /
## failed (the likelihood could not be evaluated on enough of the grid)
## (+ a non-constant flag). Exposes `$result$passport_axis` for the Passport's
## "calibrated" column.
##
## Contract: `loglik_fn` must RE-SOLVE the model at the parameter vector it is
## given (a closure over a fixed decision rule makes every structural profile
## flat) and should signal an out-of-domain vector by returning a non-finite
## value (as make_loglik_contrib() does), not by erroring.
##
## References:
##   Canova, F. (1994). Statistical inference in calibrated models. JAE, 9, S123.
##   Jorgensen, T. H. (2023). Sensitivity to calibrated parameters. REStat,
##     105(2), 474-481.
##   Hansen, L. P., & Heckman, J. J. (1996). The empirical foundations of
##     calibration. Journal of Economic Perspectives, 10(1), 87-104.
##   Mueller, U. K. (2012). Measuring prior sensitivity ... JME, 59(6), 581-597.
## --------------------------------------------------------------------------


# Default profiling range, always CENTRED on the calibrated value v (so with an
# odd n_grid v is itself a grid point and LR >= 0 by construction):
#  * v in (0,1): additive half-width min(0.3, frac*v, frac*(1-v)) -- stays
#    strictly inside (0,1) and scales down for small values (a 0.007 shock
#    std) and for values near 1 (beta = 0.99). The old clamp to [0.02, 0.98]
#    produced windows that EXCLUDED v in both cases.
#  * v == 0: +/- frac (the sign of an admissible value is unknown; points the
#    loglik rejects come back non-finite and are dropped).
#  * otherwise: multiplicative v * (1 -/+ frac).
.d36_default_range <- function(v, frac = 0.5) {
  if (length(v) != 1L || !is.finite(v)) return(NULL)
  if (v > 0 && v < 1) {
    half <- min(0.3, frac * v, frac * (1 - v))
    c(v - half, v + half)
  } else if (v == 0) {
    c(-1, 1) * frac
  } else {
    sort(c(v * (1 - frac), v * (1 + frac)))
  }
}

# Which grid edge holds the (finite) maximum: "lo", "hi" or NA (interior /
# too few finite points).
.d36_edge_max <- function(ll) {
  ok <- is.finite(ll)
  if (sum(ok) < 3) return(NA_character_)
  i <- which.max(replace(ll, !ok, -Inf))
  if (i == 1L) "lo" else if (i == length(ll)) "hi" else NA_character_
}

# Widen one side of a range by its width, never crossing the sign of v and,
# for v in (0,1), staying inside (0,1) (at most halving the distance to the
# boundary).
.d36_expand_range <- function(rng, side, v) {
  lo <- rng[1]; hi <- rng[2]; w <- hi - lo
  if (side == "lo") {
    lo <- if (v > 0) max(lo - w, lo / 2) else lo - w
  } else {
    hi <- if (v > 0 && v < 1) min(hi + w, (hi + 1) / 2)
          else if (v < 0) min(hi + w, hi / 2)
          else hi + w
  }
  c(lo, hi)
}

# Evaluate the user loglik; anything that is not a finite scalar -> NA.
.d36_ll_value <- function(loglik_fn, p) {
  v <- loglik_fn(p)
  if (is.numeric(v) && length(v) == 1L && is.finite(v)) as.numeric(v) else NA_real_
}

# Slice value at target = x (all other params held at `params`).
.d36_slice_fn <- function(loglik_fn, params, target) {
  function(x) {
    p <- params; p[[target]] <- x
    .d36_ll_value(loglik_fn, p)
  }
}

# Profile one parameter (slice): evaluate the loglik over a grid with `target`
# swept across `range`, holding all other params at `params`.
.d36_profile_one <- function(loglik_fn, params, target, range, n_grid) {
  grid <- seq(range[1], range[2], length.out = n_grid)
  f <- .d36_slice_fn(loglik_fn, params, target)
  list(grid = grid, ll = vapply(grid, f, numeric(1)))
}

# Refinement around the grid argmax -> (chat, d2 = curvature, llmax).
# Grid argmax, then (interior only) the exact 3-point parabola; when the
# one-dimensional curve `f` is supplied, a Brent search (stats::optimize) on
# the bracketing interval [x_{i-1}, x_{i+1}] and a central second difference
# at the refined maximiser (step h/4) replace the parabola, so the implied
# value and l_max are exact to `tol` rather than to the grid spacing.
.d36_refine <- function(grid, ll, f = NULL) {
  ok <- is.finite(ll)
  if (sum(ok) < 3) return(list(chat = NA_real_, d2 = NA_real_, llmax = NA_real_))
  i <- which.max(replace(ll, !ok, -Inf))
  n <- length(grid)
  if (i == 1L || i == n || !all(is.finite(ll[c(i - 1L, i + 1L)])))
    return(list(chat = grid[i], d2 = NA_real_, llmax = ll[i]))
  x0 <- grid[i - 1L]; x1 <- grid[i]; x2 <- grid[i + 1L]
  y0 <- ll[i - 1L];   y1 <- ll[i];   y2 <- ll[i + 1L]
  h  <- x1 - x0
  d2 <- (y0 - 2 * y1 + y2) / h^2
  a1 <- (y2 - y0) / (2 * h)
  chat  <- if (d2 < 0) x1 - a1 / d2 else x1
  chat  <- min(max(chat, x0), x2)
  llmax <- if (d2 < 0) y1 - a1^2 / (2 * d2) else y1
  if (is.function(f)) {
    g <- function(x) { v <- f(x); if (is.finite(v)) v else -.Machine$double.xmax }
    opt <- stats::optimize(g, c(x0, x2), maximum = TRUE,
                           tol = (x2 - x0) * 1e-6)
    if (is.finite(opt$objective) && opt$objective > -.Machine$double.xmax &&
        opt$objective >= y1) {
      chat <- opt$maximum; llmax <- opt$objective
    } else {
      chat <- x1; llmax <- y1
    }
    hh <- h / 4
    fm <- f(chat - hh); fp <- f(chat + hh)
    if (is.finite(fm) && is.finite(fp)) d2 <- (fm - 2 * llmax + fp) / hh^2
  }
  list(chat = chat, d2 = d2, llmax = max(llmax, y1))
}

# LR interval {x : l(x) >= llmax - cut/2}: bracket each crossing between the
# grid point just outside and the nearest point inside (a grid point or chat
# itself), then -- when the curve `f` is supplied -- bisect it (`n_bisect`
# steps; a non-finite midpoint counts as outside); otherwise interpolate
# linearly. NA on a side where the curve never falls below the cut-off inside
# the range (interval open beyond the range there).
.d36_lr_interval <- function(grid, ll, llmax, chat, cut, f = NULL,
                             n_bisect = 30L) {
  out <- c(lo = NA_real_, hi = NA_real_)
  if (!is.finite(llmax) || !is.finite(chat)) return(out)
  lev <- llmax - cut / 2
  z <- ll - lev                          # >= 0 inside the interval
  ok <- is.finite(z)
  crossing <- function(x_out, z_out, x_in, z_in) {
    if (!is.function(f))
      return(x_out + (x_in - x_out) * (0 - z_out) / (z_in - z_out))
    for (it in seq_len(n_bisect)) {
      xm <- (x_out + x_in) / 2
      zm <- f(xm) - lev
      if (is.finite(zm) && zm >= 0) x_in <- xm else x_out <- xm
    }
    (x_out + x_in) / 2
  }
  left <- which(ok & grid <= chat)
  below_l <- left[z[left] < 0]
  if (length(below_l)) {
    j <- max(below_l)
    k <- left[left > j]
    xk <- if (length(k)) grid[min(k)] else chat
    zk <- if (length(k)) z[min(k)] else cut / 2
    out["lo"] <- crossing(grid[j], z[j], xk, zk)
  }
  right <- which(ok & grid >= chat)
  below_r <- right[z[right] < 0]
  if (length(below_r)) {
    j <- min(below_r)
    k <- right[right < j]
    xk <- if (length(k)) grid[max(k)] else chat
    zk <- if (length(k)) z[max(k)] else cut / 2
    out["hi"] <- crossing(grid[j], z[j], xk, zk)
  }
  out
}


# Concentrated log-likelihood: fix `target` at `fixed_val`, re-optimise the
# estimated/nuisance block `est`. Returns the max loglik and the argmax (for
# warm-starting the next grid point).
.d36_concentrate <- function(loglik_fn, params, target, fixed_val, est,
                             start, lower, upper, maxit = 200) {
  p0 <- params; p0[[target]] <- fixed_val
  if (length(est) == 0) {
    return(list(value = .d36_ll_value(loglik_fn, p0), par = numeric(0)))
  }
  obj <- function(th) {
    p <- p0
    for (i in seq_along(est)) p[[est[i]]] <- th[i]
    v <- .d36_ll_value(loglik_fn, p)
    if (is.finite(v)) -v else 1e10
  }
  use_box <- !is.null(lower) && !is.null(upper) &&
             all(is.finite(lower)) && all(is.finite(upper))
  opt <- if (use_box) {
    stats::optim(start, obj, method = "L-BFGS-B", lower = lower,
                 upper = upper, control = list(maxit = maxit))
  } else if (length(est) == 1L) {
    # optim(Nelder-Mead) warns in 1-D; use it anyway (no bounds known) but
    # keep the warning local to this call.
    suppressWarnings(stats::optim(start, obj, method = "Nelder-Mead",
                                  control = list(maxit = maxit)))
  } else {
    stats::optim(start, obj, method = "Nelder-Mead",
                 control = list(maxit = maxit))
  }
  if (!is.finite(opt$value) || opt$value >= 1e10)
    list(value = NA_real_, par = start)
  else list(value = -opt$value, par = opt$par)
}

# Full profile: sweep `target` across `range`, re-optimising `est` at each
# point (warm-started). Also returns the concentrated value at `base_val` and
# the one-dimensional concentrated curve (for refinement).
.d36_profile_full <- function(loglik_fn, params, target, range, n_grid,
                              est, lower, upper, base_val) {
  grid <- seq(range[1], range[2], length.out = n_grid)
  ll   <- rep(NA_real_, n_grid)
  pars <- vector("list", n_grid)
  warm <- unlist(params[est])
  for (i in seq_along(grid)) {
    r <- .d36_concentrate(loglik_fn, params, target, grid[i], est, warm, lower, upper)
    ll[i] <- r$value
    pars[[i]] <- r$par
    if (length(r$par) && all(is.finite(r$par))) warm <- r$par
  }
  base <- .d36_concentrate(loglik_fn, params, target, base_val,
                           est, unlist(params[est]), lower, upper)$value
  # concentrated curve, warm-started from the nearest grid point's argmax
  f <- function(x) {
    j <- which.min(abs(grid - x))
    st <- pars[[j]]
    if (!length(st) || !all(is.finite(st))) st <- unlist(params[est])
    .d36_concentrate(loglik_fn, params, target, x, est, st, lower, upper)$value
  }
  list(grid = grid, ll = ll, base = base, f = f)
}


# ---------------------------------------------------------------------------
#' D36. Calibration deepness / validation
#'
#' Sweeps each calibrated deep parameter along a grid and asks whether the
#' likelihood's own maximiser is consistent with the calibrated value.
#'
#' \strong{\code{profile = "slice"} is the default, and a slice is NOT a
#' profile.} A slice holds every other parameter fixed at its current value
#' while the target moves: a CONDITIONAL likelihood. The asymptotic
#' \eqn{\chi^2_1} calibration of the reported LR belongs to the PROFILE
#' likelihood, which optimises the nuisance block out at each grid point
#' (Murphy & van der Vaart 2000). A slice's curvature and tension therefore do
#' NOT account for cross-parameter compensation and will overstate tension
#' whenever the target trades off with an estimated nuisance parameter; read
#' them as conditional-on-the-rest-of-the-calibration, and rerun with
#' \code{profile = "full"} (which needs \code{estimated}) before drawing a
#' conclusion about a parameter that matters.
#'
#' The slice default is not an ad hoc simplification: it is exactly the
#' convention of \strong{Dynare's \code{mode_check}}, which evaluates the
#' posterior kernel in a neighbourhood of the mode one parameter at a time
#' while holding all the others at their point estimates. No DSGE package in
#' wide use computes re-optimised profile likelihoods for every calibrated
#' parameter by default, because the cost is one re-optimisation per grid
#' point per parameter; published profile-likelihood exercises are one-off
#' robustness checks on one or two key parameters. For the same reason D36 is
#' not auto-built by \code{run_model_diagnostics()}: it must be called
#' directly with a \code{loglik_fn} (or reached through
#' \code{run_all_diagnostics(loglik_full_fn = )}).
#'
#' @param model     Parsed model (for the @dynhr:deep taxonomy). Optional if
#'   \code{deep_spec} is supplied.
#' @param deep_spec Optional \code{\link{build_deep_spec}}.
#' @param params    Named numeric full parameter vector (the calibration).
#' @param loglik_fn Function: a *full* named parameter vector -> scalar
#'   log-likelihood. It must RE-SOLVE the model at the vector it receives (a
#'   closure over a fixed decision rule makes every profile flat) and should
#'   return a non-finite value (not error) for an inadmissible vector.
#' @param targets   Character vector of parameters to profile. Default: the
#'   calibrated deep primitives (deep, role = primitive, not in
#'   \code{estimated}).
#' @param estimated Character vector of estimated parameter names (excluded from
#'   the default \code{targets}).
#' @param ranges    Optional named list of \code{c(lo, hi)} profiling ranges.
#'   Default: a window centred on the calibrated value (see
#'   \code{.d36_default_range}).
#' @param n_grid    Grid points per profile (default 21; odd keeps the
#'   calibrated value on the grid for the default range). The grid maximum is
#'   refined by a Brent search on its bracketing interval.
#' @param frac      Default relative half-width for ranges (default 0.5).
#' @param lr_threshold Likelihood-ratio cutoff for "tension", on the
#'   \eqn{2\Delta\ell} scale (default 3.84 = chi-square(1) 95\%, i.e. a
#'   log-likelihood drop of 1.92). Also defines the reported LR interval and
#'   the cut-off line in the plot.
#' @param min_curv_nats Minimum profile height (nats) over the range to call a
#'   parameter identifiable rather than "untestable" (default 1).
#' @param max_expand Maximum number of times a DEFAULT range is widened on
#'   the side where the maximum sits on the grid edge (default 3; 0 = never).
#' @param profile \code{"slice"} (default) holds the rest of the model fixed
#'   while sweeping the target -- a cheap \emph{conditional} likelihood, the
#'   Dynare \code{mode_check} convention, NOT a profile likelihood (see
#'   Details).
#'   \code{"full"} re-optimises the \code{estimated} nuisance block at every
#'   grid point (the concentrated / profile likelihood), which separates
#'   genuine calibration tension from a target that merely trades off with a
#'   nuisance parameter; it requires \code{estimated} and is much more
#'   expensive.
#' @param estimated_bounds Optional named list of \code{c(lo, hi)} box
#'   constraints for the \code{estimated} block (enables L-BFGS-B in the full
#'   profile; otherwise Nelder-Mead).
#' @param subsample_loglik_fns Optional named list of per-regime full-vector
#'   log-likelihood functions for the constancy check (slice profiles).
#' @param meta      Optional \code{\link{diag_meta}} provenance descriptor.
#' @references
#'   Murphy, S. A. & van der Vaart, A. W. (2000). On profile likelihood.
#'   \emph{Journal of the American Statistical Association}, 95(450), 449-465.
#' @return A \code{dynhr_diagnostic}. \code{$result$table} has one row per
#'   target: \code{calibrated}, \code{implied} (maximiser), \code{profile_se}
#'   (\eqn{1/\sqrt{-\ell''}}), \code{height_nats}, \code{LR}, \code{p_value},
#'   \code{lr_lo}/\code{lr_hi} (LR interval; \code{NA} = open beyond the
#'   range), \code{at_bound}, \code{verdict} (untestable / consistent /
#'   tension / failed), \code{non_constant}. \code{$result$passport_axis} is a
#'   named logical (TRUE = consistent/deep, FALSE = tension) over the
#'   testable calibrated parameters for the Passport's "calibrated" column.
#' @noRd
d36_calibration_deepness <- function(model = NULL,
                                     deep_spec = NULL,
                                     params = NULL,
                                     loglik_fn = NULL,
                                     targets = NULL,
                                     estimated = NULL,
                                     ranges = NULL,
                                     n_grid = 21,
                                     frac = 0.5,
                                     lr_threshold = 3.84,
                                     min_curv_nats = 1,
                                     max_expand = 3L,
                                     profile = c("slice", "full"),
                                     estimated_bounds = NULL,
                                     subsample_loglik_fns = NULL,
                                     meta = NULL) {
  profile <- match.arg(profile)
  if (is.null(loglik_fn))
    return(.make_result(pass = NA,
      summary = "D36 calibration deepness: no loglik_fn supplied."))
  if (!is.function(loglik_fn))
    .dynhr_abort("D36: `loglik_fn` must be a function.")
  if (!(is.numeric(lr_threshold) && length(lr_threshold) == 1L &&
        is.finite(lr_threshold) && lr_threshold > 0))
    .dynhr_abort("D36: `lr_threshold` must be a positive number.")
  if (!(is.numeric(n_grid) && length(n_grid) == 1L && n_grid >= 5))
    .dynhr_abort("D36: `n_grid` must be >= 5.")
  if (is.null(params)) params <- model$param_values
  if (is.null(params) || length(params) == 0)
    return(.make_result(pass = NA,
      summary = "D36 calibration deepness: no parameter values available."))
  params <- as.list(params)
  if (is.null(deep_spec)) deep_spec <- build_deep_spec(model = model)

  if (is.null(targets)) {
    cal_deep <- deep_spec$param[deep_spec$is_deep &
                                deep_spec$role == "primitive"]
    targets <- setdiff(cal_deep, estimated)
  }
  targets <- intersect(targets, names(params))
  if (length(targets) == 0)
    return(.make_result(pass = NA,
      summary = "D36 calibration deepness: no calibrated deep primitives to test."))

  base_ll <- .d36_ll_value(loglik_fn, params)
  if (!is.finite(base_ll))
    return(.make_result(pass = NA,
      summary = "D36 calibration deepness: baseline log-likelihood is not finite."))
  if (profile == "full" && (is.null(estimated) || length(estimated) == 0))
    return(.make_result(pass = NA,
      summary = "D36 (profile='full'): supply `estimated` (the nuisance block to re-optimise at each grid point)."))
  if (profile == "full" && !all(estimated %in% names(params)))
    .dynhr_abort("D36: `estimated` not in `params`: ",
                 paste(setdiff(estimated, names(params)), collapse = ", "), ".")

  .bounds_for <- function(est_tg) {
    if (is.null(estimated_bounds)) return(list(lo = NULL, hi = NULL))
    lo <- vapply(est_tg, function(e) (estimated_bounds[[e]] %||% c(NA, NA))[1], numeric(1))
    hi <- vapply(est_tg, function(e) (estimated_bounds[[e]] %||% c(NA, NA))[2], numeric(1))
    list(lo = lo, hi = hi)
  }
  .range_for <- function(tg, v) {
    if (!is.null(ranges[[tg]])) as.numeric(ranges[[tg]]) else .d36_default_range(v, frac)
  }

  rows <- list(); profiles <- list()
  for (tg in targets) {
    v   <- suppressWarnings(as.numeric(params[[tg]]))
    rng <- if (length(v) == 1L && is.finite(v)) .range_for(tg, v) else NULL
    if (is.null(rng) || length(rng) != 2L || !all(is.finite(rng)) || rng[1] >= rng[2]) {
      rows[[tg]] <- data.frame(
        param = tg, calibrated = if (length(v) == 1L) v else NA_real_,
        implied = NA_real_, profile_se = NA_real_, height_nats = NA_real_,
        LR = NA_real_, p_value = NA_real_, lr_lo = NA_real_, lr_hi = NA_real_,
        at_bound = NA, verdict = "failed", stringsAsFactors = FALSE)
      next
    }
    sweep <- function(rng) {
      if (profile == "full") {
        est_tg <- setdiff(estimated, tg)
        bd <- .bounds_for(est_tg)
        pf <- .d36_profile_full(loglik_fn, params, tg, rng, n_grid, est_tg,
                                bd$lo, bd$hi, base_val = v)
        list(grid = pf$grid, ll = pf$ll, f = pf$f, base = pf$base)
      } else {
        pr <- .d36_profile_one(loglik_fn, params, tg, rng, n_grid)
        list(grid = pr$grid, ll = pr$ll,
             f = .d36_slice_fn(loglik_fn, params, tg), base = base_ll)
      }
    }
    pr <- sweep(rng)
    # Default ranges only: while the maximum sits on a finite grid EDGE, widen
    # that side and re-sweep (a max on the edge understates LR -> a missed
    # tension); then, for a curved profile, widen a side on which the LR
    # interval is still open. User-supplied ranges are respected as given.
    if (is.null(ranges[[tg]])) {
      for (k in seq_len(max_expand)) {
        side <- .d36_edge_max(pr$ll)
        if (is.na(side) && sum(is.finite(pr$ll)) >= 3 &&
            diff(range(pr$ll, na.rm = TRUE)) >= min_curv_nats) {
          ci0 <- .d36_lr_interval(pr$grid, pr$ll, max(pr$ll, na.rm = TRUE),
                                  pr$grid[which.max(pr$ll)], lr_threshold)
          # only a side whose grid edge is evaluable can be usefully widened
          n_g <- length(pr$ll)
          side <- if (is.na(ci0[["lo"]]) && is.finite(pr$ll[1])) "lo"
                  else if (is.na(ci0[["hi"]]) && is.finite(pr$ll[n_g])) "hi"
                  else NA_character_
        }
        if (is.na(side)) break
        rng <- .d36_expand_range(rng, side, v)
        pr  <- sweep(rng)
      }
    }
    f1 <- pr$f
    base_for_lr <- pr$base
    # Default ranges only: a profile so sharp that fewer than 5 grid points
    # fall inside the LR interval is re-swept on a window around it (the
    # numbers are refined by Brent/bisection anyway; this keeps the curve,
    # the height and the plot resolved). The window always contains v.
    if (is.null(ranges[[tg]]) && sum(is.finite(pr$ll)) >= 3) {
      rf0 <- .d36_refine(pr$grid, pr$ll, f1)
      ci0 <- .d36_lr_interval(pr$grid, pr$ll, rf0$llmax, rf0$chat,
                              lr_threshold, f1)
      inside <- sum(is.finite(pr$ll) & pr$ll >= rf0$llmax - lr_threshold / 2)
      if (all(is.finite(ci0)) && inside < 5 && ci0[["hi"]] > ci0[["lo"]]) {
        w_ci <- ci0[["hi"]] - ci0[["lo"]]
        rng  <- range(c(ci0[["lo"]] - w_ci, ci0[["hi"]] + w_ci, v))
        pr   <- sweep(rng)
        f1   <- pr$f
        base_for_lr <- pr$base
      }
    }
    n_ok <- sum(is.finite(pr$ll))
    rf  <- .d36_refine(pr$grid, pr$ll, f1)
    # the calibrated point lies on the same curve: l_max can never be below it
    if (is.finite(rf$llmax) && is.finite(base_for_lr) && base_for_lr > rf$llmax) {
      rf$llmax <- base_for_lr; rf$chat <- v
    }
    failed <- n_ok < 3 || !is.finite(rf$llmax) || !is.finite(base_for_lr)
    height <- if (n_ok >= 2) diff(range(pr$ll, na.rm = TRUE)) else NA_real_
    se  <- if (is.finite(rf$d2) && rf$d2 < 0) sqrt(-1 / rf$d2) else Inf
    lr  <- if (failed) NA_real_ else max(2 * (rf$llmax - base_for_lr), 0)
    pval <- if (is.finite(lr)) stats::pchisq(lr, df = 1, lower.tail = FALSE) else NA_real_
    ci  <- .d36_lr_interval(pr$grid, pr$ll, rf$llmax, rf$chat, lr_threshold, f1)
    # "at bound" = the maximiser sits on the edge of the EVALUABLE support
    fin_x <- pr$grid[is.finite(pr$ll)]
    tol_b <- 1e-9 * max(1, abs(rng))
    at_bound <- !failed && is.finite(rf$chat) &&
      (rf$chat <= min(fin_x) + tol_b || rf$chat >= max(fin_x) - tol_b)

    identifiable <- is.finite(height) && height >= min_curv_nats
    verdict <- if (failed) "failed"
               else if (lr >= lr_threshold) "tension"
               else if (!identifiable) "untestable"
               else "consistent"

    profiles[[tg]] <- data.frame(param = tg, x = pr$grid, ll = pr$ll,
                                 stringsAsFactors = FALSE)
    rows[[tg]] <- data.frame(
      param = tg, calibrated = v, implied = rf$chat,
      profile_se = se, height_nats = height, LR = lr, p_value = pval,
      lr_lo = unname(ci["lo"]), lr_hi = unname(ci["hi"]),
      at_bound = at_bound, verdict = verdict,
      llmax = rf$llmax, base_ll = base_for_lr, stringsAsFactors = FALSE)
  }
  rows <- lapply(rows, function(r) {
    if (is.null(r$llmax)) { r$llmax <- NA_real_; r$base_ll <- NA_real_ }
    r
  })
  tab <- do.call(rbind, rows)
  rownames(tab) <- NULL

  # optional cross-regime constancy of the implied value
  constancy <- NULL
  if (!is.null(subsample_loglik_fns) && length(subsample_loglik_fns) >= 2) {
    ok_rows <- tab$param[tab$verdict != "failed"]
    cc <- lapply(names(subsample_loglik_fns), function(rn) {
      f <- subsample_loglik_fns[[rn]]
      vapply(ok_rows, function(tg) {
        v <- as.numeric(params[[tg]])
        rng <- .range_for(tg, v)
        pr <- .d36_profile_one(f, params, tg, rng, n_grid)
        .d36_refine(pr$grid, pr$ll, .d36_slice_fn(f, params, tg))$chat
      }, numeric(1))
    })
    cmat <- matrix(unlist(cc), nrow = length(ok_rows),
                   dimnames = list(ok_rows, names(subsample_loglik_fns)))
    spread <- apply(cmat, 1, function(z)
      if (sum(is.finite(z)) >= 2) diff(range(z, na.rm = TRUE)) else NA_real_)
    se_ok <- tab$profile_se[match(ok_rows, tab$param)]
    non_constant <- ok_rows[is.finite(spread) & is.finite(se_ok) &
                            spread > 2 * se_ok]
    constancy <- list(implied_by_regime = cmat, spread = spread,
                      non_constant = non_constant)
    tab$non_constant <- tab$param %in% non_constant
  } else {
    tab$non_constant <- FALSE
  }

  tension     <- tab$param[tab$verdict == "tension"]
  untestable  <- tab$param[tab$verdict == "untestable"]
  consistent  <- tab$param[tab$verdict == "consistent"]
  failed_p    <- tab$param[tab$verdict == "failed"]
  non_const   <- tab$param[tab$non_constant]

  # passport axis: consistent (+constant) -> TRUE, tension/non-constant ->
  # FALSE, untestable/failed -> omitted (stays "?")
  testable  <- tab$verdict %in% c("consistent", "tension")
  ax_params <- tab$param[testable]
  ax_vals   <- (tab$verdict[testable] == "consistent") & !tab$non_constant[testable]
  passport_axis <- stats::setNames(ax_vals, ax_params)

  # Every evaluable profile EXACTLY flat = loglik_fn ignores the targets
  # (typically a closure over a fixed decision rule). Not evidence of
  # anything -> INFO, not PASS.
  ok_h <- tab$height_nats[tab$verdict != "failed"]
  insensitive <- length(ok_h) > 0 && all(is.finite(ok_h)) &&
    all(ok_h <= 1e-10 * max(1, abs(base_ll)))
  pass <- if (length(tension) || length(non_const)) FALSE
          else if (length(failed_p) == nrow(tab) || insensitive) NA
          else TRUE

  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE) && length(profiles) > 0)
    plots$profiles <- .plot_d36_profiles(do.call(rbind, profiles), tab, meta,
                                         profile = profile,
                                         lr_threshold = lr_threshold)

  summary_txt <- sprintf(
    "D36 Calibration deepness (%s): %d calibrated deep param(s) -- %d consistent, %d tension, %d untestable%s%s.",
    if (profile == "full") "concentrated profile" else "slice",
    nrow(tab), length(consistent), length(tension), length(untestable),
    if (length(failed_p)) sprintf(", %d failed", length(failed_p)) else "",
    if (length(non_const)) sprintf(", %d non-constant", length(non_const)) else "")
  if (insensitive)
    summary_txt <- paste(summary_txt,
      "loglik_fn is insensitive to every target (does it re-solve the model?).")

  badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
  fmt_ci <- function(lo, hi) sprintf("[%s, %s]",
    if (is.finite(lo)) sprintf("%.4g", lo) else "<range",
    if (is.finite(hi)) sprintf("%.4g", hi) else ">range")
  worst <- tab[tab$verdict != "failed", , drop = FALSE]
  worst <- utils::head(worst[order(-worst$LR), , drop = FALSE], 4)
  llm <- paste(c(
    sprintf("D36 | Calibration Deepness | %s", badge),
    sprintf("  profile=%s lr_threshold=%.3g (drop %.3g nats)", profile,
            lr_threshold, lr_threshold / 2),
    sprintf("  calibrated_deep=%d consistent=%d tension=%d untestable=%d failed=%d non_constant=%d",
            nrow(tab), length(consistent), length(tension), length(untestable),
            length(failed_p), length(non_const)),
    if (nrow(worst))
      sprintf("  implied_vs_calibrated: %s",
              paste(sprintf("%s cal=%.3g data=%.3g LR=%.2f CI=%s", worst$param,
                            worst$calibrated, worst$implied, worst$LR,
                            mapply(fmt_ci, worst$lr_lo, worst$lr_hi)),
                    collapse = "; ")),
    if (length(tension))
      sprintf("  tension(data disputes the fixed value): %s",
              paste(tension, collapse = ", ")),
    if (length(non_const))
      sprintf("  non_constant(implied value drifts across regimes): %s",
              paste(non_const, collapse = ", ")),
    if (length(untestable))
      sprintf("  untestable(flat profile, data uninformative): %s",
              paste(utils::head(untestable, 8), collapse = ", ")),
    if (length(failed_p))
      sprintf("  failed(loglik non-finite on too much of the grid): %s",
              paste(utils::head(failed_p, 8), collapse = ", ")),
    sprintf("  action: %s",
            if (insensitive)
              "loglik_fn does not respond to ANY target (every profile exactly flat) -- it probably reuses a fixed decision rule; pass a closure that re-solves the model."
            else if (is.na(pass))
              "No profile could be evaluated -- check that loglik_fn re-solves the model and returns finite values near the calibration."
            else if (pass)
              "Calibrated deep parameters are either data-consistent or innocuous (flat profile); no calibration tension."
            else sprintf("%s: the data prefers a different value than the calibration -- re-calibrate, estimate it, or treat the gap as a misspecification signal.",
                         paste(utils::head(c(tension, non_const), 3), collapse = ", ")))
  ), collapse = "\n")

  .make_result(
    result = list(table = tab, profiles = profiles, base_loglik = base_ll,
                  constancy = constancy, tension = tension, profile = profile,
                  lr_threshold = lr_threshold,
                  untestable = untestable, failed = failed_p,
                  passport_axis = passport_axis),
    pass = pass, plots = plots,
    summary = summary_txt, llm_summary = llm)
}


# Faceted likelihood curves, each shifted so its maximum is 0: the dotted
# horizontal line at -lr_threshold/2 is the cut-off the verdict uses (the
# calibrated value is in "tension" exactly when its open circle lies below
# it); dashed vertical = calibrated value; filled point = maximiser.
.plot_d36_profiles <- function(prof_df, tab, meta, profile = "slice",
                               lr_threshold = 3.84) {
  vlev <- c("consistent", "tension", "untestable", "failed")
  llmax <- stats::setNames(tab$llmax, tab$param)
  prof_df$rel <- prof_df$ll - llmax[prof_df$param]
  # zoom: show each curve down to max(2 x cut-off, 1.2 x the calibrated drop)
  # so the cut-off line stays readable next to a very deep tension profile
  base_rel <- stats::setNames(tab$base_ll - tab$llmax, tab$param)
  floor_p  <- pmin(-lr_threshold, 1.2 * base_rel[prof_df$param], na.rm = TRUE)
  prof_df <- prof_df[is.finite(prof_df$rel) & prof_df$rel >= floor_p, , drop = FALSE]
  tab$verdict <- factor(tab$verdict, levels = vlev)
  lab <- stats::setNames(sprintf("%s (%s)", tab$param, tab$verdict), tab$param)
  plot_params <- tab$param[tab$param %in% unique(c(prof_df$param,
                                                   tab$param[tab$verdict == "failed"]))]
  lvl <- unname(lab[plot_params])
  prof_df$facet <- factor(lab[prof_df$param], levels = lvl)
  tb <- tab[tab$param %in% plot_params, , drop = FALSE]
  tb$facet <- factor(lab[tb$param], levels = lvl)
  tb$base_rel <- tb$base_ll - tb$llmax
  ok <- tb$verdict != "failed"
  cut_df <- data.frame(facet = tb$facet[ok], y = rep(-lr_threshold / 2, sum(ok)))

  cols <- c(consistent = dynhr_colours$teal, tension = dynhr_colours$red,
            untestable = dynhr_colours$grey, failed = dynhr_colours$orange)
  p <- ggplot2::ggplot(prof_df, ggplot2::aes(x = x, y = rel)) +
    ggplot2::geom_hline(data = cut_df, ggplot2::aes(yintercept = y),
                        linetype = "dotted", colour = dynhr_colours$red,
                        linewidth = 0.5) +
    ggplot2::geom_vline(data = tb, ggplot2::aes(xintercept = calibrated),
                        linetype = "dashed", colour = dynhr_colours$dark_blue,
                        linewidth = 0.4, na.rm = TRUE) +
    ggplot2::geom_line(colour = dynhr_colours$mid_blue, linewidth = 0.6) +
    ggplot2::geom_point(data = tb[ok, , drop = FALSE],
                        ggplot2::aes(x = calibrated, y = base_rel),
                        shape = 21, fill = "white", colour = dynhr_colours$dark_blue,
                        size = 2.2, na.rm = TRUE) +
    ggplot2::geom_point(data = tb[ok, , drop = FALSE],
                        ggplot2::aes(x = implied, y = 0, colour = verdict),
                        size = 2.6, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = cols, drop = TRUE,
                                 name = "Maximiser (verdict)")
  if (any(!ok)) {
    fl <- tb[!ok, , drop = FALSE]
    fl$x   <- ifelse(is.finite(fl$calibrated), fl$calibrated, 0)
    fl$rel <- 0
    p <- p + ggplot2::geom_text(data = fl, label = "loglik not evaluable",
                                colour = dynhr_colours$orange, size = 3)
  }
  p <- p +
    ggplot2::facet_wrap(~ facet, scales = "free") +
    ggplot2::scale_x_continuous(n.breaks = 4) +
    ggplot2::labs(
      title    = sprintf("D36: Calibration deepness -- %s per calibrated parameter",
                         if (profile == "full") "concentrated (profile) likelihood"
                         else "likelihood slice"),
      subtitle = sprintf(paste0(
        "Dashed line / open circle = calibrated value; filled point = value the data prefers.\n",
        "Dotted line = LR cut-off (-%.2f nats, LR %.2f): an open circle below it = tension."),
        lr_threshold / 2, lr_threshold),
      x = "Parameter value",
      y = if (profile == "full") "Concentrated loglik - max (nats)"
          else "Loglik slice - max (nats)")
  p <- p + theme_dynhr()
  .apply_meta(p, meta)
}
