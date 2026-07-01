## R/diag-deep-d36-calibration.R
## --------------------------------------------------------------------------
## D36. Calibration deepness / validation -- "is this *fixed* parameter
## actually deep, or is the calibration hiding misspecification?"
##
## Most DSGE deep parameters are calibrated, not estimated, so the rest of the
## deep-parameter suite (D1/D6/D33/D34/D35 -- which work off estimated draws or
## the estimated mode) reports them as "?" (unassessed). But a calibrated
## parameter still leaves a footprint in the likelihood: even though we fix it,
## we can ask the data what value IT would choose. D36 profiles the likelihood
## along each calibrated deep parameter (a one-dimensional "slice", holding the
## rest of the model fixed) and reports three things:
##
##   1. IDENTIFICATION -- is the profile curved at all? A flat profile means the
##      data cannot speak to the parameter; the calibration is then a pure
##      assumption ("untestable") and its deepness cannot be verified or refuted.
##   2. TENSION -- where does the likelihood peak relative to the calibrated
##      value? A likelihood-ratio test of H0: c = c_cal. If the data strongly
##      prefers a different value, either the parameter is not at its deep value
##      or the model is misspecified at the calibrated value ("tension").
##   3. CONSTANCY (optional) -- profiling on different sub-samples, does the
##      likelihood-implied value drift? A "deep" constant should not. This is the
##      Lucas-critique constancy test extended to *calibrated* parameters.
##
## Verdict per calibrated deep parameter: untestable / consistent / tension
## (+ a non-constant flag). Exposes `$result$passport_axis` for the Passport's
## new "calibrated" column.
##
## References:
##   Canova, F. (1994). Statistical inference in calibrated models. JAE, 9, S123.
##   Jorgensen, T. H. (2023). Sensitivity to calibrated parameters. REStat,
##     105(2), 474-481.
##   Hansen, L. P., & Heckman, J. J. (1996). The empirical foundations of
##     calibration. Journal of Economic Perspectives, 10(1), 87-104.
##   Mueller, U. K. (2012). Measuring prior sensitivity ... JME, 59(6), 581-597.
## --------------------------------------------------------------------------


# Default profiling range for a parameter value: probability-like values in
# (0,1) get a clamped additive window; others a multiplicative +/- frac window.
.d36_default_range <- function(v, frac = 0.5) {
  if (!is.finite(v)) return(NULL)
  if (v > 0 && v < 1) {
    c(max(0.02, v - 0.30), min(0.98, v + 0.30))
  } else if (v == 0) {
    c(-1, 1) * frac
  } else {
    s <- sort(c(v * (1 - frac), v * (1 + frac)))
    s
  }
}

# Profile one parameter: evaluate the (full-vector) loglik over a grid with
# `target` swept across `range`, holding all other params at `params`.
.d36_profile_one <- function(loglik_fn, params, target, range, n_grid) {
  grid <- seq(range[1], range[2], length.out = n_grid)
  ll <- vapply(grid, function(g) {
    p <- params; p[[target]] <- g
    v <- tryCatch(loglik_fn(p), error = function(e) NA_real_)
    if (is.numeric(v) && length(v) == 1 && is.finite(v)) v else NA_real_
  }, numeric(1))
  list(grid = grid, ll = ll)
}

# Quadratic refinement around the grid argmax -> (chat, curvature).
.d36_refine <- function(grid, ll) {
  ok <- is.finite(ll)
  if (sum(ok) < 3) return(list(chat = NA_real_, d2 = NA_real_, llmax = NA_real_))
  i <- which.max(replace(ll, !ok, -Inf))
  lo <- max(1, i - 1); hi <- min(length(grid), i + 1)
  idx <- lo:hi
  if (length(idx) < 3 || any(!is.finite(ll[idx])))
    return(list(chat = grid[i], d2 = NA_real_, llmax = ll[i]))
  fit <- tryCatch(stats::lm(ll[idx] ~ poly(grid[idx], 2, raw = TRUE)),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(chat = grid[i], d2 = NA_real_, llmax = ll[i]))
  b <- stats::coef(fit)
  a2 <- b[3]; a1 <- b[2]
  chat <- if (is.finite(a2) && a2 < 0) -a1 / (2 * a2) else grid[i]
  # keep within the local bracket
  chat <- min(max(chat, grid[lo]), grid[hi])
  llmax <- if (is.finite(a2)) b[1] + a1 * chat + a2 * chat^2 else ll[i]
  list(chat = chat, d2 = 2 * a2, llmax = max(llmax, ll[i]))
}


# Concentrated log-likelihood: fix `target` at `fixed_val`, re-optimise the
# estimated/nuisance block `est`. Returns the max loglik and the argmax (for
# warm-starting the next grid point).
.d36_concentrate <- function(loglik_fn, params, target, fixed_val, est,
                             start, lower, upper, maxit = 200) {
  if (length(est) == 0) {
    p <- params; p[[target]] <- fixed_val
    v <- tryCatch(loglik_fn(p), error = function(e) NA_real_)
    return(list(value = if (is.finite(v)) v else NA_real_, par = numeric(0)))
  }
  obj <- function(th) {
    p <- params; p[[target]] <- fixed_val
    for (i in seq_along(est)) p[[est[i]]] <- th[i]
    v <- tryCatch(loglik_fn(p), error = function(e) NA_real_)
    if (is.finite(v)) -v else 1e10
  }
  use_box <- !is.null(lower) && !is.null(upper) &&
             all(is.finite(lower)) && all(is.finite(upper))
  opt <- tryCatch(suppressWarnings(
    if (use_box)
      stats::optim(start, obj, method = "L-BFGS-B", lower = lower,
                   upper = upper, control = list(maxit = maxit))
    else
      stats::optim(start, obj, method = "Nelder-Mead",
                   control = list(maxit = maxit))),
    error = function(e) NULL)
  if (is.null(opt) || !is.finite(opt$value) || opt$value >= 1e10)
    list(value = NA_real_, par = start)
  else list(value = -opt$value, par = opt$par)
}

# Full profile: sweep `target` across `range`, re-optimising `est` at each
# point (warm-started). Also returns the concentrated value at `base_val`.
.d36_profile_full <- function(loglik_fn, params, target, range, n_grid,
                              est, lower, upper, base_val) {
  grid <- seq(range[1], range[2], length.out = n_grid)
  ll   <- rep(NA_real_, n_grid)
  warm <- unlist(params[est])
  for (i in seq_along(grid)) {
    r <- .d36_concentrate(loglik_fn, params, target, grid[i], est, warm, lower, upper)
    ll[i] <- r$value
    if (length(r$par) && all(is.finite(r$par))) warm <- r$par
  }
  base <- .d36_concentrate(loglik_fn, params, target, base_val,
                           est, unlist(params[est]), lower, upper)$value
  list(grid = grid, ll = ll, base = base)
}


# ---------------------------------------------------------------------------
#' D36. Calibration deepness / validation
#'
#' @param model     Parsed model (for the @dynhr:deep taxonomy). Optional if
#'   \code{deep_spec} is supplied.
#' @param deep_spec Optional \code{\link{build_deep_spec}}.
#' @param params    Named numeric full parameter vector (the calibration).
#' @param loglik_fn Function: a *full* named parameter vector -> scalar
#'   log-likelihood (re-solving the model as needed).
#' @param targets   Character vector of parameters to profile. Default: the
#'   calibrated deep primitives (deep, role = primitive, not in
#'   \code{estimated}).
#' @param estimated Character vector of estimated parameter names (excluded from
#'   the default \code{targets}).
#' @param ranges    Optional named list of \code{c(lo, hi)} profiling ranges.
#' @param n_grid    Grid points per profile (default 21).
#' @param frac      Default half-width for multiplicative ranges (default 0.5).
#' @param lr_threshold Likelihood-ratio cutoff for "tension" (default 3.84 =
#'   chi-square(1) 95\%).
#' @param min_curv_nats Minimum profile height (nats) over the range to call a
#'   parameter identifiable rather than "untestable" (default 1).
#' @param profile \code{"slice"} (default) holds the rest of the model fixed
#'   while sweeping the target -- a cheap *conditional* profile. \code{"full"}
#'   re-optimises the \code{estimated} nuisance block at every grid point (the
#'   concentrated likelihood), which separates genuine calibration tension from
#'   a target that merely trades off with a nuisance parameter; it requires
#'   \code{estimated} and is much more expensive.
#' @param estimated_bounds Optional named list of \code{c(lo, hi)} box
#'   constraints for the \code{estimated} block (enables L-BFGS-B in the full
#'   profile; otherwise Nelder-Mead).
#' @param subsample_loglik_fns Optional named list of per-regime full-vector
#'   log-likelihood functions for the constancy check.
#' @param meta      Optional \code{\link{diag_meta}} provenance descriptor.
#' @return A \code{dynhr_diagnostic}; \code{$result$passport_axis} is a named
#'   logical (TRUE = consistent/deep, FALSE = tension) over the calibrated
#'   parameters for the Passport's "calibrated" column.
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
                                     profile = c("slice", "full"),
                                     estimated_bounds = NULL,
                                     subsample_loglik_fns = NULL,
                                     meta = NULL) {
    profile <- match.arg(profile)
  tryCatch({
    if (is.null(loglik_fn))
      return(.make_result(pass = NA,
        summary = "D36 calibration deepness: no loglik_fn supplied."))
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

    base_ll <- tryCatch(loglik_fn(params), error = function(e) NA_real_)
    if (!is.finite(base_ll))
      return(.make_result(pass = NA,
        summary = "D36 calibration deepness: baseline log-likelihood is not finite."))
    if (profile == "full" && (is.null(estimated) || length(estimated) == 0))
      return(.make_result(pass = NA,
        summary = "D36 (profile='full'): supply `estimated` (the nuisance block to re-optimise at each grid point)."))

    .bounds_for <- function(est_tg) {
      if (is.null(estimated_bounds)) return(list(lo = NULL, hi = NULL))
      lo <- vapply(est_tg, function(e) (estimated_bounds[[e]] %||% c(NA, NA))[1], numeric(1))
      hi <- vapply(est_tg, function(e) (estimated_bounds[[e]] %||% c(NA, NA))[2], numeric(1))
      list(lo = lo, hi = hi)
    }

    rows <- list(); profiles <- list()
    for (tg in targets) {
      v   <- suppressWarnings(as.numeric(params[[tg]]))
      rng <- if (!is.null(ranges[[tg]])) ranges[[tg]] else .d36_default_range(v, frac)
      if (is.null(rng)) next
      if (profile == "full") {
        est_tg <- setdiff(estimated, tg)
        bd <- .bounds_for(est_tg)
        pf <- .d36_profile_full(loglik_fn, params, tg, rng, n_grid, est_tg,
                                bd$lo, bd$hi, base_val = v)
        pr <- list(grid = pf$grid, ll = pf$ll)
        base_for_lr <- pf$base
      } else {
        pr <- .d36_profile_one(loglik_fn, params, tg, rng, n_grid)
        base_for_lr <- base_ll
      }
      rf  <- .d36_refine(pr$grid, pr$ll)
      height <- diff(range(pr$ll, na.rm = TRUE))
      se  <- if (is.finite(rf$d2) && rf$d2 < 0) sqrt(-1 / rf$d2) else Inf
      lr  <- if (is.finite(rf$llmax) && is.finite(base_for_lr))
        2 * (rf$llmax - base_for_lr) else NA_real_
      lr  <- max(lr, 0, na.rm = TRUE)
      pval <- stats::pchisq(lr, df = 1, lower.tail = FALSE)
      at_bound <- is.finite(rf$chat) &&
        (abs(rf$chat - rng[1]) < 1e-9 || abs(rf$chat - rng[2]) < 1e-9)

      identifiable <- is.finite(height) && height >= min_curv_nats
      verdict <- if (!identifiable) "untestable"
                 else if (lr >= lr_threshold) "tension"
                 else "consistent"

      profiles[[tg]] <- data.frame(param = tg, x = pr$grid, ll = pr$ll,
                                   stringsAsFactors = FALSE)
      rows[[tg]] <- data.frame(
        param = tg, calibrated = v, implied = rf$chat,
        profile_se = se, height_nats = height, LR = lr, p_value = pval,
        at_bound = at_bound, verdict = verdict, stringsAsFactors = FALSE)
    }
    tab <- do.call(rbind, rows)

    # optional cross-regime constancy of the implied value
    constancy <- NULL
    if (!is.null(subsample_loglik_fns) && length(subsample_loglik_fns) >= 2) {
      cc <- lapply(names(subsample_loglik_fns), function(rn) {
        f <- subsample_loglik_fns[[rn]]
        vapply(tab$param, function(tg) {
          v <- as.numeric(params[[tg]])
          rng <- if (!is.null(ranges[[tg]])) ranges[[tg]] else .d36_default_range(v, frac)
          pr <- .d36_profile_one(f, params, tg, rng, n_grid)
          .d36_refine(pr$grid, pr$ll)$chat
        }, numeric(1))
      })
      cmat <- do.call(cbind, cc); colnames(cmat) <- names(subsample_loglik_fns)
      rownames(cmat) <- tab$param
      spread <- apply(cmat, 1, function(z) diff(range(z, na.rm = TRUE)))
      non_constant <- tab$param[is.finite(spread) & is.finite(tab$profile_se) &
                                spread > 2 * tab$profile_se]
      constancy <- list(implied_by_regime = cmat, spread = spread,
                        non_constant = non_constant)
      tab$non_constant <- tab$param %in% non_constant
    } else {
      tab$non_constant <- FALSE
    }

    tension     <- tab$param[tab$verdict == "tension"]
    untestable  <- tab$param[tab$verdict == "untestable"]
    consistent  <- tab$param[tab$verdict == "consistent"]
    non_const   <- tab$param[tab$non_constant]

    # passport axis: consistent (+constant) -> TRUE, tension/non-constant ->
    # FALSE, untestable -> omitted (stays "?")
    ax_params <- tab$param[tab$verdict != "untestable"]
    ax_vals   <- (tab$verdict[tab$verdict != "untestable"] == "consistent") &
                 !tab$non_constant[tab$verdict != "untestable"]
    passport_axis <- stats::setNames(ax_vals, ax_params)

    pass <- length(tension) == 0 && length(non_const) == 0

    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE) && length(profiles) > 0)
      plots$profiles <- .plot_d36_profiles(do.call(rbind, profiles), tab, meta)

    summary_txt <- sprintf(
      "D36 Calibration deepness (%s profile): %d calibrated deep param(s) -- %d consistent, %d tension, %d untestable%s.",
      profile, nrow(tab), length(consistent), length(tension), length(untestable),
      if (length(non_const)) sprintf(", %d non-constant", length(non_const)) else "")

    badge <- if (is.na(pass)) "INFO" else if (pass) "PASS" else "FAIL"
    llm <- paste(c(
      sprintf("D36 | Calibration Deepness | %s", badge),
      sprintf("  calibrated_deep=%d consistent=%d tension=%d untestable=%d non_constant=%d",
              nrow(tab), length(consistent), length(tension), length(untestable),
              length(non_const)),
      {
        worst <- tab[order(-tab$LR), , drop = FALSE]; worst <- utils::head(worst, 4)
        sprintf("  implied_vs_calibrated: %s",
                paste(sprintf("%s cal=%.3g data=%.3g LR=%.1f", worst$param,
                              worst$calibrated, worst$implied, worst$LR),
                      collapse = "; "))
      },
      if (length(tension))
        sprintf("  tension(data disputes the fixed value): %s",
                paste(tension, collapse = ", ")),
      if (length(non_const))
        sprintf("  non_constant(implied value drifts across regimes): %s",
                paste(non_const, collapse = ", ")),
      if (length(untestable))
        sprintf("  untestable(flat profile, data uninformative): %s",
                paste(utils::head(untestable, 8), collapse = ", ")),
      sprintf("  action: %s",
              if (isTRUE(pass))
                "Calibrated deep parameters are either data-consistent or innocuous (flat profile); no calibration tension."
              else sprintf("%s: the data prefers a different value than the calibration -- re-calibrate, estimate it, or treat the gap as a misspecification signal.",
                           paste(utils::head(c(tension, non_const), 3), collapse = ", ")))
    ), collapse = "\n")

    .make_result(
      result = list(table = tab, profiles = profiles, base_loglik = base_ll,
                    constancy = constancy, tension = tension, profile = profile,
                    untestable = untestable, passport_axis = passport_axis),
      pass = pass, plots = plots,
      summary = summary_txt, llm_summary = llm)
  }, error = function(e) {
    .make_result(pass = NA, errored = TRUE,
                 summary = paste("D36 calibration deepness: ERROR --",
                                 conditionMessage(e)))
  })
}


# Faceted profile-likelihood curves with the calibrated value (dashed) and the
# likelihood-implied value (point, coloured by verdict).
.plot_d36_profiles <- function(prof_df, tab, meta) {
  vmap <- stats::setNames(tab$verdict, tab$param)
  prof_df$verdict <- vmap[prof_df$param]
  prof_df$facet   <- sprintf("%s (%s)", prof_df$param, prof_df$verdict)
  cal <- data.frame(param = tab$param,
                    facet = sprintf("%s (%s)", tab$param, tab$verdict),
                    calibrated = tab$calibrated, implied = tab$implied,
                    verdict = tab$verdict, stringsAsFactors = FALSE)
  # implied point y-position = facet max ll
  ymax <- tapply(prof_df$ll, prof_df$facet, function(z) max(z, na.rm = TRUE))
  cal$y <- as.numeric(ymax[cal$facet])

  p <- ggplot2::ggplot(prof_df, ggplot2::aes(x = x, y = ll)) +
    ggplot2::geom_line(colour = dynhr_colours$mid_blue, linewidth = 0.6) +
    ggplot2::geom_vline(data = cal,
                        ggplot2::aes(xintercept = calibrated),
                        linetype = "dashed", colour = dynhr_colours$grey,
                        linewidth = 0.5, inherit.aes = FALSE) +
    ggplot2::geom_point(data = cal,
                        ggplot2::aes(x = implied, y = y, colour = verdict),
                        size = 2.4, inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(
      values = c(consistent = dynhr_colours$teal, tension = dynhr_colours$red,
                 untestable = dynhr_colours$grey), drop = FALSE, name = NULL) +
    ggplot2::facet_wrap(~ facet, scales = "free") +
    ggplot2::labs(
      title    = "D36: Calibration deepness -- likelihood profile per calibrated parameter",
      subtitle = "dashed = calibrated value; point = value the data prefers; a peak far from the dashed line = tension",
      x = "Parameter value", y = "Profile log-likelihood")
  p <- p + theme_dynhr()
  .apply_meta(p, meta)
}
