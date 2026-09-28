## R/diag-post-d8-irf.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D8 IRF plausibility checks; .get_irfs_long() helper
## --------------------------------------------------------------------------

.get_irfs_long <- function(irfs) {
  ## A StochSimulResult-like wrapper carries the collection in $irfs. Exact
  ## lookup: `$` partial-matches on lists, so a shock whose name merely
  ## starts with "irfs" must not be mistaken for the wrapper slot.
  inner <- if (is.list(irfs)) irfs[["irfs", exact = TRUE]] else NULL
  if (is.list(inner)) irfs <- inner
  if (!is.list(irfs) || length(irfs) == 0)
    .dynhr_abort("irfs must be a named list of matrices (IRFCollection)")
  shock_names <- names(irfs)
  rows <- list()
  for (sh in shock_names) {
    mat <- irfs[[sh]]
    if (!is.matrix(mat)) next
    var_nms <- colnames(mat)
    if (is.null(var_nms)) var_nms <- paste0("V", seq_len(ncol(mat)))
    n_periods <- nrow(mat)
    for (j in seq_along(var_nms)) {
      rows[[length(rows) + 1L]] <- data.frame(
        horizon  = seq_len(n_periods),
        variable = var_nms[j],
        shock    = sh,
        value    = mat[, j],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}


## Order-1 IRFs at a parameter draw `theta`, with shock scaling taken from
## the SAME params the likelihood uses.
##
## theta is applied by .solve_dr_for_theta() via .apply_theta_to_params(), so
## an estimated `stderr <shock>` (named after the shock) reaches
## .get_shock_cov() and each shock is scaled by ITS OWN std at theta (ghu is a
## unit-shock matrix -- the ghu unit-shock convention). The model is re-solved
## at theta through the likelihood's own solve pipeline, so ghx/ghu are also
## evaluated at theta. dr$Sigma_e is cleared: when present it is a solve-time
## snapshot of the CALIBRATED covariance, and compute_irfs() would prefer it
## over params (params stay authoritative, as on the Kalman path).
##
## Returns an IRFCollection, or NULL when theta is infeasible (no steady
## state / BK violation).
.d8_irfs_at_theta <- function(model, compiled, theta, n_periods = 40L,
                              sys_cache = NULL, state = NULL) {
  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence
  if (is.null(sys_cache)) sys_cache <- cache_system_structure(compiled)
  if (is.null(state)) {
    state <- new.env(parent = emptyenv())
    state$ss_warm <- NULL
  }
  ## lik_init = "diffuse": an IRF is well defined for a unit-root draw, so do
  ## not apply the stationary-likelihood rejection.
  sol <- .solve_dr_for_theta(model, compiled, sys_cache, theta, state,
                             lik_init = "diffuse")
  if (is.null(sol)) return(NULL)
  dr <- sol$dr
  dr$Sigma_e <- NULL
  compute_irfs(dr, model, n_periods = n_periods, params = sol$params)
}


## Default D8 benchmarks (0.9.4): SIGN ONLY, over a 1-8 period window.
##
## Neither the peak-timing window nor the peak-magnitude range is a default
## any more; both are opt-in fields of a user benchmark.
##
##  * Magnitude. The IRFs are responses to a one-standard-deviation shock in
##    the model's own units (log-deviation, percent, annualised percent ...),
##    so a fixed default range (the pre-0.9.4 defaults assumed a 25bp shock in
##    decimal units) compares incommensurable numbers and FLAGs almost every
##    model.
##  * Timing. The CEE (2005) 6-8 quarter hump is a PRODUCT of the frictions
##    they estimate (habits, investment adjustment costs, Calvo indexation),
##    not a restriction imposed on identification; Ramey (2016) finds peak
##    timing varies with the identification scheme and sample; Uhlig (2005)
##    deliberately leaves output's dynamic response unrestricted. A textbook
##    three-equation NK model without those frictions peaks on impact and is
##    not thereby misspecified. Dynare's `irf_calibration` block likewise
##    ships no default timing restrictions -- they are user-supplied.
##
## `sign_window = c(1, 8)`: the dominant response whose sign is tested is
## taken over horizons 1-8 rather than the whole 40-period path, so a
## correctly-signed short-run response is not overturned by a large, slow
## overshoot far out. The window length (8) is a PACKAGE CHOICE covering the
## 1-2 year range the monetary-shock literature reports; it has no single
## literature source.
.d8_default_benchmarks <- function() {
  mp_shocks   <- c("eps_r", "eps_mp", "eps_m", "eps_r_", "e_r", "e_m",
                   "e_mp", "er", "em")
  tech_shocks <- c("eps_a", "eps_a_", "eps_tech", "e_a", "ea")
  uip_shocks  <- c("eps_uip_", "eps_uip", "eps_rp", "e_uip", "e_rp", "es")
  list(
    monetary = list(
      shock = mp_shocks, variable = "y", direction = "negative",
      peak_range = NULL, peak_horizon = NULL, sign_window = c(1, 8),
      description = "Monetary contraction -> y falls within 1-8q"
    ),
    monetary_pi = list(
      shock = mp_shocks, variable = "pi", direction = "negative",
      peak_range = NULL, peak_horizon = NULL, sign_window = c(1, 8),
      description = "Monetary contraction -> pi falls within 1-8q"
    ),
    technology = list(
      shock = tech_shocks, variable = "y", direction = "positive",
      peak_range = NULL, peak_horizon = NULL, sign_window = c(1, 8),
      description = "Positive TFP shock -> y rises within 1-8q"
    ),
    uip = list(
      shock = uip_shocks, variable = "q", direction = "positive",
      peak_range = NULL, peak_horizon = NULL, sign_window = c(1, 8),
      description = "UIP shock -> RER depreciates within 1-8q"
    )
  )
}

#' Opt-in D8 benchmark presets
#'
#' The stylised peak-timing windows that were D8's pre-0.9.4 defaults, kept
#' as an explicit preset. Supply the result (optionally edited) as
#' \code{d8_irf_plausibility(benchmarks = )} when your model HAS the frictions
#' that generate hump-shaped responses (habit formation, investment adjustment
#' costs, price indexation); a frictionless model peaks on impact and will be
#' FLAGged by these windows, which is why they are not the default.
#'
#' The windows (output trough 4-8q, inflation trough 6-12q) follow Christiano,
#' Eichenbaum & Evans (2005) and the range surveyed by Ramey (2016); they are
#' a PACKAGE CHOICE of round numbers over that range, not a single published
#' set of bounds.
#'
#' @return A named list of benchmark specifications in
#'   \code{d8_irf_plausibility()}'s \code{benchmarks} format.
#' @references
#'   Christiano, L. J., Eichenbaum, M., & Evans, C. L. (2005). Nominal
#'   rigidities and the dynamic effects of a shock to monetary policy.
#'   \emph{Journal of Political Economy}, 113(1), 1-45.
#'   Ramey, V. A. (2016). Macroeconomic shocks and their propagation. In
#'   \emph{Handbook of Macroeconomics}, Vol. 2A, ch. 2. Elsevier.
#' @noRd
d8_hump_benchmarks <- function() {
  b <- .d8_default_benchmarks()
  b$monetary$peak_horizon <- c(4, 8)
  b$monetary$sign_window  <- NULL
  b$monetary$description  <- "Monetary contraction -> y falls, trough in 4-8q"
  b$monetary_pi$peak_horizon <- c(6, 12)
  b$monetary_pi$sign_window  <- NULL
  b$monetary_pi$description  <- "Monetary contraction -> pi falls, trough in 6-12q"
  b$technology$peak_horizon <- c(1, 6)
  b$technology$sign_window  <- NULL
  b$technology$description  <- "Positive TFP shock -> y rises, peak in 1-6q"
  b$uip$peak_horizon <- c(1, 4)
  b$uip$sign_window  <- NULL
  b$uip$description  <- "UIP shock -> RER depreciates, peak in 1-4q"
  b
}

#' D8. IRF plausibility checks
#'
#' Compares model impulse response functions against stylised economic
#' benchmarks.  Each benchmark names a (shock, variable) pair and its expected
#' direction; a peak-magnitude range and a peak-horizon window are OPTIONAL
#' fields.  The diagnostic SKIPs benchmarks whose shock or variable is absent
#' from the IRF collection and counts PASS, FLAG, and FAIL results for the
#' remaining checks.
#'
#' The peak of a response is its largest ABSOLUTE value (earliest horizon on
#' ties); horizon 1 is the impact period.  The sign check asks whether that
#' dominant response has the expected sign.  When the benchmark carries a
#' \code{sign_window} the dominant response is taken over that horizon range
#' only.
#'
#' \strong{The shipped defaults are SIGN-ONLY over horizons 1-8} (see
#' \code{.d8_default_benchmarks()}): no peak-timing window and no magnitude
#' range.  Timing is an ESTIMATED OUTPUT of a model, not an a priori
#' restriction.  The familiar 4-8 quarter output trough comes from Christiano,
#' Eichenbaum & Evans (2005), where the hump is produced by the frictions they
#' estimate (habit formation, investment adjustment costs, Calvo pricing with
#' indexation) rather than imposed; Ramey (2016) surveys the monetary-shock
#' evidence and finds peak timing (and, for some variables, peak sign) varies
#' materially with the identification scheme and sample; Uhlig (2005) imposes
#' sign restrictions only on impact through h = 5 and deliberately leaves
#' output's response unrestricted so its dynamics are estimated.  Dynare's
#' \code{irf_calibration} block, the closest published analogue, ships no
#' default restrictions at all -- they are user-supplied.  A frictionless
#' three-equation New Keynesian model peaks on impact and is not misspecified
#' for doing so, so a default timing window would FAIL a broad class of
#' correct models.  \code{d8_hump_benchmarks()} returns the pre-0.9.4 timing
#' windows for models that DO have those frictions.
#'
#' Magnitude ranges are likewise opt-in: the IRFs are responses to a
#' one-standard-deviation shock in the model's own units (log-deviation,
#' percent, annualised percent, ...), so a fixed default range compares
#' incommensurable numbers.
#'
#' The \code{sign_window} upper bound of 8 periods on the defaults, and the
#' 25\% FLAG allowance at the gate, are PACKAGE CHOICES with no single
#' literature source.
#'
#' Pass gate: \code{pass = (n_fail == 0) && (n_flag <= floor(0.25 * n_checked))},
#' and \code{pass = NA} (INFO) when no benchmark could be checked.
#' A FAIL is a wrong sign or a non-finite IRF; a FLAG is a correct sign with
#' peak magnitude or timing outside the benchmark range.  Up to 25\% of
#' checked benchmarks may be flagged without triggering a FAIL at the gate
#' level -- one benchmark in four is a lenient but non-zero bar.
#'
#' @param irf_data   Named list of matrices (one per shock, rows = horizons,
#'   columns = observable variables), or a data frame with columns
#'   \code{horizon}, \code{variable}, \code{shock}, \code{value}.
#'   A list element can also be an IRFCollection object with \code{$irfs}.
#' @param horizons   Integer: maximum horizon to evaluate (default 40).
#'   Responses beyond this horizon are dropped before benchmark checks.
#' @param benchmarks Named list of benchmark specifications.  Each element
#'   must be a list with:
#'   \describe{
#'     \item{\code{shock}}{Character vector of candidate shock names
#'       (first match in \code{irf_data} is used).}
#'     \item{\code{variable}}{Target observable name.}
#'     \item{\code{direction}}{Either \code{"positive"} or \code{"negative"}.}
#'     \item{\code{peak_range}}{Numeric vector c(lo, hi) for the expected peak
#'       value range (optional; NULL to skip magnitude check).}
#'     \item{\code{peak_horizon}}{Integer vector c(lo, hi) for the expected
#'       peak-horizon range (OPTIONAL, default NULL = no timing check).}
#'     \item{\code{sign_window}}{Integer vector c(lo, hi) restricting the
#'       horizons the dominant response is taken over (optional; NULL = the
#'       whole checked horizon).}
#'     \item{\code{description}}{Human-readable description string.}
#'   }
#'   When NULL, the sign-only defaults are used (monetary policy, technology
#'   and UIP shocks for output, inflation and the real exchange rate), each
#'   with \code{sign_window = c(1, 8)} and no peak/timing range.
#' @param metadata  Optional metadata list for plot annotation.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{Named list of per-benchmark check results, each containing
#'     \code{status} (PASS/FLAG/FAIL/SKIP), \code{peak_value},
#'     \code{peak_horizon}, and \code{direction_ok}, \code{magnitude_ok},
#'     \code{timing_ok} flags.}
#'   \item{pass}{Logical -- TRUE if no FAILs and FLAG count <= 25\% of checked;
#'     NA when every benchmark was SKIPped.}
#'   \item{plots}{ggplot2 panel of matched IRFs with peak markers.}
#'   \item{summary}{Human-readable per-benchmark summary lines.}
#'
#' @references
#'   Christiano, L. J., Eichenbaum, M., & Evans, C. L. (2005). Nominal
#'   rigidities and the dynamic effects of a shock to monetary policy.
#'   \emph{Journal of Political Economy}, 113(1), 1-45.
#'   Ramey, V. A. (2016). Macroeconomic shocks and their propagation. In
#'   \emph{Handbook of Macroeconomics}, Vol. 2A, ch. 2. Elsevier.
#'   Uhlig, H. (2005). What are the effects of monetary policy on output?
#'   Results from an agnostic identification procedure. \emph{Journal of
#'   Monetary Economics}, 52(2), 381-419.
#' @noRd
d8_irf_plausibility <- function(irf_data, horizons = 40,
                                benchmarks = NULL, metadata = NULL) {
    if (is.data.frame(irf_data) &&
        all(c("horizon","variable","shock","value") %in% names(irf_data))) {
      irf_long <- irf_data
    } else {
      irf_long <- .get_irfs_long(irf_data)
    }
    irf_long <- irf_long[irf_long$horizon <= horizons, ]
    available_shocks <- unique(irf_long$shock)
    available_vars   <- unique(irf_long$variable)

    if (is.null(benchmarks)) benchmarks <- .d8_default_benchmarks()

    check_results <- list()
    for (bname in names(benchmarks)) {
      b <- benchmarks[[bname]]
      matched_shock <- intersect(b$shock, available_shocks)
      if (length(matched_shock) == 0) {
        check_results[[bname]] <- list(
          benchmark = bname, status = "SKIP",
          message = sprintf("No matching shock (tried: %s)",
                            paste(b$shock, collapse = ", "))
        )
        next
      }
      sh <- matched_shock[1]
      target_var <- b$variable
      if (!(target_var %in% available_vars)) {
        aliases <- list(pi = "dp", q = "ds", dp = "pi", ds = "q",
                        y = "gdp", gdp = "y")
        alt <- aliases[[target_var]]
        if (!is.null(alt) && alt %in% available_vars) target_var <- alt
      }
      sub <- irf_long[irf_long$shock == sh & irf_long$variable == target_var, ]
      if (nrow(sub) == 0) {
        check_results[[bname]] <- list(
          benchmark = bname, status = "SKIP",
          message = sprintf("Variable '%s' not in IRF for shock '%s'",
                            b$variable, sh)
        )
        next
      }
      if (!(b$direction %in% c("positive", "negative")))
        .dynhr_abort(sprintf(
          "D8 benchmark '%s': direction must be \"positive\" or \"negative\", not \"%s\".",
          bname, paste(b$direction, collapse = ",")))
      if (any(!is.finite(sub$value))) {
        check_results[[bname]] <- list(
          benchmark = bname, status = "FAIL", shock = sh, variable = target_var,
          peak_value = NA_real_, peak_horizon = NA_integer_,
          direction_ok = FALSE, magnitude_ok = FALSE, timing_ok = FALSE,
          message = sprintf("shock=%s var=%s IRF has %d non-finite value(s)",
                            sh, target_var, sum(!is.finite(sub$value))),
          description = b$description)
        next
      }
      ## The PEAK is the largest absolute response; the sign check asks
      ## whether that dominant response has the expected sign. (Taking
      ## which.min for a "negative" benchmark let a large positive response
      ## with a small late undershoot pass as "falls".) Ties: earliest horizon.
      sub      <- sub[order(sub$horizon), , drop = FALSE]
      ## `sign_window` (opt-in; set on the defaults) restricts the horizons the
      ## dominant response is taken over, so a correctly-signed short-run
      ## response is not overturned by a large, slow overshoot far out.
      if (!is.null(b$sign_window)) {
        win_sub <- sub[sub$horizon >= b$sign_window[1] &
                         sub$horizon <= b$sign_window[2], , drop = FALSE]
        if (nrow(win_sub) == 0L) {
          check_results[[bname]] <- list(
            benchmark = bname, status = "SKIP",
            message = sprintf(
              "sign_window [%g, %g] contains no horizon <= %g",
              b$sign_window[1], b$sign_window[2], horizons))
          next
        }
        sub <- win_sub
      }
      peak_idx <- which.max(abs(sub$value))
      peak_val <- sub$value[peak_idx]
      peak_h   <- sub$horizon[peak_idx]
      dir_ok  <- if (b$direction == "negative") peak_val < 0 else peak_val > 0
      mag_ok  <- is.null(b$peak_range) || (peak_val >= b$peak_range[1] & peak_val <= b$peak_range[2])
      time_ok <- is.null(b$peak_horizon) || (peak_h >= b$peak_horizon[1] & peak_h <= b$peak_horizon[2])
      # --- ENHANCED: Use magnitude/timing for FLAG status, not just direction ---
      if (!dir_ok) {
        status <- "FAIL"
      } else if (!mag_ok || !time_ok) {
        status <- "FLAG"
      } else {
        status <- "PASS"
      }
      check_results[[bname]] <- list(
        benchmark = bname, status = status,
        shock = sh, variable = target_var,
        peak_value = peak_val, peak_horizon = peak_h,
        direction_ok = dir_ok, magnitude_ok = mag_ok, timing_ok = time_ok,
        expected_horizon = b$peak_horizon, expected_range = b$peak_range,
        sign_window = b$sign_window,
        message = sprintf("shock=%s var=%s peak=%.4f at h=%d dir=%s mag=%s time=%s",
                          sh, target_var, peak_val, peak_h,
                          ifelse(dir_ok,"ok","WRONG"),
                          ifelse(mag_ok,"ok",sprintf("outside [%.4f,%.4f]",
                                 b$peak_range[1], b$peak_range[2])),
                          ifelse(time_ok,"ok",sprintf("outside [%d,%d]",
                                 b$peak_horizon[1], b$peak_horizon[2]))),
        description = b$description
      )
    }

    # --- ENHANCED: Pass if direction OK AND flags <= 25% of checked ---
    # FLAG (wrong magnitude or timing) counts against the gate so that
    # models with correct sign but badly miscalibrated responses do not
    # silently PASS.  Gate: pass = (n_fail == 0 && n_flag <= floor(0.25 * n_checked))
    # where n_checked = benchmarks that were not SKIPped.
    # Reference: Christiano, Eichenbaum & Evans (1999) monetary-policy IRF
    # benchmarks; FLAG threshold 25% is one benchmark in four -- a lenient
    # but non-zero bar.
    n_fail    <- sum(vapply(check_results, function(r) identical(r$status, "FAIL"), logical(1)))
    n_flag    <- sum(vapply(check_results, function(r) identical(r$status, "FLAG"), logical(1)))
    n_pass    <- sum(vapply(check_results, function(r) identical(r$status, "PASS"), logical(1)))
    n_skip    <- sum(vapply(check_results, function(r) identical(r$status, "SKIP"), logical(1)))
    n_total   <- length(check_results)
    n_checked <- n_total - n_skip
    ## Nothing checked (no benchmark shock/variable in the IRFs) is not a
    ## PASS: it is uninformative (INFO), e.g. a model whose shocks are named
    ## outside the default alias lists.
    pass      <- if (n_checked == 0L) NA else
      (n_fail == 0L && n_flag <= floor(0.25 * n_checked))

    summary_lines <- vapply(check_results, function(x)
      sprintf("  [%s] %s: %s", x$status, x$benchmark, x$message), character(1))

    # --- Plot: matched benchmark IRFs with the peak marked, coloured by status ---
    plots <- list()
    matched <- Filter(function(r) !is.null(r$shock) && !identical(r$status, "SKIP"),
                      check_results)
    if (length(matched) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
      panel_of <- function(r) sprintf("%s: %s -> %s  [%s]", r$benchmark,
                                      r$shock, r$variable, r$status)
      rows <- do.call(rbind, lapply(matched, function(r) {
        sub <- irf_long[irf_long$shock == r$shock & irf_long$variable == r$variable &
                          is.finite(irf_long$value), ]
        if (nrow(sub) == 0) return(NULL)
        data.frame(panel = panel_of(r), horizon = sub$horizon, value = sub$value,
                   status = r$status, stringsAsFactors = FALSE)
      }))
      if (!is.null(rows) && nrow(rows) > 0) {
        pk <- do.call(rbind, lapply(matched, function(r) {
          if (!isTRUE(is.finite(r$peak_value))) return(NULL)
          data.frame(panel = panel_of(r), horizon = r$peak_horizon,
                     value = r$peak_value, status = r$status,
                     stringsAsFactors = FALSE)
        }))
        ## The expected peak-horizon window (and peak range, when given) the
        ## badge tests, shaded per panel.
        win <- do.call(rbind, lapply(matched, function(r) {
          if (is.null(r$expected_horizon) && is.null(r$expected_range) &&
              is.null(r$sign_window)) return(NULL)
          h_rng <- r$expected_horizon %||% r$sign_window %||% c(-Inf, Inf)
          v_rng <- r$expected_range %||% c(-Inf, Inf)
          data.frame(panel = panel_of(r),
                     xmin = h_rng[1] - 0.5, xmax = h_rng[2] + 0.5,
                     ymin = v_rng[1], ymax = v_rng[2], stringsAsFactors = FALSE)
        }))
        status_cols <- c(PASS = dynhr_colours$green,
                         FLAG = dynhr_colours$orange, FAIL = dynhr_colours$red)
        p_irf <- ggplot2::ggplot(rows, ggplot2::aes(x = horizon, y = value)) +
          { if (!is.null(win)) ggplot2::geom_rect(
              data = win, inherit.aes = FALSE,
              ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = dynhr_na_fill, alpha = 0.6) } +
          geom_dynhr_zero() +
          ggplot2::geom_line(ggplot2::aes(colour = status), linewidth = 0.7) +
          { if (!is.null(pk)) ggplot2::geom_point(
              data = pk, ggplot2::aes(x = horizon, y = value, colour = status),
              size = 2.4) } +
          { if (!is.null(pk)) ggplot2::geom_text(
              data = pk,
              ggplot2::aes(x = horizon, y = value, colour = status,
                           label = sprintf("peak %.3g @ h=%d", value,
                                           as.integer(horizon))),
              vjust = ifelse(pk$value >= 0, -0.8, 1.8), hjust = -0.05,
              size = 2.9, show.legend = FALSE) } +
          ggplot2::facet_wrap(~ panel, scales = "free_y") +
          ggplot2::coord_cartesian(clip = "off") +
          ggplot2::scale_colour_manual(values = status_cols, name = NULL,
                                       drop = TRUE) +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.18)) +
          theme_dynhr_compact() +
          ggplot2::theme(
            axis.text.y  = ggplot2::element_text(size = ggplot2::rel(0.75)),
            axis.title.y = ggplot2::element_text(angle = 90, size = ggplot2::rel(0.85)),
            panel.spacing.x = ggplot2::unit(1.2, "lines"),
            panel.spacing.y = ggplot2::unit(0.8, "lines"),
            legend.position = "none") +
          ggplot2::labs(
            title = "D8: Benchmark IRFs -- sign / timing / magnitude",
            subtitle = paste0("Dot = largest absolute response (its sign is the sign check); ",
                              "shaded = expected peak window"),
            x = "Horizon (1 = impact)",
            y = "Response to a 1 s.d. shock (model units)")
        plots$benchmark_irfs <- .apply_meta(p_irf, metadata)
      }
    }

    structure(
      list(
        result  = list(checks = check_results, irf_long = irf_long),
        pass    = pass,
        plots   = plots,
        summary = paste(c("D8 IRF plausibility:", summary_lines), collapse = "\n"),
        llm_summary = {
          badge <- .badge_str(list(pass = pass))
          fail_details <- vapply(
            Filter(function(r) identical(r$status, "FAIL"), check_results),
            function(r) sprintf("%s: %s", r$benchmark, r$message %||% "failed"),
            character(1)
          )
          flag_details <- vapply(
            Filter(function(r) identical(r$status, "FLAG"), check_results),
            function(r) sprintf("%s: %s", r$benchmark, r$message %||% "flagged"),
            character(1)
          )
          paste(c(
            sprintf("D8 | IRF Plausibility | %s", badge),
            sprintf("  benchmarks: total=%d pass=%d flag=%d fail=%d skip=%d",
                    n_total, n_pass, n_flag, n_fail, n_skip),
            if (length(fail_details) > 0)
              paste("  failed:", paste(fail_details, collapse = "; ")),
            if (length(flag_details) > 0)
              paste("  flagged (mag/timing):", paste(flag_details, collapse = "; ")),
            sprintf("  action: %s",
                    if (n_checked == 0L)
                      "no benchmark matched the IRF shocks/variables -- pass `benchmarks` naming this model's shocks and variables."
                    else if (n_fail > 0)
                      sprintf("%d benchmark(s) failed (wrong direction or non-finite IRF). Check shock sign conventions and Taylor rule parameters.",
                              n_fail)
                    else if (n_flag > floor(0.25 * n_checked))
                      sprintf("%d/%d benchmarks flagged for wrong magnitude/timing (threshold: <=25%% of checked). Review calibration.",
                              n_flag, n_checked)
                    else if (n_flag > 0)
                      sprintf("All directions OK; %d/%d benchmark(s) have magnitude/timing outside range (within 25%% tolerance). Review calibration.",
                              n_flag, n_checked)
                    else
                      sprintf("All %d IRF benchmarks satisfied.", n_pass))
          ), collapse = "\n")
        }
      ),
      class = "dynhr_diagnostic"
    )
}
