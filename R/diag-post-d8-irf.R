## R/diag-post-d8-irf.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D8 IRF plausibility checks; .get_irfs_long() helper
## --------------------------------------------------------------------------

.get_irfs_long <- function(irfs) {
  if (!is.null(irfs$irfs) && is.list(irfs$irfs)) irfs <- irfs$irfs
  if (!is.list(irfs) || length(irfs) == 0)
    stop("irfs must be a named list of matrices (IRFCollection)")
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


#' D8. IRF plausibility checks
#'
#' Compares model impulse response functions against stylised economic
#' benchmarks (sign, peak magnitude, and peak timing).  Each benchmark
#' specifies the expected direction, peak range, and peak-horizon range for
#' a (shock, variable) pair.  The diagnostic SKIPs benchmarks whose shock or
#' variable is absent from the IRF collection and counts PASS, FLAG, and FAIL
#' results for the remaining checks.
#'
#' Pass gate: \code{pass = (n_fail == 0) && (n_flag <= floor(0.25 * n_checked))}.
#' A FAIL is a wrong sign; a FLAG is a correct sign with peak magnitude or
#' timing outside the benchmark range.  Up to 25\% of checked benchmarks may
#' be flagged without triggering a FAIL at the gate level -- one benchmark in
#' four is a lenient but non-zero bar.
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
#'       peak-horizon range (optional; NULL to skip timing check).}
#'     \item{\code{description}}{Human-readable description string.}
#'   }
#'   When NULL, a set of standard macroeconomic benchmarks is used (monetary
#'   policy, technology, and UIP shocks for output, inflation, and the real
#'   exchange rate).
#' @param metadata  Optional metadata list for plot annotation.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{Named list of per-benchmark check results, each containing
#'     \code{status} (PASS/FLAG/FAIL/SKIP), \code{peak_value},
#'     \code{peak_horizon}, and \code{direction_ok}, \code{magnitude_ok},
#'     \code{timing_ok} flags.}
#'   \item{pass}{Logical -- TRUE if no FAILs and FLAG count <= 25\% of checked.}
#'   \item{plots}{ggplot2 panel of matched IRFs with peak markers.}
#'   \item{summary}{Human-readable per-benchmark summary lines.}
#'
#' @references Christiano, L. J., Eichenbaum, M., & Evans, C. L. (1999). Monetary
#'   policy shocks: What have we learned and to what end? In J. B. Taylor & M.
#'   Woodford (eds), \emph{Handbook of Macroeconomics}, Vol. 1A. North-Holland.
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

    if (is.null(benchmarks)) {
      benchmarks <- list(
        monetary = list(
          shock = c("eps_r","eps_mp","eps_r_","er"),
          variable = "y", direction = "negative",
          peak_range = c(-0.005, -0.002), peak_horizon = c(4, 8),
          description = "25bp contraction -> y falls 0.2-0.5% in 4-8q"
        ),
        monetary_pi = list(
          shock = c("eps_r","eps_mp","eps_r_","er"),
          variable = "pi", direction = "negative",
          peak_range = c(-0.003, -0.001), peak_horizon = c(6, 12),
          description = "25bp contraction -> pi falls in 6-12q"
        ),
        technology = list(
          shock = c("eps_a","eps_a_","eps_tech"),
          variable = "y", direction = "positive",
          peak_range = c(0.001, 0.02), peak_horizon = c(1, 6),
          description = "Positive TFP shock -> y rises"
        ),
        uip = list(
          shock = c("eps_uip_","eps_uip","eps_rp","es"),
          variable = "q", direction = "positive",
          peak_range = c(0.001, 0.05), peak_horizon = c(1, 4),
          description = "UIP shock -> RER depreciates"
        )
      )
    }

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
      peak_idx <- if (b$direction == "negative") which.min(sub$value) else which.max(sub$value)
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
    pass      <- (n_fail == 0L && n_flag <= floor(0.25 * n_checked))

    summary_lines <- vapply(check_results, function(x)
      sprintf("  [%s] %s: %s", x$status, x$benchmark, x$message), character(1))

    # --- Plot: matched benchmark IRFs with the peak marked, coloured by status ---
    plots <- list()
    matched <- Filter(function(r) !is.null(r$shock) && !identical(r$status, "SKIP"),
                      check_results)
    if (length(matched) > 0 && requireNamespace("ggplot2", quietly = TRUE)) {
      rows <- do.call(rbind, lapply(matched, function(r) {
        sub <- irf_long[irf_long$shock == r$shock & irf_long$variable == r$variable, ]
        if (nrow(sub) == 0) return(NULL)
        data.frame(panel = sprintf("%s: %s -> %s", r$benchmark, r$shock, r$variable),
                   horizon = sub$horizon, value = sub$value, status = r$status,
                   stringsAsFactors = FALSE)
      }))
      if (!is.null(rows) && nrow(rows) > 0) {
        pk <- do.call(rbind, lapply(matched, function(r) {
          if (is.null(r$peak_horizon)) return(NULL)
          data.frame(panel = sprintf("%s: %s -> %s", r$benchmark, r$shock, r$variable),
                     horizon = r$peak_horizon, value = r$peak_value, status = r$status,
                     stringsAsFactors = FALSE)
        }))
        status_cols <- c(PASS = dynhr_colours$green %||% "#117733",
                         FLAG = dynhr_colours$orange, FAIL = dynhr_colours$red)
        p_irf <- ggplot2::ggplot(rows, ggplot2::aes(x = horizon, y = value)) +
          geom_dynhr_zero() +
          ggplot2::geom_line(ggplot2::aes(colour = status), linewidth = 0.7) +
          { if (!is.null(pk)) ggplot2::geom_point(
              data = pk, ggplot2::aes(x = horizon, y = value, colour = status),
              size = 2.4) } +
          { if (!is.null(pk)) ggplot2::geom_text(
              data = pk,
              ggplot2::aes(x = horizon, y = value, colour = status,
                           label = sprintf("peak %.3g @h%d", value, horizon)),
              vjust = -0.8, hjust = 0.0, size = 2.7, show.legend = FALSE) } +
          ggplot2::facet_wrap(~ panel, scales = "free_y") +
          ggplot2::coord_cartesian(clip = "off") +
          ggplot2::scale_colour_manual(values = status_cols, name = NULL) +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.18)) +
          theme_dynhr_compact() +
          ggplot2::labs(
            title = "D8: Benchmark IRFs (peak marked)",
            subtitle = "Matched textbook priors -- sign / magnitude / timing",
            x = "Horizon", y = NULL)
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
          badge <- if (pass) "PASS" else "FAIL"
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
                    if (n_fail > 0)
                      sprintf("%d benchmark(s) failed (wrong direction). Check shock sign conventions and Taylor rule parameters.",
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
