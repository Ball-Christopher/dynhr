## ======================================================================
## dynhr_transform.R (v3) -- Gap construction for NZSIM estimation
## ======================================================================
## Clean single-definition rewrite of dynhr_transform.
## Key changes from v2:
##   - Single definition (no duplicates)
##   - Potential model integration (Kalman path for y_trend)
##   - Removed ghost column deletions
##   - Gap formulas verified against batch_2_old.R
##   - Government share trend as residual (matching batch_2)
##
## KNOWN FIX (v3.1): attr<- in R creates a shallow copy of data.tables,
## breaking by-reference modification in downstream functions. All attr
## assignments on dt are replaced with setattr() or deferred to gaps.
##
## Dependencies: data.table
##   For Kalman path: dynhr_backend.R, dynhr_potential.R,
##                    dynhr_potential_data.R
## ======================================================================

## Phase-0 packaging: top-level library() call removed. data.table is
## listed in Suggests; each function that uses data.table calls
## `requireNamespace("data.table", quietly = TRUE)` and emits a clear
## error if missing. See .ensure_dt() helper below.

## data.table's `[` checks the calling package for this flag (cedta());
## without it, every dt[...] call from this namespace is evaluated with
## data.frame semantics and NSE like dt[order(date)] fails.
.datatable.aware <- TRUE

## Helper: ensure data.table is available, return the namespace.
.ensure_dt <- function() {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required for gap construction. ",
         "Install it with install.packages('data.table').")
  }
  getNamespace("data.table")
}


## ====================================================================
## HP Filter (self-contained, pentadiagonal dense solve)
## ====================================================================

.hp_trend <- function(x, lambda = 1600) {
  if (all(is.na(x))) return(x)
  n_full <- length(x)
  
  ## Locate valid range
  valid <- !is.na(x)
  first <- which(valid)[1]
  last  <- tail(which(valid), 1)
  y <- x[first:last]
  
  ## Interpolate interior NAs
  if (any(is.na(y))) {
    idx <- seq_along(y)
    y <- approx(idx[!is.na(y)], y[!is.na(y)], idx, rule = 2)$y
  }
  
  n <- length(y)
  if (n < 4) {
    out <- rep(NA_real_, n_full)
    out[first:last] <- y
    return(out)
  }
  
  ## Pentadiagonal system: (I + lambda * D2'D2) * tau = y
  d0 <- rep(1 + 6 * lambda, n)
  d0[1] <- 1 + lambda; d0[2] <- 1 + 5 * lambda
  d0[n] <- 1 + lambda; d0[n - 1] <- 1 + 5 * lambda
  
  d1 <- rep(-4 * lambda, n - 1)
  d1[1] <- -2 * lambda; d1[n - 1] <- -2 * lambda
  
  d2 <- rep(lambda, n - 2)
  
  A <- diag(d0)
  for (i in seq_len(n - 1)) { A[i, i + 1] <- d1[i]; A[i + 1, i] <- d1[i] }
  for (i in seq_len(n - 2)) { A[i, i + 2] <- d2[i]; A[i + 2, i] <- d2[i] }
  
  tau <- as.numeric(solve(A, y))
  
  out <- rep(NA_real_, n_full)
  out[first:last] <- tau
  out
}


## Log-domain HP: filter log(x), return exp(trend)
.hp_trend_log <- function(x, lambda = 1600) {
  lx <- log(x)
  lx[!is.finite(lx)] <- NA
  exp(.hp_trend(lx, lambda))
}


## ====================================================================
## Anchored HP filter
## ====================================================================

.hp_anchored <- function(anchors, lambda = 1600) {
  n <- length(anchors)
  if (all(is.na(anchors))) return(rep(NA_real_, n))
  
  idx  <- which(!is.na(anchors))
  vals <- anchors[idx]
  
  ## Linear interpolation between anchor points, constant extrapolation
  interp <- approx(idx, vals, xout = 1:n, rule = 2)$y
  
  ## HP-smooth the interpolated series (removes kinks at anchor points)
  .hp_trend(interp, lambda)
}


## ====================================================================
## Neutral rate builders (from batch_2_old.R)
## ====================================================================

## -- Neutral OCR --
.build_neutral_r <- function(dt) {
  dt[, r_anchor := NA_real_]
  dt[date >= data.table::as.IDate("1982-01-01") & date <= data.table::as.IDate("2008-01-01"),
     r_anchor := 5.75 / 400]
  dt[date >= data.table::as.IDate("2010-01-01") & date <= data.table::as.IDate("2010-04-01"),
     r_anchor := 4.25 / 400]
  dt[date == data.table::as.IDate("2016-04-01"), r_anchor := 3.50 / 400]
  dt[date == data.table::as.IDate("2017-10-01"), r_anchor := 3.25 / 400]
  dt[date >= data.table::as.IDate("2020-07-01") & date <= data.table::as.IDate("2021-07-01"),
     r_anchor := 2.00 / 400]
  dt[date >= data.table::as.IDate("2024-01-01"), r_anchor := 3.00 / 400]
  
  dt[, r_trend := .hp_anchored(r_anchor, lambda = 1600)]
  invisible(dt)
}


## -- Neutral mortgage rate --
.build_neutral_rh <- function(dt) {
  preGFCspread  <- 1.197 / 400
  postGFCspread <- 2.117 / 400
  
  dt[, rh_anchor := NA_real_]
  dt[date >= data.table::as.IDate("1982-01-01") & date <= data.table::as.IDate("2008-01-01"),
     rh_anchor := r_trend + preGFCspread]
  dt[date >= data.table::as.IDate("2010-01-01"),
     rh_anchor := r_trend + postGFCspread]
  
  dt[, rh_trend := .hp_anchored(rh_anchor, lambda = 1600)]
  invisible(dt)
}


## -- Neutral world interest rate --
.build_neutral_rstar <- function(dt) {
  dt[, rstar_anchor := NA_real_]
  dt[date >= data.table::as.IDate("1982-01-01") & date <= data.table::as.IDate("2007-07-01"),
     rstar_anchor := 3.75 / 400]
  dt[date >= data.table::as.IDate("2010-01-01"),
     rstar_anchor := 3.25 / 400]
  
  dt[, rstar_trend := .hp_anchored(rstar_anchor, lambda = 1600)]
  invisible(dt)
}


## ====================================================================
## Main transform function
## ====================================================================

### ======================================================================
### dynhr_transform v3 -- Gap construction with Kalman potential model
### (SUPERSEDED by v4 below; renamed to avoid shadowing v4)
### ======================================================================
### Changes from v2:
###   - potential_params: NOT IMPLEMENTED. The Kalman potential-output
###     path was never built (its builder functions do not exist); supplying
###     potential_params errors immediately with a fail-loud stop(). Use
###     potential_y_trend_override= or the HP-filter fallback instead.
###   - potential_obs_override: pre-built (extended) obs table for filter
###   - potential_y_trend_override: directly inject y_trend (bypass filter)
###   - u_k=0 override for Kalman path (capital stock is deterministic)
###   - HP fallback for dates outside Kalman coverage
### ======================================================================

## NOTE: potential_params is NOT IMPLEMENTED. Passing a non-NULL value
## selects the Kalman potential-output path (PATH B below), which errors
## immediately via stop() because its three builder functions (the
## potential-obs constructor, its endpoint extender, and the state-space
## spec builder) do not exist anywhere in the package. Use
## potential_y_trend_override= or the HP-filter fallback instead.
.dynhr_transform_v3 <- function(
    est,
    nzsim_ref_path   = NULL,
    est_start        = "1993-01-01",
    est_end          = NULL,
    potential_params  = NULL,
    potential_alpha   = 2/3,
    potential_delta   = 0.025,
    potential_start   = "1982-01-01",
    potential_n_extend = 60L,
    potential_y_trend_override = NULL
) {

  .ensure_dt()
  full <- data.table::copy(attr(est, "full_data"))
  dt <- full[order(date)]
  
  message("=== dynhr_transform v3: Building gap dataset ===")
  
  ## -- Calibrated NZSIM share parameters ---------------------------
  S <- list(cy = 0.577793, xy = 0.286209, my = 0.279613,
            ihi = 0.060565, iki = 0.125050, gy = 0.229997)
  
  ## -- HP lambdas --------------------------------------------------
  L <- list(
    y = 200000, r = 1600, rh = 1600, rstar = 1600,
    rs = 80000, p = 50000, pn_p = 25000, ph_p = 5000,
    wrlci = 25000, pstar = 32000, pm_ps = 80000, px_ps = 80000,
    ln = 50000, b_ngdp = 56000,
    c_gdp = 20000, ik_gdp = 80000, ih_gdp = 100000,
    x_gdp = 40000, m_gdp = 80000, g_gdp = 20000
  )
  
  ## ----------------------------------------------------------------
  ## STEP 1: Construct intermediate variables
  ## ----------------------------------------------------------------
  
  message("-- Step 1: Intermediate variables --")
  
  dt[is.na(r_constructed) & !is.na(r90d), r_constructed := r90d]
  
  dt[, r_q := r_constructed / 400]
  dt[, rh_q := rh_constructed / 400]
  dt[, rstar_q := rshortw_ocr / 400]
  
  dt[, rs_real := rtwi * pcpis / wcpi]
  
  dt[, pn_rel := pnt / pcpis]
  dt[, ph_rel := pqhpiz / pcpis]
  dt[, w_real := llisai / pcpis]
  
  dt[, b_ratio := -tiin / ngdpz / 4]
  
  dt[, c_sh := ncp_z / ngdpp_z]
  dt[, ik_sh := nik_z / ngdpp_z]
  dt[, ih_sh := nitd_z / ngdpp_z]
  dt[, x_sh := nx_z / ngdpp_z]
  dt[, m_sh := nm_z / ngdpp_z]
  dt[, g_sh := ncg_z / ngdpp_z]
  
  dt[, pm_ps := pmstar_constructed / wcpi]
  dt[, px_ps := pxstar_constructed / wcpi]
  
  dt[, ln_norm := lmig_z / (lhpwa_z * 1000)]
  
  n_vars <- sum(!is.na(dt$r_q) & !is.na(dt$rh_q) & !is.na(dt$c_sh))
  message(sprintf("  %d quarters with core variables available", n_vars))
  
  ## ----------------------------------------------------------------
  ## STEP 2: Trends
  ## ----------------------------------------------------------------
  
  message("-- Step 2: Trends --")
  
  ## ==============================================================
  ## y_trend: three paths (in priority order)
  ## ==============================================================
  
  if (!is.null(potential_y_trend_override)) {
    ## -- PATH A: Direct y_trend injection ----------------------
    ## Pre-computed y_trend (e.g. from standalone extended Kalman run)
    yt <- data.table::copy(potential_y_trend_override)
    yt[, date := data.table::as.IDate(date)]
    dt[yt, y_trend := i.y_trend, on = .(date)]
    n_filled <- dt[!is.na(y_trend), .N]
    n_total  <- dt[!is.na(ngdpp_z), .N]
    message(sprintf("  y_trend: using override (%d/%d quarters)", n_filled, n_total))
    
    ## Fill gaps with HP fallback
    if (n_filled < n_total) {
      dt[is.na(y_trend), y_trend := .hp_trend_log(ngdpp_z, L$y)]
      n_hp <- n_total - n_filled
      message(sprintf("  Filling %d y_trend gaps with HP fallback", n_hp))
    }
    
  } else if (!is.null(potential_params)) {
    ## -- PATH B: Kalman filter potential model -----------------
    ## NOT IMPLEMENTED: this path used to call three builder functions (the
    ## potential-obs constructor, its endpoint extender, and the state-space
    ## spec builder), none of which are defined anywhere in the package.
    ## Fail loud instead of letting a "could not find function" error
    ## surface deep in the call stack. See
    ## .claude/orchestration/track-p-hardening/brief-A2-pathb-failloud.md.
    stop("dynhr_transform: the Kalman potential-output path ",
         "(potential_params=) was never implemented (its builder ",
         "functions do not exist). Use potential_y_trend_override= or ",
         "the HP-filter fallback instead.", call. = FALSE)
  } else {
    ## -- PATH C: HP filter (default) --------------------------
    message("  y_trend: using HP filter (lambda=200000)")
    dt[, y_trend := .hp_trend_log(ngdpp_z, L$y)]
  }
  
  ## -- Interest rates: judged neutrals -------------------------
  message("  Building judged neutral rates (batch_2c anchors)")
  .build_neutral_r(dt)
  .build_neutral_rh(dt)
  .build_neutral_rstar(dt)
  
  ## -- Other trends (same for all paths) -----------------------
  dt[, rs_trend := .hp_trend_log(rs_real, L$rs)]
  dt[, p_trend := .hp_trend_log(pcpis, L$p)]
  dt[, pn_rel_trend := .hp_trend_log(pn_rel, L$pn_p)]
  dt[, ph_rel_trend := .hp_trend_log(ph_rel, L$ph_p)]
  dt[, w_real_trend := .hp_trend_log(w_real, L$wrlci)]
  dt[, pstar_trend := .hp_trend_log(wcpi, L$pstar)]
  
  dt[, pm_ps_trend := .hp_trend_log(pm_ps, L$pm_ps)]
  dt[, pmstar_trend := pm_ps_trend * pstar_trend]
  dt[, px_ps_trend := .hp_trend_log(px_ps, L$px_ps)]
  dt[, pxstar_trend := px_ps_trend * pstar_trend]
  
  if ("iwgdp_pt" %in% names(dt)) {
    dt[, ystar_trend := iwgdp_pt]
    dt[is.na(ystar_trend), ystar_trend := .hp_trend_log(iwgdp_z, L$y)]
    message("  ystar: using iwgdp_pt as trend (batch_2 equivalent)")
  } else {
    dt[, ystar_trend := .hp_trend_log(iwgdp_z, L$y)]
  }
  
  dt[, ln_trend := .hp_trend(ln_norm, L$ln)]
  dt[, b_ratio_trend := .hp_trend_log(b_ratio, L$b_ngdp)]
  
  share_map <- list(
    c = list(col = "c_sh", lambda = L$c_gdp),
    ik = list(col = "ik_sh", lambda = L$ik_gdp),
    ih = list(col = "ih_sh", lambda = L$ih_gdp),
    x = list(col = "x_sh", lambda = L$x_gdp),
    m = list(col = "m_sh", lambda = L$m_gdp),
    g = list(col = "g_sh", lambda = L$g_gdp)
  )
  for (nm in names(share_map)) {
    scol <- share_map[[nm]]$col
    tcol <- paste0(nm, "_sh_trend")
    dt[, (tcol) := .hp_trend_log(get(scol), share_map[[nm]]$lambda)]
  }
  
  dt[, c_trend := c_sh_trend * y_trend]
  dt[, ik_trend := ik_sh_trend * y_trend]
  dt[, ih_trend := ih_sh_trend * y_trend]
  dt[, x_trend := x_sh_trend * y_trend]
  dt[, m_trend := m_sh_trend * y_trend]
  dt[, g_trend := g_sh_trend * y_trend]
  
  message("  Trends computed for all 20 varobs")
  
  ## ----------------------------------------------------------------
  ## STEP 3: Compute gaps
  ## ----------------------------------------------------------------
  
  message("-- Step 3: Gaps --")
  
  dt[, r_ := r_q - r_trend]
  dt[, rh_ := rh_q - rh_trend]
  dt[, rstar_ := rstar_q - rstar_trend]
  dt[, rs_ := log(rs_real / rs_trend)]
  
  dt[, c_ := (ncp_z - c_trend) / y_trend / S$cy]
  dt[, ik_ := (nik_z - ik_trend) / y_trend / S$iki]
  dt[, ih_ := (nitd_z - ih_trend) / y_trend / S$ihi]
  dt[, x_ := (nx_z - x_trend) / y_trend / S$xy]
  dt[, m_ := (nm_z - m_trend) / y_trend / S$my]
  dt[, g_ := (ncg_z - g_trend) / y_trend / S$gy]
  
  dt[, y_ := g_ * S$gy + S$cy * c_ + S$iki * ik_ +
       S$ihi * ih_ + S$xy * x_ - S$my * m_]
  
  dt[, dp_actual := pcpis / data.table::shift(pcpis) - 1]
  dt[, dp_trend := p_trend / data.table::shift(p_trend) - 1]
  dt[, dp_ := dp_actual - dp_trend]
  
  dt[, pn_p_ := log(pn_rel / pn_rel_trend)]
  dt[, ph_p_ := log(ph_rel / ph_rel_trend)]
  dt[, wrlci_ := log(w_real / w_real_trend)]
  
  dt[, pstar_ := log(wcpi / pstar_trend)]
  dt[, pmstar_ := log(pmstar_constructed / pmstar_trend)]
  dt[, pxstar_ := log(pxstar_constructed / pxstar_trend)]
  dt[, ystar_ := log(iwgdp_z / ystar_trend)]
  
  dt[, ln_ := ln_norm - ln_trend]
  dt[, b_ := b_ratio - b_ratio_trend]
  
  ## ----------------------------------------------------------------
  ## STEP 4: Assemble estimation matrix
  ## ----------------------------------------------------------------
  
  message("-- Step 4: Assemble --")
  
  gap_cols <- c("r_", "dp_", "rs_", "y_", "c_", "x_", "m_",
                "ik_", "ih_", "ln_", "b_", "ph_p_", "pn_p_",
                "rh_", "pmstar_", "pxstar_", "ystar_", "pstar_",
                "rstar_", "wrlci_")
  
  gaps <- dt[, c("date", gap_cols), with = FALSE]
  
  if (!is.null(est_start)) gaps <- gaps[date >= data.table::as.IDate(est_start)]
  
  if (!is.null(est_end)) {
    gaps <- gaps[date <= data.table::as.IDate(est_end)]
  } else {
    complete <- complete.cases(gaps[, -"date"])
    if (any(complete)) gaps <- gaps[1:max(which(complete))]
  }
  
  message(sprintf("  Final: %d quarters (%s to %s)",
                  nrow(gaps), min(gaps$date), max(gaps$date)))
  
  na_check <- sapply(gaps[, -"date"], function(x) sum(is.na(x)))
  if (any(na_check > 0)) {
    message("  WARNING: Missing data:")
    print(na_check[na_check > 0])
  } else {
    message("  [OK] No missing data in estimation window")
  }
  
  ## ----------------------------------------------------------------
  ## STEP 5: Compare to nzsim_data.csv reference
  ## ----------------------------------------------------------------
  
  if (!is.null(nzsim_ref_path) && file.exists(nzsim_ref_path)) {
    message("-- Step 5: Comparison to nzsim_data.csv --")
    comp <- .compare_gaps(gaps, nzsim_ref_path)
    attr(gaps, "comparison") <- comp
  }
  
  attr(gaps, "trends") <- dt
  attr(gaps, "lambdas") <- L
  attr(gaps, "shares") <- S
  message("=== Done ===")
  return(gaps)
}

## ====================================================================
## Comparison to nzsim_data.csv
## ====================================================================

.compare_gaps <- function(gaps, ref_path) {
  .ensure_dt()
  ref    <- data.table::fread(ref_path)
  n_ref  <- nrow(ref)
  message(sprintf("  Reference: %d rows, %d columns", n_ref, ncol(ref)))
  
  gap_cols <- intersect(names(gaps), names(ref))
  message(sprintf("  Matching columns: %d/%d", length(gap_cols), ncol(ref)))
  
  ## ---- Date alignment via r_ correlation scan ----
  candidate_starts <- gaps[, unique(date)]
  best_corr <- -1; best_start <- candidate_starts[1]
  
  for (sd in candidate_starts) {
    idx <- which(gaps$date >= sd)
    if (length(idx) < n_ref) next
    our <- gaps[idx[1:n_ref], r_]
    if (any(is.na(our)) || any(is.na(ref$r_))) next
    cc <- cor(our, ref$r_, use = "complete.obs")
    if (!is.na(cc) && cc > best_corr) {
      best_corr <- cc; best_start <- sd
    }
  }
  message(sprintf("  Best alignment: start = %s (r_ corr = %.4f)",
                  best_start, best_corr))
  
  ## ---- Extract aligned window ----
  idx    <- which(gaps$date >= best_start)
  n_comp <- min(length(idx), n_ref)
  our_aligned <- gaps[idx[1:n_comp]]
  
  ## ---- Per-variable comparison ----
  results <- list()
  for (col in gap_cols) {
    ours <- our_aligned[[col]]
    refs <- ref[[col]][1:n_comp]
    valid <- !is.na(ours) & !is.na(refs)
    if (sum(valid) < 10) {
      results[[col]] <- data.table::data.table(
        varobs = col, n = sum(valid), corr = NA_real_,
        sd_ours = NA_real_, sd_ref = NA_real_,
        sd_ratio = NA_real_, max_abs = NA_real_, mean_abs = NA_real_
      )
      next
    }
    o <- ours[valid]; r <- refs[valid]
    results[[col]] <- data.table::data.table(
      varobs   = col,
      n        = sum(valid),
      corr     = cor(o, r),
      sd_ours  = sd(o),
      sd_ref   = sd(r),
      sd_ratio = sd(o) / sd(r),
      max_abs  = max(abs(o - r)),
      mean_abs = mean(abs(o - r))
    )
  }
  
  summary_dt <- data.table::rbindlist(results)
  
  ## ---- Print summary ----
  message("\n  == Gap comparison summary ==")
  message(sprintf("  Aligned: %s to %s (%d quarters)\n",
                  our_aligned$date[1], our_aligned$date[n_comp], n_comp))
  
  print(summary_dt[order(corr)])
  
  good    <- summary_dt[!is.na(corr) & corr > 0.95]
  ok      <- summary_dt[!is.na(corr) & corr > 0.80 & corr <= 0.95]
  bad     <- summary_dt[!is.na(corr) & corr <= 0.80]
  missing <- summary_dt[is.na(corr)]
  
  if (nrow(good) > 0) message(sprintf("\n  GOOD (>0.95): %s",
                                      paste(good$varobs, collapse = ", ")))
  if (nrow(ok) > 0) message(sprintf("  OK (0.80-0.95): %s",
                                    paste(ok$varobs, collapse = ", ")))
  if (nrow(bad) > 0) message(sprintf("  POOR (<0.80): %s",
                                     paste(bad$varobs, collapse = ", ")))
  if (nrow(missing) > 0) message(sprintf("  MISSING: %s",
                                         paste(missing$varobs, collapse = ", ")))
  
  ## Worst-3 dates for flagged series
  for (col in c(bad$varobs, ok$varobs)) {
    ours <- our_aligned[[col]]
    refs <- ref[[col]][1:n_comp]
    valid <- !is.na(ours) & !is.na(refs)
    diffs <- abs(ours[valid] - refs[valid])
    worst_idx <- head(order(-diffs), 3)
    worst_dates <- our_aligned$date[which(valid)[worst_idx]]
    message(sprintf("  %s: worst at %s (diff=%.4f, %.4f, %.4f)",
                    col, paste(worst_dates, collapse = ", "),
                    diffs[worst_idx[1]],
                    ifelse(length(worst_idx) > 1, diffs[worst_idx[2]], NA),
                    ifelse(length(worst_idx) > 2, diffs[worst_idx[3]], NA)))
  }
  
  return(list(summary = summary_dt, aligned_start = best_start,
              n_compared = n_comp, our = our_aligned, ref = ref))
}


## ====================================================================
## Diagnostic plot
## ====================================================================

dynhr_plot_comparison <- function(gaps, ncol = 4) {
  comp <- attr(gaps, "comparison")
  if (is.null(comp)) stop("No comparison data. Run dynhr_transform with nzsim_ref_path.")
  
  our <- comp$our
  ref <- comp$ref
  n   <- comp$n_compared
  cols <- intersect(names(our), names(ref))
  cols <- setdiff(cols, "date")
  
  nrow_plot <- ceiling(length(cols) / ncol)
  par(mfrow = c(nrow_plot, ncol), mar = c(2, 3, 2, 1), cex = 0.7)
  
  for (col in cols) {
    o <- our[[col]][1:n]
    r <- ref[[col]][1:n]
    cc <- cor(o, r, use = "complete.obs")
    ylim <- range(c(o, r), na.rm = TRUE)
    
    plot(seq_len(n), r, type = "l", col = "black", lwd = 1.5,
         ylim = ylim, main = sprintf("%s (r=%.3f)", col, cc),
         xlab = "", ylab = "")
    lines(seq_len(n), o, col = "steelblue", lwd = 1.5, lty = 2)
    abline(h = 0, col = "grey60", lty = 3)
  }
}

## ======================================================================
## COVID lockdown smoothing -- new helper + dynhr_transform v3 patch
## ======================================================================
## Replicates batch_2's ypsmoothed approach:
##   ypsmoothed = y; ypsmoothed[lockdown_periods] = NA
##   ypsmoothed = interp_iris(ypsmoothed)
##   lockdown = y / ypsmoothed
##   temp_t_ngdp[>=2020Q1] = ngdp / lockdown
##
## Effects:
##   1. GDP shares use smoothed GDP denominator -> smoother through lockdowns
##   2. Debt/GDP ratio uses adjusted nominal GDP -> no COVID spike
##   3. HP fallback y_trend filters smoothed GDP -> no lockdown pull-down
## ======================================================================

## -- Helper: lockdown interpolation -----------------------------------
## Matches batch_2: set lockdown periods to NA, cubic-spline interpolate,
## return smoothed series and lockdown ratio.
##
## x:       numeric vector (levels, e.g. ngdpp_z)
## dates:   IDate vector aligned with x
## lockdown_periods: list of c(start_date, end_date) character pairs
## method:  "spline" (natural cubic, matches interp_iris) or "linear"

.lockdown_smooth <- function(x, dates,
                             lockdown_periods = list(
                               c("2020-01-01", "2020-04-01"),   # 2020Q1-Q2
                               c("2021-07-01", "2022-01-01")    # 2021Q3-2022Q1
                             ),
                             method = c("spline", "linear")) {
  method <- match.arg(method)
  xs <- x
  
  ## 1. Set lockdown periods to NA
  n_removed <- 0L
  for (lp in lockdown_periods) {
    mask <- dates >= data.table::as.IDate(lp[1]) & dates <= data.table::as.IDate(lp[2]) & !is.na(x)
    n_removed <- n_removed + sum(mask)
    xs[mask] <- NA
  }
  
  ## 2. Interpolate through NAs (on levels, matching interp_iris default)
  valid <- which(!is.na(xs) & is.finite(xs))
  if (length(valid) < 4) {
    warning("lockdown_smooth: too few valid obs after exclusion")
    return(list(smoothed = x, ratio = rep(1, length(x)),
                n_excluded = 0L))
  }
  
  all_idx <- seq_along(x)
  if (method == "spline") {
    xs_interp <- spline(valid, xs[valid], xout = all_idx,
                        method = "natural")$y
  } else {
    xs_interp <- approx(valid, xs[valid], xout = all_idx, rule = 2)$y
  }
  
  ## 3. Lockdown ratio = actual / smoothed
  ##    < 1 during lockdowns (GDP dropped below counterfactual)
  ##    > 1 during recovery bounce
  ##    -> 1 as COVID effects fade
  ratio <- x / xs_interp
  
  ## Only apply from first lockdown start onwards
  first_ld <- data.table::as.IDate(lockdown_periods[[1]][1])
  pre_ld <- which(dates < first_ld)
  ratio[pre_ld] <- 1
  xs_interp[pre_ld] <- x[pre_ld]
  
  list(smoothed = xs_interp, ratio = ratio, n_excluded = n_removed)
}


## ======================================================================
## dynhr_transform v3 -- full function with lockdown smoothing
## ======================================================================
## Changes from v2:
##   - New param: lockdown_periods (default = batch_2 periods)
##   - Step 1b: compute lockdown ratio from real GDP
##   - Shares computed with smoothed GDP denominator
##   - Debt/GDP computed with adjusted nominal GDP
##   - y_trend HP fallback uses smoothed GDP
##   - Removed stale cleanup of non-existent columns (L673 bug)
##   - Three-path y_trend (override / Kalman / HP) preserved
##
## NOTE: potential_params is NOT IMPLEMENTED. Passing a non-NULL value
## selects the Kalman potential-output path (Path B below), which errors
## immediately via stop() because its three builder functions (the
## potential-obs constructor, its endpoint extender, and the state-space
## spec builder) do not exist anywhere in the package. Use
## potential_y_trend_override= or the HP-filter fallback instead.

dynhr_transform <- function(est,
                            nzsim_ref_path = NULL,
                            est_start = "1993-01-01",
                            est_end = NULL,
                            potential_params = NULL,
                            potential_n_extend = 60L,
                            potential_y_trend_override = NULL,
                            lockdown_periods = list(
                              c("2020-01-01", "2020-04-01"),
                              c("2021-07-01", "2022-01-01")
                            )) {

  .ensure_dt()
  full <- data.table::copy(attr(est, "full_data"))
  dt <- full[order(date)]
  
  message("=== dynhr_transform v3: Building gap dataset ===")
  
  ## -- Calibrated NZSIM share parameters ---------------------------
  S <- list(cy = 0.577793, xy = 0.286209, my = 0.279613,
            ihi = 0.060565, iki = 0.125050, gy = 0.229997)
  
  ## -- HP lambdas -------------------------------------------------
  L <- list(
    y = 200000,
    r = 1600, rh = 1600, rstar = 1600,
    rs = 80000,
    p = 50000,
    pn_p = 25000, ph_p = 5000,
    wrlci = 25000,
    pstar = 32000,
    pm_ps = 80000, px_ps = 80000,
    ln = 50000,
    b_ngdp = 56000,
    c_gdp = 20000, ik_gdp = 80000, ih_gdp = 100000,
    x_gdp = 40000, m_gdp = 80000, g_gdp = 20000
  )
  
  ## ----------------------------------------------------------------
  ## STEP 1: Construct intermediate variables
  ## ----------------------------------------------------------------
  message("-- Step 1: Intermediate variables --")
  
  dt[is.na(r_constructed) & !is.na(r90d), r_constructed := r90d]
  
  ## Quarterly interest rates
  dt[, r_q := r_constructed / 400]
  dt[, rh_q := rh_constructed / 400]
  dt[, rstar_q := rshortw_ocr / 400]
  
  ## Real exchange rate
  dt[, rs_real := rtwi * pcpis / wcpi]
  
  ## Relative prices
  dt[, pn_rel := pnt / pcpis]
  dt[, ph_rel := pqhpiz / pcpis]
  dt[, w_real := llisai / pcpis]
  
  ## World price ratios
  dt[, pm_ps := pmstar_constructed / wcpi]
  dt[, px_ps := pxstar_constructed / wcpi]
  
  ## -- Step 1b: COVID lockdown smoothing --------------------------
  ## batch_2: ypsmoothed = interp(y with lockdowns set to NA)
  ##          lockdown   = y / ypsmoothed
  ##          temp_t_ngdp[>=2020Q1] = ngdp / lockdown
  if (!is.null(lockdown_periods) && length(lockdown_periods) > 0) {
    message("-- Step 1b: COVID lockdown smoothing --")
    
    ## Smooth real GDP (production side)
    ld <- .lockdown_smooth(dt$ngdpp_z, dt$date, lockdown_periods)
    dt[, ngdpp_z_smooth := ld$smoothed]
    dt[, lockdown_ratio  := ld$ratio]
    
    ## Adjust nominal GDP: batch_2 temp_t_ngdp = ngdp / lockdown
    ## lockdown = y / ypsmoothed, so ngdp / lockdown = ngdp * ypsmoothed / y
    dt[, ngdpz_smooth := ngdpz / lockdown_ratio]
    
    message(sprintf("  %d lockdown obs excluded, ratio range [%.4f, %.4f]",
                    ld$n_excluded,
                    min(dt$lockdown_ratio, na.rm = TRUE),
                    max(dt$lockdown_ratio, na.rm = TRUE)))
  } else {
    dt[, ngdpp_z_smooth := ngdpp_z]
    dt[, ngdpz_smooth   := ngdpz]
    dt[, lockdown_ratio  := 1]
    message("  Lockdown smoothing: disabled")
  }
  
  ## GDP expenditure shares -- smoothed GDP denominator
  dt[, c_sh  := ncp_z  / ngdpp_z_smooth]
  dt[, ik_sh := nik_z  / ngdpp_z_smooth]
  dt[, ih_sh := nitd_z / ngdpp_z_smooth]
  dt[, x_sh  := nx_z   / ngdpp_z_smooth]
  dt[, m_sh  := nm_z   / ngdpp_z_smooth]
  dt[, g_sh  := ncg_z  / ngdpp_z_smooth]
  
  ## Debt/GDP ratio -- adjusted nominal GDP
  dt[, b_ratio := -tiin / ngdpz_smooth / 4]
  
  ## Normalized migration
  dt[, ln_norm := lmig_z / (lhpwa_z * 1000)]
  
  n_vars <- sum(!is.na(dt$r_q) & !is.na(dt$rh_q) & !is.na(dt$c_sh))
  message(sprintf("  %d quarters with core variables available", n_vars))
  
  ## ----------------------------------------------------------------
  ## STEP 2: Trends
  ## ----------------------------------------------------------------
  message("-- Step 2: Trends --")
  
  ## -- Potential output: three paths (A/B/C) --
  if (!is.null(potential_y_trend_override)) {
    ## Path A: injected y_trend
    dt <- merge(dt, potential_y_trend_override[, .(date, y_trend)],
                by = "date", all.x = TRUE)
    ## HP fallback for dates outside override
    dt[is.na(y_trend), y_trend := .hp_trend_log(ngdpp_z_smooth, L$y)]
    message("  y_trend: Path A (override + HP fallback for gaps)")
    
  } else if (!is.null(potential_params)) {
    ## Path B: Kalman potential model pipeline
    ## NOT IMPLEMENTED: this path used to call three builder functions (the
    ## potential-obs constructor, its endpoint extender, and the state-space
    ## spec builder), none of which are defined anywhere in the package.
    ## Fail loud instead of letting a "could not find function" error
    ## surface deep in the call stack. See
    ## .claude/orchestration/track-p-hardening/brief-A2-pathb-failloud.md.
    stop("dynhr_transform: the Kalman potential-output path ",
         "(potential_params=) was never implemented (its builder ",
         "functions do not exist). Use potential_y_trend_override= or ",
         "the HP-filter fallback instead.", call. = FALSE)
  } else {
    ## Path C: HP filter on lockdown-smoothed GDP
    dt[, y_trend := .hp_trend_log(ngdpp_z_smooth, L$y)]
    message("  y_trend: Path C (HP filter on lockdown-smoothed GDP)")
  }
  
  ## -- Judged neutral interest rates --
  .build_neutral_r(dt)
  .build_neutral_rh(dt)
  .build_neutral_rstar(dt)
  message("  Interest rates: judged neutrals (batch_2c anchors)")
  
  ## -- Real exchange rate --
  dt[, rs_trend := .hp_trend_log(rs_real, L$rs)]
  
  ## -- CPI level --
  dt[, p_trend := .hp_trend_log(pcpis, L$p)]
  
  ## -- Relative prices --
  dt[, pn_rel_trend := .hp_trend_log(pn_rel, L$pn_p)]
  dt[, ph_rel_trend := .hp_trend_log(ph_rel, L$ph_p)]
  
  ## -- Real wage --
  dt[, w_real_trend := .hp_trend_log(w_real, L$wrlci)]
  
  ## -- World CPI --
  dt[, pstar_trend := .hp_trend_log(wcpi, L$pstar)]
  
  ## -- World trade price ratios -> reconstruct levels --
  dt[, pm_ps_trend  := .hp_trend_log(pm_ps, L$pm_ps)]
  dt[, pmstar_trend := pm_ps_trend * pstar_trend]
  dt[, px_ps_trend  := .hp_trend_log(px_ps, L$px_ps)]
  dt[, pxstar_trend := px_ps_trend * pstar_trend]
  
  ## -- World output --
  if ("iwgdp_pt" %in% names(dt)) {
    dt[, ystar_trend := iwgdp_pt]
    dt[is.na(ystar_trend), ystar_trend := .hp_trend_log(iwgdp_z, L$y)]
    message("  ystar: using iwgdp_pt as trend")
  } else {
    dt[, ystar_trend := .hp_trend_log(iwgdp_z, L$y)]
  }
  
  ## -- Migration (normalized) --
  dt[, ln_trend := .hp_trend(ln_norm, L$ln)]
  
  ## -- Debt/GDP ratio --
  dt[, b_ratio_trend := .hp_trend_log(b_ratio, L$b_ngdp)]
  
  ## -- GDP shares --
  share_map <- list(
    c  = list(col = "c_sh",  lambda = L$c_gdp),
    ik = list(col = "ik_sh", lambda = L$ik_gdp),
    ih = list(col = "ih_sh", lambda = L$ih_gdp),
    x  = list(col = "x_sh",  lambda = L$x_gdp),
    m  = list(col = "m_sh",  lambda = L$m_gdp),
    g  = list(col = "g_sh",  lambda = L$g_gdp)
  )
  for (nm in names(share_map)) {
    scol <- share_map[[nm]]$col
    tcol <- paste0(nm, "_sh_trend")
    dt[, (tcol) := .hp_trend_log(get(scol), share_map[[nm]]$lambda)]
  }
  
  ## Component level trends = share_trend x potential GDP
  dt[, c_trend  := c_sh_trend  * y_trend]
  dt[, ik_trend := ik_sh_trend * y_trend]
  dt[, ih_trend := ih_sh_trend * y_trend]
  dt[, x_trend  := x_sh_trend  * y_trend]
  dt[, m_trend  := m_sh_trend  * y_trend]
  dt[, g_trend  := g_sh_trend  * y_trend]
  
  message("  Trends computed for all 20 varobs")
  
  ## ----------------------------------------------------------------
  ## STEP 3: Compute gaps
  ## ----------------------------------------------------------------
  message("-- Step 3: Gaps --")
  
  ## Interest rate gaps
  dt[, r_    := r_q - r_trend]
  dt[, rh_   := rh_q - rh_trend]
  dt[, rstar_ := rstar_q - rstar_trend]
  
  ## Real exchange rate
  dt[, rs_ := log(rs_real / rs_trend)]
  
  ## Expenditure gaps (actual data, not smoothed -- captures lockdown dip)
  dt[, c_  := (ncp_z  - c_trend)  / y_trend / S$cy]
  dt[, ik_ := (nik_z  - ik_trend) / y_trend / S$iki]
  dt[, ih_ := (nitd_z - ih_trend) / y_trend / S$ihi]
  dt[, x_  := (nx_z   - x_trend)  / y_trend / S$xy]
  dt[, m_  := (nm_z   - m_trend)  / y_trend / S$my]
  dt[, g_  := (ncg_z  - g_trend)  / y_trend / S$gy]
  
  ## Output gap: expenditure identity
  dt[, y_ := g_ * S$gy + S$cy * c_ + S$iki * ik_ +
       S$ihi * ih_ + S$xy * x_ - S$my * m_]
  
  ## Inflation gap
  dt[, dp_actual := pcpis / data.table::shift(pcpis) - 1]
  dt[, dp_trend  := p_trend / data.table::shift(p_trend) - 1]
  dt[, dp_ := dp_actual - dp_trend]
  
  ## Relative price gaps
  dt[, pn_p_ := log(pn_rel / pn_rel_trend)]
  dt[, ph_p_ := log(ph_rel / ph_rel_trend)]
  
  ## Real wage gap
  dt[, wrlci_ := log(w_real / w_real_trend)]
  
  ## World variable gaps
  dt[, pstar_  := log(wcpi / pstar_trend)]
  dt[, pmstar_ := log(pmstar_constructed / pmstar_trend)]
  dt[, pxstar_ := log(pxstar_constructed / pxstar_trend)]
  dt[, ystar_  := log(iwgdp_z / ystar_trend)]
  
  ## Migration gap (normalized)
  dt[, ln_ := ln_norm - ln_trend]
  
  ## Debt/GDP gap
  dt[, b_ := b_ratio - b_ratio_trend]
  
  ## ----------------------------------------------------------------
  ## STEP 4: Assemble estimation matrix
  ## ----------------------------------------------------------------
  message("-- Step 4: Assemble --")
  
  gap_cols <- c("r_", "dp_", "rs_", "y_", "c_", "x_", "m_",
                "ik_", "ih_", "ln_", "b_", "ph_p_", "pn_p_",
                "rh_", "pmstar_", "pxstar_", "ystar_", "pstar_",
                "rstar_", "wrlci_")
  
  gaps <- dt[, c("date", gap_cols), with = FALSE]
  
  if (!is.null(est_start)) gaps <- gaps[date >= data.table::as.IDate(est_start)]
  if (!is.null(est_end)) {
    gaps <- gaps[date <= data.table::as.IDate(est_end)]
  } else {
    complete <- complete.cases(gaps[, -"date"])
    if (any(complete)) gaps <- gaps[1:max(which(complete))]
  }
  
  message(sprintf("  Final: %d quarters (%s to %s)",
                  nrow(gaps), min(gaps$date), max(gaps$date)))
  
  na_check <- sapply(gaps[, -"date"], function(x) sum(is.na(x)))
  if (any(na_check > 0)) {
    message("  [!] Missing data:")
    print(na_check[na_check > 0])
  } else {
    message("  [OK] No missing data in estimation window")
  }
  
  ## ----------------------------------------------------------------
  ## STEP 5: Compare to nzsim_data.csv reference
  ## ----------------------------------------------------------------
  if (!is.null(nzsim_ref_path) && file.exists(nzsim_ref_path)) {
    message("-- Step 5: Comparison to nzsim_data.csv --")
    comp <- .compare_gaps(gaps, nzsim_ref_path)
    attr(gaps, "comparison") <- comp
  }
  
  attr(gaps, "trends") <- dt
  attr(gaps, "lambdas") <- L
  attr(gaps, "shares") <- S
  attr(gaps, "lockdown_ratio") <- dt[, .(date, lockdown_ratio)]
  message("=== Done ===")
  return(gaps)
}
