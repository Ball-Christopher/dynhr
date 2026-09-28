## R/diag-helpers.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## Shared statistical helpers: the package's ONLY ESS / split-R-hat
## estimators (.d5_*), .ljung_box(), .safe_sym_inv(),
## .nz_events(), .add_nz_event_markers(), .numerical_jacobian()
## --------------------------------------------------------------------------

#' Ljung-Box test for serial correlation
#'
#' @param x       Numeric vector (e.g. smoothed shocks).
#' @param max_lag Integer -- maximum lag to test. Default 8.
#' @return A list with $statistic, $df, $p_value, $pass (p > 0.05)
#' @noRd
.ljung_box <- function(x, max_lag = 8) {
  n <- length(x)
  if (n < max_lag + 2) {
    return(list(statistic = NA, df = NA, p_value = NA, pass = NA,
                message = "Series too short for Ljung-Box test"))
  }

  acf_vals <- acf(x, lag.max = max_lag, plot = FALSE)$acf[-1, , 1]  # drop lag-0
  k_seq <- seq_len(max_lag)
  Q <- n * (n + 2) * sum(acf_vals^2 / (n - k_seq))
  df <- max_lag
  p_val <- 1 - pchisq(Q, df = df)

  list(
    statistic = Q,
    df        = df,
    p_value   = p_val,
    pass      = p_val > 0.05,
    message   = sprintf("Ljung-Box Q(%d) = %.2f, p = %.4f", max_lag, Q, p_val)
  )
}


#' Overlay an estimated-parameter vector onto a full calibration
#'
#' Returns \code{params} with the entries named in \code{theta} overwritten by
#' \code{theta}'s values; names in \code{theta} that are absent from
#' \code{params} are ignored.  Used by the \code{model_solve_fn} closures that
#' re-solve the model at a perturbed \code{theta}.
#'
#' @param params Named list or numeric vector -- the baseline calibration.
#' @param theta  Named list or numeric vector -- values to overlay.
#' @return \code{params} with overlapping entries replaced.
#' @noRd
.update_params <- function(params, theta) {
  for (nm in names(theta))
    if (nm %in% names(params)) params[[nm]] <- theta[[nm]]
  params
}


#' Robust symmetric-matrix inverse with pseudo-inverse fallback
#'
#' Tries \code{qr.solve}; on a singular matrix falls back to
#' \code{MASS::ginv} (or a manual SVD pseudo-inverse when MASS is absent).
#' This is the exact fallback chain previously inlined in D27 and D28.
#' Note: callers that need an adaptive ridge (D29) or a rank-deficiency
#' pre-check (D1) keep their bespoke logic -- this helper deliberately
#' covers only the genuinely identical D27/D28 pattern.
#'
#' @param M A square (symmetric) matrix to invert.
#' @return The inverse (or Moore-Penrose pseudo-inverse) of \code{M}.
#' @noRd
.safe_sym_inv <- function(M) {
  k <- nrow(M)
  tryCatch(qr.solve(M, diag(k)), error = function(e) {
    if (requireNamespace("MASS", quietly = TRUE)) {
      MASS::ginv(M)
    } else {
      sv <- svd(M); d <- sv$d
      d[d < max(d) * 1e-12] <- Inf
      sv$v %*% diag(1 / d, k) %*% t(sv$u)
    }
  })
}


#' NZ macroeconomic event dates for annotation in plots
#'
#' @return A data.frame with columns: date, label, type
#' @noRd
.nz_events <- function() {
  data.frame(
    date = as.Date(c(
      "1997-07-01",   # Asian Financial Crisis onset
      "2001-03-01",   # Dot-com recession
      "2007-12-01",   # GFC onset
      "2009-03-01",   # GFC trough
      "2010-09-04",   # Canterbury earthquake
      "2011-02-22",   # Christchurch earthquake
      "2016-11-14",   # Kaikoura earthquake
      "2020-03-25",   # COVID-19 NZ lockdown (Alert Level 4)
      "2021-08-17",   # Delta lockdown
      "2022-06-01",   # Inflation surge peak
      "2023-05-24"    # Cyclone Gabrielle aftermath / rate peak
    )),
    label = c(
      "Asian Crisis",
      "Dot-com",
      "GFC onset",
      "GFC trough",
      "Canterbury EQ",
      "Christchurch EQ",
      "Kaikoura EQ",
      "COVID lockdown",
      "Delta lockdown",
      "Inflation peak",
      "Rate peak"
    ),
    type = c(
      "international", "international", "international", "international",
      "domestic", "domestic", "domestic",
      "domestic", "domestic",
      "domestic", "domestic"
    ),
    stringsAsFactors = FALSE
  )
}


#' Add NZ event markers to a ggplot with a date x-axis
#'
#' @param p       An existing ggplot object with a Date x-axis
#' @param events  Data frame from .nz_events() (or subset thereof)
#' @param y_pos   Vertical position for labels (default: Inf = top)
#' @return Modified ggplot with vertical lines and labels
#' @noRd
.add_nz_event_markers <- function(p, events = NULL, y_pos = Inf) {
  if (is.null(events)) events <- .nz_events()
  .ensure_ggplot2()

  # inherit.aes = FALSE: these annotation layers must NOT pick up the host
  # plot's global aesthetics (e.g. fill = shock), which the events data lacks.
  #
  # clip = "off" + expanded top margin prevent label text from being clipped at
  # the panel boundary when labels are at y = Inf (rotated 90 degrees at top).
  # ggrepel is not used here to keep the dependency footprint minimal; instead
  # we use a slightly negative hjust offset so the text sits just above the
  # vline and within the expanded margin.
  p +
    ggplot2::geom_vline(data = events,
                        ggplot2::aes(xintercept = date),
                        inherit.aes = FALSE,
                        linetype = "dashed", colour = dynhr_colours$grey,
                        linewidth = 0.3, alpha = 0.7) +
    ggplot2::geom_text(data = events,
                       ggplot2::aes(x = date, y = y_pos, label = label),
                       inherit.aes = FALSE,
                       angle = 90, vjust = 0.3, hjust = 1.05,
                       size = 2.2, colour = dynhr_colours$grey,
                       family = "sans") +
    ggplot2::theme(
      # Expand top margin so rotated labels at y=Inf do not get clipped by the
      # panel boundary. The coord clip is turned off so the text can overflow.
      plot.margin = ggplot2::margin(t = 55, r = 5.5, b = 5.5, l = 5.5)
    ) +
    ggplot2::coord_cartesian(clip = "off")
}


# ---------------------------------------------------------------------------
# Rank-normalised R-hat and Bulk/Tail-ESS  (Vehtari et al. 2021)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# CONVERGENCE ESTIMATORS: split-R-hat and ESS (Vehtari et al. 2021)
# ---------------------------------------------------------------------------
# THE package's ESS / R-hat implementation -- there is exactly one, and this is
# it. It follows the paper (and the reference implementation in the `posterior`
# package, verified against posterior 1.7.0): split chains, rank-normalise with
# the Blom offset 3/8, R-hat = max(bulk, folded) split-R-hat, and multi-chain
# ESS with Geyer's initial positive + monotone sequence using the
# between-chain-aware autocorrelation 1 - (W - mean acov_t) / var_plus.
#
# 0.9.4 (ledger A5): these replace a SECOND, WRONG set of shared helpers that
# used to live here -- `.effective_sample_size()`, `.rank_normalise()`,
# `.rhat_rank_norm()`, `.rhat_classic_multi()`, `.ess_bulk_multi()`,
# `.ess_tail_multi()`, `.convergence_summary()`. Three defects, all silent:
#
#   1. the Geyer sum started at `acf_vals[1]`, which is LAG 0 (= 1), so the
#      first pair was 1 + rho_1 instead of rho_1 + rho_2. For iid draws that
#      gives ESS = n / (1 + 2*1) = n/3 -- every ESS the package reported was
#      about 3x too low, and "bulk ESS >= 1000" was really "ESS >= ~3000";
#   2. multi-chain ESS CONCATENATED the chains (`unlist`) instead of using the
#      between-chain variance, so between-chain drift inflated it;
#   3. R-hat omitted the FOLDED (scale) component, so chains agreeing in mean
#      but differing in variance passed with R-hat ~ 1.
#
# The `.d5_` prefix is historical (the correct estimators were first written
# for D5); it is kept so that there is one name per function and no aliases.

#' Split every chain of an n x m matrix in half (odd n drops the middle draw)
#' @noRd
.d5_split <- function(x) {
  n <- nrow(x)
  half <- floor(n / 2)
  cbind(x[seq_len(half), , drop = FALSE],
        x[(n - half + 1L):n, , drop = FALSE])
}

#' Rank-normalised z-scores over all draws (ties averaged)
#' @noRd
.d5_zscale <- function(x) {
  r <- rank(as.vector(x), ties.method = "average")
  z <- stats::qnorm((r - 3 / 8) / (length(r) + 1 / 4))
  matrix(z, nrow = nrow(x), ncol = ncol(x))
}

#' TRUE if the n x m draws cannot carry an R-hat / ESS
#' @noRd
.d5_degenerate <- function(x) {
  nrow(x) < 4L || any(!is.finite(x)) ||
    abs(max(x) - min(x)) < .Machine$double.eps
}

#' Basic split-free R-hat of an n x m matrix
#' @noRd
.d5_rhat_basic <- function(x) {
  if (.d5_degenerate(x) || ncol(x) < 2L) return(NA_real_)
  n  <- nrow(x)
  B  <- n * stats::var(colMeans(x))
  W  <- mean(apply(x, 2, stats::var))
  sqrt((B / W + n - 1) / n)
}

#' Autocovariance (denominator n) via FFT, lags 0..n-1
#' @noRd
.d5_autocov <- function(y) {
  n  <- length(y)
  # A constant (split) chain -- common for tail indicators -- has zero
  # autocovariance at every lag, not 0/0.
  if (stats::var(y) == 0) return(numeric(n))
  M  <- stats::nextn(n)
  yc <- c(y - mean(y), rep.int(0, 2L * M - n))
  f  <- stats::fft(yc)
  ac <- Re(stats::fft(Conj(f) * f, inverse = TRUE))[seq_len(n)]
  ac / ac[1] * (sum((y - mean(y))^2) / n)
}

#' Basic multi-chain ESS of an n x m matrix (Vehtari et al. 2021, eqs 10-11)
#' @noRd
.d5_ess_basic <- function(x) {
  if (.d5_degenerate(x)) return(NA_real_)
  n <- nrow(x); m <- ncol(x)
  acov <- vapply(seq_len(m), function(j) .d5_autocov(x[, j]), numeric(n))
  acov <- matrix(acov, nrow = n)
  mean_var <- mean(acov[1, ]) * n / (n - 1)
  var_plus <- mean_var * (n - 1) / n
  if (m > 1L) var_plus <- var_plus + stats::var(colMeans(x))
  rho <- function(t) 1 - (mean_var - mean(acov[t + 1L, ])) / var_plus
  rho_t <- numeric(n)
  t <- 0L
  r_even <- 1
  r_odd  <- rho(1L)
  rho_t[1:2] <- c(r_even, r_odd)
  while (t < n - 5L && is.finite(r_even + r_odd) && r_even + r_odd > 0) {
    t <- t + 2L
    r_even <- rho(t)
    r_odd  <- rho(t + 1L)
    if (r_even + r_odd >= 0) {
      rho_t[t + 1L] <- r_even
      rho_t[t + 2L] <- r_odd
    }
  }
  max_t <- t
  if (r_even > 0) rho_t[max_t + 1L] <- r_even
  # Geyer initial monotone sequence
  t <- 0L
  while (t <= max_t - 4L) {
    t <- t + 2L
    if (rho_t[t + 1L] + rho_t[t + 2L] > rho_t[t - 1L] + rho_t[t]) {
      rho_t[t + 1L] <- (rho_t[t - 1L] + rho_t[t]) / 2
      rho_t[t + 2L] <- rho_t[t + 1L]
    }
  }
  S   <- n * m
  tau <- -1 + 2 * sum(rho_t[seq_len(max_t)]) + rho_t[max_t + 1L]
  tau <- max(tau, 1 / log10(S))
  S / tau
}

#' Per-parameter rank-normalised split R-hat, Bulk-ESS and Tail-ESS
#'
#' @param chains_list List of equal-shape n x p draw matrices (one per chain;
#'   a single chain is allowed and is split in two for R-hat).
#' @return data.frame(param, rhat, ess_bulk, ess_tail).
#' @noRd
.d5_convergence <- function(chains_list, param_names) {
  do.call(rbind, lapply(seq_along(param_names), function(j) {
    x  <- vapply(chains_list, function(mm) as.numeric(mm[, j]),
                 numeric(nrow(chains_list[[1]])))
    x  <- matrix(x, nrow = nrow(chains_list[[1]]))
    xs <- .d5_split(x)
    if (.d5_degenerate(xs)) {
      rhat <- eb <- et <- NA_real_
    } else {
      rhat <- max(.d5_rhat_basic(.d5_zscale(xs)),
                  .d5_rhat_basic(.d5_zscale(.d5_split(abs(x - stats::median(x))))))
      eb   <- .d5_ess_basic(.d5_zscale(xs))
      et   <- min(vapply(c(0.05, 0.95), function(p) {
        q <- stats::quantile(x, p, names = FALSE)
        .d5_ess_basic(.d5_split(x <= q) + 0)
      }, numeric(1)))
    }
    data.frame(param = param_names[j], rhat = rhat, ess_bulk = eb,
               ess_tail = et, stringsAsFactors = FALSE)
  }))
}


#' Bayesian Fraction of Missing Information (BFMI) for one chain
#'
#' BFMI = Var(E[m+1] - E[m]) / Var(E[m])  where E[m] = -H[m] (negative Hamiltonian).
#' Values < 0.3 indicate the sampler is not exploring heavy tails adequately.
#'
#' @param energy Numeric vector of Hamiltonian values per post-warmup iteration.
#'   These are the negated log joint (logpost - kinetic_energy) stored by
#'   the NUTS sampler.
#' @return Numeric scalar BFMI.
#' @noRd
.bfmi <- function(energy) {
  if (length(energy) < 3 || !all(is.finite(energy))) return(NA_real_)
  E        <- -energy   # E = -H = negative Hamiltonian
  diffs    <- diff(E)
  var(diffs) / var(E)
}


#' Compute numerical Jacobian via central finite differences
#'
#' @param fn    Function mapping parameter vector theta -> output vector
#' @param theta Parameter vector at which to evaluate
#' @param eps   Step size for finite differences (default 1e-5)
#' @return Jacobian matrix (length(fn(theta)) x length(theta))
#' @noRd
.numerical_jacobian <- function(fn, theta, eps = 1e-5) {
  # Baseline evaluation
  f0 <- fn(theta)
  n_out <- length(f0)
  n_par <- length(theta)

  J <- matrix(NA_real_, nrow = n_out, ncol = n_par)

  for (j in seq_len(n_par)) {
    theta_plus  <- theta
    theta_minus <- theta
    theta_plus[j]  <- theta[j] + eps
    theta_minus[j] <- theta[j] - eps

    f_plus  <- fn(theta_plus)
    f_minus <- fn(theta_minus)

    # Guard: if fn returns wrong-length output or NAs, use NA column
    if (length(f_plus) != n_out || length(f_minus) != n_out) {
      J[, j] <- NA_real_
    } else {
      J[, j] <- (f_plus - f_minus) / (2 * eps)
    }
  }

  J
}


## ---------------------------------------------------------------------------
## Shared identification-diagnostic helpers (D27 / D28)
##
## Both diagnostics used to carry byte-identical private copies
## (.stoch_simul_internal_d27 / _d28 and .moments_from_dr_d27 /
## .compute_moments_from_dr).  One implementation each, 2026-09.
## ---------------------------------------------------------------------------

#' Compile + solve a model at a parameter draw (identification diagnostics)
#'
#' @param dr_order Perturbation order (default 1).
#' @return DecisionRules, or NULL if compilation / steady state / perturbation
#'   failed.
#' @noRd
.stoch_simul_internal_diag <- function(model, params, dr_order = 1L) {
  compiled <- compile_model(model, verbose = FALSE)
  if (is.null(compiled)) return(NULL)
  ss <- solve_steady(compiled, params = params, verbose = FALSE)
  if (is.null(ss) || !isTRUE(ss$converged)) return(NULL)
  dr <- solve_perturbation(model, compiled, ss$values, params, order = dr_order, verbose = FALSE)
  dr
}

#' Model-implied moments from decision rules (identification diagnostics)
#'
#' @param dr DecisionRules object.
#' @return List with $sigma_y, $acf_y, or NULL.
#' @noRd
.moments_from_dr <- function(dr, model = NULL, params = NULL) {
  if (!is.null(model)) {
    moments <- compute_moments(dr, model, params = params)
  } else {
    moments <- compute_moments(dr)
  }
  if (!is.list(moments) || is.null(moments$var_cov)) return(NULL)
  list(
    sigma_y = moments$var_cov,
    acf_y = if (!is.null(moments$autocorr)) moments$autocorr else NULL
  )
}
