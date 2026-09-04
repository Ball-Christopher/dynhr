## R/diag-helpers.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## Shared statistical helpers: .effective_sample_size(),
## .ljung_box(), .svd_rank(), .safe_sym_inv(),
## .nz_events(), .add_nz_event_markers(), .numerical_jacobian()
## --------------------------------------------------------------------------

#' Effective sample size (ESS) from a univariate MCMC chain
#'
#' Uses the initial positive-sequence estimator (Geyer 1992): sum
#' autocorrelation pairs until the first negative pair, then apply
#' the standard ESS formula.
#'
#' @param x Numeric vector -- a single MCMC chain for one parameter.
#' @return Numeric scalar -- the effective sample size.
#' @noRd
.effective_sample_size <- function(x) {
  n <- length(x)
  if (n < 4) return(n)

  # Cap lags at 5000: Geyer estimator exits in <200 lags for typical DSGE chains;
  # computing floor(n/2) lags on 400K pooled draws takes minutes for no benefit.
  max_lag <- min(n - 1, floor(n / 2), 5000L)
  acf_vals <- acf(x, lag.max = max_lag, plot = FALSE)$acf[, , 1]

  # Initial positive-sequence estimator: sum pairs (rho_{2k}, rho_{2k+1})
  # until the sum is negative
  sum_rho <- 0
  k <- 1
  while (k + 1 <= length(acf_vals)) {
    pair_sum <- acf_vals[k] + acf_vals[k + 1]
    if (pair_sum < 0) break
    sum_rho <- sum_rho + pair_sum
    k <- k + 2
  }

  # ESS = n / (1 + 2 * sum of autocorrelations)
  # acf_vals[1] is lag-0 = 1.0; we start summing from lag-1
  ess <- n / (1 + 2 * sum_rho)
  return(max(1, ess))
}


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


#' Numerical rank of a matrix from its singular values
#'
#' Applies the standard LAPACK-style tolerance
#' \code{max(dim) * max(sv) * .Machine$double.eps} and counts the singular
#' values above it.  Extracted from the (previously duplicated) identical
#' rank-check used by D1/D25/D27/D28/D30/D37.
#'
#' @param d    Numeric vector of singular values (e.g. \code{svd(M)$d}).
#' @param dims Integer vector of matrix dimensions (e.g. \code{dim(M)}).
#' @return Integer scalar -- the numerical rank.
#' @noRd
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


.svd_rank <- function(d, dims) {
  if (length(d) == 0L) return(0L)
  tol <- max(dims) * max(d) * .Machine$double.eps
  sum(d > tol)
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

#' Rank-normalise a set of draws across chains
#'
#' Combines all chains, ranks all draws, then applies the normal-scores
#' transformation z = qnorm((r - 3/8) / (N + 1/4)).
#' Returns a list of per-chain normalised vectors (same structure as input).
#'
#' @param chains_list List of numeric vectors (one per chain), equal length.
#' @return List of same structure with rank-normalised values.
#' @noRd
.rank_normalise <- function(chains_list) {
  n_chains <- length(chains_list)
  n        <- length(chains_list[[1]])
  all_x    <- unlist(chains_list)
  N        <- length(all_x)
  r        <- rank(all_x, ties.method = "average")
  z        <- qnorm((r - 3/8) / (N + 1/4))
  split(z, rep(seq_len(n_chains), each = n))
}


#' Rank-normalised R-hat (Vehtari et al. 2021, split version)
#'
#' More robust than the classic R-hat for heavy-tailed posteriors.  Each
#' chain is also split in half so that the test is sensitive to within-chain
#' non-stationarity.
#'
#' @param chains_list List of numeric vectors (one per chain).
#' @return Named list: $rhat (scalar), $pass (rhat < 1.01), $pass_loose (< 1.05)
#' @noRd
.rhat_rank_norm <- function(chains_list) {
  m <- length(chains_list)
  if (m < 1) return(list(rhat = NA_real_, pass = NA, pass_loose = NA))
  n <- length(chains_list[[1]])
  if (n < 4) return(list(rhat = NA_real_, pass = NA, pass_loose = NA))

  # Split each chain in half -> 2m chains of length floor(n/2)
  half <- floor(n / 2L)
  split_chains <- unlist(lapply(chains_list, function(ch) {
    list(ch[1:half], ch[(half + 1):(2 * half)])
  }), recursive = FALSE)

  # Rank-normalise across all 2m split chains
  z_chains <- .rank_normalise(split_chains)

  # Classic R-hat on normalised draws
  rhat_val <- .rhat_classic_multi(z_chains)

  list(
    rhat       = rhat_val,
    pass       = !is.na(rhat_val) && rhat_val < 1.01,
    pass_loose = !is.na(rhat_val) && rhat_val < 1.05
  )
}


#' Classic multi-chain R-hat (internal, called after rank normalisation)
#' @noRd
.rhat_classic_multi <- function(chains_list) {
  m <- length(chains_list)
  n <- length(chains_list[[1]])
  if (m < 2 || n < 2) return(NA_real_)

  chain_means <- vapply(chains_list, mean, numeric(1))
  grand_mean  <- mean(chain_means)
  B           <- n / (m - 1) * sum((chain_means - grand_mean)^2)
  W           <- mean(vapply(chains_list, var, numeric(1)))
  if (W < .Machine$double.eps) return(1.0)   # constant chain

  V_hat <- (n - 1) / n * W + B / n
  sqrt(max(V_hat / W, 0))
}


#' Bulk-ESS for a single parameter across multiple chains
#'
#' ESS computed on rank-normalised draws; sensitive to location and spread.
#'
#' @param chains_list List of numeric vectors (one per chain).
#' @return Numeric scalar.
#' @noRd
.ess_bulk_multi <- function(chains_list) {
  z_chains <- .rank_normalise(chains_list)
  n_chains <- length(z_chains)
  n        <- length(z_chains[[1]])

  # Pool all chains then use the Geyer estimator
  .effective_sample_size(unlist(z_chains)) * min(1, n_chains)
}


#' Tail-ESS for a single parameter across multiple chains
#'
#' ESS computed on indicators I(x <= q0.05) and I(x <= q0.95);
#' sensitive to tail behaviour.
#'
#' @param chains_list List of numeric vectors.
#' @return Numeric scalar (min of lower- and upper-tail ESS).
#' @noRd
.ess_tail_multi <- function(chains_list) {
  all_x <- unlist(chains_list)
  q05   <- quantile(all_x, 0.05)
  q95   <- quantile(all_x, 0.95)

  ind_lo <- lapply(chains_list, function(ch) as.numeric(ch <= q05))
  ind_hi <- lapply(chains_list, function(ch) as.numeric(ch <= q95))

  ess_lo <- .effective_sample_size(unlist(ind_lo))
  ess_hi <- .effective_sample_size(unlist(ind_hi))
  min(ess_lo, ess_hi)
}


#' Per-parameter multi-chain convergence summary
#'
#' Computes rank-normalised R-hat, Bulk-ESS, and Tail-ESS for every
#' parameter given a list of draw matrices (one per chain).
#'
#' @param chains_list List of matrices (n_draws x n_params), same dims.
#' @return data.frame with columns: param, rhat, ess_bulk, ess_tail
#' @noRd
.convergence_summary <- function(chains_list) {
  if (length(chains_list) < 1) return(NULL)
  p_names <- colnames(chains_list[[1]])
  n_par   <- ncol(chains_list[[1]])
  if (is.null(p_names)) p_names <- paste0("theta_", seq_len(n_par))

  do.call(rbind, lapply(seq_len(n_par), function(j) {
    per_chain <- lapply(chains_list, function(m) m[, j])
    rh  <- .rhat_rank_norm(per_chain)
    eb  <- .ess_bulk_multi(per_chain)
    et  <- .ess_tail_multi(per_chain)
    data.frame(param    = p_names[j],
               rhat     = rh$rhat,
               ess_bulk = eb,
               ess_tail = et,
               stringsAsFactors = FALSE)
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
