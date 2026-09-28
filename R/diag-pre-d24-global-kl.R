## R/diag-pre-d24-global-kl.R
## --------------------------------------------------------------------------
## Phase C: D24 Global Identification via KL Minimisation.
##
## Implements the Kullback-Leibler divergence-based global identification
## diagnostic of Qu & Tkachenko (2017, RES): for a Gaussian model
## the per-observation KL divergence between the observable processes at
## theta0 and theta is a frequency-domain integral of the two spectral
## density matrices, and theta0 is globally identified (at distance c) when
## the MINIMUM of that KL over all theta at distance >= c from theta0 is
## bounded away from zero.
##
## References:
##   Qu, Z., & Tkachenko, D. (2017). Global identification in DSGE models
##     allowing for indeterminacy. Review of Economic Studies 84(3),
##     1306-1345 (KL criterion, Section 3).
##   Qu, Z., & Tkachenko, D. (2012). Identification and frequency domain
##     analysis of DSGE models. Quantitative Economics 3(1), 95-132.
## --------------------------------------------------------------------------

## Per-observation KL rate KL(f0 || f1) between two zero-mean stationary
## Gaussian processes from their (complex Hermitian) spectral density
## matrices on a uniform grid over [0, pi]:
##   KL = 1/(4 pi) int_{-pi}^{pi} [tr(f1^-1 f0) - log det(f1^-1 f0) - n] dw
##      = 1/(2 pi) int_0^pi  [...] dw        (integrand is even in w)
## evaluated by the trapezoid rule (spectrally accurate for a smooth periodic
## integrand).  The integrand is sum(lam - log(lam) - 1) over the eigenvalues
## lam of f1^{-1/2} f0 f1^{-1/2}, which is >= 0 term by term.  Returns NA when
## f1 is singular at any grid frequency (KL undefined / infinite).
## @noRd
.d24_kl_rate <- function(S0_list, S1_list, freq_grid) {
  n_f <- length(freq_grid)
  vals <- numeric(n_f)
  for (k in seq_len(n_f)) {
    # .eigen_hermitian_safe(): complex Hermitian eigen() (zheev) segfaults
    # under Apple vecLib, see R/whittle-likelihood.R.
    e1 <- .eigen_hermitian_safe(S1_list[[k]])
    d1 <- e1$values
    if (!all(is.finite(d1)) || min(d1) <= 1e-12 * max(abs(d1), 1e-300))
      return(NA_real_)
    W <- e1$vectors %*% (Conj(t(e1$vectors)) / sqrt(d1))
    M <- W %*% S0_list[[k]] %*% Conj(t(W))
    lam <- .eigen_hermitian_safe(M, only_values = TRUE)$values
    if (!all(is.finite(lam)) || min(lam) <= 0) return(NA_real_)
    vals[k] <- sum(lam - log(lam) - 1)
  }
  h <- freq_grid[2] - freq_grid[1]
  (sum(vals) - (vals[1] + vals[n_f]) / 2) * h / (2 * pi)
}

#' D24. Global identification via KL divergence (Qu & Tkachenko 2017)
#'
#' Assesses global parameter identification from the Kullback-Leibler
#' divergence between the Gaussian processes of observables implied by the
#' model at the baseline \eqn{\theta_0} and at other parameter values.
#'
#' The per-observation KL divergence is computed in the frequency domain from
#' the (complex Hermitian) spectral density matrices \eqn{f_\theta(\omega)}:
#'   \deqn{KL(\theta_0, \theta) = \frac{1}{4\pi} \int_{-\pi}^{\pi} \left[
#'     \mathrm{tr}(f_\theta^{-1} f_{\theta_0}) -
#'     \log\det(f_\theta^{-1} f_{\theta_0}) - n_y \right] d\omega}
#' (trapezoid rule on \code{n_freq} points over \eqn{[0,\pi]}).  \eqn{T \cdot KL}
#' is the expected log-likelihood ratio of a sample of length \eqn{T}.
#'
#' For each parameter \eqn{k} the diagnostic reports
#' \eqn{\min KL(\theta_0,\theta)} over the shell
#' \eqn{\max_i |\theta_i - \theta_{0,i}| / \delta_i = 1} with
#' \eqn{\theta_k = \theta_{0,k} \pm \delta_k}, where
#' \eqn{\delta_i = \mathrm{perturb\_scale}_i \cdot |\theta_{0,i}|}
#' (\code{perturb_scale_i} itself when \eqn{\theta_{0,i} = 0}).  With
#' \code{search = "joint"} the other parameters are free to COMPENSATE within
#' the box (L-BFGS-B), so directions of observational equivalence
#' (e.g. only a product of two parameters identified) are found; with
#' \code{search = "axis"} only one-at-a-time moves are evaluated, which cannot
#' detect compensating directions.  This is the Qu-Tkachenko neighbourhood
#' exclusion restricted to the shell at distance \eqn{c = 1} (in units of
#' \eqn{\delta}); the one-at-a-time profile on \code{n_perturb} grid points is
#' also returned for plotting.
#'
#' \strong{What gates the badge: a three-tier band on the per-observation KL,
#' not a single \eqn{T \cdot KL} cutoff.} Qu & Tkachenko (2017, Review of
#' Economic Studies 84(3):1306-1345, "Global Identification in DSGE Models
#' Allowing for Indeterminacy") show that global identification fails
#' \emph{iff} the minimised frequency-domain KL divergence is exactly zero.
#' Their criterion is therefore a zero/non-zero statement; they specify no
#' numeric per-observation cutoff, and in their applications they report the
#' KL values and plots and leave the detectability judgement to the reader.
#' Nothing in the literature supports a single bright line, so the diagnostic
#' keeps three bands -- \code{min_kl > 10 * kl_threshold} "strong",
#' \code{> kl_threshold} "adequate", otherwise "weak / unidentified" -- and
#' only the "weak" band FAILs. \code{kl_threshold = 0.01} per observation, and
#' the factor 10 between bands, are package choices with no literature source.
#'
#' \eqn{T \cdot KL} is reported (when \code{T_obs} is given) as the
#' \emph{power} interpretation of the same number, and never gates: it is, up
#' to convention, the expected log-likelihood ratio a sample of length \eqn{T}
#' delivers against the nearest observationally-close alternative. One nat is
#' barely detectable in principle; a conventional 80\% power at the 5\% level
#' for a single free direction corresponds to a noncentrality near 7.85, i.e.
#' roughly 4 nats of expected log-LR. (That benchmark is a standard
#' noncentral-chi-square power calculation applied to Qu-Tkachenko's own
#' \eqn{T \cdot KL} quantity, not a number stated in their paper.)
#'
#' A re-solve closure (\code{dr_solve_fn} or \code{abcd_solve_fn}) is required:
#' without one the spectral density cannot change with \eqn{\theta} and the
#' diagnostic is skipped.  The baseline spectral density is computed from the
#' SAME closure at \code{params}; a supplied \code{dr} (or Ramsey DR) must be
#' reproduced by \code{dr_solve_fn(params)} or the call aborts.
#'
#' @param dr             First-order decision rules (\code{ghx}, \code{ghu};
#'   \code{ghu} excludes \eqn{\Sigma_e}).  Used for dimensions and
#'   \code{Z}/\code{Sigma_e} defaults.
#' @param model          dynhr_mod (only \code{model$obs_mat} is used).
#' @param params         Named parameter vector at the baseline.
#' @param ramsey_result  Optional \code{dynhr_ramsey_result2}; its DR replaces
#'   \code{dr}, and \code{dr_solve_fn} must then return Ramsey DRs.
#' @param obs_mat        Observation matrix (n_obs x n_state), current-state
#'   convention; default \code{dr$Z}, \code{model$obs_mat}, else identity.
#' @param D_mat          Direct shock loading (n_obs x n_shock); default 0.
#'   Used by \code{spectral = "exact"} only.
#' @param Sigma_e        Shock covariance (n_shock x n_shock).  Default
#'   \code{dr$Sigma_e} (or the \code{abcd_solve_fn} baseline's), else identity
#'   with a warning.  A closure result carrying \code{$Sigma_e} overrides it at
#'   that \eqn{\theta}, so shock-std parameters are identified.
#' @param param_names    Optional parameter names.
#' @param n_freq         Frequency grid points on \eqn{[0,\pi]} (default 256).
#' @param spectral       \code{"exact"} (default; lagged-state spectral density
#'   including \code{D_mat}) or \code{"companion"} (\eqn{Z(I-Tz)^{-1}R}, ignores
#'   D).
#' @param perturb_scale  Relative distance \eqn{c} (scalar or named vector;
#'   default 0.10 = 10\% of \eqn{|\theta_0|}).
#' @param n_perturb      Grid points of the one-at-a-time profile (>= 3).
#' @param search         \code{"auto"} (joint when \eqn{n_{par} \le 10}),
#'   \code{"joint"} or \code{"axis"}.
#' @param maxit          Max L-BFGS-B iterations per face in the joint search.
#' @param T_obs          Optional sample size; adds \eqn{T \cdot KL} (expected
#'   log-likelihood ratio) to the output. Reported, never gated -- see
#'   Details.
#' @param kl_threshold   Per-observation KL defining the three-tier band
#'   (default 0.01): \code{> 10 * kl_threshold} is "strong",
#'   \code{> kl_threshold} is "adequate", at or below it is
#'   "weak / unidentified" and drives the FAIL. This is a package choice;
#'   Qu & Tkachenko (2017) specify no numeric cutoff (see Details).
#' @param verbose        Print progress messages.
#' @param dr_solve_fn    Function theta -> DR list (\code{ghx}, \code{ghu},
#'   optional \code{state_idx}, \code{D_mat}, \code{Sigma_e}); return
#'   \code{NULL} where the model cannot be solved.
#' @param abcd_solve_fn  Function theta -> \code{list(A, B, C, D[, Sigma_e])}
#'   in the lagged convention \eqn{s_t = A s_{t-1} + B\epsilon_t},
#'   \eqn{y_t = C s_{t-1} + D\epsilon_t} (as for D37).  Takes precedence over
#'   \code{dr}/\code{dr_solve_fn}.
#' @param meta           Plot metadata.
#'
#' @return A \code{dynhr_diagnostic}; \code{result} holds \code{kl_matrix}
#'   (n_par x n_perturb, one-at-a-time KL), \code{global_strength} (per
#'   parameter: \code{axis_min_kl}, \code{min_kl}, \code{T_kl},
#'   \code{globally_identified}, \code{ident_quality}), \code{global_min}
#'   (minimum KL over the whole shell and its argmin \eqn{\theta}),
#'   \code{kl_threshold}, \code{search}, \code{spectral}, \code{Sigma_e}.
#'
#' @noRd
d24_global_kl_identification <- function(dr = NULL,
                                          model = NULL,
                                          params = NULL,
                                          ramsey_result = NULL,
                                          obs_mat = NULL,
                                          D_mat = NULL,
                                          Sigma_e = NULL,
                                          param_names = NULL,
                                          n_freq = 256L,
                                          spectral = c("exact", "companion"),
                                          perturb_scale = 0.10,
                                          n_perturb = 5L,
                                          search = c("auto", "joint", "axis"),
                                          maxit = 50L,
                                          T_obs = NULL,
                                          verbose = FALSE,
                                          dr_solve_fn = NULL,
                                          abcd_solve_fn = NULL,
                                          kl_threshold = 0.01,
                                          meta = NULL) {
  spectral <- match.arg(spectral)
  search   <- match.arg(search)

  .skip <- function(msg, reason) {
    .make_result(pass = NA,
                 summary = paste("D24 Global KL identification:", msg),
                 llm_summary = sprintf("[INFO] D24 | status=skipped reason=%s", reason))
  }

  # ---- 1. Resolve decision rules ----
  use_ramsey <- FALSE
  if (!is.null(ramsey_result) && inherits(ramsey_result, "dynhr_ramsey_result2")) {
    ramsey_dr <- ramsey_result$ramsey_dr$ramsey_dr
    if (!is.null(ramsey_dr) && !is.null(ramsey_dr$ghx)) {
      dr <- ramsey_dr
      use_ramsey <- TRUE
      if (verbose) .dynhr_cat("[d24] Using Ramsey-optimal decision rules.\n")
    }
  }
  use_abcd <- is.function(abcd_solve_fn)

  if (!use_abcd && (is.null(dr) || is.null(dr$ghx) || is.null(dr$ghu)))
    return(.skip("decision rules missing ghx/ghu.", "no_dr"))
  if (is.null(params)) return(.skip("params required.", "no_params"))
  if (!use_abcd && !is.function(dr_solve_fn))
    return(.skip(paste("no re-solve closure (dr_solve_fn / abcd_solve_fn); the",
                       "spectral density cannot vary with theta."), "no_resolve_fn"))

  n_par <- length(params)
  if (is.null(param_names))
    param_names <- names(params) %||% paste0("theta_", seq_len(n_par))
  if (length(param_names) != n_par)
    .dynhr_abort(sprintf("D24: length(param_names) (%d) != length(params) (%d).",
                         length(param_names), n_par))
  params <- stats::setNames(as.numeric(params), param_names)
  if (!all(is.finite(params))) .dynhr_abort("D24: params must be finite.")

  # ---- 2. State-space at theta (NULL when unsolvable / nonstationary) ----
  sigma_default <- NULL
  if (use_abcd) {
    ss0_raw <- abcd_solve_fn(params)
    if (is.null(ss0_raw) || is.null(ss0_raw$A) || is.null(ss0_raw$B) ||
        is.null(ss0_raw$C) || is.null(ss0_raw$D))
      .dynhr_abort("D24: abcd_solve_fn(params) must return list(A, B, C, D[, Sigma_e]).")
    n_state <- nrow(as.matrix(ss0_raw$A))
    n_shock <- ncol(as.matrix(ss0_raw$B))
    n_obs   <- nrow(as.matrix(ss0_raw$C))
    sigma_default <- ss0_raw$Sigma_e
  } else {
    ghx <- as.matrix(dr$ghx); ghu <- as.matrix(dr$ghu)
    if (nrow(ghu) != nrow(ghx))
      .dynhr_abort(sprintf("D24: ghu rows (%d) != ghx rows (%d).", nrow(ghu), nrow(ghx)))
    n_state <- ncol(ghx)
    n_shock <- ncol(ghu)
    sidx0 <- dr$state_idx %||% seq_len(n_state)
    if (length(sidx0) != n_state)
      .dynhr_abort(sprintf("D24: length(state_idx) (%d) != ncol(ghx) (%d).",
                           length(sidx0), n_state))
    if (is.null(obs_mat)) {
      if (!is.null(dr$Z)) {
        obs_mat <- dr$Z
      } else if (!is.null(model$obs_mat)) {
        obs_mat <- model$obs_mat
      } else if (!is.null(ramsey_result$augmented_model$obs_mat)) {
        obs_mat <- ramsey_result$augmented_model$obs_mat
      } else if (use_ramsey && !is.null(ramsey_result$meta$n_orig_vars)) {
        obs_mat <- diag(n_state)[seq_len(min(ramsey_result$meta$n_orig_vars, n_state)), , drop = FALSE]
      } else {
        obs_mat <- diag(n_state)
      }
    }
    obs_mat <- as.matrix(obs_mat)
    if (ncol(obs_mat) != n_state)
      .dynhr_abort(sprintf("D24: obs_mat columns (%d) != n_state (%d).", ncol(obs_mat), n_state))
    n_obs <- nrow(obs_mat)
    if (is.null(D_mat)) D_mat <- matrix(0, n_obs, n_shock)
    D_mat <- as.matrix(D_mat)
    if (nrow(D_mat) != n_obs || ncol(D_mat) != n_shock)
      .dynhr_abort(sprintf("D24: D_mat dims (%dx%d) != expected (%dx%d).",
                           nrow(D_mat), ncol(D_mat), n_obs, n_shock))
    sigma_default <- dr$Sigma_e
  }

  sigma_source <- "argument"
  if (is.null(Sigma_e)) {
    if (!is.null(sigma_default)) {
      Sigma_e <- sigma_default; sigma_source <- "model"
    } else {
      Sigma_e <- diag(n_shock); sigma_source <- "identity_default"
      .dynhr_warn("D24: Sigma_e not supplied; using the identity shock covariance. ",
                  "ghu/B exclude Sigma_e, so pass sigma_e = shock_cov(model) unless ",
                  "the shocks really have unit variance.")
    }
  }
  Sigma_e <- as.matrix(Sigma_e)
  if (nrow(Sigma_e) != n_shock || ncol(Sigma_e) != n_shock)
    .dynhr_abort(sprintf("D24: Sigma_e dims (%dx%d) != n_shock=%d.",
                         nrow(Sigma_e), ncol(Sigma_e), n_shock))

  .ss_at <- function(theta) {
    if (use_abcd) {
      o <- abcd_solve_fn(theta)
      if (is.null(o) || is.null(o$A) || is.null(o$B) || is.null(o$C) || is.null(o$D))
        return(NULL)
      Tm <- as.matrix(o$A); Rm <- as.matrix(o$B); Zm <- as.matrix(o$C); Dm <- as.matrix(o$D)
      timing <- "lagged"
    } else {
      o <- dr_solve_fn(theta)
      if (is.null(o) || is.null(o$ghx) || is.null(o$ghu)) return(NULL)
      gx <- as.matrix(o$ghx); gu <- as.matrix(o$ghu)
      sidx <- o$state_idx %||% seq_len(ncol(gx))
      if (ncol(gx) != n_state || length(sidx) != n_state || ncol(gu) != n_shock ||
          nrow(gu) != nrow(gx))
        return(NULL)
      Tm <- gx[sidx, , drop = FALSE]; Rm <- gu[sidx, , drop = FALSE]
      Zm <- obs_mat
      Dm <- if (!is.null(o$D_mat)) as.matrix(o$D_mat) else D_mat
      timing <- "current"
    }
    Sg <- if (!is.null(o$Sigma_e)) as.matrix(o$Sigma_e) else Sigma_e
    if (!all(dim(Tm) == c(n_state, n_state)) || !all(dim(Rm) == c(n_state, n_shock)) ||
        !all(dim(Zm) == c(n_obs, n_state)) || !all(dim(Dm) == c(n_obs, n_shock)) ||
        !all(dim(Sg) == c(n_shock, n_shock)))
      return(NULL)
    if (!all(is.finite(Tm)) || !all(is.finite(Rm)) || !all(is.finite(Sg)))
      return(NULL)
    if (max(Mod(eigen(Tm, only.values = TRUE)$values)) >= 1 - 1e-8)
      return(NULL)                                  # nonstationary: no spectrum
    list(T = Tm, R = Rm, Z = Zm, D = Dm, Sigma_e = Sg, timing = timing,
         ghx = if (use_abcd) NULL else as.matrix(o$ghx))
  }

  n_freq <- max(as.integer(n_freq), 8L)
  freq_grid <- seq(0, pi, length.out = n_freq)

  .spec <- function(ss) {
    lapply(freq_grid, function(om) {
      if (spectral == "exact") {
        # current-state y_t = Z s_t + D e_t  ==  lagged Z T s_{t-1} + (Z R + D) e_t
        lagged <- ss$timing == "lagged"
        .spectral_density_core(om, TT = ss$T, RR = ss$R,
                               ZZ = if (lagged) ss$Z else ss$Z %*% ss$T,
                               DD = if (lagged) ss$D else ss$Z %*% ss$R + ss$D,
                               Sigma_e = ss$Sigma_e)
      } else {
        .spectral_density_core_current_no_d(om, T_mat = ss$T, R_mat = ss$R,
                                            obs_mat = ss$Z, Sigma_e = ss$Sigma_e)
      }
    })
  }

  ss_base <- .ss_at(params)
  if (is.null(ss_base))
    .dynhr_abort("D24: the re-solve closure failed (or is nonstationary) at the baseline params.")
  if (!use_abcd) {
    same <- all(dim(ss_base$ghx) == dim(as.matrix(dr$ghx))) &&
      max(abs(ss_base$ghx - as.matrix(dr$ghx))) <= 1e-6 * max(1, max(abs(dr$ghx)))
    if (!same)
      .dynhr_abort("D24: dr_solve_fn(params)$ghx does not reproduce the supplied ",
                   if (use_ramsey) "Ramsey " else "", "decision rules; the baseline ",
                   "and perturbed spectra would come from different models.")
  }
  S0 <- .spec(ss_base)

  # Stochastic singularity (n_obs > rank of f): KL is undefined.
  S_mid <- S0[[max(2L, round(n_freq / 3))]]
  ev_mid <- .eigen_hermitian_safe(S_mid, only_values = TRUE)$values
  if (min(ev_mid) <= 1e-10 * max(ev_mid))
    return(.skip(sprintf(paste("baseline spectral density is singular (%d observables,",
                               "rank < %d; stochastic singularity), so the KL divergence",
                               "is undefined. Supply obs_mat with n_obs <= n_shock or add",
                               "measurement error."), n_obs, n_obs), "singular_spectrum"))

  .kl_at <- function(theta) {
    ss <- .ss_at(theta)
    if (is.null(ss)) return(NA_real_)
    .d24_kl_rate(S0, .spec(ss), freq_grid)
  }
  kl_self <- .d24_kl_rate(S0, S0, freq_grid)

  # ---- 3. Perturbation scales ----
  if (length(perturb_scale) == 1L) {
    perturb_scale <- stats::setNames(rep(perturb_scale, n_par), param_names)
  } else if (!is.null(names(perturb_scale)) && all(param_names %in% names(perturb_scale))) {
    perturb_scale <- perturb_scale[param_names]
  } else if (length(perturb_scale) == n_par) {
    names(perturb_scale) <- param_names
  } else {
    .dynhr_abort("D24: perturb_scale must be a scalar, length n_par, or named for every parameter.")
  }
  if (!all(is.finite(perturb_scale)) || any(perturb_scale <= 0))
    .dynhr_abort("D24: perturb_scale must be positive and finite.")
  delta <- perturb_scale * ifelse(params == 0, 1, abs(params))

  # ---- 4. One-at-a-time profile ----
  n_perturb <- max(as.integer(n_perturb), 3L)
  grid_frac <- seq(-1, 1, length.out = n_perturb)
  grid_frac[abs(grid_frac) < 1e-12] <- 0
  kl_matrix <- matrix(NA_real_, nrow = n_par, ncol = n_perturb,
                      dimnames = list(param_names, ifelse(grid_frac == 0, "0", sprintf("%+.2gc", grid_frac))))
  for (i in seq_len(n_par)) {
    if (verbose) .dynhr_cat(sprintf("[d24] Parameter %d/%d: %s\n", i, n_par, param_names[i]))
    for (j in seq_len(n_perturb)) {
      if (grid_frac[j] == 0) { kl_matrix[i, j] <- 0; next }
      th <- params; th[i] <- params[i] + grid_frac[j] * delta[i]
      kl_matrix[i, j] <- .kl_at(th)
    }
  }
  axis_min <- apply(kl_matrix[, c(1L, n_perturb), drop = FALSE], 1, function(v)
    if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE))

  # ---- 5. Joint (compensated) search over each face of the shell ----
  if (search == "auto") search <- if (n_par <= 10L) "joint" else "axis"
  face_min <- axis_min
  face_arg <- lapply(seq_len(n_par), function(i) NULL)
  penalty <- 1e6
  if (search == "joint" && n_par > 1L) {
    for (i in seq_len(n_par)) {
      best <- axis_min[i]; best_u <- NULL
      for (sgn in c(-1, 1)) {
        obj <- function(u_free) {
          u <- numeric(n_par); u[i] <- sgn; u[-i] <- u_free
          v <- .kl_at(params + u * delta)
          if (is.finite(v)) v else penalty
        }
        opt <- stats::optim(rep(0, n_par - 1L), obj, method = "L-BFGS-B",
                            lower = rep(-1, n_par - 1L), upper = rep(1, n_par - 1L),
                            control = list(maxit = as.integer(maxit)))
        if (is.finite(opt$value) && opt$value < penalty &&
            (is.na(best) || opt$value < best)) {
          best <- opt$value
          best_u <- numeric(n_par); best_u[i] <- sgn; best_u[-i] <- opt$par
        }
      }
      face_min[i] <- best
      if (!is.null(best_u)) face_arg[[i]] <- params + best_u * delta
    }
  }
  if (search == "axis" || n_par == 1L) {
    for (i in seq_len(n_par)) {
      j <- if (is.na(kl_matrix[i, 1L]) ||
               (!is.na(kl_matrix[i, n_perturb]) && kl_matrix[i, n_perturb] < kl_matrix[i, 1L]))
        n_perturb else 1L
      th <- params; th[i] <- params[i] + grid_frac[j] * delta[i]
      face_arg[[i]] <- th
    }
  }

  global_strength <- data.frame(
    parameter   = param_names,
    baseline    = unname(params),
    delta       = unname(delta),
    axis_min_kl = unname(axis_min),
    min_kl      = unname(face_min),
    stringsAsFactors = FALSE
  )
  global_strength$T_kl <- if (!is.null(T_obs)) T_obs * global_strength$min_kl else NA_real_
  global_strength$globally_identified <- global_strength$min_kl > kl_threshold
  global_strength$ident_quality <- ifelse(
    is.na(global_strength$min_kl), "not assessed",
    ifelse(global_strength$min_kl > 10 * kl_threshold, "strong",
           ifelse(global_strength$min_kl > kl_threshold, "adequate", "weak / unidentified")))

  gm_idx <- if (all(is.na(face_min))) NA_integer_ else which.min(face_min)
  global_min <- list(
    min_kl = if (is.na(gm_idx)) NA_real_ else face_min[gm_idx],
    face   = if (is.na(gm_idx)) NA_character_ else param_names[gm_idx],
    theta  = if (is.na(gm_idx)) NULL else face_arg[[gm_idx]]
  )

  n_strong   <- sum(global_strength$ident_quality == "strong")
  n_adequate <- sum(global_strength$ident_quality == "adequate")
  n_weak     <- sum(global_strength$ident_quality == "weak / unidentified")
  n_na       <- sum(global_strength$ident_quality == "not assessed")
  pass <- if (n_weak > 0) FALSE else if (n_na > 0) NA else TRUE
  weak_params_str <- if (n_weak > 0)
    paste(param_names[global_strength$ident_quality == "weak / unidentified"], collapse = ", ")
  else "none"
  na_params_str <- if (n_na > 0)
    paste(param_names[global_strength$ident_quality == "not assessed"], collapse = ", ")
  else "none"
  flat <- isTRUE(max(kl_matrix, na.rm = TRUE) < 1e-12)

  # ---- 6. Plots ----
  plots <- list()
  kl_profile_df <- do.call(rbind, lapply(seq_len(n_par), function(i) {
    data.frame(parameter = param_names[i], perturbation = grid_frac,
               side = ifelse(grid_frac < 0, "minus", "plus"),
               kl_divergence = kl_matrix[i, ], stringsAsFactors = FALSE)
  }))
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    floor_kl <- min(kl_threshold * 1e-3, 1e-8)
    pdf_ <- kl_profile_df[kl_profile_df$perturbation != 0 &
                            is.finite(kl_profile_df$kl_divergence), ]
    pdf_$kl_plot <- pmax(pdf_$kl_divergence, floor_kl)
    pdf_$line_id <- paste(pdf_$parameter, pdf_$side)
    n_hl <- min(6L, n_par)
    ord <- order(global_strength$min_kl, na.last = TRUE)   # weakest first
    hl <- param_names[ord][seq_len(n_hl)]
    pdf_$group <- factor(ifelse(pdf_$parameter %in% hl, pdf_$parameter, "other"),
                         levels = c(hl, "other"))
    col_map <- c(stats::setNames(dynhr_palette_vibrant[seq_len(n_hl)], hl),
                 other = unname(tol_vibrant["grey"]))
    p_prof <- ggplot2::ggplot(pdf_, ggplot2::aes(x = .data$perturbation, y = .data$kl_plot,
                                                 group = .data$line_id,
                                                 colour = .data$group)) +
      ggplot2::geom_line(linewidth = 0.7) +
      ggplot2::geom_point(size = 1.4) +
      ggplot2::geom_hline(yintercept = kl_threshold, linetype = "dashed",
                          colour = tol_vibrant[["red"]], linewidth = 0.5) +
      ggplot2::scale_y_log10() +
      ggplot2::scale_colour_manual(values = col_map, name = NULL,
                                   breaks = if (n_par > n_hl) c(hl, "other") else hl) +
      theme_dynhr_diagnostic() +
      ggplot2::labs(
        title = paste0("D24: KL divergence, one parameter moved at a time",
                       if (use_ramsey) " (Ramsey)" else ""),
        subtitle = sprintf(paste0("KL(theta0, theta) per observation; dashed = threshold %.3g.",
                                  "%s"), kl_threshold,
                           if (n_par > n_hl) sprintf(" Weakest %d coloured.", n_hl) else ""),
        x = if (length(unique(perturb_scale)) == 1L)
          sprintf("Parameter move in units of c (c = %g%% of |baseline|)", 100 * perturb_scale[[1]])
        else "Parameter move in units of c (c = perturb_scale x |baseline|)",
        y = "KL per observation (log scale)") +
      ggplot2::theme(legend.position = "bottom")
    plots$kl_profile <- .apply_meta(p_prof, meta)

    gs <- global_strength
    gs$status <- factor(ifelse(is.na(gs$min_kl), "not assessed",
                               ifelse(gs$globally_identified, "identified", "weak")),
                        levels = c("identified", "weak", "not assessed"))
    gs$min_kl_plot <- ifelse(is.na(gs$min_kl), floor_kl, pmax(gs$min_kl, floor_kl))
    gs$axis_plot <- pmax(gs$axis_min_kl, floor_kl)
    gs$param_label <- factor(gs$parameter, levels = gs$parameter[order(gs$min_kl_plot)])
    p_gs <- ggplot2::ggplot(gs, ggplot2::aes(x = .data$param_label)) +
      ggplot2::geom_segment(ggplot2::aes(xend = .data$param_label, y = floor_kl,
                                         yend = .data$min_kl_plot, colour = .data$status),
                            linewidth = 1.2) +
      ggplot2::geom_point(ggplot2::aes(y = .data$min_kl_plot, colour = .data$status),
                          size = 3) +
      ggplot2::geom_hline(yintercept = kl_threshold, linetype = "dashed",
                          colour = tol_vibrant[["red"]], linewidth = 0.5) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_log10() +
      ggplot2::scale_colour_manual(
        values = c(identified = tol_vibrant[["blue"]], weak = tol_vibrant[["red"]],
                   `not assessed` = tol_vibrant[["grey"]]),
        name = NULL, drop = TRUE) +
      theme_dynhr_diagnostic() +
      ggplot2::labs(
        title = "D24: Global identification strength",
        subtitle = if (search == "joint")
          sprintf("Min KL at distance c, others compensating (x: one at a time). Dashed = %.3g", kl_threshold)
        else sprintf("Min KL at distance c, one parameter at a time. Dashed = %.3g", kl_threshold),
        x = NULL, y = "Min KL per observation (log scale)")
    gs_floor <- gs[!is.na(gs$min_kl) & gs$min_kl <= floor_kl, , drop = FALSE]
    if (nrow(gs_floor) > 0)
      p_gs <- p_gs + ggplot2::geom_text(
        data = gs_floor, ggplot2::aes(y = .data$min_kl_plot),
        label = sprintf("KL <= %.0e (unidentified)", floor_kl),
        hjust = -0.08, size = 3.2, colour = "grey20")
    if (search == "joint")
      p_gs <- p_gs + ggplot2::geom_point(ggplot2::aes(y = .data$axis_plot), shape = 4,
                                         size = 2.5, colour = "grey30")
    plots$global_strength <- .apply_meta(p_gs, meta)
  }

  # ---- 7. Summary ----
  summary_text <- sprintf(
    "D24 Global identification via KL (%s search, c = %s): %d strong, %d adequate, %d weak%s. %s%s%s",
    search, paste(unique(sprintf("%g", perturb_scale)), collapse = "/"),
    n_strong, n_adequate, n_weak,
    if (n_na > 0) sprintf(", %d not assessed (%s)", n_na, na_params_str) else "",
    if (n_weak == 0 && n_na == 0) "All parameters globally identified."
    else if (n_weak > 0) sprintf("Globally weak: %s.", weak_params_str) else "",
    if (!is.na(global_min$min_kl))
      sprintf(" Smallest KL %.3g (moving %s)%s.", global_min$min_kl, global_min$face,
              if (!is.null(T_obs)) sprintf(", T*KL = %.3g at T = %d", T_obs * global_min$min_kl,
                                           as.integer(T_obs)) else "")
    else "",
    if (flat) " KL is zero everywhere: check that the re-solve closure uses theta." else "")

  llm_summary <- paste(c(
    sprintf("D24 | Global KL Identification (Qu-Tkachenko) | %s",
            if (isTRUE(pass)) "PASS" else if (identical(pass, FALSE)) "FAIL" else "INFO"),
    sprintf("  mode=%s spectral=%s search=%s n_par=%d n_freq=%d sigma_e=%s",
            if (use_ramsey) "ramsey" else "competitive", spectral, search, n_par,
            n_freq, sigma_source),
    sprintf("  strong=%d adequate=%d weak=%d not_assessed=%d bands=(>%.4g strong, >%.4g adequate) per obs",
            n_strong, n_adequate, n_weak, n_na, 10 * kl_threshold, kl_threshold),
    "  note: the three-tier KL band gates; T*KL is reported for power only. Qu-Tkachenko give a zero/non-zero criterion and no numeric cutoff, so kl_threshold=0.01 and the 10x band factor are package choices.",
    if (!is.null(T_obs) && is.finite(global_min$min_kl))
      sprintf("  power: T*KL=%.3g nats at T=%d (~1 nat = barely detectable in principle; ~4 nats ~ 80%% power at the 5%% level for one free direction)",
              T_obs * global_min$min_kl, as.integer(T_obs)),
    sprintf("  weak_global_params: %s", weak_params_str),
    sprintf("  global_min_kl=%.3g face=%s%s kl_self=%.1e",
            global_min$min_kl, global_min$face,
            if (!is.null(T_obs)) sprintf(" T_kl=%.3g", T_obs * global_min$min_kl) else "",
            kl_self),
    sprintf("  action: %s",
            if (n_weak > 0) sprintf(paste("Parameters %s are nearly observationally",
                                          "equivalent to a point at distance c; inspect",
                                          "result$global_min$theta."), weak_params_str)
            else "none")
  ), collapse = "\n")

  .make_result(
    result = list(
      kl_matrix = kl_matrix,
      kl_profile = kl_profile_df,
      global_strength = global_strength,
      global_min = global_min,
      face_argmin = stats::setNames(face_arg, param_names),
      frequency_grid = freq_grid,
      baseline_spectral = S0,
      kl_self = kl_self,
      kl_threshold = kl_threshold,
      T_obs = T_obs,
      use_ramsey = use_ramsey,
      spectral = spectral,
      search = search,
      Sigma_e = Sigma_e,
      sigma_e_source = sigma_source
    ),
    pass = pass,
    plots = plots,
    summary = summary_text,
    llm_summary = llm_summary
  )
}
