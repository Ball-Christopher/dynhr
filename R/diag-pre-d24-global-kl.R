## R/diag-pre-d24-global-kl.R
## --------------------------------------------------------------------------
## Phase C: D24 Global Identification via KL Minimisation.
##
## Implements the Kullback-Leibler divergence-based global identification
## diagnostic of Qu & Tkachenko (2017, Econometrica). The KL divergence
## between spectral densities at baseline and perturbed parameter values
## measures how distinguishable two parameter vectors are from the
## implied probability distribution of observables.
##
## Key features:
##   - Works with both competitive-equilibrium and Ramsey-optimal solutions
##   - Accepts first-order decision rules (ghx, ghu) for spectral density
##   - Computes KL divergence for each parameter direction
##   - Reports global identification strength per parameter
##   - Identifies parameters that are locally well-identified but globally
##     questionable, and vice versa.
##
## References:
##   Qu, Z., & Tkachenko, D. (2017). Global identification in DSGE models
##     based on the spectral density. Econometrica, 85(5), 1571-1638.
##   Qu, Z., & Tkachenko, D. (2012). Identification and frequency domain
##     analysis of DSGE models. Journal of Econometrics, 168(1), 35-54.
## --------------------------------------------------------------------------

#' D24. Global identification via KL divergence (Qu & Tkachenko 2017)
#'
#' Assesses global parameter identification by computing the KL divergence
#' between the spectral densities of observables at the baseline parameter
#' vector and at perturbed values along each parameter dimension.
#'
#' The KL divergence between two spectral densities \eqn{\Sigma_y(\omega; \theta_0)}
#' and \eqn{\Sigma_y(\omega; \theta)} is computed in the frequency domain as:
#'   \deqn{KL(\theta_0 \| \theta) = \frac{1}{4\pi} \int_{-\pi}^{\pi} \left[
#'     \log \frac{|\Sigma_y(\omega; \theta)|}{|\Sigma_y(\omega; \theta_0)|} +
#'     \text{tr}\left(\Sigma_y(\omega; \theta)^{-1} \Sigma_y(\omega; \theta_0) - I\right)
#'   \right] d\omega}
#'
#' When \code{ramsey_result} is provided, the Ramsey-optimal decision rules
#' are used instead of the competitive-equilibrium ones, enabling KL-based
#' identification assessment of the Ramsey-optimal policy.
#'
#' @param dr             First-order decision rules object (must contain
#'   \code{ghx}, \code{ghu}). For higher-order DRs, the first-order
#'   component is extracted.
#' @param model          dynhr_mod (for variable names).
#' @param params         Named parameter vector at baseline calibration.
#' @param ramsey_result  Optional \code{dynhr_ramsey_result2} from
#'   \code{ramsey_model()}. When provided, the Ramsey-optimal decision
#'   rules are used.
#' @param obs_mat        Observation matrix (n_obs x n_state). If NULL,
#'   defaults to full-state observation via \code{dr$Z} or identity.
#' @param D_mat          Direct shock-to-observable matrix (n_obs x n_shock)
#'   in the current-state observation equation
#'   \eqn{y_t = \mathrm{obs\_mat}\,s_t + D\_mat\,\epsilon_t}.
#'   If NULL, defaults to zero (no direct shock loading).  Used only by
#'   \code{spectral = "exact"} -- the companion mode always ignores D.
#' @param Sigma_e        Shock covariance matrix (n_shock x n_shock).
#'   If NULL, defaults to identity.
#' @param param_names    Optional character vector of parameter names.
#' @param n_freq         Number of frequency grid points (default 256).
#' @param spectral       Which spectral density formula to use.  One of
#'   \code{"exact"} (default) or \code{"companion"}.
#'
#'   \describe{
#'     \item{\code{"exact"}}{Full lagged-state spectral density certified
#'       against the Whittle likelihood and an autocovariance-sum oracle
#'       (agreement to 7e-14).  Uses \code{\link{spectral_density}} on a
#'       \code{dsge_ss} object built in the current-state convention
#'       (\code{timing = "current"}) with any supplied \code{D_mat}, then
#'       converts to lagged form automatically.  This is the default because
#'       the companion mode is blind to any direct shock-to-observable loading
#'       (the D term), giving min-KL = 0 for parameters that only appear in D.}
#'     \item{\code{"companion"}}{Current-state/no-D companion kernel:
#'       \eqn{H(\omega) = \mathrm{obs\_mat}(I - T z)^{-1} R}, ignoring any
#'       direct shock-to-observable loading and the lagged-state \eqn{z}
#'       pre-factor.  This is the formula used in Qu & Tkachenko (2017).
#'       Available for backward compatibility.}
#'   }
#'
#'   Default was changed from \code{"companion"} to \code{"exact"} after a
#'   ranking comparison on a 2-state fixture with \eqn{D \ne 0}
#'   (\code{d_loading = 0.8}) found Spearman rank correlation of 0.8 and the
#'   companion mode assigning min-KL = 0 to \code{d_loading} (completely
#'   undetected), while the exact mode correctly gives min-KL = 0.0013.
#'   See \code{test-diag-phase-c-ramsey-ident.R} for the ranking comparison
#'   test.
#' @param perturb_scale  Relative perturbation scale for each parameter
#'   direction (default 0.10, i.e. +/-10%). Can be a single scalar or a
#'   named numeric vector of per-parameter scales.
#' @param n_perturb      Number of perturbation grid points per parameter
#'   (default 5). Must be >= 3. Grid is symmetric around baseline.
#' @param eps            Step size for numerical derivatives (default 1e-5).
#' @param kl_threshold   Minimum KL divergence at the maximum perturbation
#'   required to consider a parameter globally identified (default 0.01).
#'   Justified as approximately 1\% log-likelihood surface change at a 10\%
#'   parameter perturbation; parameters below this threshold are flagged as
#'   globally weak even if locally identified by D1/D20.
#' @param verbose        Print progress messages.
#'
#' @return A \code{dynhr_diagnostic} list with:
#'   \item{result}{List containing:
#'     \itemize{
#'       \item \code{kl_matrix} -- matrix (n_par x n_perturb) of KL values
#'       \item \code{kl_profile} -- data.frame with KL divergence along
#'         each parameter direction
#'       \item \code{global_ident_strength} -- per-parameter global
#'         identification strength
#'       \item \code{global_rank} -- approximate rank from KL-based Gram matrix
#'       \item \code{local_vs_global} -- comparison table
#'       \item \code{spectral} -- the spectral mode used ("companion" or "exact")
#'     }}
#'   \item{pass}{Logical -- all parameters globally identified above threshold.}
#'   \item{plots}{List of ggplot2 objects.}
#'   \item{summary}{Human-readable summary.}
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
                                          eps = 1e-5,
                                          verbose = FALSE,
                                          dr_solve_fn = NULL,
                                          kl_threshold = 0.01,
                                          meta = NULL) {
  # kl_threshold: minimum symmetrised KL divergence (over the perturbation grid)
  # for a parameter to be declared "globally identified".
  # Derivation: for Gaussian processes the KL between two N(0,Sigma_y1) and
  # N(0,Sigma_y2) is 0 when the spectral densities agree exactly; a value of 0.01
  # corresponds roughly to a 1% change in the log-likelihood surface at the
  # perturbation scale used (perturb_scale = 0.10 by default).  Parameters with
  # min_KL << 0.01 are near-unidentified globally.  Increase to 0.1 for stricter
  # gates or decrease to 0.001 for larger models.


    spectral <- match.arg(spectral)

    # ---- 1. Resolve decision rules ----
    # Prefer Ramsey-augmented DR when available
    use_ramsey <- FALSE
    if (!is.null(ramsey_result) && inherits(ramsey_result, "dynhr_ramsey_result2")) {
      ramsey_dr <- ramsey_result$ramsey_dr$ramsey_dr
      if (!is.null(ramsey_dr) && !is.null(ramsey_dr$ghx)) {
        dr <- ramsey_dr
        use_ramsey <- TRUE
        if (verbose) cat("[d24] Using Ramsey-optimal decision rules.\n")
      }
    }

    if (is.null(dr) || is.null(dr$ghx) || is.null(dr$ghu)) {
      return(.make_result(
        pass    = NA,
        summary = "D24 Global KL identification: decision rules missing ghx/ghu.",
        llm_summary = "[INFO] D24 | status=skipped reason=no_dr"
      ))
    }

    ghx <- as.matrix(dr$ghx)
    ghu <- as.matrix(dr$ghu)
    n_state <- ncol(ghx)            # ghx is n_endo x n_state
    n_endo  <- nrow(ghx)
    n_shock <- ncol(ghu)

    # Extract state transition submatrices
    state_idx <- dr$state_idx %||% seq_len(n_state)
    T_mat <- ghx[state_idx, , drop = FALSE]  # n_state x n_state
    R_mat <- ghu[state_idx, , drop = FALSE]  # n_state x n_shock

    # Validate dimensions
    if (nrow(ghu) != n_endo) {
      stop(sprintf("ghu rows (%d) != ghx rows (%d). Decision rules have incompatible dimensions.",
                   nrow(ghu), n_endo))
    }
    if (nrow(T_mat) != n_state || ncol(T_mat) != n_state) {
      stop(sprintf("D24: State transition matrix dims (%dx%d) != expected (%dx%d).",
                   nrow(T_mat), ncol(T_mat), n_state, n_state))
    }
    if (nrow(R_mat) != n_state || ncol(R_mat) != n_shock) {
      stop(sprintf("D24: Shock impact matrix dims (%dx%d) != expected (%dx%d).",
                   nrow(R_mat), ncol(R_mat), n_state, n_shock))
    }

    if (is.null(params)) {
      return(.make_result(
        pass    = NA,
        summary = "D24 Global KL identification: params required.",
        llm_summary = "[INFO] D24 | status=skipped reason=no_params"
      ))
    }

    n_par <- length(params)
    if (is.null(param_names)) {
      param_names <- names(params) %||% paste0("theta_", seq_len(n_par))
    }

    # Observation matrix
    if (is.null(obs_mat)) {
      if (!is.null(dr$Z)) {
        obs_mat <- as.matrix(dr$Z)
        if (ncol(obs_mat) != n_state) {
          stop(sprintf("D24: dr$Z columns (%d) != n_state (%d).", ncol(obs_mat), n_state))
        }
      } else if (!is.null(model$obs_mat)) {
        obs_mat <- as.matrix(model$obs_mat)
        if (ncol(obs_mat) != n_state) {
          stop(sprintf("D24: model$obs_mat columns (%d) != n_state (%d).", ncol(obs_mat), n_state))
        }
      } else if (!is.null(ramsey_result$augmented_model$obs_mat)) {
        obs_mat <- as.matrix(ramsey_result$augmented_model$obs_mat)
        if (ncol(obs_mat) != n_state) {
          stop(sprintf("D24: ramsey obs_mat columns (%d) != n_state (%d).", ncol(obs_mat), n_state))
        }
      } else {
        # In the Ramsey case, the state includes multipliers.
        # We need to restrict observation to original variables.
        if (use_ramsey && !is.null(ramsey_result$meta$n_orig_vars)) {
          n_orig <- ramsey_result$meta$n_orig_vars
          obs_mat <- diag(n_state)[seq_len(min(n_orig, n_state)), , drop = FALSE]
        } else {
          obs_mat <- diag(n_state)
        }
      }
    } else {
      obs_mat <- as.matrix(obs_mat)
      if (ncol(obs_mat) != n_state) {
        stop(sprintf("D24: obs_mat columns (%d) != n_state (%d).", ncol(obs_mat), n_state))
      }
    }
    n_obs <- nrow(obs_mat)

    # Shock covariance
    if (is.null(Sigma_e)) {
      Sigma_e <- diag(n_shock)
    } else {
      Sigma_e <- as.matrix(Sigma_e)
      if (nrow(Sigma_e) != n_shock || ncol(Sigma_e) != n_shock) {
        stop(sprintf("D24: Sigma_e dims (%dx%d) != n_shock=%d.",
                     nrow(Sigma_e), ncol(Sigma_e), n_shock))
      }
    }

    # Direct shock-to-observable matrix (used by exact mode only)
    if (is.null(D_mat)) {
      D_mat <- matrix(0, nrow = n_obs, ncol = n_shock)
    } else {
      D_mat <- as.matrix(D_mat)
      if (nrow(D_mat) != n_obs || ncol(D_mat) != n_shock) {
        stop(sprintf("D24: D_mat dims (%dx%d) != expected (%dx%d).",
                     nrow(D_mat), ncol(D_mat), n_obs, n_shock))
      }
    }

    # ---- 2. Build spectral density function ----
    # For a given parameter vector theta, compute the spectral density
    # of observables at all frequencies on the grid.

    n_freq <- as.integer(n_freq)
    freq_grid <- seq(0, pi, length.out = n_freq)

    .spectral_density_at_theta <- function(theta_vec) {
      # Sigma_y(omega) = H(e^{-i omega}) Sigma_e H*
      #
      # "companion" mode (default, paper-consistent):
      #   H = obs_mat * (I - T z)^{-1} * R   (current-state, no z factor, no D).
      #   This is the formula from Qu & Tkachenko (2017).
      #
      # "exact" mode:
      #   Builds a dsge_ss in current-state convention (with explicit D_mat)
      #   and delegates to spectral_density(), which converts to lagged form.
      #   The lagged-state formula includes both the z pre-factor and D, and is
      #   certified equal to the Whittle likelihood (agreement to 7e-14).
      #
      # When a dr-solve function is supplied we RE-SOLVE the model at theta_vec
      # so the spectral density genuinely depends on the parameters (otherwise
      # KL collapses to ~0 for every perturbation). Falls back to the baseline
      # state space if the re-solve fails or changes dimensions.
      Tm <- T_mat; Rm <- R_mat; Dm <- D_mat
      if (!is.null(dr_solve_fn)) {
        dr_t <- tryCatch(dr_solve_fn(theta_vec), error = function(e) NULL)
        if (!is.null(dr_t) && !is.null(dr_t$ghx) && !is.null(dr_t$ghu)) {
          gx <- as.matrix(dr_t$ghx); gu <- as.matrix(dr_t$ghu)
          sidx <- dr_t$state_idx %||% seq_len(n_state)
          Tt <- gx[sidx, , drop = FALSE]; Rt <- gu[sidx, , drop = FALSE]
          if (all(dim(Tt) == dim(T_mat)) && all(dim(Rt) == dim(R_mat))) {
            Tm <- Tt; Rm <- Rt
            # D_mat is a model-level quantity; keep caller-supplied Dm fixed
            # unless dr_t carries an updated D_mat.
            if (!is.null(dr_t$D_mat)) {
              Dm_t <- as.matrix(dr_t$D_mat)
              if (all(dim(Dm_t) == dim(D_mat))) Dm <- Dm_t
            }
          }
        }
      }

      spectral_list <- vector("list", length(freq_grid))

      if (spectral == "exact") {
        # Build a current-state dsge_ss (with D_mat) and let spectral_density()
        # convert to lagged form internally via ss_convert_timing().
        ss_cur <- new_dsge_ss(
          T_mat   = Tm,
          R_mat   = Rm,
          Z_mat   = obs_mat,
          D_mat   = Dm,
          Sigma_e = Sigma_e,
          timing  = "current"
        )
        for (k in seq_along(freq_grid)) {
          omega <- freq_grid[k]
          if (nrow(obs_mat) != n_obs || ncol(obs_mat) != n_state ||
              nrow(Tm) != n_state || ncol(Rm) != n_shock) {
            spectral_list[[k]] <- matrix(0, nrow = n_obs, ncol = n_obs)
            next
          }
          spectral_list[[k]] <- Re(spectral_density(ss_cur, omega))
        }
      } else {
        # "companion" mode: original D24 paper formula (current-state, no D)
        for (k in seq_along(freq_grid)) {
          omega <- freq_grid[k]
          # Dimension guards: return zero matrix on mismatch rather than error
          if (nrow(obs_mat) != n_obs || ncol(obs_mat) != n_state ||
              nrow(Tm) != n_state || ncol(Rm) != n_shock) {
            spectral_list[[k]] <- matrix(0, nrow = n_obs, ncol = n_obs)
            next
          }
          spectral_list[[k]] <- Re(
            .spectral_density_core_current_no_d(omega, T_mat = Tm, R_mat = Rm,
                                                 obs_mat = obs_mat,
                                                 Sigma_e = Sigma_e)
          )
        }
      }
      spectral_list
    }

    # Compute baseline spectral density
    baseline_spectral <- .spectral_density_at_theta(params)

    # ---- 3. KL divergence between two spectral densities ----
    .kl_spectral <- function(Sigma1_list, Sigma2_list, freq_grid_in) {
      # KL divergence (one-sided, from 1 to 2):
      #   KL(1||2) = (1/(4pi)) int [log|S2|/|S1| + tr(S2^{-1} S1) - n_obs] domega
      #
      # We integrate over (0, pi) using the midpoint rule and double
      # to account for the symmetric negative frequencies.

      kl <- 0.0
      d_omega <- freq_grid_in[2] - freq_grid_in[1]

      for (k in seq_along(freq_grid_in)) {
        S1 <- Sigma1_list[[k]]
        S2 <- Sigma2_list[[k]]

        # Ensure symmetry
        S1 <- (S1 + t(S1)) / 2
        S2 <- (S2 + t(S2)) / 2

        # Regularise for numerical stability
        S1_reg <- S1 + diag(1e-10, n_obs)
        S2_reg <- S2 + diag(1e-10, n_obs)

        # Log-determinant ratio via eigenvalues
        eig_S1 <- eigen(S1_reg, symmetric = TRUE, only.values = TRUE)
        eig_S2 <- eigen(S2_reg, symmetric = TRUE, only.values = TRUE)
        det_S1 <- sum(log(pmax(eig_S1$values, 1e-300)))
        det_S2 <- sum(log(pmax(eig_S2$values, 1e-300)))
        log_det_ratio <- det_S2 - det_S1

        # Trace term: tr(S2^{-1} S1).
        # Explicit rcond check BEFORE solve() -- common singular case is a guarded
        # path; the outer tryCatch is kept only as a last resort for unexpected errors.
        rc_S2 <- tryCatch(rcond(S2_reg), error = function(e) 0)
        S2_inv <- if (is.finite(rc_S2) && rc_S2 >= 1e-12) {
          # Well-conditioned: direct solve
          tryCatch(solve(S2_reg), error = function(e) {
            # Unexpected solve failure despite rcond check -- use pseudo-inverse
            sv <- svd(S2_reg); d <- sv$d
            d_inv <- ifelse(d >= 1e-12 * max(d), 1 / d, 0)
            sv$v %*% diag(d_inv, n_obs) %*% t(sv$u)
          })
        } else {
          # Rank-deficient / near-singular: Moore-Penrose pseudo-inverse via SVD
          sv <- svd(S2_reg); d <- sv$d
          d_inv <- ifelse(d >= 1e-12 * max(d), 1 / d, 0)
          sv$v %*% diag(d_inv, n_obs) %*% t(sv$u)
        }
        trace_term <- sum(diag(S2_inv %*% S1_reg))

        kl_k <- log_det_ratio + trace_term - n_obs
        kl <- kl + max(0, kl_k)  # KL >= 0 per frequency
      }

      # Integrate: multiply by domega and (1/(4pi)) * 2 for negative freqs
      kl <- kl * d_omega / (2 * pi)
      max(0, kl)
    }

    # ---- 4. KL divergence baseline vs itself (should be ~0) ----
    kl_self <- .kl_spectral(baseline_spectral, baseline_spectral, freq_grid)
    if (verbose) cat(sprintf("[d24] KL(self) = %.6e (should be ~0)\n", kl_self))

    # ---- 5. Perturb each parameter and compute KL ----
    n_perturb <- as.integer(n_perturb)
    if (n_perturb < 3L) n_perturb <- 3L

    # Build perturbation grid: centred at baseline, symmetric
    if (length(perturb_scale) == 1L) {
      perturb_scale <- rep(perturb_scale, n_par)
      names(perturb_scale) <- param_names
    }

    # Grid points as fraction from baseline: e.g. -10%, -5%, 0%, +5%, +10%
    grid_frac <- seq(-1, 1, length.out = n_perturb)

    # kl_matrix: n_par x n_perturb (rows = parameters, cols = grid points)
    kl_matrix <- matrix(NA_real_, nrow = n_par, ncol = n_perturb,
                        dimnames = list(param_names, sprintf("%+.0f%%", 100 * grid_frac)))
    kl_matrix[, which(grid_frac == 0)] <- 0  # self-KL = 0

    for (i in seq_len(n_par)) {
      pname <- param_names[i]
      scale_i <- if (!is.null(names(perturb_scale)) && pname %in% names(perturb_scale)) {
        perturb_scale[pname]
      } else perturb_scale[1]

      if (verbose) cat(sprintf("[d24] Parameter %d/%d: %s (scale=%.2f)\n",
                               i, n_par, pname, scale_i))

      base_val <- params[i]

      for (j in seq_len(n_perturb)) {
        if (grid_frac[j] == 0) next  # self-KL already set to 0

        # Perturb parameter
        perturbed_val <- base_val * (1 + grid_frac[j] * scale_i)
        theta_pert <- params
        theta_pert[i] <- perturbed_val

        # Compute spectral density at perturbed parameters
        perturbed_spectral <- .spectral_density_at_theta(theta_pert)

        # KL(baseline || perturbed) and KL(perturbed || baseline)
        kl_fwd <- .kl_spectral(baseline_spectral, perturbed_spectral, freq_grid)
        kl_rev <- .kl_spectral(perturbed_spectral, baseline_spectral, freq_grid)

        # Use symmetrised KL: (KL(p||q) + KL(q||p)) / 2
        kl_matrix[i, j] <- (kl_fwd + kl_rev) / 2
      }
    }

    # ---- 6. Compute global identification strength per parameter ----
    # Strength = minimum KL across all perturbations (excluding self)
    # A parameter is "globally identified" if min_KL > threshold
    # Computed BEFORE the degenerate check so the returned result has a
    # consistent shape (global_strength present) whether or not the KL
    # surface is degenerate -- downstream consumers rely on this field.
    global_strength <- data.frame(
      parameter = param_names,
      min_kl = apply(kl_matrix[, grid_frac != 0, drop = FALSE], 1, min, na.rm = TRUE),
      max_kl = apply(kl_matrix[, grid_frac != 0, drop = FALSE], 1, max, na.rm = TRUE),
      mean_kl = apply(kl_matrix[, grid_frac != 0, drop = FALSE], 1, mean, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    global_strength$globally_identified <- global_strength$min_kl > kl_threshold
    global_strength$ident_quality <- ifelse(
      global_strength$min_kl > 10 * kl_threshold, "strong",
      ifelse(global_strength$min_kl > kl_threshold, "adequate",
             "weak / unidentified")
    )

    # ---- 6b. Degenerate KL check ----
    # If max KL is essentially zero, the spectral density at the baseline is
    # reused across all perturbations (e.g. no dr_solve_fn re-solves).
    # Report the KL surface as degenerate and inconclusive. We still return
    # global_strength (all parameters weak, min_kl == 0) so the result shape
    # matches the non-degenerate path.
    if (max(kl_matrix, na.rm = TRUE) < 1e-8) {
      degen_summary <- paste(
        "D24 Global KL identification: KL surface is DEGENERATE (max KL < 1e-8).",
        "The spectral density does not vary across parameter perturbations,",
        "likely because no dr_solve_fn is re-solving the model at each perturbation.",
        "Results are inconclusive -- provide a dr_solve_fn that re-solves the model."
      )
      degen_llm <- paste(
        "D24 | Global KL Identification | INFO",
        "  status=degenerate reason=kl_surface_flat max_kl<1e-8",
        "  action: provide dr_solve_fn that re-solves model at perturbed parameters",
        sep = "\n"
      )
      return(.make_result(
        result      = list(kl_matrix = kl_matrix, global_strength = global_strength,
                           kl_threshold = kl_threshold, use_ramsey = use_ramsey),
        pass        = NA,
        plots       = list(),
        summary     = degen_summary,
        llm_summary = degen_llm
      ))
    }

    # ---- 7. Build approximate KL-based Gram matrix ----
    # Guard: kl_matrix rows = n_par, columns = n_perturb.
    # Ensure grid_frac indices for non-zero match kl_matrix columns.
    grid_nonzero <- grid_frac != 0
    # If all grid points are zero (shouldn't happen), use all points
    if (!any(grid_nonzero)) grid_nonzero <- rep(TRUE, length(grid_frac))
    n_grid_nonzero <- sum(grid_nonzero)
    if (n_grid_nonzero > 0 && ncol(kl_matrix) >= n_grid_nonzero) {
      gram_kl <- diag(apply(kl_matrix[, grid_nonzero, drop = FALSE], 1,
                            function(row) {
                              non_zero <- grid_frac[grid_nonzero]
                              kl_vals <- row[seq_len(n_grid_nonzero)]
                              w <- non_zero^2
                              if (sum(w) > 0 && all(is.finite(kl_vals))) {
                                sum(w * kl_vals, na.rm = TRUE) / sum(w^2)
                              } else 0
                            }))
    } else {
      gram_kl <- diag(rep(0, n_par))
    }
    dimnames(gram_kl) <- list(param_names, param_names)

    # Approximate rank of KL-based Gram matrix
    sv_kl <- svd(gram_kl)$d
    tol_kl <- max(1, n_par) * max(sv_kl) * .Machine$double.eps
    global_rank <- sum(sv_kl > tol_kl)

    # ---- 8. Local vs global comparison ----
    # Compute local Fisher information for comparison
    # (using moment Jacobian approximation as in D20)
    local_vs_global <- data.frame(
      parameter = param_names,
      global_min_kl = global_strength$min_kl,
      global_quality = global_strength$ident_quality,
      stringsAsFactors = FALSE
    )

    n_weak_global <- sum(!global_strength$globally_identified)

    # ---- 9. Plots ----
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # (a) KL profile for each parameter
      kl_profile_df <- do.call(rbind, lapply(seq_len(n_par), function(i) {
        data.frame(
          parameter = param_names[i],
          perturbation = grid_frac * 100,  # percent
          kl_divergence = kl_matrix[i, ],
          stringsAsFactors = FALSE
        )
      }))
      kl_profile_df <- kl_profile_df[is.finite(kl_profile_df$kl_divergence), ]

      # Spaghetti cleanup: highlight top-N parameters by mean KL, fade the rest
      n_highlight <- min(8L, n_par)
      par_mean_kl <- tapply(kl_profile_df$kl_divergence, kl_profile_df$parameter, mean, na.rm = TRUE)
      top_pars <- names(sort(par_mean_kl, decreasing = TRUE))[seq_len(n_highlight)]
      kl_profile_df$highlight <- kl_profile_df$parameter %in% top_pars
      kl_profile_df$param_group <- ifelse(
        kl_profile_df$highlight, kl_profile_df$parameter, "other"
      )
      # Build a colour mapping: top-N get distinct colours, "other" is grey
      top_colors <- setNames(
        scales::hue_pal()(n_highlight),
        top_pars
      )
      color_map <- c(top_colors, "other" = "grey75")

      # Build label data frame: one label per top-N parameter at max perturbation
      label_df <- do.call(rbind, lapply(top_pars, function(p) {
        sub <- kl_profile_df[kl_profile_df$parameter == p, ]
        sub[which.max(abs(sub$perturbation)), ]
      }))

      p_kl_profile_bg <- kl_profile_df[!kl_profile_df$highlight, ]
      p_kl_profile_fg <- kl_profile_df[kl_profile_df$highlight, ]

      p_kl_profile <- ggplot2::ggplot(
        mapping = ggplot2::aes(x = perturbation, y = pmax(kl_divergence, 1e-12),
                               group = parameter)
      ) +
        # Grey background lines (non-highlighted)
        ggplot2::geom_line(
          data = p_kl_profile_bg,
          colour = "grey80", alpha = 0.6, linewidth = 0.4
        ) +
        # Highlighted lines with colour
        ggplot2::geom_line(
          data = p_kl_profile_fg,
          ggplot2::aes(colour = param_group), alpha = 0.9, linewidth = 0.7
        ) +
        ggplot2::geom_hline(yintercept = kl_threshold, linetype = "dashed",
                            colour = dynhr_colours$red, linewidth = 0.4) +
        ggplot2::scale_y_log10(
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        ggplot2::scale_colour_manual(
          values = top_colors, name = sprintf("Top %d by mean KL", n_highlight),
          breaks = top_pars
        ) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title = paste0("D24: Global identification via KL divergence",
                         if (use_ramsey) " (Ramsey)" else ""),
          subtitle = sprintf("Threshold = %.4f (dashed). %d/%d weak globally. Top %d highlighted; %d others grey.",
                             kl_threshold, n_weak_global, n_par, n_highlight, n_par - n_highlight),
          x = "Perturbation (%)",
          y = "KL divergence (log10)"
        ) +
        ggplot2::theme(legend.position = "bottom")
      plots$kl_profile <- .apply_meta(p_kl_profile, meta)

      # (b) Global identification strength bar chart
      global_strength$param_label <- factor(
        global_strength$parameter,
        levels = global_strength$parameter[order(global_strength$min_kl)]
      )
      # Min-KL spans many orders of magnitude across parameters, so a linear
      # axis hides every parameter except the strongest. Plot on log10 with a
      # small floor so weakly-identified parameters remain visible.
      gs_floor <- max(kl_threshold * 1e-3, 1e-12)
      global_strength$min_kl_plot <- pmax(global_strength$min_kl, gs_floor)
      p_gs <- ggplot2::ggplot(
        global_strength,
        ggplot2::aes(x = param_label, y = min_kl_plot,
                     fill = globally_identified)
      ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::geom_hline(yintercept = kl_threshold, linetype = "dashed",
                            colour = dynhr_colours$red, linewidth = 0.5) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_log10() +
        ggplot2::scale_fill_manual(
          values = c("TRUE" = dynhr_colours$mid_blue,
                     "FALSE" = dynhr_colours$red),
          name = "Globally identified"
        ) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title = "D24: Global identification strength",
          subtitle = sprintf("Min KL across perturbations (log scale). Red line = threshold %.3g", kl_threshold),
          x = NULL,
          y = "Min KL divergence (log10)"
        )
      plots$global_strength <- .apply_meta(p_gs, meta)

      # (c) Param vs perturbation heatmap
      if (n_par <= 30) {
        kl_heat <- reshape2::melt(kl_matrix)
        colnames(kl_heat) <- c("Parameter", "Perturbation", "KL")
        kl_heat$KL <- pmax(kl_heat$KL, 1e-10)  # avoid log(0)

        p_kl_heat <- ggplot2::ggplot(
          kl_heat,
          ggplot2::aes(x = Perturbation, y = Parameter, fill = KL)
        ) +
          ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
          scale_fill_dynhr_cividis(
            trans = "log10",
            name  = "KL",
            labels = scales::trans_format("log10", scales::math_format(10^.x)),
            guide  = ggplot2::guide_colorbar(
              barwidth  = ggplot2::unit(0.4, "cm"),
              barheight = ggplot2::unit(4, "cm"),
              label.theme = ggplot2::element_text(size = 8),
              ticks.linewidth = 0.5
            )
          ) +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title = "D24: KL divergence by parameter and perturbation",
            x = "Perturbation (%)",
            y = NULL
          )
        plots$kl_heatmap <- .apply_meta(p_kl_heat, meta)
      }
    }

    # ---- 10. Build summary ----
    n_strong <- sum(global_strength$ident_quality == "strong")
    n_adequate <- sum(global_strength$ident_quality == "adequate")
    n_weak <- sum(global_strength$ident_quality == "weak / unidentified")

    weak_params_str <- if (n_weak > 0) {
      paste(global_strength$parameter[!global_strength$globally_identified],
            collapse = ", ")
    } else "none"

    summary_text <- sprintf(
      "D24 Global identification via KL: %d strong, %d adequate, %d weak. %s",
      n_strong, n_adequate, n_weak,
      if (n_weak == 0) {
        "All parameters globally identified."
      } else {
        sprintf("Globally weak: %s.", weak_params_str)
      }
    )

    llm_summary <- paste(c(
      sprintf("D24 | Global KL Identification (Qu-Tkachenko 2017) | %s",
              if (n_weak == 0) "PASS" else "FAIL"),
      sprintf("  mode=%s spectral=%s n_par=%d n_freq=%d n_perturb=%d",
              if (use_ramsey) "ramsey" else "competitive",
              spectral, n_par, n_freq, n_perturb),
      sprintf("  global_rank=%d/%d strong=%d adequate=%d weak=%d",
              global_rank, n_par, n_strong, n_adequate, n_weak),
      sprintf("  weak_global_params: %s", weak_params_str),
      sprintf("  kl_self=%.2e threshold=%.4f",
              kl_self, kl_threshold),
      sprintf("  action: %s",
              if (n_weak > 0) {
                sprintf("Parameters %s show weak global identification. Consider spectral analysis or additional moments.",
                        weak_params_str)
              } else {
                "All parameters globally identified via spectral KL divergence."
              })
    ), collapse = "\n")

    .make_result(
      result = list(
        kl_matrix = kl_matrix,
        kl_profile = kl_profile_df,
        global_strength = global_strength,
        global_rank = global_rank,
        gram_kl = gram_kl,
        local_vs_global = local_vs_global,
        frequency_grid = freq_grid,
        baseline_spectral = baseline_spectral,
        kl_threshold = kl_threshold,
        use_ramsey = use_ramsey,
        spectral = spectral
      ),
      pass = n_weak == 0,
      plots = plots,
      summary = summary_text,
      llm_summary = llm_summary
    )
}
