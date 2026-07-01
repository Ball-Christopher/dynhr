## R/diag-pre-d23-spectral.R
## --------------------------------------------------------------------------
## Phase A implementation: D23 Spectral identification check.
##
## Implements the Qu & Tkachenko (2012, QE) spectral identification test.
##
## Algorithm:
##   1. Given first-order solution (ghx, ghu, C, D, Sigma_u), compute the
##      transfer function H(e^{iω}) = C (I - ghx e^{-iω})^{-1} ghu + D
##   2. Compute the spectral density Σ_y(ω) = H(e^{iω}) Σ_u H(e^{iω})^*
##      on a frequency grid ω ∈ (0, π)
##   3. Numerically differentiate Σ_y(ω) w.r.t. θ (re-use .numerical_jacobian)
##   4. Assemble Gram matrix G(θ) by integrating the Jacobian of the spectral
##      density over the frequency domain
##   5. Check rank via SVD; report frequency bands driving identification
## --------------------------------------------------------------------------

#' D23. Spectral identification check (Qu & Tkachenko 2012)
#'
#' Computes the Gram matrix of the spectral density derivatives with respect
#' to model parameters. The spectral density is derived from the first-order
#' state-space representation:
#'   s_t = ghx * s_{t-1} + ghu * ε_t
#'   y_t = C * s_t + D * ε_t   (observation equation)
#'
#' Identification is assessed by checking whether the Gram matrix is
#' full rank. The method reports which frequency bands contribute most
#' to parameter identification, which helps guide data filtering choices.
#'
#' \strong{REQUIRES} a \code{model_solve_fn} that returns a decision-rules
#' list (containing \code{ghx}).  When this is not supplied, the Gram matrix
#' cannot be computed from spectral derivatives and the diagnostic returns
#' \code{rank = 0} (all parameters flagged weak) with a warning.  This is
#' a data-availability failure, not a model identification failure.
#'
#' @param dr Decision rules object (must contain ghx, ghu, and optionally
#'           C/obs_mat and D). If C and D are not provided, defaults to
#'           observing all states with no direct impact.
#' @param model_solve_fn Function: theta -> decision-rules list with
#'   \code{ghx} (must contain the state-space matrices). If it returns
#'   generic moments instead of a DR list, the spectral Gram falls back to
#'   zero-rank -- supply a DR-solve closure (e.g.
#'   \code{function(th) solve_perturbation(model, th)}) for a meaningful
#'   result.
#' @param theta Numeric vector of parameter values at the calibration point.
#' @param param_names Optional character vector of parameter names.
#' @param Sigma_e Shock covariance matrix (n_exo x n_exo). If NULL, defaults
#'                to identity (unit-variance structural shocks).
#' @param n_freq Number of frequency grid points (default 256). More points
#'               gives more accurate Gram matrix integration.
#' @param eps Step size for finite differences (default 1e-5).
#' @return dynhr_diagnostic list with results including:
#'   - Gram matrix
#'   - Singular values with frequency band decomposition
#'   - Frequency bands driving identification per parameter
#' @references
#'   Qu, Z., & Tkachenko, D. (2012). Identification and frequency domain
#'     analysis of DSGE models. *Journal of Econometrics*, 168(1), 35-54.
#'   Qu, Z., & Tkachenko, D. (2017). Global identification in DSGE models
#'     based on the spectral density. *Econometrica*, 85(5), 1571-1638.
#' @noRd
d23_spectral_identification_impl <- function(dr,
                                        model_solve_fn = NULL,
                                        theta          = NULL,
                                        param_names    = NULL,
                                        Sigma_e        = NULL,
                                        n_freq         = 256L,
                                        eps            = 1e-5,
                                        meta           = NULL) {

    # ---------------------------------------------------------------
    # 1. Extract state-space matrices from decision rules
    # ---------------------------------------------------------------
    ghx <- as.matrix(dr$ghx)
    ghu <- as.matrix(dr$ghu)
    n_state <- ncol(ghx)            # ghx is n_endo x n_state
    n_endo  <- nrow(ghx)            # number of endogenous variables
    n_shock <- ncol(ghu)

    # Extract state transition submatrices
    # ghx/ghu are n_endo x n_state / n_endo x n_shock.
    # The state-to-state transition and shock-to-state impact
    # are the rows corresponding to state variables.
    state_idx <- dr$state_idx %||% seq_len(n_state)
    T_mat <- ghx[state_idx, , drop = FALSE]  # n_state x n_state
    R_mat <- ghu[state_idx, , drop = FALSE]  # n_state x n_shock

    # Validate dimensions
    if (nrow(ghu) != n_endo) {
      stop(sprintf("ghu rows (%d) != ghx rows (%d). Decision rules have incompatible dimensions.",
                   nrow(ghu), n_endo))
    }
    if (nrow(T_mat) != n_state || ncol(T_mat) != n_state) {
      stop(sprintf("State transition matrix dims (%dx%d) != expected (%dx%d). Check state_idx.",
                   nrow(T_mat), ncol(T_mat), n_state, n_state))
    }
    if (nrow(R_mat) != n_state || ncol(R_mat) != n_shock) {
      stop(sprintf("Shock impact matrix dims (%dx%d) != expected (%dx%d). Check state_idx.",
                   nrow(R_mat), ncol(R_mat), n_state, n_shock))
    }

    # Observation matrix C and direct impact D
    if (!is.null(dr$obs_mat)) {
      C <- as.matrix(dr$obs_mat)
      if (ncol(C) != n_state) {
        stop(sprintf("obs_mat columns (%d) != n_state (%d). Observation matrix incompatible with state vector.",
                     ncol(C), n_state))
      }
    } else if (!is.null(dr$Z)) {
      C <- as.matrix(dr$Z)
      if (ncol(C) != n_state) {
        stop(sprintf("Z matrix columns (%d) != n_state (%d). Observation matrix incompatible.",
                     ncol(C), n_state))
      }
    } else {
      # Default: observe all states
      C <- diag(n_state)
    }
    n_obs <- nrow(C)

    if (!is.null(dr$D_mat)) {
      D <- as.matrix(dr$D_mat)
      if (nrow(D) != n_obs || ncol(D) != n_shock) {
        warning(sprintf("D_mat dims (%dx%d) != expected (%dx%d). Resizing to zero matrix.",
                        nrow(D), ncol(D), n_obs, n_shock))
        D <- matrix(0, nrow = n_obs, ncol = n_shock)
      }
    } else {
      D <- matrix(0, nrow = n_obs, ncol = n_shock)
    }

    # Shock covariance (L5 fix: prefer model-derived covariance when available)
    # When model is attached to dr (dr$model), use .get_shock_cov so that the
    # spectral density matches the Whittle likelihood and Kalman filter which
    # always call .get_shock_cov.  Fall back to identity only as a last resort,
    # with a warning that the Gram matrix may be wrong.
    if (is.null(Sigma_e)) {
      if (!is.null(dr$model)) {
        Sigma_e <- tryCatch(
          .get_shock_cov(dr$model, dr$model$varexo_names, dr$model$param_values),
          error = function(e) NULL
        )
      }
      if (is.null(Sigma_e)) {
        warning("d23: Sigma_e is NULL; defaulting to identity. ",
                "Pass Sigma_e explicitly or attach the model object as dr$model ",
                "so that D23 uses the same shock covariance as Whittle/KF. ",
                "When shocks have unequal variances or are correlated the Gram ",
                "matrix will be incorrect.",
                call. = FALSE)
        Sigma_e <- diag(n_shock)
      }
    } else {
      Sigma_e <- as.matrix(Sigma_e)
      if (nrow(Sigma_e) != n_shock || ncol(Sigma_e) != n_shock) {
        stop(sprintf("Sigma_e dims (%dx%d) != n_shock=%d.", nrow(Sigma_e), ncol(Sigma_e), n_shock))
      }
    }

    n_par <- length(theta)
    if (is.null(param_names)) {
      param_names <- names(theta) %||% paste0("theta_", seq_len(n_par))
    }

    # ---------------------------------------------------------------
    # 2. Construct frequency grid (0, π)  — half spectrum
    # ---------------------------------------------------------------
    n_freq <- as.integer(n_freq)
    freq_grid <- seq(0, pi, length.out = n_freq)

    # ---------------------------------------------------------------
    # 3. Compute spectral density on frequency grid at baseline theta
    # ---------------------------------------------------------------
    .spectral_density <- function(omega, T_mat, R_mat, C_mat, D_mat, Sigma)
      .spectral_density_core(omega, TT = T_mat, RR = R_mat, ZZ = C_mat,
                              DD = D_mat, Sigma_e = Sigma)

    # Compute baseline spectral density at each frequency
    # Store as a vectorized list
    spectral_list <- lapply(freq_grid, function(omega) {
      .spectral_density(omega, T_mat, R_mat, C, D, Sigma_e)
    })

    # ---------------------------------------------------------------
    # 4. Build spectral Gram matrix by numerically differentiating the
    #    frequency-domain spectral density directly from the state-space
    #    matrices already extracted above (T_mat, R_mat, C, D, Sigma_e).
    #
    #    This path is ALWAYS preferred over differentiating model_solve_fn
    #    because it operates in the frequency domain: the resulting Gram
    #    matrix, rank, and SV plot are genuinely different from D1's
    #    time-domain moment Gram.
    #
    #    Strategy (Qu & Tkachenko 2012):
    #      - Vectorise the unique upper-triangle elements of Σ_y(ω) at each
    #        frequency on a sub-sampled grid.
    #      - Stack into a single moment vector m(θ) across frequencies.
    #      - Finite-difference m(θ) w.r.t. θ to get J_spec (n_mom × n_par).
    #      - Gram = J_spec' J_spec  (approximates the spectral information
    #        matrix integrated over the frequency domain).
    #
    #    When theta is non-NULL, finite-difference the spectral moments
    #    directly via the state-space transfer function.  When theta is
    #    NULL we fall back to the dr-based spectral density without
    #    parameter derivatives.
    # ---------------------------------------------------------------

    # Frequency bands (in cycles per period, relative to π):
    #   Low:        ω < π/6   (~ > 3-year cycles for quarterly data)
    #   Business:   π/6 ≤ ω ≤ π/1.5  (~ 6-18 quarter cycles)
    #   High:       ω > π/1.5  (~ < 6 quarter cycles)
    freq_bands <- list(
      low      = which(freq_grid < pi / 6),
      business = which(freq_grid >= pi / 6 & freq_grid <= pi / 1.5),
      high     = which(freq_grid > pi / 1.5)
    )
    band_names <- c("low", "business", "high")

    if (!is.null(theta) && n_par > 0) {
      # ------------------------------------------------------------------
      # Build a closure that maps theta -> spectral-moment vector.
      # The closure re-computes the transfer function from the CURRENT
      # state-space extracted at the top (T_mat, R_mat, C, D, Sigma_e).
      # When model_solve_fn is also supplied we use it ONLY to re-extract
      # updated ghx/ghu at perturbed theta; otherwise we perturb theta
      # symbolically (useful when theta controls T_mat/R_mat directly).
      #
      # Practical note: if neither model_solve_fn nor a dr-re-solve path
      # is wired in, the spectral moments at perturbed theta are identical
      # to baseline (constant function), yielding a zero Jacobian.  In that
      # case we report the spectral density but set rank_J = 0 and warn.
      # ------------------------------------------------------------------

      # Sub-sample frequencies for tractability: take at most 64 evenly
      # spaced points from the full grid (avoids n_freq*n_unique rows).
      n_freq_sub  <- min(64L, n_freq)
      freq_sub_idx <- round(seq(1, n_freq, length.out = n_freq_sub))
      freq_sub     <- freq_grid[freq_sub_idx]

      # Real coordinate dimension of the complex Hermitian Sigma_y:
      #   - n_obs diagonal elements contribute Re only (Im == 0 by Hermitian symmetry)
      #   - n_obs*(n_obs-1)/2 strictly-upper-triangle elements each contribute
      #     Re AND Im as separate real coordinates
      # Total = n_obs + 2 * n_obs*(n_obs-1)/2 = n_obs^2
      n_unique     <- as.integer(n_obs^2L)
      n_spec_mom   <- n_freq_sub * n_unique

      # Helper: vectorise the Hermitian spectral density into real coordinates.
      # Order: diagonal Re values first, then for each strictly-upper-triangle
      # (i < j) pair: Re(Sy[i,j]) then Im(Sy[i,j]).
      .spec_vec_at <- function(T_loc, R_loc, C_loc, D_loc, Sig_loc, freqs) {
        out <- numeric(length(freqs) * n_unique)
        for (fi in seq_along(freqs)) {
          Sy  <- .spectral_density(freqs[fi], T_loc, R_loc, C_loc, D_loc, Sig_loc)
          uv  <- numeric(n_unique)
          idx <- 1L
          # Diagonal (real only)
          for (k in seq_len(n_obs)) {
            uv[idx] <- Re(Sy[k, k])
            idx <- idx + 1L
          }
          # Strictly upper triangle: Re then Im
          if (n_obs > 1L) {
            for (j in seq(2L, n_obs)) {
              for (i in seq_len(j - 1L)) {
                uv[idx]     <- Re(Sy[i, j])
                uv[idx + 1L] <- Im(Sy[i, j])
                idx <- idx + 2L
              }
            }
          }
          out[((fi - 1L) * n_unique + 1L):(fi * n_unique)] <- uv
        }
        out
      }

      if (!is.null(model_solve_fn)) {
        # model_solve_fn: theta -> moment vector.  We check whether it
        # returns decision-rule matrices or generic moments.
        # For spectral D23 we MUST have it return moments that vary with
        # theta; use it to rebuild the state-space at perturbed theta.
        # Attempt to detect if model_solve_fn returns a dr list.
        test_out <- tryCatch(model_solve_fn(theta), error = function(e) NULL)
        is_dr_fn <- is.list(test_out) && ("ghx" %in% names(test_out))

        if (is_dr_fn) {
          # model_solve_fn is a dr-solve function: rebuild spectral density
          # at each perturbed theta.
          .spec_mom_fn <- function(th) {
            dr_th <- tryCatch(model_solve_fn(th), error = function(e) NULL)
            if (is.null(dr_th)) return(rep(NA_real_, n_spec_mom))
            ghx_th  <- as.matrix(dr_th$ghx)
            ghu_th  <- as.matrix(dr_th$ghu)
            n_st_th <- nrow(ghx_th)
            si_th   <- dr_th$state_idx %||% seq_len(ncol(ghx_th))
            T_th    <- ghx_th[si_th, , drop = FALSE]
            R_th    <- ghu_th[si_th, , drop = FALSE]
            C_th    <- if (!is.null(dr_th$obs_mat)) as.matrix(dr_th$obs_mat)
                       else if (!is.null(dr_th$Z))   as.matrix(dr_th$Z)
                       else diag(nrow(T_th))
            D_th    <- if (!is.null(dr_th$D_mat)) as.matrix(dr_th$D_mat)
                       else matrix(0, nrow = nrow(C_th), ncol = ncol(R_th))
            .spec_vec_at(T_th, R_th, C_th, D_th, Sigma_e, freq_sub)
          }
        } else {
          # model_solve_fn returns generic moments (not a dr list).
          # Fall through to the baseline-only path below.
          is_dr_fn <- FALSE
          .spec_mom_fn <- NULL
        }
      } else {
        .spec_mom_fn <- NULL
      }

      if (is.null(.spec_mom_fn)) {
        # No re-solve available: spectral moments are constant w.r.t. theta.
        # Report spectral density but set rank to 0 with a note.
        warning("d23: cannot finite-difference spectral density w.r.t. theta ",
                "because no dr-solve function was detected. ",
                "Gram matrix will be zero. Provide a dr-solve function as ",
                "model_solve_fn to obtain the spectral identification rank.")
        J_mom          <- matrix(0, nrow = n_spec_mom, ncol = n_par)
        colnames(J_mom) <- param_names
        singular_values <- rep(0, n_par)
        names(singular_values) <- paste0("sv_", seq_len(n_par))
        rank_J         <- 0L
        full_rank      <- FALSE
        gram_matrix    <- matrix(0, n_par, n_par)
        weak_params    <- param_names
        weak_threshold <- 0
        param_band_contrib <- NULL
      } else {
        # Compute spectral Jacobian via finite differences of spectral moments
        J_mom <- .numerical_jacobian(.spec_mom_fn, theta, eps = eps)

        if (is.null(J_mom) || !all(is.finite(J_mom))) {
          warning("d23: spectral Jacobian contains non-finite values. ",
                  "Some parameters may not affect the spectral density at theta.")
          if (!is.null(J_mom)) J_mom[!is.finite(J_mom)] <- 0
          else J_mom <- matrix(0, nrow = n_spec_mom, ncol = n_par)
        }

        # Guard: ensure dimnames match Jacobian dimensions
        if (ncol(J_mom) != length(param_names)) {
          warning(sprintf("d23: param_names length (%d) != Jacobian columns (%d). Using generic labels.",
                          length(param_names), ncol(J_mom)))
          param_names <- paste0("theta_", seq_len(ncol(J_mom)))
        }
        colnames(J_mom) <- param_names

        # Spectral Gram matrix (approximates ∫ J_spec(ω)' J_spec(ω) dω)
        gram_matrix     <- crossprod(J_mom)
        sv              <- svd(J_mom)
        singular_values <- sv$d
        names(singular_values) <- paste0("sv_", seq_along(singular_values))

        tol        <- max(dim(J_mom)) * max(singular_values) * .Machine$double.eps
        rank_J     <- sum(singular_values > tol)
        full_rank  <- (rank_J == n_par)

        weak_threshold <- max(singular_values) * 1e-3
        n_weak         <- sum(singular_values < weak_threshold)
        weak_params    <- character(0)
        if (n_weak > 0 && n_par > 0) {
          n_v_cols  <- ncol(sv$v)
          start_col <- max(1L, n_v_cols - n_weak + 1L)
          if (start_col <= n_v_cols) {
            V_weak <- sv$v[, start_col:n_v_cols, drop = FALSE]
            for (k in seq_len(ncol(V_weak))) {
              idx <- which.max(abs(V_weak[, k]))
              weak_params <- c(weak_params, param_names[idx])
            }
            weak_params <- unique(weak_params)
          }
        }

        # ---------------------------------------------------------------
        # Frequency band decomposition: attribute Gram contribution to
        # each of the three spectral bands (low / business / high).
        # For each band b we take the rows of J_mom corresponding to
        # frequencies in that band, compute their sub-Gram, and report
        # the diagonal (per-parameter column norm) as a % of total.
        # ---------------------------------------------------------------
        n_par_bands <- ncol(J_mom)
        band_param_names <- head(param_names, n_par_bands)
        param_band_contrib <- matrix(NA_real_, nrow = n_par_bands, ncol = 3L,
                                     dimnames = list(band_param_names, band_names))
        for (b in seq_along(band_names)) {
          b_freq_idx <- freq_bands[[band_names[b]]]
          # Map from full freq_grid indices to sub-sampled freq_sub_idx positions
          b_sub_positions <- which(freq_sub_idx %in% b_freq_idx)
          if (length(b_sub_positions) == 0L) {
            param_band_contrib[, b] <- 0
            next
          }
          # Row indices in J_mom for this band (each freq contributes n_unique rows)
          row_idx <- unlist(lapply(b_sub_positions, function(fi)
            ((fi - 1L) * n_unique + 1L):(fi * n_unique)))
          row_idx <- row_idx[row_idx <= nrow(J_mom)]
          J_band  <- J_mom[row_idx, , drop = FALSE]
          param_band_contrib[, b] <- colSums(J_band^2)
        }
        # Normalise rows to percentages
        row_sums <- rowSums(param_band_contrib, na.rm = TRUE)
        for (i in seq_len(n_par_bands)) {
          if (!is.na(row_sums[i]) && row_sums[i] > 0)
            param_band_contrib[i, ] <- param_band_contrib[i, ] / row_sums[i] * 100
        }
      }

    } else {
      # No theta available: report spectral density only; no Gram / rank info.
      singular_values <- numeric(0)
      rank_J          <- 0L
      full_rank       <- NA
      weak_params     <- character(0)
      gram_matrix     <- matrix(NA_real_, 0, 0)
      param_band_contrib <- NULL
      J_mom           <- NULL
      weak_threshold  <- 0
    }

    pass <- if (is.na(full_rank)) NA else full_rank

    # ---------------------------------------------------------------
    # 7. Build plots
    # ---------------------------------------------------------------
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE) &&
        requireNamespace("reshape2", quietly = TRUE)) {

      # (a) Spectral density heatmap across frequencies
      # For plotting purposes use real-valued summaries:
      #   - diagonal elements: Re (auto-spectra, real by Hermitian symmetry)
      #   - off-diagonal elements: Mod (cross-spectral amplitude)
      spec_mat <- do.call(cbind, lapply(spectral_list, function(S) {
        S_real <- matrix(0, nrow = n_obs, ncol = n_obs)
        for (k in seq_len(n_obs)) S_real[k, k] <- Re(S[k, k])
        if (n_obs > 1L) {
          for (j in seq(2L, n_obs)) {
            for (i in seq_len(j - 1L)) {
              S_real[i, j] <- Mod(S[i, j])
              S_real[j, i] <- S_real[i, j]
            }
          }
        }
        S_vec <- S_real[lower.tri(S_real, diag = TRUE)]
        S_vec[is.na(S_vec) | is.nan(S_vec) | is.infinite(S_vec)] <- 0
        S_vec
      }))
      spec_long <- reshape2::melt(spec_mat)
      colnames(spec_long) <- c("Element", "FreqIdx", "Value")
      spec_long$Frequency <- freq_grid[spec_long$FreqIdx]

      # Only plot first 20 elements to avoid overcrowding
      n_spec_elements <- min(20, nrow(spec_mat))
      spec_plot <- spec_long[spec_long$Element <= n_spec_elements, ]

      # Build observable labels. Element index maps to the lower-triangle
      # (including diagonal) of the n_obs x n_obs spectral matrix in
      # column-major order.  Here we just label element k as the k-th
      # diagonal observable (auto-spectra) or "cross-k" for off-diagonal.
      obs_labels <- if (!is.null(dr$obs_names)) dr$obs_names
                    else if (!is.null(dr$model$obs_names)) dr$model$obs_names
                    else paste0("obs", seq_len(n_obs))
      diag_idx  <- cumsum(c(1L, seq_len(n_obs - 1L) + cumsum(seq_len(n_obs - 1L))))
      elem_labs <- vapply(seq_len(n_spec_elements), function(k) {
        # Is this element a diagonal entry? If so use observable name.
        di <- which(diag_idx == k)
        if (length(di) > 0 && di[1] <= length(obs_labels)) {
          obs_labels[di[1]]
        } else {
          sprintf("cross-%d", k)
        }
      }, character(1))
      spec_plot$ElemLabel <- factor(
        elem_labs[spec_plot$Element],
        levels = unique(elem_labs[seq_len(n_spec_elements)])
      )

      # Floor small values so log scale doesn't crash on zeroes or NA.
      # Then drop elements that are entirely at the floor (e.g., zero cross-spectra
      # for uncorrelated observables) since they add nothing but noise lines.
      spec_plot$Value <- pmax(spec_plot$Value, 1e-10)
      spec_plot <- spec_plot[is.finite(spec_plot$Value) & spec_plot$Value > 0, ]
      # Filter: keep only elements with max value > 1e-8 (at least some signal)
      elem_maxval <- tapply(spec_plot$Value, spec_plot$Element, max, na.rm = TRUE)
      active_elements <- as.integer(names(elem_maxval[elem_maxval > 1e-8]))
      spec_plot <- spec_plot[spec_plot$Element %in% active_elements, ]

      p_spec <- ggplot2::ggplot(
        spec_plot,
        ggplot2::aes(x = Frequency, y = Value, colour = ElemLabel)
      ) +
        ggplot2::geom_line(alpha = 0.8, linewidth = 0.6) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_log10(
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        scale_colour_dynhr_vibrant() +
        ggplot2::guides(colour = ggplot2::guide_legend(
          title = "Observable", ncol = 3L, override.aes = list(linewidth = 1.2)
        )) +
        theme_dynhr_diagnostic() +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::labs(
          title = "D23: Spectral density of observables (log scale)",
          subtitle = sprintf("n_obs = %d, n_freq = %d. Auto-spectra shown with observable labels.", n_obs, n_freq),
          x = "Frequency (radians)", y = "Spectral density (log10)"
        )
      plots$spectral_density <- .apply_meta(p_spec, meta)

      # (b) Singular values if available
      if (length(singular_values) > 0) {
        sv_df <- data.frame(
          index = seq_along(singular_values),
          value = singular_values
        )
        sv_df$identified <- ifelse(
          sv_df$value >= weak_threshold, "Identified", "Weak / unidentified"
        )

        # Detect SVs that are floored to 1e-300 (truly zero / unreachable)
        # vs genuinely small-but-nonzero values.  Floored SVs get an annotation.
        floor_val <- .Machine$double.xmin  # ~2.2e-308
        sv_df$floored <- sv_df$value < floor_val * 1e10  # below ~2e-298

        p_sv <- ggplot2::ggplot(
          sv_df,
          ggplot2::aes(x = factor(index), y = pmax(value, 1e-300), fill = identified)
        ) +
          ggplot2::geom_col(width = 0.7) +
          ggplot2::geom_hline(
            yintercept = weak_threshold,
            linetype = "dashed",
            colour = dynhr_colours$red,
            linewidth = 0.5
          ) +
          ggplot2::scale_fill_manual(
            values = c(
              "Identified" = dynhr_colours$mid_blue,
              "Weak / unidentified" = dynhr_colours$red
            ),
            name = NULL
          ) +
          ggplot2::scale_y_log10(
            labels = scales::trans_format("log10", scales::math_format(10^.x))
          ) +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title = "D23: Singular values of frequency-domain spectral Gram",
            subtitle = sprintf(
              "Spectral Gram J'J (freq-domain), threshold = %.2e (0.1%% of max SV). Bars at floor = structurally zero.",
              weak_threshold),
            x = "Singular value index",
            y = "Value (log10)"
          )

        # Annotate floored SVs so the reader isn't misled by the log-axis floor
        if (any(sv_df$floored)) {
          p_sv <- p_sv +
            ggplot2::annotate("text",
              x     = which(sv_df$floored),
              y     = 1e-250,
              label = "~0",
              size  = 2.5,
              colour = dynhr_colours$red,
              vjust  = 0)
        }
        plots$singular_values <- .apply_meta(p_sv, meta)
      }

      # (c) Frequency band contribution heatmap
      if (!is.null(param_band_contrib) && nrow(param_band_contrib) > 0) {
        band_df <- reshape2::melt(param_band_contrib)
        colnames(band_df) <- c("Parameter", "Band", "Contribution")
        band_df$Band <- factor(band_df$Band, levels = c("low", "business", "high"))

        p_band <- ggplot2::ggplot(
          band_df,
          ggplot2::aes(x = Band, y = Parameter, fill = Contribution)
        ) +
          ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
          scale_fill_dynhr_cividis(name = "% contribution") +
          theme_dynhr_diagnostic() +
          ggplot2::labs(
            title = "D23: Frequency band contribution to identification",
            subtitle = "Low (< pi/6), Business cycle, High (> pi/1.5)",
            x = "Frequency band",
            y = NULL
          )
        plots$band_contribution <- .apply_meta(p_band, meta)
      }
    }

    # ---------------------------------------------------------------
    # 8. Build summary and result
    # ---------------------------------------------------------------
    # Detect whether Sigma_e was explicitly provided or defaulted to identity.
    # The warning path above sets it to identity; we surface this in the summary
    # so users know the Gram matrix may be wrong for models with unequal shock variances.
    sigma_e_note <- if (!is.null(dr$model) || !is.null(match.call()$Sigma_e)) ""
                    else " [NOTE: Sigma_e defaulted to identity -- pass Sigma_e explicitly for correct Gram matrix.]"

    summary_text <- sprintf(
      "D23 Spectral identification (Qu-Tkachenko 2012): %s. %s%s Gram matrix built from frequency-domain spectral Jacobian (J_spec'J_spec, %d sub-sampled frequencies); band decomposition uses column norms of J_spec by frequency band.%s",
      if (isTRUE(pass)) "PASS -- full rank" else if (identical(pass, FALSE))
        "FAIL -- rank deficient" else "INFO -- spectral density computed",
      if (isTRUE(full_rank))
        sprintf("Gram matrix rank = %d / %d parameters.", rank_J, n_par)
      else if (identical(full_rank, FALSE))
        sprintf("Gram matrix rank = %d / %d. RANK DEFICIENT.", rank_J, n_par)
      else "",
      if (length(weak_params) > 0)
        sprintf(" Weakly identified: %s.", paste(weak_params, collapse = ", "))
      else "",
      if (exists("n_freq_sub", inherits = FALSE)) n_freq_sub else n_freq,
      sigma_e_note
    )

    .make_result(
      result = list(
        gram_matrix = gram_matrix,
        singular_values = singular_values,
        rank = rank_J,
        full_rank = full_rank,
        weak_params = weak_params,
        frequency_grid = freq_grid,
        spectral_density = spectral_list,
        band_contribution = param_band_contrib,
        state_space = list(T = ghx, R = ghu, C = C, D = D)
      ),
      pass = pass,
      plots = plots,
      summary = summary_text,
      llm_summary = {
        badge <- if (isTRUE(pass)) "PASS"
        else if (identical(pass, FALSE)) "FAIL"
        else "INFO"
        sv_str <- if (length(singular_values) > 0) {
          paste(sprintf("%.3e", head(sort(singular_values), 5)), collapse = ", ")
        } else "N/A"
        paste(c(
          sprintf("D23 | Spectral Identification (Qu-Tkachenko 2012) | %s", badge),
          sprintf("  n_freq=%d n_par=%d gram_rank=%d",
                  n_freq, n_par, rank_J),
          sprintf("  smallest_sv: %s", sv_str),
          if (length(weak_params) > 0)
            sprintf("  weak_params: %s", paste(weak_params, collapse = ", ")),
          sprintf("  action: %s",
                  if (isTRUE(pass))
                    "Model is spectrally identified. Frequency band decomposition available."
                  else if (identical(pass, FALSE))
                    sprintf("Rank deficient (rank=%d, expected=%d). Check for collinear parameter effects on spectral density. Consider data filtering or frequency-domain restrictions.",
                            rank_J, n_par)
                  else
                    "Spectral density computed. Provide model_solve_fn for full identification rank check.")
        ), collapse = "\n")
      }
    )
}


#' Build a spectral-moment function for D23 numerical differentiation
#'
#' Creates a closure that, given a first-order solution function, computes
#' the unique elements of the spectral density matrix at each frequency on
#' a grid, vectorized across frequencies. This is used as the
#' \code{model_solve_fn} for D23 to compute the spectral Jacobian.
#'
#' @param dr_solve_fn Function: theta -> decision rules list (must contain
#'                     ghx, ghu, and optionally obs_mat/D_mat).
#' @param theta0 Numeric vector of baseline parameter values.
#' @param param_names Character vector of parameter names.
#' @param Sigma_e Shock covariance matrix.
#' @param n_freq Number of frequency grid points (default 256).
#' @param n_mom_max Maximum number of moment rows to return (default 500).
#'                   Spectral density can produce many moments
#'                   (n_obs^2 * n_freq), so we cap for tractability.
#' @return A function: theta -> numeric vector of spectral moments.
#' @noRd
.build_spectral_moment_fn <- function(dr_solve_fn,
                                      theta0,
                                      param_names,
                                      Sigma_e = NULL,
                                      n_freq = 256L,
                                      n_mom_max = 500L) {
  # Get baseline decision rules
  dr0 <- dr_solve_fn(theta0)

  ghx <- as.matrix(dr0$ghx)
  ghu <- as.matrix(dr0$ghu)
  n_state <- nrow(ghx)
  n_shock <- ncol(ghu)

  C <- if (!is.null(dr0$obs_mat)) as.matrix(dr0$obs_mat)
  else if (!is.null(dr0$Z)) as.matrix(dr0$Z)
  else diag(n_state)
  n_obs <- nrow(C)

  D <- if (!is.null(dr0$D_mat)) as.matrix(dr0$D_mat)
  else matrix(0, nrow = n_obs, ncol = n_shock)

  if (is.null(Sigma_e)) Sigma_e <- diag(n_shock)

  n_freq <- as.integer(n_freq)
  freq_grid <- seq(0, pi, length.out = n_freq)
  # Real coordinate dimension: n_obs diagonal (Re only) + 2*n_obs*(n_obs-1)/2
  # strictly-upper-triangle (Re + Im each) = n_obs^2
  n_unique <- as.integer(n_obs^2L)
  n_total <- n_unique * n_freq

  # Select a subset of frequency-moment pairs if too large
  if (n_total > n_mom_max) {
    # Stratified sampling: pick indices evenly spaced across frequencies
    n_per_freq <- max(1, floor(n_mom_max / n_freq))
    select_idx <- unlist(lapply(seq_len(n_freq), function(f) {
      base <- (f - 1) * n_unique
      freq_idx <- base + sort(sample.int(n_unique, min(n_per_freq, n_unique)))
      freq_idx
    }))
    select_idx <- head(select_idx, n_mom_max)
  } else {
    select_idx <- seq_len(n_total)
  }

  function(theta) {
    dr <- dr_solve_fn(theta)

    ghx_new <- as.matrix(dr$ghx)
    ghu_new <- as.matrix(dr$ghu)

    # Build spectral densities at each frequency
    spec_vec <- numeric(0)
    for (omega in freq_grid) {
      Sigma_y <- .spectral_density_core(omega, TT = ghx_new, RR = ghu_new,
                                         ZZ = C, DD = D, Sigma_e = Sigma_e)

      # Extract real coordinates from Hermitian Sigma_y:
      #   - diagonal: Re only
      #   - strictly upper triangle: Re then Im (as separate coordinates)
      uv <- numeric(n_unique)
      idx <- 1L
      for (k in seq_len(n_obs)) {
        uv[idx] <- Re(Sigma_y[k, k])
        idx <- idx + 1L
      }
      if (n_obs > 1L) {
        for (j in seq(2L, n_obs)) {
          for (i in seq_len(j - 1L)) {
            uv[idx]     <- Re(Sigma_y[i, j])
            uv[idx + 1L] <- Im(Sigma_y[i, j])
            idx <- idx + 2L
          }
        }
      }
      spec_vec <- c(spec_vec, uv)
    }

    # Return selected subset
    spec_vec[select_idx]
  }
}
