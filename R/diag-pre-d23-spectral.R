## R/diag-pre-d23-spectral.R
## --------------------------------------------------------------------------
## D23 Spectral identification check (Qu & Tkachenko 2012) -- implementation.
## The public-facing entry point `d23_spectral_identification()` lives in
## diag-pre-d2-spectral.R and delegates here.
##
## Definition (QT 2012, Theorem 1):
##   f_theta(w) = (1 / 2 pi) H(e^{-iw}) Sigma_e H(e^{-iw})^*
##   G(theta)   = int_{-pi}^{pi} (d vec f_theta(w) / d theta')^*
##                               (d vec f_theta(w) / d theta') dw
## theta is locally identified from the second-order properties of the
## observables iff G(theta) has full rank.
##
## Implementation:
##   * H is the dynhr lagged-state transfer function (.spectral_density_core):
##       s_t = T s_{t-1} + R e_t,  y_t = Z s_{t-1} + D e_t.
##   * The integral is the trapezoid rule on [0, pi] doubled -- identical to
##     the Fourier-frequency Riemann sum on [-pi, pi) with N = 2 (n_freq - 1)
##     points, because the integrand is even (f(-w) = conj f(w)) and
##     2 pi-periodic (so the rule is spectrally accurate).
##   * Each frequency contributes the real coordinates of vec(df) scaled by the
##     square-root quadrature weight: Re of the diagonal, and sqrt(2) x (Re, Im)
##     of the strict upper triangle (the lower triangle is its conjugate).
##     With that scaling crossprod(J) == G exactly (up to quadrature).
##   * d/dtheta by central differences through a theta -> state-space closure.
##     Sigma_e is re-read at every theta when the closure supplies it, so
##     estimated shock standard deviations get their true derivative.
##   * Rank uses the column-equilibrated J and an FD-noise-aware tolerance
##     (same device as D1): a second Jacobian at step 2*eps estimates the
##     numerical error; singular values within 10x of it are not separable
##     from zero.
## --------------------------------------------------------------------------

#' D23. Spectral identification check (Qu & Tkachenko 2012)
#'
#' Computes the Qu-Tkachenko Gram matrix
#' \deqn{G(\theta) = \int_{-\pi}^{\pi} \left(\frac{\partial\, vec f_\theta(\omega)}
#'   {\partial \theta'}\right)^{*} \frac{\partial\, vec f_\theta(\omega)}
#'   {\partial \theta'}\, d\omega,\qquad
#'   f_\theta(\omega) = \frac{1}{2\pi} H(e^{-i\omega}) \Sigma_e H(e^{-i\omega})^{*}}
#' for the first-order state space in dynhr's lagged-state convention
#' \eqn{s_t = T s_{t-1} + R e_t}, \eqn{y_t = Z s_{t-1} + D e_t}
#' (as returned by \code{build_dsge_state_space()}), and checks its rank.
#'
#' The derivative is taken by central finite differences, so a
#' \code{model_solve_fn} that re-solves the state space at a perturbed
#' \code{theta} is REQUIRED for the rank check. Without one (or when it
#' returns something that is not a state-space list) the result is
#' \code{pass = NA}: that is a missing input, not an identification failure.
#'
#' @param dr State space at \code{theta}: a list with \code{ghx} (state
#'   transition, or the full \code{n_endo x n_state} decision rule together
#'   with \code{state_idx}), \code{ghu}, and optionally \code{obs_mat}
#'   (or \code{Z}; lagged-state loading, default: all states), \code{D_mat}
#'   (default 0), \code{Sigma_e} and \code{obs_names}.
#' @param model_solve_fn Function \code{theta -> } list of the same form as
#'   \code{dr} (\code{NULL} for an infeasible \code{theta}). A
#'   \code{Sigma_e} element in its output is used at that \code{theta}, which
#'   is how estimated shock standard deviations enter the derivative.
#' @param theta Named numeric parameter vector.
#' @param param_names Optional parameter names (default \code{names(theta)}).
#' @param Sigma_e Shock covariance used when the state-space lists carry
#'   none. Resolution order: the closure's \code{Sigma_e}; this argument;
#'   \code{dr$Sigma_e}; \code{.get_shock_cov(dr$model)}. There is no identity
#'   default -- with none of these the function aborts.
#' @param n_freq Number of frequencies on \eqn{[0, \pi]} (trapezoid rule).
#' @param eps Central finite-difference step.
#' @param tol_rank Optional relative rank tolerance (fraction of the largest
#'   singular value); \code{NULL} uses the finite-difference noise estimate.
#' @param weak_rel Singular values below \code{weak_rel} times the largest
#'   (but above the rank tolerance) are reported as weak; informational only.
#' @param meta Provenance descriptor from \code{diag_meta()}.
#' @return dynhr_diagnostic list; \code{result} holds \code{gram_matrix}
#'   (the QT \eqn{G(\theta)}, parameter units), \code{singular_values} (of the
#'   column-equilibrated Jacobian, which the rank uses), \code{rank},
#'   \code{rank_tolerance}, \code{fd_error}, \code{unidentified_params},
#'   \code{weak_params}, \code{band_contribution} (\% of each diagonal
#'   \eqn{G_{kk}} from the low / business-cycle / high bands),
#'   \code{frequency_grid} and \code{spectral_density} (the \eqn{2\pi f}
#'   matrices, as before).
#' @references
#'   Qu, Z., & Tkachenko, D. (2012). Identification and frequency domain
#'     quasi-maximum likelihood estimation of linearized dynamic stochastic
#'     general equilibrium models. \emph{Quantitative Economics}, 3(1), 95-132.
#' @noRd
d23_spectral_identification_impl <- function(dr,
                                             model_solve_fn = NULL,
                                             theta          = NULL,
                                             param_names    = NULL,
                                             Sigma_e        = NULL,
                                             n_freq         = 256L,
                                             eps            = 1e-5,
                                             tol_rank       = NULL,
                                             weak_rel       = 1e-3,
                                             meta           = NULL) {

  ss <- .d23_state_space(dr)
  n_obs   <- nrow(ss$ZZ)
  n_shock <- ncol(ss$RR)

  n_par <- length(theta)
  test_out <- NULL
  if (n_par > 0L && is.function(model_solve_fn)) test_out <- model_solve_fn(theta)
  solve_ok <- is.list(test_out) && !is.null(test_out$ghx) && !is.null(test_out$ghu)

  ## ---- Shock covariance: no identity default (CLAUDE.md: Q=I is a bug) ----
  if (!is.null(Sigma_e)) {
    Sigma_e <- as.matrix(Sigma_e)
  } else if (!is.null(ss$Sigma_e)) {
    Sigma_e <- ss$Sigma_e
  } else if (!is.null(dr$model)) {
    Sigma_e <- .get_shock_cov(dr$model, dr$model$varexo_names,
                              dr$model$param_values)
  } else if (solve_ok && !is.null(test_out$Sigma_e)) {
    Sigma_e <- as.matrix(test_out$Sigma_e)
  } else {
    .dynhr_abort(paste0(
      "D23: no shock covariance. Pass Sigma_e (e.g. shock_cov(model)), or ",
      "supply it as dr$Sigma_e / in the model_solve_fn output."))
  }
  if (!identical(dim(Sigma_e), c(n_shock, n_shock)))
    .dynhr_abort(sprintf("D23: Sigma_e is %s, expected %d x %d.",
                         paste(dim(Sigma_e), collapse = " x "), n_shock, n_shock))

  obs_labels <- dr$obs_names %||% dr$model$obs_names %||% paste0("obs", seq_len(n_obs))
  if (length(obs_labels) != n_obs) obs_labels <- paste0("obs", seq_len(n_obs))

  ## ---- Frequency grid and trapezoid weights for int_{-pi}^{pi} ----
  n_freq <- as.integer(n_freq)
  if (n_freq < 3L) .dynhr_abort("D23: n_freq must be >= 3.")
  freq_grid <- seq(0, pi, length.out = n_freq)
  h <- pi / (n_freq - 1L)
  w_freq <- c(h, rep(2 * h, n_freq - 2L), h)

  spectral_list <- lapply(freq_grid, function(om)
    .spectral_density_core(om, TT = ss$TT, RR = ss$RR, ZZ = ss$ZZ,
                           DD = ss$DD, Sigma_e = Sigma_e))

  band_breaks <- c(2 * pi / 32, 2 * pi / 6)   # 32- and 6-period cycles
  band_names  <- c("low", "business", "high")
  band_of_freq <- ifelse(freq_grid < band_breaks[1], "low",
                  ifelse(freq_grid <= band_breaks[2], "business", "high"))

  if (is.null(param_names)) param_names <- names(theta) %||% paste0("theta_", seq_len(n_par))
  if (n_par > 0L && length(param_names) != n_par)
    .dynhr_abort(sprintf("D23: %d param_names for %d parameters.",
                         length(param_names), n_par))

  ## ---- Status of the rank computation ----
  status <- if (n_par == 0L) "no_theta" else if (!solve_ok) "no_solve_fn" else "ok"

  if (status == "ok" && !is.null(test_out$Sigma_e) &&
      !identical(dim(as.matrix(test_out$Sigma_e)), dim(Sigma_e)))
    .dynhr_abort("D23: model_solve_fn Sigma_e has the wrong dimension.")
  if (status == "ok" && !is.null(test_out$Sigma_e) &&
      max(abs(as.matrix(test_out$Sigma_e) - Sigma_e)) > 1e-10 * max(1, abs(Sigma_e))) {
    .dynhr_warn("D23: the Sigma_e returned by model_solve_fn(theta) differs from ",
                "the one supplied; the Gram matrix uses model_solve_fn's.")
  }

  ## rows per frequency: n_obs (Re diag) + n_obs (n_obs - 1) (Re, Im upper)
  n_row_f <- n_obs * n_obs
  up <- which(upper.tri(diag(n_obs)), arr.ind = TRUE)
  row_freq <- rep(seq_len(n_freq), each = n_row_f)

  spec_vec <- function(th) {
    d_th <- model_solve_fn(th)
    if (!is.list(d_th) || is.null(d_th$ghx)) return(rep(NA_real_, n_freq * n_row_f))
    s_th <- .d23_state_space(d_th)
    Sig  <- if (!is.null(s_th$Sigma_e)) s_th$Sigma_e else Sigma_e
    if (nrow(s_th$ZZ) != n_obs || ncol(s_th$RR) != n_shock ||
        !identical(dim(Sig), c(n_shock, n_shock)))
      .dynhr_abort("D23: model_solve_fn changed the state-space dimensions.")
    out <- numeric(n_freq * n_row_f)
    for (fi in seq_len(n_freq)) {
      S <- .spectral_density_core(freq_grid[fi], TT = s_th$TT, RR = s_th$RR,
                                  ZZ = s_th$ZZ, DD = s_th$DD, Sigma_e = Sig)
      S <- S / (2 * pi)
      offd <- S[up]
      out[(fi - 1L) * n_row_f + seq_len(n_row_f)] <-
        sqrt(w_freq[fi]) * c(Re(diag(S)), sqrt(2) * Re(offd), sqrt(2) * Im(offd))
    }
    out
  }

  gram_matrix <- NULL; singular_values <- numeric(0); V <- NULL
  rank_J <- NA_integer_; full_rank <- NA; tol <- NA_real_; fd_noise <- NA_real_
  tol_source <- NA_character_; weak_threshold <- NA_real_
  unidentified_params <- character(0); weak_params <- character(0)
  sv_class <- character(0); param_band_contrib <- NULL

  if (status == "ok") {
    J <- .numerical_jacobian(spec_vec, theta, eps = eps)
    if (!all(is.finite(J))) status <- "non_finite"
  }

  if (status == "ok") {
    colnames(J) <- param_names
    gram_matrix <- crossprod(J)
    dimnames(gram_matrix) <- list(param_names, param_names)

    ## FD error estimate: same Jacobian at twice the step.
    J2 <- NULL
    if (is.null(tol_rank)) {
      J2 <- .numerical_jacobian(spec_vec, theta, eps = 2 * eps)
      if (!all(is.finite(J2))) J2 <- NULL
    }
    col_sig <- sqrt(colSums(J^2))
    col_err <- if (!is.null(J2)) sqrt(colSums((J - J2)^2)) else rep(0, n_par)
    ## A column no larger than 10x its own FD error is numerically zero:
    ## the parameter does not move the spectral density at all.
    zero_col <- col_sig <= pmax(10 * col_err, .Machine$double.eps * max(col_sig, 1e-300))
    scl <- ifelse(zero_col, 1, col_sig)
    Je <- sweep(J, 2, scl, "/")
    Je[, zero_col] <- 0

    sv <- svd(Je, nu = 0, nv = n_par)
    singular_values <- c(sv$d, rep(0, n_par - length(sv$d)))[seq_len(n_par)]
    names(singular_values) <- paste0("sv_", seq_len(n_par))
    V <- sv$v
    dimnames(V) <- list(param_names, names(singular_values))
    sv_max <- max(singular_values)

    tol_machine <- max(dim(Je)) * sv_max * .Machine$double.eps
    if (!is.null(tol_rank)) {
      tol <- max(tol_machine, tol_rank * sv_max); tol_source <- "user"
    } else if (!is.null(J2)) {
      E <- sweep(J - J2, 2, scl, "/"); E[, zero_col] <- 0
      fd_noise <- max(svd(E, nu = 0, nv = 0)$d)
      tol <- max(tol_machine, 10 * fd_noise); tol_source <- "finite-difference"
    } else {
      tol <- tol_machine; tol_source <- "machine"
    }
    weak_threshold <- max(weak_rel * sv_max, tol)
    rank_J    <- sum(singular_values > tol)
    full_rank <- rank_J == n_par
    sv_class  <- ifelse(singular_values <= tol, "Unidentified",
                 ifelse(singular_values < weak_threshold, "Weak", "Identified"))

    loaders <- function(k) {
      v <- abs(V[, k])
      if (max(v) <= 0) character(0) else param_names[v >= 0.5 * max(v)]
    }
    unidentified_params <- unique(c(param_names[zero_col],
      unlist(lapply(which(sv_class == "Unidentified"), loaders))))
    weak_params <- setdiff(unique(unlist(lapply(which(sv_class == "Weak"), loaders))),
                           unidentified_params)
    unidentified_params <- as.character(unidentified_params)
    weak_params <- as.character(weak_params)

    ## Share of each G_kk from each band.
    band_rows <- band_of_freq[row_freq]
    param_band_contrib <- matrix(NA_real_, n_par, 3L,
                                 dimnames = list(param_names, band_names))
    for (b in band_names)
      param_band_contrib[, b] <- colSums(J[band_rows == b, , drop = FALSE]^2)
    tot <- rowSums(param_band_contrib)
    param_band_contrib <- param_band_contrib / ifelse(tot > 0, tot, NA) * 100
    param_band_contrib[zero_col, ] <- NA
  }

  pass <- if (status == "ok") full_rank else NA

  ## ---- Plots ----
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    auto <- vapply(spectral_list, function(S) Re(diag(S)) / (2 * pi), numeric(n_obs))
    auto <- matrix(auto, nrow = n_obs)
    spec_df <- data.frame(
      Frequency  = rep(freq_grid, each = n_obs),
      Value      = as.vector(auto),
      Observable = factor(rep(obs_labels, times = n_freq), levels = obs_labels)
    )
    ## Drop points more than 8 decades below the peak (e.g. the exact zero of
    ## a growth-rate observable at w = 0) so they do not flatten the axis.
    spec_df <- spec_df[is.finite(spec_df$Value) & spec_df$Value > 0, ]
    if (nrow(spec_df) > 0L)
      spec_df <- spec_df[spec_df$Value >= 1e-8 * max(spec_df$Value), ]
    brk_df <- data.frame(x = band_breaks, lab = c("32-period cycle", "6-period cycle"))
    p_spec <- ggplot2::ggplot(spec_df,
                              ggplot2::aes(x = Frequency, y = Value, colour = Observable)) +
      ggplot2::geom_vline(data = brk_df, ggplot2::aes(xintercept = x),
                          linetype = "dashed", colour = dynhr_colours$grey,
                          linewidth = 0.4) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::scale_x_continuous(
        limits = c(0, pi), expand = c(0, 0),
        breaks = c(0, pi / 4, pi / 2, 3 * pi / 4, pi),
        labels = c("0", "pi/4", "pi/2", "3pi/4", "pi"),
        sec.axis = ggplot2::sec_axis(~ ., name = "Cycle length (periods)",
          breaks = c(2 * pi / 32, 2 * pi / 16, 2 * pi / 8, 2 * pi / 6, 2 * pi / 4, pi),
          labels = c("32", "16", "8", "6", "4", "2"))) +
      ggplot2::scale_y_log10() +
      scale_colour_dynhr_vibrant(name = NULL) +
      theme_dynhr_diagnostic() +
      ggplot2::labs(
        title = "D23: Auto-spectra of the observables",
        subtitle = "f(w) = H Sigma_e H' / (2 pi) at theta (values > 8 decades below the peak omitted). Dashed: band edges",
        x = "Frequency w (radians per period)", y = "Spectral density f(w) (log scale)"
      )
    plots$spectral_density <- .apply_meta(p_spec, meta)

    if (length(singular_values) > 0L) {
      pos  <- singular_values[singular_values > 0]
      cand <- c(pos, tol, weak_threshold)
      cand <- cand[is.finite(cand) & cand > 0]
      floor_val <- 10^(floor(log10(if (length(cand)) min(cand) else 1e-16)) - 1)
      lead <- vapply(seq_len(n_par), function(k) {
        v <- abs(V[, k]); if (max(v) <= 0) "" else param_names[which.max(v)]
      }, character(1))
      sv_df <- data.frame(
        index = factor(seq_len(n_par), labels = sprintf("%d\n%s", seq_len(n_par), lead)),
        value = pmax(singular_values, floor_val),
        class = factor(sv_class, levels = c("Identified", "Weak", "Unidentified")),
        zero  = singular_values <= 0
      )
      ref_df <- data.frame(
        y = c(tol, weak_threshold),
        what = factor(c(sprintf("rank tolerance (%s) = %.1e", tol_source, tol),
                        sprintf("weak threshold = %.0e x largest", weak_rel)))
      )
      p_sv <- ggplot2::ggplot(sv_df, ggplot2::aes(x = index, y = value, colour = class)) +
        ggplot2::geom_segment(ggplot2::aes(xend = index, y = floor_val, yend = value),
                              linewidth = 0.8) +
        ggplot2::geom_point(ggplot2::aes(shape = zero), size = 2.6) +
        ggplot2::geom_hline(data = ref_df,
                            ggplot2::aes(yintercept = y, linetype = what),
                            colour = dynhr_colours$grey, linewidth = 0.5) +
        ggplot2::scale_colour_manual(values = c("Identified" = dynhr_colours$mid_blue,
                                                "Weak" = dynhr_colours$orange,
                                                "Unidentified" = dynhr_colours$red),
                                     name = NULL) +
        ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                                    breaks = "TRUE",
                                    labels = c(`TRUE` = "exactly 0 (drawn at floor)"),
                                    name = NULL) +
        ggplot2::scale_linetype_manual(values = c("dashed", "dotted"), name = NULL) +
        ggplot2::scale_y_log10() +
        theme_dynhr_diagnostic() +
        ggplot2::theme(legend.box = "vertical") +
        ggplot2::labs(
          title = "D23: Singular values of the spectral Jacobian (Qu-Tkachenko G)",
          subtitle = sprintf("Columns scaled to unit norm. Rank = %d / %d. Axis label: parameter with the largest loading",
                             rank_J, n_par),
          x = "Singular value (index / dominant parameter)",
          y = "Singular value (log scale)"
        )
      plots$singular_values <- .apply_meta(p_sv, meta)
    }

    if (!is.null(param_band_contrib)) {
      band_df <- data.frame(
        Parameter = factor(rep(param_names, 3L), levels = rev(param_names)),
        Band = factor(rep(c("low\n(> 32 periods)", "business cycle\n(6-32 periods)",
                            "high\n(< 6 periods)"), each = n_par),
                      levels = c("low\n(> 32 periods)", "business cycle\n(6-32 periods)",
                                 "high\n(< 6 periods)")),
        Contribution = as.vector(param_band_contrib)
      )
      band_df$label <- ifelse(is.na(band_df$Contribution), "n/a",
                              sprintf("%.0f%%", band_df$Contribution))
      p_band <- ggplot2::ggplot(band_df,
                                ggplot2::aes(x = Band, y = Parameter, fill = Contribution)) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        ggplot2::geom_text(ggplot2::aes(label = label,
                                        colour = !is.na(Contribution) & Contribution > 50),
                           size = 3.2, show.legend = FALSE) +
        ggplot2::scale_colour_manual(values = c(`TRUE` = "#000000", `FALSE` = "#FFFFFF")) +
        scale_fill_dynhr_cividis(name = "% of G_kk", limits = c(0, 100)) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title = "D23: Where in the spectrum each parameter is identified",
          subtitle = "Share of the diagonal Gram entry G_kk contributed by each frequency band (rows sum to 100%)",
          x = "Frequency band", y = NULL
        )
      plots$band_contribution <- .apply_meta(p_band, meta)
    }
  }

  ## ---- Summary ----
  badge <- if (isTRUE(pass)) "PASS" else if (identical(pass, FALSE)) "FAIL" else "INFO"
  detail <- switch(status,
    ok = sprintf("G(theta) rank = %d / %d (tolerance %.1e, %s)%s%s.",
                 rank_J, n_par, tol, tol_source,
                 if (length(unidentified_params))
                   paste0("; not identified: ", paste(unidentified_params, collapse = ", "))
                 else "",
                 if (length(weak_params))
                   paste0("; weak: ", paste(weak_params, collapse = ", "))
                 else ""),
    no_theta = "No theta supplied: spectral density only, no rank check.",
    no_solve_fn = paste0("Rank NOT checked: needs a model_solve_fn returning the ",
                         "state space at theta (ghx/ghu[/obs_mat/D_mat/Sigma_e])."),
    non_finite = paste0("Rank NOT checked: the spectral density is non-finite at ",
                        "a perturbed theta (unit root or infeasible solve)."))
  action <- switch(status,
    ok = if (isTRUE(full_rank))
      "Parameters are locally identified from the spectral density of the observables."
    else "Rank deficient: some parameter combinations leave the observables' spectrum unchanged; fix or reparameterise them, or add observables.",
    no_theta = "Supply theta and a state-space model_solve_fn for the rank check.",
    no_solve_fn = "Supply a state-space model_solve_fn (run_all_diagnostics builds one from model + compiled + obs_names).",
    non_finite = "Evaluate at a stationary theta.")

  summary_text <- sprintf(
    "D23 Spectral identification (Qu-Tkachenko 2012): %s. %s Gram over %d frequencies on [0, pi] (trapezoid, doubled for [-pi, pi]); %d observable(s).",
    badge, detail, n_freq, n_obs)

  .make_result(
    result = list(
      gram_matrix         = gram_matrix,
      singular_values     = singular_values,
      rank                = rank_J,
      full_rank           = full_rank,
      rank_tolerance      = tol,
      tolerance_source    = tol_source,
      fd_error            = fd_noise,
      unidentified_params = unidentified_params,
      weak_params         = weak_params,
      status              = status,
      frequency_grid      = freq_grid,
      spectral_density    = spectral_list,
      band_contribution   = param_band_contrib,
      Sigma_e             = Sigma_e,
      state_space         = list(T = ss$TT, R = ss$RR, Z = ss$ZZ, D = ss$DD)
    ),
    pass    = pass,
    plots   = plots,
    summary = summary_text,
    llm_summary = paste(c(
      sprintf("D23 | Spectral Identification (Qu-Tkachenko 2012) | %s", badge),
      sprintf("  n_freq=%d n_obs=%d n_par=%d rank=%s", n_freq, n_obs, n_par,
              if (is.na(rank_J)) "NA" else as.character(rank_J)),
      if (length(singular_values))
        sprintf("  smallest_sv (col-scaled): %s", paste(sprintf("%.3e",
                head(sort(singular_values), 5)), collapse = ", ")),
      if (length(unidentified_params))
        sprintf("  unidentified: %s", paste(unidentified_params, collapse = ", ")),
      if (length(weak_params))
        sprintf("  weak: %s", paste(weak_params, collapse = ", ")),
      sprintf("  action: %s", action)
    ), collapse = "\n")
  )
}


## Normalise a D23 state-space list to lagged-convention TT/RR/ZZ/DD.
## `ghx` may be the square state transition, or the n_endo x n_state decision
## rule with `state_idx` selecting the state rows.
.d23_state_space <- function(dr) {
  ghx <- as.matrix(dr$ghx)
  ghu <- as.matrix(dr$ghu)
  if (nrow(ghu) != nrow(ghx))
    .dynhr_abort(sprintf("D23: ghu has %d rows, ghx has %d.", nrow(ghu), nrow(ghx)))
  n_state <- ncol(ghx)
  state_idx <- dr$state_idx %||% seq_len(n_state)
  if (length(state_idx) != n_state || any(state_idx > nrow(ghx)))
    .dynhr_abort("D23: state_idx does not select ncol(ghx) rows of ghx.")
  TT <- ghx[state_idx, , drop = FALSE]
  RR <- ghu[state_idx, , drop = FALSE]
  ZZ <- if (!is.null(dr$obs_mat)) as.matrix(dr$obs_mat)
        else if (!is.null(dr$Z)) as.matrix(dr$Z)
        else diag(n_state)
  if (ncol(ZZ) != n_state)
    .dynhr_abort(sprintf("D23: observation matrix has %d columns, expected %d states.",
                         ncol(ZZ), n_state))
  DD <- if (!is.null(dr$D_mat)) as.matrix(dr$D_mat)
        else matrix(0, nrow(ZZ), ncol(RR))
  if (!identical(dim(DD), c(nrow(ZZ), ncol(RR))))
    .dynhr_abort(sprintf("D23: D_mat is %s, expected %d x %d.",
                         paste(dim(DD), collapse = " x "), nrow(ZZ), ncol(RR)))
  Sig <- if (!is.null(dr$Sigma_e)) as.matrix(dr$Sigma_e) else NULL
  list(TT = TT, RR = RR, ZZ = ZZ, DD = DD, Sigma_e = Sig)
}


## theta -> observables' state space (the D23 model_solve_fn) for a parsed +
## compiled model. Re-solves the steady state and the first-order solution at
## theta and re-derives Sigma_e from params (so estimated `stderr` entries
## move the spectrum). Returns NULL where the model does not solve.
.d23_state_space_solve_fn <- function(model, compiled, obs_names) {
  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence
  sys_cache <- cache_system_structure(compiled)
  state <- new.env(parent = emptyenv())
  state$ss_warm <- NULL
  function(theta) {
    sol <- .solve_dr_for_theta(model, compiled, sys_cache, theta, state,
                               lik_init = "diffuse")
    if (is.null(sol)) return(NULL)
    d <- sol$dr
    d$Sigma_e <- NULL
    sp <- build_dsge_state_space(model, d, obs_names, verbose = FALSE,
                                 params = sol$params)
    list(ghx = sp$T_mat, ghu = sp$R_mat, obs_mat = sp$Z_mat,
         D_mat = sp$D_mat, Sigma_e = sp$Sigma_e, obs_names = obs_names)
  }
}
