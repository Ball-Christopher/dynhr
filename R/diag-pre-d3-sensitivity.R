## R/diag-pre-d3-sensitivity.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D3 Morris sensitivity screening
## --------------------------------------------------------------------------

#' D3. Sensitivity analysis -- Morris screening (Morris 1991)
#'
#' Implements the Morris method of elementary effects to identify which
#' parameters most influence which moments. Returns mu* (absolute mean
#' effect) and sigma (standard deviation of effect) for each
#' parameter-moment pair.
#'
#' @param model_solve_fn  Function: theta -> named numeric vector of moments
#' @param theta           Numeric vector -- baseline parameter values
#' @param param_bounds    Matrix (n_par x 2) -- lower/upper bounds for screening
#' @param param_names     Character vector
#' @param moment_names    Character vector
#' @param n_paths         Number of Morris trajectories (default 10).  Each
#'   trajectory costs \code{n_par + 1} model solves, so D3 is the dominant
#'   solver cost in the suite; 10 trajectories is a standard screening default
#'   (Campolongo et al. 2007).  Raise it (e.g. 20-50) for less noisy
#'   \eqn{\mu^*} estimates when solves are cheap.
#' @param n_levels        Number of grid levels (default 4)
#' @return dynhr_diagnostic list with heatmap of mu*/sigma
#' @references Morris, M. D. (1991). Factorial sampling plans for preliminary
#'   computational experiments. \emph{Technometrics}, 33(2), 161-174.
#'   Saltelli, A., Ratto, M., Andres, T., Campolongo, F., Cariboni, J., Gatelli, D.,
#'   Saisana, M., & Tarantola, S. (2008). \emph{Global Sensitivity Analysis: The
#'   Primer}. John Wiley & Sons.
#' @noRd
d3_sensitivity_morris <- function(model_solve_fn,
                                  theta,
                                  param_bounds,
                                  param_names  = NULL,
                                  moment_names = NULL,
                                  n_paths  = 10,
                                  n_levels = 4,
                                  meta     = NULL) {

    n_par <- length(theta)
    if (is.null(param_names))  param_names  <- paste0("theta_", seq_len(n_par))

    lb <- param_bounds[, 1]
    ub <- param_bounds[, 2]

    # Grid step size
    delta <- n_levels / (2 * (n_levels - 1))

    # Normalise theta to [0, 1]
    theta_norm <- (theta - lb) / (ub - lb)

    # Storage for elementary effects
    f0 <- model_solve_fn(theta)
    n_mom <- length(f0)
    if (is.null(moment_names)) {
      moment_names <- if (!is.null(names(f0))) names(f0) else paste0("m_", seq_len(n_mom))
    }
    EE <- array(NA_real_, dim = c(n_paths, n_par, n_mom))

    for (r in seq_len(n_paths)) {
      # Random starting point on the grid
      x0 <- runif(n_par)
      x0 <- round(x0 * (n_levels - 1)) / (n_levels - 1)

      # Random permutation of parameter indices
      perm <- sample(n_par)

      x_current <- x0
      # Evaluate at starting position (unnormalise). IMPORTANT: model_solve_fn
      # maps parameters BY NAME, so the theta vectors handed to it must carry
      # param_names -- otherwise every perturbation silently reuses the baseline
      # calibration and all elementary effects collapse to zero.
      theta_curr <- stats::setNames(lb + x_current * (ub - lb), param_names)
      f_current  <- model_solve_fn(theta_curr)

      for (j in seq_along(perm)) {
        k <- perm[j]

        # Perturb parameter k by +/- delta
        x_next <- x_current
        if (x_current[k] + delta <= 1) {
          x_next[k] <- x_current[k] + delta
        } else {
          x_next[k] <- x_current[k] - delta
        }

        theta_next <- stats::setNames(lb + x_next * (ub - lb), param_names)
        f_next <- model_solve_fn(theta_next)

        # Elementary effect in NORMALISED [0,1] parameter space: dividing by the
        # normalised step `delta` (not the original-unit range) keeps mu* directly
        # comparable across parameters with different prior widths.
        ee <- (f_next - f_current) / delta
        EE[r, k, ] <- ee

        x_current <- x_next
        f_current <- f_next
      }
    }

    # Compute mu* (mean of |EE|) and sigma (sd of EE) per parameter-moment pair
    mu_star <- apply(abs(EE), c(2, 3), mean, na.rm = TRUE)
    sigma_ee <- apply(EE, c(2, 3), sd, na.rm = TRUE)

    colnames(mu_star) <- moment_names
    rownames(mu_star) <- param_names
    colnames(sigma_ee) <- moment_names
    rownames(sigma_ee) <- param_names

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    # Heatmap of mu*
    mu_long <- reshape2::melt(mu_star)
    colnames(mu_long) <- c("Parameter", "Moment", "mu_star")

    # Heatmap: log-scale fill so mid-range values are visible
    # (linear scale crushes everything except the maximum cell)
    mu_long$mu_star_plot <- pmax(mu_long$mu_star, 1e-6)  # floor for log scale
    p_mu <- ggplot2::ggplot(
      mu_long, ggplot2::aes(x = Moment, y = Parameter, fill = mu_star_plot)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
      scale_fill_dynhr_cividis(
        name   = expression(mu * "*"),
        trans  = "log10",
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        guide  = ggplot2::guide_colorbar(
          barwidth     = ggplot2::unit(0.4, "cm"),
          barheight    = ggplot2::unit(4, "cm"),
          label.theme  = ggplot2::element_text(size = 8)
        )
      ) +
      theme_dynhr_diagnostic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(
        title = "D3: Morris screening -- mean absolute elementary effect (log scale)",
        subtitle = sprintf("%d paths, %d levels. Colour shows log10(mu*) for contrast.", n_paths, n_levels),
        x = NULL, y = NULL
      )
    plots$morris_mu_star <- .apply_meta(p_mu, meta)

    # mu* vs sigma scatter for each moment (faceted)
    scatter_df <- data.frame(
      Parameter = rep(param_names, n_mom),
      Moment    = rep(moment_names, each = n_par),
      mu_star   = as.vector(mu_star),
      sigma     = as.vector(sigma_ee)
    )

    # Label only the top-3 parameters by |mu*| per facet to reduce clutter.
    # If ggrepel is available use repel labels; otherwise fall back to geom_text.
    top3_df <- do.call(rbind, lapply(unique(scatter_df$Moment), function(m) {
      sub <- scatter_df[scatter_df$Moment == m, , drop = FALSE]
      sub[order(sub$mu_star, decreasing = TRUE)[seq_len(min(3L, nrow(sub)))], ]
    }))

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      label_layer <- ggrepel::geom_text_repel(
        data       = top3_df,
        ggplot2::aes(label = Parameter),
        size       = 2.2,
        colour     = dynhr_colours$dark_blue,
        max.overlaps = 10,
        seed       = 42
      )
    } else {
      label_layer <- ggplot2::geom_text(
        data   = top3_df,
        ggplot2::aes(label = Parameter),
        size   = 2.2,
        vjust  = -0.5,
        colour = dynhr_colours$dark_blue
      )
    }

    p_sc <- ggplot2::ggplot(
      scatter_df, ggplot2::aes(x = mu_star, y = sigma, label = Parameter)
    ) +
      ggplot2::geom_point(colour = dynhr_colours$mid_blue, size = 1.5) +
      ggplot2::geom_abline(slope = 0.5, intercept = 0, linetype = "dashed",
                           colour = dynhr_colours$grey) +
      label_layer +
      ggplot2::facet_wrap(~ Moment, scales = "free", ncol = 4L) +
      ggplot2::scale_x_continuous(n.breaks = 3L) +
      ggplot2::scale_y_continuous(n.breaks = 3L) +
      theme_dynhr_compact() +
      ggplot2::theme(
        axis.title.y = ggplot2::element_text(
          size   = ggplot2::rel(0.85),
          angle  = 90,
          margin = ggplot2::margin(r = 6)
        ),
        axis.text.y  = ggplot2::element_text(size = ggplot2::rel(0.7)),
        axis.ticks.y = ggplot2::element_line(colour = "#222222", linewidth = 0.2),
        axis.line.y  = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        title = "D3: Morris screening -- parameter sensitivity",
        subtitle = "Points above the diagonal indicate non-linear / interaction effects. Top-3 by |mu*| labelled per facet.",
        x = "mu* (mean absolute elementary effect)",
        y = "sigma (sd of elementary effects)"
      )
    plots$morris_scatter <- .apply_meta(p_sc, meta)

    }  # end requireNamespace guard

    .make_result(
      result  = list(mu_star = mu_star, sigma = sigma_ee, EE = EE),
      pass    = NA,  # Informational diagnostic
      plots   = plots,
      summary = sprintf(
        "D3 Morris screening: %d paths, %d levels, %d params, %d moments. Most influential parameters per moment: %s",
        n_paths, n_levels, n_par, n_mom,
        paste(apply(mu_star, 2, function(col) param_names[which.max(col)]), collapse = ", ")
      ),
      llm_summary = {
        # mu_star is a matrix (n_par x n_mom) or a vector -- handle both
        if (is.matrix(mu_star)) {
          overall_mu <- rowMeans(mu_star)
        } else {
          overall_mu <- mu_star
        }
        top5    <- head(sort(overall_mu, decreasing = TRUE), 5)
        low5    <- head(sort(overall_mu), 5)
        negligi <- names(low5[low5 < 0.05 * max(overall_mu)])
        paste(c(
          "D3 | Morris Sensitivity | INFO",
          sprintf("  params=%d moments=%d n_paths=%d",
                  n_par, n_mom, n_paths),
          sprintf("  high_sensitivity: %s",
                  paste(sprintf("%s=%.2f", names(top5), top5), collapse = ", ")),
          if (length(negligi) > 0)
            sprintf("  negligible_sensitivity: %s (may not be identifiable from these moments)",
                    paste(negligi, collapse = ", ")),
          sprintf("  action: %s",
                  if (length(negligi) > 0)
                    sprintf("%s have negligible sensitivity. Consider calibrating or removing from estimation.",
                            paste(head(negligi, 3), collapse = ", "))
                  else
                    "All parameters have meaningful sensitivity to chosen moments.")
        ), collapse = "\n")
      }
    )
}


#' D22. Observable informativeness ranking (Iskrev 2019 inspired)
#'
#' Uses the moment Jacobian to apportion parameter sensitivity across
#' observables. Moments are mapped to observables using prefix matching on
#' `moment_names` (e.g. "y_sd", "y_acf1" -> observable "y").
#'
#' @param model_solve_fn Function: theta -> named numeric vector of moments.
#' @param theta Numeric parameter vector.
#' @param obs_names Character vector of observable names.
#' @param param_names Optional parameter names.
#' @param moment_names Optional moment names.
#' @param eps Step size for numerical derivatives.
#' @param meta Optional metadata list (passed to \code{.apply_meta} for plot
#'   annotation).
#' @return dynhr_diagnostic object.
#' @noRd
d22_observable_informativeness <- function(model_solve_fn,
                                           theta,
                                           obs_names,
                                           param_names  = NULL,
                                           moment_names = NULL,
                                           eps = 1e-5,
                                           meta = NULL) {
  
    n_par <- length(theta)
    if (is.null(param_names)) param_names <- names(theta) %||% paste0("theta_", seq_len(n_par))
    f0 <- model_solve_fn(theta)
    if (is.null(moment_names)) moment_names <- names(f0) %||% paste0("m_", seq_len(length(f0)))
    if (is.null(obs_names) || length(obs_names) == 0L) {
      stop("obs_names must be supplied for D22 observable informativeness.")
    }

    J <- .numerical_jacobian(model_solve_fn, theta, eps = eps)
    if (nrow(J) != length(moment_names)) {
      warning(sprintf(
        "d22: moment_names length (%d) != Jacobian rows (%d). Using generic labels.",
        length(moment_names), nrow(J)))
      moment_names <- paste0("m_", seq_len(nrow(J)))
    }
    # Guard: ensure param_names length matches Jacobian column count
    if (ncol(J) != length(param_names)) {
      warning(sprintf(
        "d22: param_names length (%d) != Jacobian columns (%d). Using generic param labels.",
        length(param_names), ncol(J)))
      param_names <- paste0("theta_", seq_len(ncol(J)))
    }
    colnames(J) <- param_names
    rownames(J) <- moment_names

    moment_to_obs <- vapply(moment_names, function(mn) {
      hit <- obs_names[startsWith(mn, paste0(obs_names, "_"))]
      if (length(hit) > 0) return(hit[1])
      hit2 <- obs_names[vapply(obs_names, function(v) grepl(paste0("(^|_)", v, "($|_)"), mn), logical(1))]
      if (length(hit2) > 0) return(hit2[1])
      NA_character_
    }, character(1))

    valid <- !is.na(moment_to_obs)
    if (!any(valid)) {
      stop("Could not map moment_names to observables. Use names like 'obs_sd' or 'obs_acf1'.")
    }

    obs_set <- unique(moment_to_obs[valid])
    info_mat <- matrix(0, nrow = length(obs_set), ncol = n_par,
                       dimnames = list(obs_set, param_names))
    for (ob in obs_set) {
      idx <- which(moment_to_obs == ob)
      if (length(idx) == 1) {
        info_mat[ob, ] <- J[idx, ]^2
      } else {
        info_mat[ob, ] <- colSums(J[idx, , drop = FALSE]^2)
      }
    }
    shares <- sweep(info_mat, 2, colSums(info_mat) + 1e-16, "/")
    dominant_obs <- apply(shares, 2, function(x) names(which.max(x)))

    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      long <- reshape2::melt(shares)
      colnames(long) <- c("Observable", "Parameter", "Share")
      plots$informativeness_heatmap <- ggplot2::ggplot(
        long, ggplot2::aes(x = Parameter, y = Observable, fill = Share)
      ) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        scale_fill_dynhr_cividis(name = "Share") +
        theme_dynhr_diagnostic() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            angle = 45, hjust = 1, vjust = 1, size = ggplot2::rel(0.8)
          )
        ) +
        ggplot2::labs(
          title = "D22: Observable informativeness by parameter",
          subtitle = "Shares of Jacobian sensitivity attributed to each observable",
          x = NULL, y = NULL
        )
    }

    .make_result(
      result = list(
        jacobian = J,
        moment_to_observable = moment_to_obs,
        informativeness = shares,
        dominant_observable = dominant_obs
      ),
      pass = NA,
      plots = plots,
      summary = sprintf(
        "D22 Observable informativeness: mapped %d/%d moments to %d observables. Dominant observable by parameter: %s",
        sum(valid), length(valid), nrow(shares),
        paste(sprintf("%s->%s", names(dominant_obs), dominant_obs), collapse = ", ")
      ),
      llm_summary = {
        # Shorten dominant_observable: show top-3 by share + count of the rest
        dom_sorted <- sort(table(dominant_obs), decreasing = TRUE)
        top3_dom   <- head(names(dom_sorted), 3L)
        top3_str   <- paste(
          sprintf("%s(%d)", top3_dom, as.integer(dom_sorted[top3_dom])),
          collapse = ", "
        )
        n_rest <- length(dom_sorted) - length(top3_dom)
        dom_str <- if (n_rest > 0)
          sprintf("%s + %d more", top3_str, n_rest)
        else
          top3_str
        paste(c(
          "D22 | Observable Informativeness | INFO",
          sprintf("  params=%d moments_mapped=%d/%d observables=%d",
                  n_par, sum(valid), length(valid), nrow(shares)),
          sprintf("  dominant_observable (top-3): %s", dom_str),
          "  action: Use low-share observables/parameters map to redesign measurement set before re-estimation."
        ), collapse = "\n")
      }
    )
}
