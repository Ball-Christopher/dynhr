## R/diag-pre-d3-sensitivity.R
## --------------------------------------------------------------------------
## Phase-3 split from diagnostics-monolith.R.
##
## D3 Morris sensitivity screening
## --------------------------------------------------------------------------

#' D3. Sensitivity analysis -- Morris screening (Morris 1991)
#'
#' Implements the Morris method of elementary effects to identify which
#' parameters most influence which moments.  Each parameter is mapped to
#' \eqn{[0, 1]} via \code{param_bounds}; on a \code{n_levels}-level grid with
#' step \eqn{\Delta = p / (2(p - 1))} the elementary effect of parameter
#' \eqn{k} on moment \eqn{m} is
#' \deqn{EE_{km} = [f_m(x + s\Delta e_k) - f_m(x)] / (s\Delta),}
#' with \eqn{s = \pm 1} the direction of the step (a step DOWN is divided by
#' \eqn{-\Delta}, so a linear moment has identical EEs and \eqn{\sigma = 0}).
#' Returns \eqn{\mu^*} (mean \eqn{|EE|}, Campolongo et al. 2007) and
#' \eqn{\sigma} (sd of EE) per parameter-moment pair, in units of the moment
#' per full parameter range.
#'
#' Because moments carry different units, cross-moment rankings use the
#' RELATIVE \eqn{\mu^*}: each moment's column divided by its largest entry.
#' A parameter is flagged negligible only when its relative \eqn{\mu^*} is
#' below \code{negligible_tol} for EVERY moment.
#'
#' @param model_solve_fn  Function: named theta -> numeric vector of moments.
#'   A failed solve should RETURN non-finite values (e.g. \code{NA}) rather
#'   than error; those yield \code{NA} elementary effects (counted in
#'   \code{$result$n_failed}).  An error propagates.
#' @param theta           Numeric vector -- baseline parameter values (only
#'   used for names and the moment count; Morris screening is global over
#'   \code{param_bounds}).
#' @param param_bounds    Matrix (n_par x 2) -- finite lower/upper bounds.  If
#'   it has row names they are matched to \code{param_names}.
#' @param param_names     Character vector (default \code{names(theta)}).
#' @param moment_names    Character vector (default \code{names(f(theta))}).
#' @param n_paths         Number of Morris trajectories (default 10, >= 2).
#'   Each trajectory costs \code{n_par + 1} model solves, so D3 is the
#'   dominant solver cost in the suite; 10 trajectories is a standard
#'   screening default (Campolongo et al. 2007).  Raise it (e.g. 20-50) for
#'   less noisy \eqn{\mu^*} estimates when solves are cheap.
#' @param n_levels        Number of grid levels (default 4; must be even so
#'   that every step stays on the grid).
#' @param negligible_tol  Relative-\eqn{\mu^*} threshold for the negligible
#'   flag (default 0.05).
#' @param meta            Optional plot-provenance list.
#' @return dynhr_diagnostic list with \code{result$mu_star},
#'   \code{result$sigma}, \code{result$mu_star_rel}, \code{result$EE}
#'   (n_paths x n_par x n_mom), \code{result$negligible},
#'   \code{result$n_failed}; plots \code{morris_mu_star} (relative
#'   \eqn{\mu^*} heatmap) and \code{morris_scatter}.
#' @references Morris, M. D. (1991). Factorial sampling plans for preliminary
#'   computational experiments. \emph{Technometrics}, 33(2), 161-174.
#'   Campolongo, F., Cariboni, J., & Saltelli, A. (2007). An effective
#'   screening design for sensitivity analysis of large models.
#'   \emph{Environmental Modelling & Software}, 22(10), 1509-1518.
#'   Saltelli, A., Ratto, M., Andres, T., Campolongo, F., Cariboni, J., Gatelli, D.,
#'   Saisana, M., & Tarantola, S. (2008). \emph{Global Sensitivity Analysis: The
#'   Primer}. John Wiley & Sons.
#' @noRd
d3_sensitivity_morris <- function(model_solve_fn,
                                  theta,
                                  param_bounds,
                                  param_names    = NULL,
                                  moment_names   = NULL,
                                  n_paths        = 10,
                                  n_levels       = 4,
                                  negligible_tol = 0.05,
                                  meta           = NULL) {

    n_par <- length(theta)
    if (is.null(param_names)) {
      param_names <- names(theta) %||% paste0("theta_", seq_len(n_par))
    }
    if (length(param_names) != n_par) {
      .dynhr_abort(sprintf("D3: param_names has length %d but theta has length %d.",
                           length(param_names), n_par))
    }
    n_paths <- as.integer(n_paths)
    if (length(n_paths) != 1L || is.na(n_paths) || n_paths < 2L) {
      .dynhr_abort("D3: n_paths must be an integer >= 2 (sigma needs >= 2 effects).")
    }
    n_levels <- as.integer(n_levels)
    if (length(n_levels) != 1L || is.na(n_levels) || n_levels < 2L || n_levels %% 2L != 0L) {
      .dynhr_abort("D3: n_levels must be an even integer >= 2 (odd p puts steps off the grid).")
    }

    param_bounds <- as.matrix(param_bounds)
    if (ncol(param_bounds) != 2L || nrow(param_bounds) != n_par) {
      .dynhr_abort(sprintf("D3: param_bounds must be a %d x 2 matrix (got %d x %d).",
                           n_par, nrow(param_bounds), ncol(param_bounds)))
    }
    if (!is.null(rownames(param_bounds))) {
      if (!setequal(rownames(param_bounds), param_names)) {
        .dynhr_abort("D3: rownames(param_bounds) do not match param_names.")
      }
      param_bounds <- param_bounds[param_names, , drop = FALSE]
    }
    lb <- as.numeric(param_bounds[, 1])
    ub <- as.numeric(param_bounds[, 2])
    if (any(!is.finite(lb)) || any(!is.finite(ub)) || any(ub <= lb)) {
      bad <- param_names[!is.finite(lb) | !is.finite(ub) | !(ub > lb)]
      .dynhr_abort(sprintf(
        "D3: param_bounds must be finite with lower < upper; offending: %s.",
        paste(bad, collapse = ", ")))
    }

    # Grid step size (Morris 1991): delta = p / (2 (p - 1)).
    delta <- n_levels / (2 * (n_levels - 1))

    # model_solve_fn maps parameters BY NAME, so every theta handed to it must
    # carry param_names -- otherwise each perturbation silently reuses the
    # baseline calibration and all elementary effects collapse to zero.
    f0 <- model_solve_fn(stats::setNames(as.numeric(theta), param_names))
    n_mom <- length(f0)
    if (n_mom == 0L) .dynhr_abort("D3: model_solve_fn(theta) returned no moments.")
    if (is.null(moment_names)) {
      moment_names <- names(f0) %||% paste0("m_", seq_len(n_mom))
    }
    if (length(moment_names) != n_mom) {
      .dynhr_abort(sprintf("D3: moment_names has length %d but model_solve_fn returns %d moments.",
                           length(moment_names), n_mom))
    }

    n_failed <- 0L
    checked_solve <- function(x) {
      out <- as.numeric(model_solve_fn(stats::setNames(lb + x * (ub - lb), param_names)))
      if (length(out) != n_mom || any(!is.finite(out))) {
        return(rep(NA_real_, n_mom))
      }
      out
    }

    EE <- array(NA_real_, dim = c(n_paths, n_par, n_mom),
                dimnames = list(NULL, param_names, moment_names))

    for (r in seq_len(n_paths)) {
      # Random starting point, uniform over the p grid levels {0, 1/(p-1), ..., 1}.
      x_current <- (sample.int(n_levels, n_par, replace = TRUE) - 1L) / (n_levels - 1L)
      f_current <- checked_solve(x_current)
      if (anyNA(f_current)) n_failed <- n_failed + 1L

      for (k in sample.int(n_par)) {
        # Step +delta if it stays inside [0, 1], otherwise -delta.
        step <- if (x_current[k] + delta <= 1 + 1e-12) delta else -delta
        x_next <- x_current
        x_next[k] <- x_current[k] + step
        f_next <- checked_solve(x_next)
        if (anyNA(f_next)) n_failed <- n_failed + 1L

        # Divide by the SIGNED step: a downward step is a forward difference
        # taken from the other end.  Normalised [0,1] units keep mu* comparable
        # across parameters with different bound widths.
        EE[r, k, ] <- (f_next - f_current) / step

        x_current <- x_next
        f_current <- f_next
      }
    }

    # mu* (mean |EE|) and sigma (sd EE) per parameter-moment pair.
    mu_star  <- apply(abs(EE), c(2, 3), mean, na.rm = TRUE)
    sigma_ee <- apply(EE, c(2, 3), stats::sd, na.rm = TRUE)
    mu_star[is.nan(mu_star)] <- NA_real_
    # Snap round-off: a linear moment gives identical EEs up to ~1e-16 relative
    # noise; report sigma = 0 rather than a meaningless 1e-16.
    snap <- !is.na(sigma_ee) & !is.na(mu_star) & sigma_ee <= 1e-8 * mu_star
    sigma_ee[snap] <- 0
    dim(mu_star) <- dim(sigma_ee) <- c(n_par, n_mom)
    dimnames(mu_star) <- dimnames(sigma_ee) <- list(param_names, moment_names)

    # Relative mu*: each moment scaled by its most influential parameter, so
    # moments with different units can be compared / aggregated.
    col_max <- apply(mu_star, 2, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
    mu_star_rel <- sweep(mu_star, 2, ifelse(is.finite(col_max) & col_max > 0, col_max, NA_real_), "/")
    max_rel <- apply(mu_star_rel, 1, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE))
    negligible <- param_names[!is.na(max_rel) & max_rel < negligible_tol]
    flat_moments <- moment_names[!is.na(col_max) & col_max == 0]

    top_by_moment <- vapply(seq_len(n_mom), function(j) {
      if (is.na(col_max[j]) || col_max[j] == 0) return(NA_character_)
      param_names[which.max(mu_star[, j])]
    }, character(1))

    # --- Plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {

    heat_df <- data.frame(
      Parameter = factor(rep(param_names, n_mom), levels = rev(param_names)),
      Moment    = factor(rep(moment_names, each = n_par), levels = moment_names),
      rel       = as.vector(mu_star_rel),
      mu_star   = as.vector(mu_star)
    )
    heat_df$label <- ifelse(is.na(heat_df$mu_star), "NA",
                     ifelse(is.na(heat_df$rel), "flat",
                            formatC(heat_df$mu_star, digits = 2, format = "g")))
    heat_df$txt_col <- ifelse(!is.na(heat_df$rel) & heat_df$rel > 0.6, "black", "white")
    show_txt <- n_par * n_mom <= 150L
    p_mu <- ggplot2::ggplot(
      heat_df, ggplot2::aes(x = Moment, y = Parameter, fill = rel)
    ) +
      ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
      scale_fill_dynhr_cividis(
        name   = "relative mu*",
        limits = c(0, 1),
        guide  = ggplot2::guide_colourbar(barwidth = ggplot2::unit(5, "cm"),
                                          barheight = ggplot2::unit(0.35, "cm"))
      ) +
      theme_dynhr_diagnostic() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(
        title = "D3: Morris screening -- relative mean absolute elementary effect",
        subtitle = sprintf(paste0(
          "Colour: mu* / max mu* within each moment (moments have different units); ",
          "text: raw mu*.\n%d paths, %d levels; negligible if relative mu* < %.2g for every moment."),
          n_paths, n_levels, negligible_tol),
        x = NULL, y = NULL
      )
    if (show_txt) {
      p_mu <- p_mu + ggplot2::geom_text(ggplot2::aes(label = label, colour = txt_col),
                                        size = 2.6, show.legend = FALSE) +
        ggplot2::scale_colour_identity()
    }
    plots$morris_mu_star <- .apply_meta(p_mu, meta)

    # mu* vs sigma scatter for each moment (faceted, in moment order)
    scatter_df <- data.frame(
      Parameter = rep(param_names, n_mom),
      Moment    = factor(rep(moment_names, each = n_par), levels = moment_names),
      mu_star   = as.vector(mu_star),
      sigma     = as.vector(sigma_ee)
    )
    # Drop non-finite pairs and flat moments (nothing to show: all points at 0).
    scatter_df <- scatter_df[is.finite(scatter_df$mu_star) & is.finite(scatter_df$sigma) &
                               !(scatter_df$Moment %in% flat_moments), , drop = FALSE]
    scatter_df$Moment <- droplevels(scatter_df$Moment)

    # Label only the top-3 parameters by mu* per facet to reduce clutter.
    top3_df <- do.call(rbind, lapply(split(scatter_df, scatter_df$Moment, drop = TRUE), function(sub) {
      sub <- sub[sub$mu_star > 0, , drop = FALSE]
      sub[order(sub$mu_star, decreasing = TRUE)[seq_len(min(3L, nrow(sub)))], , drop = FALSE]
    }))
    if (is.null(top3_df)) top3_df <- scatter_df[0, , drop = FALSE]

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      label_layer <- ggrepel::geom_text_repel(
        data       = top3_df,
        ggplot2::aes(label = Parameter),
        size       = 2.6,
        colour     = dynhr_colours$dark_blue,
        max.overlaps = 10,
        seed       = 42
      )
    } else {
      label_layer <- ggplot2::geom_text(
        data   = top3_df,
        ggplot2::aes(label = Parameter),
        size   = 2.6,
        vjust  = -0.5,
        colour = dynhr_colours$dark_blue
      )
    }

    # Invisible per-facet anchor so an all-linear facet (sigma == 0 for every
    # parameter) still gets a y-range that shows the sigma = 0.1 mu* line.
    blank_df <- do.call(rbind, lapply(split(scatter_df, scatter_df$Moment, drop = TRUE), function(sub)
      data.frame(Moment = sub$Moment[1], mu_star = 0,
                 sigma = max(sub$sigma, 0.15 * max(sub$mu_star)))))
    if (is.null(blank_df)) blank_df <- scatter_df[0, c("Moment", "mu_star", "sigma")]

    sc_sub <- paste0(
      "Dashed: sigma = mu* (above: non-linear / interaction effects); ",
      "dotted: sigma = 0.1 mu* (below: ~linear).\nTop-3 parameters by mu* labelled per facet.",
      if (length(flat_moments) > 0)
        sprintf(" Flat moments omitted: %s.", paste(flat_moments, collapse = ", ")) else "")

    if (nrow(scatter_df) == 0L) {
      p_sc <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0,
                          label = "No finite, non-flat elementary effects to plot.") +
        ggplot2::theme_void() +
        ggplot2::labs(title = "D3: Morris screening -- mu* vs sigma per moment",
                      subtitle = sc_sub)
    } else {
    p_sc <- ggplot2::ggplot(
      scatter_df, ggplot2::aes(x = mu_star, y = sigma)
    ) +
      ggplot2::geom_blank(data = blank_df) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                           colour = dynhr_colours$grey) +
      ggplot2::geom_abline(slope = 0.1, intercept = 0, linetype = "dotted",
                           colour = dynhr_colours$grey) +
      ggplot2::geom_point(colour = dynhr_colours$mid_blue, size = 1.8) +
      label_layer +
      ggplot2::facet_wrap(~ Moment, scales = "free", ncol = 4L) +
      ggplot2::scale_x_continuous(n.breaks = 3L, limits = c(0, NA)) +
      ggplot2::scale_y_continuous(n.breaks = 3L, limits = c(0, NA)) +
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
        title = "D3: Morris screening -- mu* vs sigma per moment",
        subtitle = sc_sub,
        x = "mu* (mean absolute elementary effect; moment units per full parameter range)",
        y = "sigma (sd of elementary effects)"
      )
    }
    plots$morris_scatter <- .apply_meta(p_sc, meta)

    }  # end requireNamespace guard

    overall_rel <- sort(max_rel, decreasing = TRUE)
    fail_txt <- if (n_failed > 0L)
      sprintf(" %d of %d solves failed (NA effects).", n_failed, n_paths * (n_par + 1L)) else ""

    .make_result(
      result  = list(mu_star = mu_star, sigma = sigma_ee, mu_star_rel = mu_star_rel,
                     EE = EE, negligible = negligible, n_failed = n_failed,
                     n_paths = n_paths, n_levels = n_levels),
      pass    = NA,  # Informational diagnostic
      plots   = plots,
      summary = sprintf(
        "D3 Morris screening: %d paths, %d levels, %d params, %d moments. Most influential parameter per moment: %s.%s%s",
        n_paths, n_levels, n_par, n_mom,
        paste(sprintf("%s->%s", moment_names, ifelse(is.na(top_by_moment), "none", top_by_moment)),
              collapse = ", "),
        if (length(negligible) > 0)
          sprintf(" Negligible for every moment: %s.", paste(negligible, collapse = ", ")) else "",
        fail_txt
      ),
      llm_summary = paste(c(
        "D3 | Morris Sensitivity | INFO",
        sprintf("  params=%d moments=%d n_paths=%d failed_solves=%d",
                n_par, n_mom, n_paths, n_failed),
        sprintf("  max_relative_mu* (top 5): %s",
                paste(sprintf("%s=%.2f", names(head(overall_rel, 5)), head(overall_rel, 5)),
                      collapse = ", ")),
        if (length(flat_moments) > 0)
          sprintf("  flat_moments (no parameter moves them): %s",
                  paste(flat_moments, collapse = ", ")),
        if (length(negligible) > 0)
          sprintf("  negligible_sensitivity: %s (relative mu* < %.2g for every moment; may not be identifiable from these moments)",
                  paste(negligible, collapse = ", "), negligible_tol),
        sprintf("  action: %s",
                if (length(negligible) > 0)
                  sprintf("%s have negligible sensitivity. Consider calibrating or removing from estimation.",
                          paste(head(negligible, 3), collapse = ", "))
                else
                  "All parameters have meaningful sensitivity to at least one chosen moment.")
      ), collapse = "\n")
    )
}


#' D22. Observable informativeness ranking (Iskrev 2019 inspired)
#'
#' Apportions each parameter's local moment sensitivity across observables.
#' With \eqn{J_{mk} = \partial f_m / \partial\theta_k} (central differences)
#' and \eqn{M_o} the moments attributed to observable \eqn{o}, the
#' informativeness share of \eqn{o} for parameter \eqn{k} is
#' \deqn{s_{ok} = \sum_{m \in M_o} J_{mk}^2 / \sum_{o'} \sum_{m \in M_{o'}} J_{mk}^2,}
#' i.e. observable \eqn{o}'s part of the diagonal of the unweighted
#' information \eqn{J'J}.  Each column sums to 1.  This is a MARGINAL
#' measure: it ignores collinearity between parameters and moment weighting,
#' so it depends on the units of the moments (rescaling an observable's
#' moments rescales its share).  Use D1/D20 for joint identification.
#'
#' \strong{The shares never gate a badge.}  D22 always reports \code{INFO}.
#' No source in the identification literature treats a sensitivity share as
#' identification-decisive, precisely because a share is scale- and
#' collinearity-blind: two nearly collinear observables can each show a large
#' share while jointly adding no rank.  Iskrev (2010) answers "is this
#' observable needed?" with a \emph{drop-one} experiment -- remove the
#' observable/moment block, recompute the rank and the collinearity indices --
#' and Andrle (2010) likewise frames observable informativeness as a
#' collinearity/identification-pattern question, not a variance-share number.
#' Dynare's \code{identification} command follows the same route (drop-one
#' rank and collinearity output), and ships no share statistic.  So the
#' PASS/FAIL verdict on identification belongs to the D1/D20 rank test; D22
#' is an exploratory ranking aid to read alongside it.
#'
#' Moments are attributed to observables BY NAME: \code{moment_names} are
#' split on \code{"_"} and an observable matches when its own
#' \code{"_"}-token sequence appears in the moment name (\code{"y_sd"},
#' \code{"var_y"}, \code{"acv1_y"} -> \code{"y"}).  The longest match wins
#' (\code{"y_gap_sd"} -> \code{"y_gap"}, not \code{"y"}).  A moment naming two
#' different observables (e.g. \code{"cov_y_pi"}) is a cross-moment: it is
#' excluded from the shares and listed in \code{result$cross_moments}.
#' Jacobian rows are matched to \code{moment_names} by \code{names(f(theta))}
#' and columns to \code{param_names} by \code{names(theta)} when those exist.
#'
#' A parameter that moves no mapped moment (or whose derivative is
#' non-finite) gets \code{NA} shares and an \code{NA} dominant observable
#' rather than an arbitrary one.
#'
#' @param model_solve_fn Function: theta -> named numeric vector of moments.
#' @param theta Numeric parameter vector (named, or labelled by
#'   \code{param_names} positionally).
#' @param obs_names Character vector of observable names.
#' @param param_names Optional parameter names (subset/order of
#'   \code{names(theta)}).
#' @param moment_names Optional moment names (subset/order of
#'   \code{names(model_solve_fn(theta))}; positional only when the moments
#'   are unnamed).
#' @param eps Step size for numerical derivatives.
#' @param meta Optional metadata list (passed to \code{.apply_meta} for plot
#'   annotation).
#' @return dynhr_diagnostic object; \code{result} holds \code{jacobian},
#'   \code{moment_to_observable}, \code{information} (obs x param sums of
#'   squared sensitivities), \code{informativeness} (shares),
#'   \code{dominant_observable}, \code{uninformed_params},
#'   \code{failed_params}, \code{cross_moments}, \code{unmapped_moments},
#'   \code{unmapped_observables}.
#' @references Iskrev, N. (2019). What to observe: Understanding the
#'   information content of observables in DSGE models. \emph{Working paper}.
#'
#'   Iskrev, N. (2010). Local identification in DSGE models.
#'   \emph{Journal of Monetary Economics} 57(2), 189-202. (Rank/collinearity
#'   and drop-one analysis -- the rank-based route this diagnostic defers to.)
#'
#'   Andrle, M. (2010). A note on identification patterns in DSGE models.
#'   \emph{ECB Working Paper} 1235. (Identification patterns as a
#'   collinearity question rather than a share decomposition.)
#' @noRd
d22_observable_informativeness <- function(model_solve_fn,
                                           theta,
                                           obs_names,
                                           param_names  = NULL,
                                           moment_names = NULL,
                                           eps = 1e-5,
                                           meta = NULL) {
    if (!is.function(model_solve_fn))
      .dynhr_abort("D22: model_solve_fn must be a function.")
    if (!is.numeric(theta) || length(theta) == 0L)
      .dynhr_abort("D22: theta must be a non-empty numeric vector.")
    if (is.null(obs_names) || length(obs_names) == 0L)
      .dynhr_abort("D22: obs_names must be supplied for observable informativeness.")
    obs_names <- unique(as.character(obs_names))

    # --- parameter labels (by name when theta is named) ---
    th_names <- names(theta)
    if (is.null(th_names)) {
      if (is.null(param_names)) param_names <- paste0("theta_", seq_along(theta))
      if (length(param_names) != length(theta))
        .dynhr_abort(sprintf("D22: theta is unnamed and param_names has length %d, not %d.",
                             length(param_names), length(theta)))
      th_names <- param_names
    } else if (is.null(param_names)) {
      param_names <- th_names
    } else if (!all(param_names %in% th_names)) {
      if (length(param_names) != length(theta))
        .dynhr_abort(sprintf("D22: param_names not in names(theta): %s",
                             paste(setdiff(param_names, th_names), collapse = ", ")))
      .dynhr_warn("D22: param_names do not match names(theta); labelling parameters by position.")
      th_names <- param_names
    }
    if (anyDuplicated(th_names))
      .dynhr_abort("D22: parameter names must be unique.")

    # --- moment labels (by name when the moments are named) ---
    f0 <- model_solve_fn(theta)
    n_out <- length(f0)
    if (n_out == 0L) .dynhr_abort("D22: model_solve_fn(theta) returned no moments.")
    f_names <- names(f0)
    if (is.null(f_names)) {
      if (is.null(moment_names))
        .dynhr_abort("D22: model_solve_fn returns unnamed moments; supply moment_names.")
      if (length(moment_names) != n_out)
        .dynhr_abort(sprintf("D22: moment_names has length %d but model_solve_fn returns %d moments.",
                             length(moment_names), n_out))
      f_names <- moment_names
    } else if (is.null(moment_names)) {
      moment_names <- f_names
    } else if (!all(moment_names %in% f_names)) {
      if (length(moment_names) != n_out)
        .dynhr_abort(sprintf("D22: moment_names not in names(model_solve_fn(theta)): %s",
                             paste(setdiff(moment_names, f_names), collapse = ", ")))
      .dynhr_warn(sprintf(
        "D22: moment_names do not match names(model_solve_fn(theta)) (e.g. '%s'); labelling moments by position.",
        f_names[1L]))
      f_names <- moment_names
    }
    if (anyDuplicated(f_names))
      .dynhr_abort("D22: moment names must be unique.")

    J <- .numerical_jacobian(model_solve_fn, theta, eps = eps)
    dimnames(J) <- list(f_names, th_names)
    J <- J[moment_names, param_names, drop = FALSE]
    n_par <- length(param_names)

    moment_to_obs <- .d22_moment_to_observable(moment_names, obs_names)
    is_cross <- !is.na(moment_to_obs) & moment_to_obs == "<cross>"
    cross_moments <- moment_names[is_cross]
    moment_to_obs[is_cross] <- NA_character_
    unmapped_moments <- setdiff(moment_names[is.na(moment_to_obs)], cross_moments)
    valid <- !is.na(moment_to_obs)
    if (!any(valid)) {
      .dynhr_abort(paste0(
        "D22: could not map any moment_names to a single observable. ",
        "Use names like 'obs_sd', 'var_obs' or 'acv1_obs'."))
    }

    obs_set <- obs_names[obs_names %in% moment_to_obs]
    unmapped_obs <- setdiff(obs_names, obs_set)
    info_mat <- matrix(0, nrow = length(obs_set), ncol = n_par,
                       dimnames = list(obs_set, param_names))
    for (ob in obs_set) {
      idx <- which(moment_to_obs == ob)
      info_mat[ob, ] <- colSums(J[idx, , drop = FALSE]^2)
    }
    tot <- colSums(info_mat)
    failed <- !is.finite(tot)
    # Zero-information floor: central-difference round-off is ~ eps_mach*|f|/eps.
    f_scale <- max(1, abs(f0[is.finite(f0)]))
    uninformed <- !failed & sqrt(tot) <= 1e-9 * f_scale
    shares <- sweep(info_mat, 2, tot, "/")
    shares[, failed | uninformed] <- NA_real_
    dominant_obs <- vapply(seq_len(n_par), function(k) {
      x <- shares[, k]
      if (anyNA(x)) NA_character_ else obs_set[which.max(x)]
    }, character(1))
    names(dominant_obs) <- param_names

    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      long <- data.frame(
        Observable = factor(rep(obs_set, times = n_par), levels = rev(obs_set)),
        Parameter  = factor(rep(param_names, each = length(obs_set)), levels = param_names),
        Share      = as.vector(shares)
      )
      long$label <- ifelse(is.na(long$Share), "NA", sprintf("%.0f%%", 100 * long$Share))
      long$txt_col <- ifelse(is.na(long$Share) | long$Share > 0.6, "black", "white")
      p <- ggplot2::ggplot(
        long, ggplot2::aes(x = Parameter, y = Observable, fill = Share)
      ) +
        ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
        scale_fill_dynhr_cividis(
          name = "share of sensitivity", limits = c(0, 1),
          labels = function(x) sprintf("%.0f%%", 100 * x),
          guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(5, "cm"),
                                           barheight = ggplot2::unit(0.35, "cm"))
        ) +
        theme_dynhr_diagnostic() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            angle = 45, hjust = 1, vjust = 1, size = ggplot2::rel(0.8)
          )
        ) +
        ggplot2::labs(
          title = "D22: Observable informativeness by parameter",
          subtitle = paste0(
            "Share of sum of (d moment / d theta)^2 over each observable's moments; ",
            "columns sum to 100%.\n",
            "Depends on moment units; ignores parameter collinearity (see D1/D20).\n",
            "NA = parameter moves no mapped moment",
            if (length(cross_moments)) sprintf("; %d cross-moment(s) excluded", length(cross_moments)) else "",
            "."),
          x = "Parameter", y = "Observable"
        )
      if (length(obs_set) * n_par <= 150L) {
        p <- p + ggplot2::geom_text(ggplot2::aes(label = label, colour = txt_col),
                                    size = 2.8, show.legend = FALSE) +
          ggplot2::scale_colour_identity()
      }
      plots$informativeness_heatmap <- .apply_meta(p, meta)
    }

    uninformed_params <- param_names[uninformed]
    failed_params <- param_names[failed]
    ok <- !is.na(dominant_obs)
    dom_pairs <- if (any(ok))
      paste(sprintf("%s->%s(%.0f%%)", param_names[ok], dominant_obs[ok],
                    100 * apply(shares[, ok, drop = FALSE], 2, max)), collapse = ", ")
    else "none"

    .make_result(
      result = list(
        jacobian = J,
        moment_to_observable = moment_to_obs,
        information = info_mat,
        informativeness = shares,
        dominant_observable = dominant_obs,
        uninformed_params = uninformed_params,
        failed_params = failed_params,
        cross_moments = cross_moments,
        unmapped_moments = unmapped_moments,
        unmapped_observables = unmapped_obs
      ),
      pass = NA,
      plots = plots,
      summary = paste0(
        sprintf("D22 Observable informativeness: mapped %d/%d moments to %d observables. Dominant observable by parameter: %s",
                sum(valid), length(valid), length(obs_set), dom_pairs),
        if (length(uninformed_params))
          sprintf(". No mapped moment moves: %s", paste(uninformed_params, collapse = ", ")) else "",
        if (length(failed_params))
          sprintf(". Non-finite sensitivity: %s", paste(failed_params, collapse = ", ")) else "",
        if (length(cross_moments))
          sprintf(". Cross-moments excluded: %s", paste(cross_moments, collapse = ", ")) else "",
        if (length(unmapped_moments))
          sprintf(". Unmapped moments: %s", paste(unmapped_moments, collapse = ", ")) else "",
        if (length(unmapped_obs))
          sprintf(". Observables with no moments: %s", paste(unmapped_obs, collapse = ", ")) else ""
      ),
      llm_summary = {
        dom_tab <- sort(table(dominant_obs[ok]), decreasing = TRUE)
        top3 <- head(names(dom_tab), 3L)
        dom_str <- if (length(top3))
          paste(sprintf("%s(%d)", top3, as.integer(dom_tab[top3])), collapse = ", ")
        else "none"
        if (length(dom_tab) > 3L) dom_str <- sprintf("%s + %d more", dom_str, length(dom_tab) - 3L)
        paste(c(
          "D22 | Observable Informativeness | INFO",
          sprintf("  params=%d moments_mapped=%d/%d observables=%d cross_moments=%d",
                  n_par, sum(valid), length(valid), length(obs_set), length(cross_moments)),
          sprintf("  dominant_observable counts (top-3): %s", dom_str),
          if (length(uninformed_params))
            sprintf("  uninformed_params (move no mapped moment): %s",
                    paste(uninformed_params, collapse = ", ")),
          if (length(failed_params))
            sprintf("  failed_params (non-finite sensitivity): %s",
                    paste(failed_params, collapse = ", ")),
          if (length(unmapped_obs))
            sprintf("  observables_without_moments: %s", paste(unmapped_obs, collapse = ", ")),
          sprintf("  action: %s",
                  if (length(uninformed_params))
                    sprintf("%s are not moved by any mapped moment; add observables/moments or calibrate them.",
                            paste(head(uninformed_params, 3), collapse = ", "))
                  else
                    "Observables that dominate no parameter add little marginal sensitivity; check D1/D20 before dropping them.")
        ), collapse = "\n")
      }
    )
}

# Attribute each moment name to one observable by "_"-token matching.
# Returns the observable name, NA (no match) or "<cross>" (>= 2 distinct
# observables after dropping matches nested inside a longer match).
.d22_moment_to_observable <- function(moment_names, obs_names) {
  obs_tok <- strsplit(obs_names, "_", fixed = TRUE)
  out <- vapply(moment_names, function(mn) {
    tok <- strsplit(mn, "_", fixed = TRUE)[[1]]
    spans <- NULL
    for (j in seq_along(obs_names)) {
      ot <- obs_tok[[j]]
      k <- length(ot)
      if (k == 0L || k > length(tok)) next
      for (i in seq_len(length(tok) - k + 1L)) {
        if (all(tok[i:(i + k - 1L)] == ot))
          spans <- rbind(spans, data.frame(obs = obs_names[j], s = i, e = i + k - 1L))
      }
    }
    if (is.null(spans)) return(NA_character_)
    len <- spans$e - spans$s
    nested <- vapply(seq_len(nrow(spans)), function(r)
      any(spans$s <= spans$s[r] & spans$e >= spans$e[r] & len > len[r]), logical(1))
    hits <- unique(spans$obs[!nested])
    if (length(hits) == 1L) hits else "<cross>"
  }, character(1))
  names(out) <- moment_names
  out
}
