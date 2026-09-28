## R/diag-pre-d38-sloppiness.R
## --------------------------------------------------------------------------
## D38. Posterior-curvature "sloppiness" spectrum
##
## Diagnoses ridge/flat-direction difficulty by the eigen-spectrum of the
## negative log-posterior Hessian at a supplied or found mode, taken with
## respect to LOG parameters (Gutenkunst et al. 2007) so that the spectrum is
## unit-free: rescaling a parameter (theta_i -> c * theta_i) is a shift of
## log(theta_i) and leaves every eigenvalue unchanged.
##
## References:
##   Gutenkunst et al. (2007) PLoS Comput Biol 3(10):e189 — sloppy models
##   Transtrum et al. (2011) J Chem Phys 135:014102 — manifold boundary
## --------------------------------------------------------------------------


#' D38. Posterior-curvature sloppiness spectrum (Gutenkunst et al. 2007)
#'
#' Eigen-spectrum of the negative log-posterior Hessian
#' \eqn{H = -\partial^2 \log p(\theta \mid Y)} at \code{theta}. Small
#' eigenvalues are "sloppy" (flat/ridged) directions, large ones "stiff".
#'
#' @section Coordinates:
#' With \code{scale = "log"} (default) the Hessian is taken with respect to
#' \eqn{u_i = \log|\theta_i|}, i.e. \eqn{H_u = D H_\theta D} with
#' \eqn{D = \mathrm{diag}(\theta)} at a stationary point. This is the
#' Gutenkunst et al. convention and makes eigenvalues unit-free. A parameter
#' that is exactly zero has no log coordinate; it is kept on
#' its raw scale, listed in \code{$result$raw_scale_params}, and a warning is
#' raised because its row/column of \eqn{H_u} then carries units.
#' With \code{scale = "raw"} the raw-parameter Hessian is decomposed (its
#' spectrum depends on the units of each parameter).
#'
#' @section Hessian sources:
#' \enumerate{
#'   \item \code{hessian} supplied: the RAW-parameter negative log-posterior
#'     Hessian; for \code{scale = "log"} it is transformed as \eqn{D H D}
#'     (exact at a mode; off-mode the gradient term is not available).
#'   \item \code{log_post_fn} supplied: central differences of
#'     \code{log_post_fn} directly in the chosen coordinates (step \code{h}
#'     in log units for \code{scale = "log"}), repeated at step \code{2h};
#'     the spectral norm of the difference is the finite-difference noise
#'     estimate (the idea of \code{.ident_equilibrated_rank()}).
#' }
#'
#' @section Eigenvalue classes:
#' The noise floor is \code{max(n * eps * max|lambda|, 10 * fd_noise)}.
#' Each eigenvalue is classed \code{"negative"} (\eqn{\lambda < -}floor: not
#' a local maximum of the posterior), \code{"unresolved"}
#' (\eqn{|\lambda| \le} floor: indistinguishable from zero / FD noise),
#' \code{"sloppy"} (floor \eqn{< \lambda <} \code{sloppy_cutoff}), or
#' \code{"stiff"}. \code{sloppy_cutoff = lambda_max * sloppy_threshold}.
#'
#' @section Where the six decades come from, and what they are worth:
#' Gutenkunst et al. (2007) established the empirical regularity that
#' sensitivity/Fisher eigenvalues in nonlinear models are roughly
#' log-uniformly spaced over many orders of magnitude -- commonly six or
#' more. Transtrum, Machta & Sethna (2011) give the information-geometric
#' account of the same picture. Neither fixes a universal numeric cutoff
#' separating "sloppy" from "stiff": sloppiness is a property of the
#' spectrum's spread, not of one eigenvalue. \code{sloppy_threshold = 1e-6}
#' simply encodes the six-decade Gutenkunst convention and is a package
#' choice, not a value taken from a paper -- which is why this diagnostic is
#' INFO and reports \code{spread_decades}, \code{condition_number} and
#' \code{participation_ratio} alongside the counts rather than a verdict.
#'
#' This is a transplant from the physics / systems-biology sloppiness
#' literature: no DSGE paper or toolbox (Dynare, IRIS, MacroModelling.jl,
#' RISE) uses the sloppy/stiff framework. Applied Bayesian DSGE work
#' diagnoses the same underlying problem through the mode Hessian's condition
#' number or through poor MCMC mixing / low effective sample size.
#'
#' \strong{A meaningful "sloppy" verdict needs an ANALYTIC Hessian.} With a
#' finite-difference Hessian the noise floor
#' \code{max(n * eps * max|lambda|, 10 * fd_noise)} is typically far above
#' the eigenvalues that make a model sloppy, so genuinely sloppy directions
#' are returned as \code{"unresolved"} rather than \code{"sloppy"} and only
#' the most extreme sloppiness survives. That is the correct behaviour --
#' the floor exists precisely so the diagnostic cannot mistake FD noise for
#' flatness -- but it means a high \code{n_unresolved} under
#' \code{log_post_fn} is a statement about the differencing, not about the
#' posterior. Supply \code{hessian} from an analytic/AD source when the
#' spread itself is the question.
#'
#' @param theta         Named numeric vector — the evaluation point.
#' @param log_post_fn   Function \code{theta -> scalar} log-posterior (or
#'   list with \code{$logpost}), evaluated on the RAW parameter scale.
#' @param hessian       Optional pre-computed RAW-scale negative
#'   log-posterior Hessian (n x n). When supplied, \code{log_post_fn} is
#'   ignored.
#' @param param_names   Character vector of parameter names (length n).
#'   Defaults to \code{names(theta)}.
#' @param n_flat        Integer. Number of flattest (smallest |lambda|)
#'   eigendirections, and of stiffest ones, to tabulate (default 3L).
#' @param sloppy_threshold  Numeric in (0, 1). Eigenvalues below
#'   \code{lambda_max * sloppy_threshold} are sloppy (default 1e-6, i.e. a
#'   spread of six decades, the Gutenkunst et al. scale). A package choice
#'   following that convention; the sloppiness literature specifies no
#'   universal cutoff. It never gates -- D38 is INFO.
#' @param n_top_loadings Integer. Number of top parameter loadings kept per
#'   tabulated direction (default 5L). Full eigenvectors are always in
#'   \code{$result$eigenvectors}.
#' @param h             Finite-difference step (default 1e-4; in log units
#'   for \code{scale = "log"}; each step is \code{h * max(1, |coord|)}).
#' @param scale         \code{"log"} (default) or \code{"raw"}.
#' @param meta          Optional metadata list (passed to \code{.apply_meta}).
#'
#' @return A \code{dynhr_diagnostic} (\code{pass = NA}, informational) with
#'   \code{$result} containing \code{eigenvalues} (descending),
#'   \code{eigenvectors} (n x n, rows named by parameter), \code{class},
#'   \code{condition_number} (\code{Inf} if any eigenvalue is unresolved or
#'   negative, else lambda_max / lambda_min), \code{spread_decades}
#'   (log10 of lambda_max over the smallest resolved positive eigenvalue),
#'   \code{participation_ratio} (on eigenvalues clipped at 0; in [1, n]),
#'   \code{spectral_gap}, \code{n_sloppy}, \code{n_unresolved},
#'   \code{n_negative}, \code{sloppy_cutoff}, \code{noise_floor},
#'   \code{fd_noise}, \code{scale}, \code{raw_scale_params},
#'   \code{flat_directions} / \code{stiff_directions} (lists of loading
#'   data.frames with attributes \code{eigenvalue}, \code{class_label}),
#'   \code{hessian} (the matrix decomposed) and \code{hessian_raw}
#'   (supplied path only).
#'
#' @references
#'   Gutenkunst, R. N., Waterfall, J. J., Casey, F. P., Brown, K. S.,
#'   Myers, C. R., & Sethna, J. P. (2007). Universally sloppy parameter
#'   sensitivities in systems biology models. \emph{PLoS Computational
#'   Biology}, 3(10), e189.
#'
#'   Transtrum, M. K., Machta, B. B., & Sethna, J. P. (2011). Geometry of
#'   nonlinear least squares with applications to sloppy models and
#'   optimization. \emph{Journal of Chemical Physics}, 135(1), 014102.
#' @noRd
d38_sloppiness <- function(theta,
                           log_post_fn      = NULL,
                           hessian          = NULL,
                           param_names      = NULL,
                           n_flat           = 3L,
                           sloppy_threshold = 1e-6,
                           n_top_loadings   = 5L,
                           h                = 1e-4,
                           scale            = c("log", "raw"),
                           meta             = NULL) {

  scale <- match.arg(scale)

  ## ------------------------------------------------------------------
  ## 0. Setup / validation
  ## ------------------------------------------------------------------
  if (!is.numeric(theta) || length(theta) < 1L || !all(is.finite(theta)))
    .dynhr_abort("d38: `theta` must be a non-empty, finite numeric vector")
  n <- length(theta)
  if (is.null(param_names))
    param_names <- if (!is.null(names(theta))) names(theta) else paste0("theta_", seq_len(n))
  if (length(param_names) != n)
    .dynhr_abort(sprintf("d38: param_names length (%d) != length(theta) (%d)",
                         length(param_names), n))
  if (!is.numeric(sloppy_threshold) || length(sloppy_threshold) != 1L ||
      !is.finite(sloppy_threshold) || sloppy_threshold <= 0 || sloppy_threshold >= 1)
    .dynhr_abort("d38: `sloppy_threshold` must be a single number in (0, 1)")
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0)
    .dynhr_abort("d38: `h` must be a single positive number")

  n_flat <- max(0L, min(as.integer(n_flat), n))
  n_top_loadings <- max(1L, as.integer(n_top_loadings))

  ## Log coordinates: u_i = log|theta_i|, theta_i = sign_i * exp(u_i),
  ## d theta_i / d u_i = theta_i. Zero / non-finite params stay raw.
  use_log <- if (scale == "log") theta != 0 else rep(FALSE, n)
  raw_scale_params <- if (scale == "log") param_names[!use_log] else character(0)
  if (length(raw_scale_params) > 0L)
    .dynhr_warn(sprintf(paste0(
      "d38: parameter(s) %s are exactly zero and have no log ",
      "coordinate; kept on the raw scale, so the spectrum is not unit-free ",
      "in those directions."), paste(raw_scale_params, collapse = ", ")))

  ## ------------------------------------------------------------------
  ## 1. Obtain the negative log-posterior Hessian in the chosen coordinates
  ## ------------------------------------------------------------------
  H_raw    <- NULL
  fd_noise <- NA_real_

  if (!is.null(hessian)) {
    if (!is.matrix(hessian) || nrow(hessian) != n || ncol(hessian) != n)
      .dynhr_abort("d38: `hessian` must be an n x n matrix where n = length(theta)")
    if (!all(is.finite(hessian)))
      .dynhr_abort("d38: `hessian` contains non-finite entries")
    H_raw <- 0.5 * (hessian + t(hessian))
    dimnames(H_raw) <- list(param_names, param_names)
    d <- ifelse(use_log, theta, 1)
    H <- H_raw * outer(d, d)
  } else if (!is.null(log_post_fn)) {
    if (!is.function(log_post_fn))
      .dynhr_abort("d38: `log_post_fn` must be a function")
    sgn <- sign(theta)
    u0  <- ifelse(use_log, log(abs(theta)), theta)
    names(u0) <- names(theta)
    to_theta <- function(u) {
      th <- ifelse(use_log, sgn * exp(u), u)
      names(th) <- names(theta)
      th
    }
    g <- function(u) log_post_fn(to_theta(u))
    H1 <- -num_hessian(g, u0, h = h)
    H2 <- -num_hessian(g, u0, h = 2 * h)
    if (!all(is.finite(H1)))
      .dynhr_abort(paste0(
        "d38: the finite-difference Hessian is non-finite; log_post_fn ",
        "returned a non-finite value near theta (reduce `h` or move theta ",
        "away from a boundary)"))
    H1 <- 0.5 * (H1 + t(H1))
    if (all(is.finite(H2))) {
      H2 <- 0.5 * (H2 + t(H2))
      fd_noise <- max(abs(eigen(H1 - H2, symmetric = TRUE, only.values = TRUE)$values))
    }
    H <- H1
  } else {
    .dynhr_abort("d38: supply either `hessian` (pre-computed) or `log_post_fn`")
  }

  dimnames(H) <- list(param_names, param_names)

  ## ------------------------------------------------------------------
  ## 2. Eigen-decomposition of the symmetric Hessian
  ## ------------------------------------------------------------------
  ev       <- eigen(H, symmetric = TRUE)
  lambdas  <- ev$values                       # descending
  evectors <- ev$vectors
  dimnames(evectors) <- list(param_names, paste0("lambda_", seq_len(n)))

  ## ------------------------------------------------------------------
  ## 3. Noise floor, classes and scalar summaries
  ## ------------------------------------------------------------------
  abs_max     <- max(abs(lambdas))
  noise_floor <- max(n * .Machine$double.eps * abs_max,
                     if (is.finite(fd_noise)) 10 * fd_noise else 0)
  lambda_max  <- lambdas[1]
  sloppy_cutoff <- if (lambda_max > noise_floor) lambda_max * sloppy_threshold else NA_real_

  cls <- ifelse(lambdas < -noise_floor, "negative",
         ifelse(abs(lambdas) <= noise_floor, "unresolved",
         ifelse(!is.na(sloppy_cutoff) & lambdas < sloppy_cutoff, "sloppy", "stiff")))
  names(cls) <- colnames(evectors)

  n_sloppy     <- sum(cls == "sloppy")
  n_unresolved <- sum(cls == "unresolved")
  n_negative   <- sum(cls == "negative")

  resolved_pos <- lambdas[lambdas > noise_floor]
  spread_decades <- if (length(resolved_pos) >= 1L)
    log10(max(resolved_pos) / min(resolved_pos)) else NA_real_
  condition_number <- if (length(resolved_pos) == 0L || n_unresolved + n_negative > 0L)
    Inf else max(resolved_pos) / min(resolved_pos)

  lam_clip <- pmax(lambdas, 0)
  participation_ratio <- if (sum(lam_clip^2) > 0)
    sum(lam_clip)^2 / sum(lam_clip^2) else NA_real_

  spectral_gap <- if (n >= 2L && lambdas[2] > noise_floor)
    lambdas[1] / lambdas[2] else NA_real_

  ## ------------------------------------------------------------------
  ## 4. Attribution — loading tables (flattest = smallest |lambda|)
  ## ------------------------------------------------------------------
  loading_table <- function(col_idx) {
    evec <- evectors[, col_idx]
    df <- data.frame(parameter = param_names, loading = unname(evec),
                     abs_loading = unname(abs(evec)), stringsAsFactors = FALSE)
    df <- df[order(df$abs_loading, decreasing = TRUE), , drop = FALSE]
    df <- df[seq_len(min(n_top_loadings, nrow(df))), , drop = FALSE]
    rownames(df) <- NULL
    attr(df, "eigenvalue") <- lambdas[col_idx]
    attr(df, "class_label") <- cls[[col_idx]]
    attr(df, "index") <- col_idx
    df
  }
  flat_order  <- order(abs(lambdas))[seq_len(n_flat)]
  stiff_order <- seq_len(n_flat)
  flat_directions  <- lapply(flat_order, loading_table)
  stiff_directions <- lapply(stiff_order, loading_table)
  for (k in seq_along(flat_directions)) attr(flat_directions[[k]], "direction") <- k
  names(flat_directions)  <- sprintf("flat_%d", seq_len(n_flat))
  names(stiff_directions) <- sprintf("stiff_%d", seq_len(n_flat))

  ## ------------------------------------------------------------------
  ## 5. Plots
  ## ------------------------------------------------------------------
  coord_lab <- if (scale == "log") "log-parameters" else "raw parameters"
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plots <- .d38_plots(lambdas, cls, evectors, flat_order, n_flat,
                        noise_floor, sloppy_cutoff, sloppy_threshold,
                        coord_lab, n, n_sloppy, n_unresolved, n_negative,
                        spread_decades, meta)
  }

  ## ------------------------------------------------------------------
  ## 6. Assemble result
  ## ------------------------------------------------------------------
  fmt_dir <- function(fd) {
    top <- head(fd, 3L)
    paste(sprintf("%s%s(%.2f)", ifelse(top$loading >= 0, "+", "-"),
                  top$parameter, top$abs_loading), collapse = " ")
  }
  flagged <- n_sloppy + n_unresolved + n_negative > 0L
  flat_txt <- if (flagged && n_flat >= 1L)
    sprintf(" | flattest (lambda=%.3g, %s): %s",
            attr(flat_directions[[1]], "eigenvalue"),
            attr(flat_directions[[1]], "class_label"),
            fmt_dir(flat_directions[[1]]))
  else ""
  raw_txt <- if (length(raw_scale_params))
    sprintf(" | raw-scale: %s", paste(raw_scale_params, collapse = ", ")) else ""

  summary_text <- sprintf(paste0(
    "D38 Sloppiness spectrum (%s): n=%d | spread=%.1f decades | cond=%.2e | ",
    "PR=%.2f | sloppy=%d unresolved=%d negative=%d%s%s"),
    coord_lab, n, spread_decades, condition_number, participation_ratio,
    n_sloppy, n_unresolved, n_negative, flat_txt, raw_txt)

  action <- if (n_negative > 0L)
    "Hessian has negative curvature: theta is not a posterior mode; re-run the mode finder before reading the spectrum."
  else if (n_unresolved > 0L)
    "Some curvature is indistinguishable from zero/FD noise: the listed combination is (locally) unidentified; add information or fix a parameter."
  else if (n_sloppy > 0L)
    sprintf("Sloppy model: %d direction(s) below lambda_max*%.0e. Consider tightening priors, reparameterising, or adding observables for the flat combination.",
            n_sloppy, sloppy_threshold)
  else
    "No eigenvalue below the sloppy cutoff; posterior curvature spread is moderate."

  llm_summary <- paste(c(
    "D38 | Sloppiness Spectrum | INFO",
    sprintf("  coords=%s n_params=%d spread_decades=%.2f cond=%.2e PR=%.2f",
            coord_lab, n, spread_decades, condition_number, participation_ratio),
    sprintf("  n_sloppy=%d n_unresolved=%d n_negative=%d cutoff=%.3e noise_floor=%.3e",
            n_sloppy, n_unresolved, n_negative, sloppy_cutoff, noise_floor),
    sprintf("  eigenvalues: %s",
            paste(sprintf("%.3e", head(lambdas, 10L)), collapse = ", ")),
    if (n_flat >= 1L)
      sprintf("  flattest direction (lambda=%.3e, %s): %s",
              attr(flat_directions[[1]], "eigenvalue"),
              attr(flat_directions[[1]], "class_label"),
              fmt_dir(flat_directions[[1]])),
    if (n_flat >= 1L)
      sprintf("  stiffest direction (lambda=%.3e): %s",
              attr(stiff_directions[[1]], "eigenvalue"),
              fmt_dir(stiff_directions[[1]])),
    if (length(raw_scale_params))
      sprintf("  raw_scale_params (zero-valued, not unit-free): %s",
              paste(raw_scale_params, collapse = ", ")),
    sprintf("  action: %s", action)
  ), collapse = "\n")

  .make_result(
    result  = list(
      eigenvalues         = lambdas,
      eigenvectors        = evectors,
      class               = cls,
      condition_number    = condition_number,
      spread_decades      = spread_decades,
      participation_ratio = participation_ratio,
      spectral_gap        = spectral_gap,
      n_sloppy            = n_sloppy,
      n_unresolved        = n_unresolved,
      n_negative          = n_negative,
      sloppy_cutoff       = sloppy_cutoff,
      noise_floor         = noise_floor,
      fd_noise            = fd_noise,
      scale               = scale,
      raw_scale_params    = raw_scale_params,
      flat_directions     = flat_directions,
      stiff_directions    = stiff_directions,
      hessian             = H,
      hessian_raw         = H_raw
    ),
    pass        = NA,   # informational: sloppiness is a spectrum
    plots       = plots,
    summary     = summary_text,
    llm_summary = llm_summary
  )
}


## Plot builder for D38: eigenvalue ladder + flat-direction loadings.
#' @noRd
.d38_plots <- function(lambdas, cls, evectors, flat_order, n_flat,
                       noise_floor, sloppy_cutoff, sloppy_threshold,
                       coord_lab, n, n_sloppy, n_unresolved, n_negative,
                       spread_decades, meta) {
  cls_cols <- c(stiff      = unname(tol_vibrant["blue"]),
                sloppy     = unname(tol_vibrant["orange"]),
                unresolved = unname(tol_vibrant["grey"]),
                negative   = unname(tol_vibrant["red"]))
  cls_shapes <- c(stiff = 16, sloppy = 16, unresolved = 1, negative = 4)
  cls_labs <- c(stiff = "stiff", sloppy = "sloppy",
                unresolved = "unresolved (<= noise floor)",
                negative = "negative (magnitude shown)")
  present <- intersect(names(cls_cols), unique(cls))

  abs_l <- abs(lambdas)
  pos   <- abs_l[abs_l > 0]
  y_min <- min(c(pos, noise_floor[noise_floor > 0], sloppy_cutoff[is.finite(sloppy_cutoff)],
                 if (!length(pos)) 1)) / 10
  eig_df <- data.frame(index = seq_along(lambdas),
                       y     = pmax(abs_l, y_min),
                       class = factor(cls, levels = names(cls_cols)))

  ## Reference lines: cutoff labelled above-right, noise floor below-left so
  ## the two labels never collide when the lines are close.
  ref <- data.frame(y = numeric(0), lab = character(0), lt = character(0),
                    x = numeric(0), hj = numeric(0), vj = numeric(0))
  if (is.finite(sloppy_cutoff))
    ref <- rbind(ref, data.frame(
      y = sloppy_cutoff, lt = "dashed", x = n + 0.45, hj = 1, vj = -0.4,
      lab = sprintf("sloppy cutoff = lambda_max x %.0e", sloppy_threshold)))
  if (noise_floor > 0)
    ref <- rbind(ref, data.frame(y = noise_floor, lt = "dotted",
                                 x = 0.55, hj = 0, vj = 1.4,
                                 lab = "noise floor"))

  p_spec <- ggplot2::ggplot(eig_df, ggplot2::aes(x = index, y = y)) +
    ggplot2::geom_segment(ggplot2::aes(x = index - 0.35, xend = index + 0.35,
                                       yend = y, colour = class),
                          linewidth = 1.1) +
    ggplot2::geom_point(ggplot2::aes(colour = class, shape = class), size = 2.4)
  if (nrow(ref) > 0L) {
    p_spec <- p_spec +
      ggplot2::geom_hline(data = ref, ggplot2::aes(yintercept = y),
                          linetype = ref$lt, colour = dynhr_colours$grey) +
      ggplot2::annotate("text", x = ref$x, y = ref$y, label = ref$lab,
                        hjust = ref$hj, vjust = ref$vj, size = 3.2,
                        colour = "grey30")
  }
  p_spec <- p_spec +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(breaks = seq_len(n),
                                limits = c(0.5, n + 0.5)) +
    ggplot2::scale_colour_manual(values = cls_cols[present],
                                 labels = cls_labs[present],
                                 breaks = present, name = NULL, drop = TRUE) +
    ggplot2::scale_shape_manual(values = cls_shapes[present],
                                labels = cls_labs[present],
                                breaks = present, name = NULL, drop = TRUE) +
    theme_dynhr_diagnostic() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      title    = sprintf("D38: Hessian eigenvalue ladder (%s)", coord_lab),
      subtitle = sprintf(
        "spread = %.1f decades | sloppy = %d, unresolved = %d, negative = %d",
        spread_decades, n_sloppy, n_unresolved, n_negative),
      x = "Eigen-direction (1 = stiffest)",
      y = "Eigenvalue magnitude of -Hessian (log scale)")
  plots <- list(spectrum = .apply_meta(p_spec, meta))

  if (n_flat >= 1L) {
    pn <- rownames(evectors)
    ld <- do.call(rbind, lapply(seq_len(n_flat), function(k) {
      j <- flat_order[k]
      v <- evectors[, j]
      keep <- order(abs(v), decreasing = TRUE)[seq_len(min(10L, length(v)))]
      data.frame(
        panel = sprintf("flat %d (%s)\nlambda = %.3g", k, cls[[j]], lambdas[j]),
        k = k, parameter = pn[keep], loading = unname(v[keep]),
        key = sprintf("%s___%d", pn[keep], k),
        stringsAsFactors = FALSE)
    }))
    ld$panel <- factor(ld$panel, levels = unique(ld$panel[order(ld$k)]))
    ld$key   <- factor(ld$key, levels = ld$key[order(ld$k, abs(ld$loading))])
    ld$sign  <- factor(ifelse(ld$loading >= 0, "positive", "negative"),
                       levels = c("positive", "negative"))
    p_load <- ggplot2::ggplot(ld, ggplot2::aes(x = loading, y = key, fill = sign)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_vline(xintercept = 0, colour = dynhr_colours$grey,
                          linewidth = 0.3) +
      ggplot2::facet_wrap(~panel, scales = "free_y",
                          ncol = min(3L, n_flat)) +
      ggplot2::scale_y_discrete(labels = function(x) sub("___[0-9]+$", "", x)) +
      ggplot2::scale_x_continuous(limits = c(-1, 1)) +
      ggplot2::scale_fill_manual(
        values = c(positive = unname(tol_vibrant["blue"]),
                   negative = unname(tol_vibrant["orange"])),
        drop = FALSE, name = "loading sign") +
      theme_dynhr_diagnostic() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::labs(
        title    = sprintf("D38: Flattest eigen-directions (%s)", coord_lab),
        subtitle = "Unit eigenvector loadings; the overall sign of each direction is arbitrary",
        x = "Eigenvector loading", y = NULL)
    plots$flat_directions <- .apply_meta(p_load, meta)
  }
  plots
}


## ---------------------------------------------------------------------------
## Public entry point
## ---------------------------------------------------------------------------

#' D38. Posterior-curvature sloppiness spectrum
#'
#' Compute the eigen-spectrum of the negative log-posterior Hessian at
#' \code{theta} to diagnose parameter-space ridge / flat-direction difficulty
#' ("sloppy" directions following Gutenkunst et al. 2007). By default the
#' Hessian is taken with respect to log-parameters, so the eigenvalues are
#' unit-free; zero-valued parameters stay on their raw scale (with a warning).
#' Eigenvalues are classed stiff / sloppy / unresolved (below the
#' finite-difference noise floor) / negative (theta is not a mode).
#'
#' @param theta       Named numeric parameter vector at which to assess
#'   sloppiness (raw scale).
#' @param log_post_fn Optional log-posterior function of raw \code{theta}; if
#'   supplied (or built from \code{model}), the Hessian is formed numerically
#'   in the coordinates given by \code{scale}, at steps \code{h} and
#'   \code{2h} (their difference sets the noise floor).
#' @param hessian     Optional pre-computed RAW-scale Hessian of the negative
#'   log-posterior at \code{theta}; supplied instead of \code{log_post_fn}.
#'   With \code{scale = "log"} it is transformed to
#'   \code{diag(theta) \%*\% hessian \%*\% diag(theta)} (exact at a mode).
#' @param param_names Optional character vector of parameter names (for
#'   labelling the eigenvector loadings).
#' @param n_flat      Integer number of flattest (smallest |eigenvalue|) and
#'   stiffest eigendirections to tabulate (default 3).
#' @param sloppy_threshold Relative threshold in (0, 1): eigenvalues below
#'   \code{lambda_max * sloppy_threshold} are "sloppy" (default 1e-6).
#' @param n_top_loadings Integer number of top parameter loadings kept per
#'   tabulated direction (default 5).
#' @param h           Finite-difference step (default 1e-4; log units when
#'   \code{scale = "log"}).
#' @param scale       \code{"log"} (default; unit-free) or \code{"raw"}.
#' @param meta        Optional \code{diag_meta()} list for plot captions.
#' @param model       Optional parsed dynhr model (used to build
#'   \code{log_post_fn} when neither \code{hessian} nor \code{log_post_fn}
#'   are supplied).
#' @param data        n_obs x T data matrix (observables x time), required
#'   when \code{model} is supplied.
#' @param prior_spec  Prior specification list, as from
#'   \code{make_log_posterior}.
#' @param obs_vars    Character vector of observable variable names, required
#'   when \code{model} is supplied.
#' @param compiled    Compiled model object (\code{compile_model} output),
#'   required when \code{model} is supplied.
#' @param me_variance Scalar measurement-error variance (default 0).
#' @param \dots       Passed to \code{make_log_posterior} when constructing
#'   the log-posterior internally.
#'
#' @return A \code{dynhr_diagnostic} with \code{pass = NA} (informational).
#'   \code{$result} holds \code{eigenvalues}, \code{eigenvectors} (rows named
#'   by parameter), \code{class}, \code{condition_number},
#'   \code{spread_decades}, \code{participation_ratio}, \code{spectral_gap},
#'   \code{n_sloppy}, \code{n_unresolved}, \code{n_negative},
#'   \code{sloppy_cutoff}, \code{noise_floor}, \code{fd_noise}, \code{scale},
#'   \code{raw_scale_params}, \code{flat_directions},
#'   \code{stiff_directions}, \code{hessian} and \code{hessian_raw}.
#'   Plots: \code{spectrum} (eigenvalue ladder with the cutoff and noise
#'   floor drawn) and \code{flat_directions} (loadings).
#'
#' @seealso \code{d1_local_identification},
#'   \code{d20_fisher_identification_strength}
#' @export
diag_sloppiness <- function(theta,
                            log_post_fn      = NULL,
                            hessian          = NULL,
                            model            = NULL,
                            data             = NULL,
                            prior_spec       = NULL,
                            obs_vars         = NULL,
                            compiled         = NULL,
                            me_variance      = 0,
                            param_names      = NULL,
                            n_flat           = 3L,
                            sloppy_threshold = 1e-6,
                            n_top_loadings   = 5L,
                            h                = 1e-4,
                            scale            = c("log", "raw"),
                            meta             = NULL,
                            ...) {

  if (is.null(hessian) && is.null(log_post_fn)) {
    if (is.null(model) || is.null(data) || is.null(obs_vars) || is.null(compiled))
      .dynhr_abort(paste0(
        "diag_sloppiness: supply either `hessian`, `log_post_fn`, ",
        "or (model + data + obs_vars + compiled)"
      ))
    log_post_fn <- make_log_posterior(
      model        = model,
      data         = data,
      prior_spec   = prior_spec,
      obs_vars     = obs_vars,
      compiled     = compiled,
      me_variance  = me_variance,
      ...
    )
  }

  d38_sloppiness(
    theta            = theta,
    log_post_fn      = log_post_fn,
    hessian          = hessian,
    param_names      = param_names,
    n_flat           = n_flat,
    sloppy_threshold = sloppy_threshold,
    n_top_loadings   = n_top_loadings,
    h                = h,
    scale            = match.arg(scale),
    meta             = meta
  )
}
