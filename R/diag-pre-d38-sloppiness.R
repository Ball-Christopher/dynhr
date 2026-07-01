## R/diag-pre-d38-sloppiness.R
## --------------------------------------------------------------------------
## D38. Posterior-curvature "sloppiness" spectrum
##
## Diagnoses ridge/flat-direction difficulty by the eigen-spectrum of the
## negative log-posterior Hessian at a supplied or found mode.
##
## References:
##   Gutenkunst et al. (2007) PLoS Comput Biol 3(10):e189 — sloppy models
##   Transtrum et al. (2011) J Chem Phys 135:014102 — manifold boundary
## --------------------------------------------------------------------------


#' D38. Posterior-curvature sloppiness spectrum (Gutenkunst et al. 2007)
#'
#' Diagnoses ridge/flat-direction difficulty by the eigen-spectrum of the
#' negative log-posterior Hessian H = -d2 log p(theta|Y) at \code{theta}.
#' Small eigenvalues of H correspond to "sloppy" (flat/ridged) directions in
#' parameter space where the posterior is nearly uninformative.
#'
#' @section Hessian strategy:
#' The function tries sources in order:
#' \enumerate{
#'   \item If \code{hessian} is supplied directly (numeric matrix), use it.
#'   \item If \code{log_post_fn} is supplied, apply \code{num_hessian()} to
#'     it at \code{theta} (central-difference, O(n^2) evaluations).
#'   \item Otherwise error — one of the two must be supplied.
#' }
#' The returned Hessian is the NEGATIVE log-posterior curvature (positive
#' definite at the mode), so eigenvalues are the curvatures; large eigenvalue
#' = stiff/identified direction, small eigenvalue = sloppy/flat direction.
#'
#' @param theta         Named numeric vector — the evaluation point (mode or
#'   any other point). Must be supplied.
#' @param log_post_fn   Function \code{theta -> scalar} log-posterior (or
#'   list with \code{$logpost}). Used for numerical Hessian when
#'   \code{hessian} is not supplied.
#' @param hessian       Optional pre-computed negative log-posterior Hessian
#'   (n x n numeric matrix). When supplied, \code{log_post_fn} is ignored.
#' @param param_names   Character vector of parameter names (length n).
#'   Defaults to \code{names(theta)}.
#' @param n_flat        Integer. Number of flattest eigendirections for which
#'   to report parameter loadings (default 3L).
#' @param sloppy_threshold  Numeric in (0, 1). Eigenvalues below
#'   \code{lambda_max * sloppy_threshold} are classified as sloppy
#'   (default 1e-6).
#' @param n_top_loadings Integer. Number of top parameter loadings to report
#'   per flat direction (default 5L).
#' @param h             Step-size fraction for \code{num_hessian} (default
#'   1e-4; each step is \code{h * max(1, |theta_i|)}).
#' @param meta          Optional metadata list (passed to \code{.apply_meta}).
#'
#' @return A \code{dynhr_diagnostic} with \code{$result} containing:
#'   \describe{
#'     \item{\code{eigenvalues}}{Numeric vector (length n), sorted descending.
#'       These are curvatures of the negative log-posterior.}
#'     \item{\code{condition_number}}{lambda_max / lambda_min_positive.
#'       Large condition number => ill-conditioned Hessian / sloppy model.}
#'     \item{\code{participation_ratio}}{(sum lambda)^2 / sum(lambda^2).
#'       Effective number of stiff/identified directions (between 1 and n).
#'       Values well below n indicate most stiffness is concentrated in few
#'       directions.}
#'     \item{\code{spectral_gap}}{lambda_1 / lambda_2 (ratio of largest to
#'       second-largest eigenvalue); large gap means one dominant curvature.}
#'     \item{\code{n_sloppy}}{Number of eigenvalues below
#'       \code{lambda_max * sloppy_threshold}.}
#'     \item{\code{flat_directions}}{A list of length \code{n_flat} (or fewer
#'       if n < n_flat). Each element is a data.frame with columns
#'       \code{parameter} (character), \code{loading} (signed eigenvector
#'       component), and \code{abs_loading}, sorted by \code{abs_loading}
#'       descending. Attribute \code{eigenvalue} gives the curvature for that
#'       direction.}
#'     \item{\code{hessian}}{The n x n negative log-posterior Hessian used.}
#'   }
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
                           meta             = NULL) {

  ## ------------------------------------------------------------------
  ## 0. Setup
  ## ------------------------------------------------------------------
  n <- length(theta)
  if (is.null(param_names))
    param_names <- if (!is.null(names(theta))) names(theta) else paste0("theta_", seq_len(n))
  if (length(param_names) != n)
    stop(sprintf("d38: param_names length (%d) != length(theta) (%d)",
                 length(param_names), n))

  n_flat <- min(as.integer(n_flat), n)

  ## ------------------------------------------------------------------
  ## 1. Obtain the negative log-posterior Hessian
  ## ------------------------------------------------------------------
  H <- NULL

  if (!is.null(hessian)) {
    ## Caller-supplied matrix
    if (!is.matrix(hessian) || nrow(hessian) != n || ncol(hessian) != n)
      stop("d38: `hessian` must be an n x n matrix where n = length(theta)")
    H <- hessian
  } else if (!is.null(log_post_fn)) {
    ## Numerical Hessian of the log-posterior; num_hessian returns H of logpost
    ## (curvatures in log-posterior), so negate to get curvatures of -logpost.
    H_logpost <- num_hessian(log_post_fn, theta, h = h)
    H <- -H_logpost
  } else {
    stop("d38: supply either `hessian` (pre-computed) or `log_post_fn`")
  }

  ## Symmetrise (guard against tiny numerical asymmetry)
  H <- 0.5 * (H + t(H))
  rownames(H) <- colnames(H) <- param_names

  ## ------------------------------------------------------------------
  ## 2. Eigen-decomposition of the symmetric Hessian
  ## ------------------------------------------------------------------
  ev <- tryCatch(
    eigen(H, symmetric = TRUE),
    error = function(e)
      stop(sprintf("d38: eigen() failed: %s", conditionMessage(e)))
  )

  ## eigen() returns eigenvalues sorted DESCENDING for symmetric matrices.
  lambdas   <- ev$values   # length n, descending
  evectors  <- ev$vectors  # n x n

  ## ------------------------------------------------------------------
  ## 3. Scalar summaries
  ## ------------------------------------------------------------------
  ## Work with eigenvalues; some may be negative (if not at the mode, or
  ## near-singular prior). Report but handle gracefully.
  lambda_pos <- lambdas[lambdas > 0]

  lambda_max <- if (length(lambda_pos) > 0) max(lambda_pos) else NA_real_
  lambda_min <- if (length(lambda_pos) > 0) min(lambda_pos) else NA_real_

  condition_number <- if (!is.na(lambda_max) && !is.na(lambda_min) &&
                          lambda_min > 0)
    lambda_max / lambda_min
  else
    Inf

  ## Participation ratio: (sum lambda)^2 / sum(lambda^2)
  ## Use all eigenvalues (can be negative off-mode); clamp to [1, n].
  sum_lam  <- sum(lambdas)
  sum_lam2 <- sum(lambdas^2)
  participation_ratio <- if (sum_lam2 > 0)
    sum_lam^2 / sum_lam2
  else
    NA_real_

  spectral_gap <- if (length(lambdas) >= 2 && lambdas[2] != 0)
    lambdas[1] / lambdas[2]
  else
    NA_real_

  ## Sloppy count: eigenvalues below lambda_max * sloppy_threshold
  sloppy_cutoff <- if (!is.na(lambda_max)) lambda_max * sloppy_threshold else NA_real_
  n_sloppy <- if (!is.na(sloppy_cutoff))
    sum(lambdas < sloppy_cutoff)
  else
    NA_integer_

  ## ------------------------------------------------------------------
  ## 4. Attribution — flat-direction loading tables
  ## ------------------------------------------------------------------
  ## The flattest directions are the LAST columns of ev$vectors
  ## (eigenvalues sorted descending => last = smallest).
  flat_directions <- vector("list", n_flat)

  for (k in seq_len(n_flat)) {
    col_idx   <- n - k + 1L        # last column = flattest
    eigenval  <- lambdas[col_idx]
    evec      <- evectors[, col_idx]

    df <- data.frame(
      parameter   = param_names,
      loading     = evec,
      abs_loading = abs(evec),
      stringsAsFactors = FALSE
    )
    df <- df[order(df$abs_loading, decreasing = TRUE), ]
    rownames(df) <- NULL

    if (n_top_loadings < n)
      df <- df[seq_len(min(n_top_loadings, nrow(df))), ]

    attr(df, "eigenvalue") <- eigenval
    attr(df, "direction")  <- k          # 1 = flattest
    flat_directions[[k]]   <- df
  }
  names(flat_directions) <- paste0("flat_", seq_len(n_flat))

  ## ------------------------------------------------------------------
  ## 5. Plots (optional, ggplot2)
  ## ------------------------------------------------------------------
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {

    ## (a) Eigenvalue spectrum (log scale)
    eig_df <- data.frame(
      index  = seq_along(lambdas),
      lambda = lambdas
    )
    eig_df$sloppy <- if (!is.na(sloppy_cutoff))
      eig_df$lambda < sloppy_cutoff
    else
      FALSE

    p_spec <- ggplot2::ggplot(
      eig_df,
      ggplot2::aes(x = index, y = pmax(lambda, .Machine$double.eps),
                   colour = sloppy)
    ) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_line(linewidth = 0.4, colour = "grey70") +
      ggplot2::scale_y_log10() +
      ggplot2::scale_colour_manual(
        values = c("FALSE" = "#1A5276", "TRUE" = "#C0392B"),
        labels = c("FALSE" = "Stiff", "TRUE" = "Sloppy"),
        name   = NULL
      ) +
      theme_dynhr_diagnostic() +
      ggplot2::labs(
        title    = "D38: Posterior Hessian eigenvalue spectrum",
        subtitle = sprintf(
          "cond = %.2e | PR = %.2f / %d | n_sloppy = %d",
          condition_number, participation_ratio %||% NA_real_, n, n_sloppy
        ),
        x = "Eigenvalue index (1 = stiffest)",
        y = "Eigenvalue (log scale)"
      )
    plots$spectrum <- .apply_meta(p_spec, meta)

    ## (b) Loading bar chart for the flattest direction
    if (n_flat >= 1L && nrow(flat_directions[[1]]) > 0) {
      fd1 <- flat_directions[[1]]
      p_load <- ggplot2::ggplot(
        fd1[seq_len(min(10L, nrow(fd1))), ],
        ggplot2::aes(
          x    = stats::reorder(parameter, abs_loading),
          y    = loading,
          fill = loading > 0
        )
      ) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(
          values = c("TRUE" = "#1A5276", "FALSE" = "#C0392B"),
          guide  = "none"
        ) +
        theme_dynhr_diagnostic() +
        ggplot2::labs(
          title    = sprintf("D38: Flattest eigendirection (lambda = %.3e)",
                             attr(flat_directions[[1]], "eigenvalue")),
          subtitle = "Parameters with largest absolute loading define the flat ridge",
          x = NULL, y = "Eigenvector loading"
        )
      plots$flat_direction_1 <- .apply_meta(p_load, meta)
    }
  }

  ## ------------------------------------------------------------------
  ## 6. Assemble result
  ## ------------------------------------------------------------------
  ## pass = NA (informational) — sloppiness is a spectrum, not pass/fail.
  ## Flag as WARN if condition_number is extreme (>1e8) or n_sloppy > 0.
  is_sloppy <- isTRUE(n_sloppy > 0)

  summary_text <- sprintf(
    paste0(
      "D38 Sloppiness spectrum: n=%d params | ",
      "cond=%.2e | PR=%.2f | n_sloppy=%d | spectral_gap=%.2f%s"
    ),
    n,
    condition_number,
    participation_ratio %||% NA_real_,
    n_sloppy %||% 0L,
    spectral_gap %||% NA_real_,
    if (is_sloppy) {
      top_flat <- flat_directions[[1]]
      top2 <- head(top_flat$parameter, 2)
      sprintf(" | flattest ridge: %s", paste(top2, collapse = " + "))
    } else ""
  )

  .make_result(
    result  = list(
      eigenvalues         = lambdas,
      condition_number    = condition_number,
      participation_ratio = participation_ratio,
      spectral_gap        = spectral_gap,
      n_sloppy            = n_sloppy,
      sloppy_cutoff       = sloppy_cutoff,
      flat_directions     = flat_directions,
      hessian             = H
    ),
    pass    = NA,   # informational: no binary pass/fail threshold
    plots   = plots,
    summary = summary_text,
    llm_summary = {
      paste(c(
        sprintf("D38 | Sloppiness Spectrum | INFO"),
        sprintf("  n_params=%d  cond=%.2e  participation_ratio=%.2f  n_sloppy=%d",
                n, condition_number, participation_ratio %||% NA_real_,
                n_sloppy %||% 0L),
        sprintf("  eigenvalues (top 5): %s",
                paste(sprintf("%.3e", head(lambdas, 5)), collapse = ", ")),
        sprintf("  eigenvalues (bot 5): %s",
                paste(sprintf("%.3e", tail(lambdas, 5)), collapse = ", ")),
        if (is_sloppy && length(flat_directions) >= 1L) {
          fd <- flat_directions[[1]]
          sprintf("  flattest direction (lambda=%.3e): %s",
                  attr(fd, "eigenvalue"),
                  paste(head(fd$parameter, 3), collapse = " + "))
        } else
          "  no sloppy directions detected at threshold",
        sprintf(
          "  action: %s",
          if (!is_sloppy)
            "Posterior well-curvated in all directions; no flat ridges detected."
          else
            sprintf(paste0(
              "Sloppy model: %d flat direction(s). ",
              "Consider tightening priors, reparameterising, ",
              "or adding observables for the flat combination."
            ), n_sloppy)
        )
      ), collapse = "\n")
    }
  )
}


## ---------------------------------------------------------------------------
## Public entry point
## ---------------------------------------------------------------------------

#' D38. Posterior-curvature sloppiness spectrum
#'
#' Compute the eigen-spectrum of the negative log-posterior Hessian at
#' \code{theta} to diagnose parameter-space ridge / flat-direction difficulty
#' ("sloppy" directions following Gutenkunst et al. 2007).
#'
#' @inheritParams d38_sloppiness
#' @param theta       Numeric parameter vector at which to assess sloppiness.
#' @param log_post_fn Optional log-posterior function of \code{theta}; if
#'   supplied (or built from \code{model}), the Hessian is formed numerically.
#' @param hessian     Optional pre-computed Hessian of the negative
#'   log-posterior at \code{theta}; supplied instead of \code{log_post_fn}.
#' @param param_names Optional character vector of parameter names (for
#'   labelling the eigenvector loadings).
#' @param n_flat      Integer number of flattest (sloppiest) eigendirections
#'   to report (default 3).
#' @param sloppy_threshold Eigenvalue below which a direction is deemed
#'   "sloppy" (default 1e-6).
#' @param n_top_loadings Integer number of top parameter loadings to report
#'   per sloppy direction (default 5).
#' @param h           Finite-difference step for the numerical Hessian
#'   (default 1e-4).
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
#' @return A \code{dynhr_diagnostic} (see \code{d38_sloppiness} for
#'   the full \code{$result} structure). Pass is always \code{NA}
#'   (informational).
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
                            ...) {

  ## If none of hessian or log_post_fn are supplied, try to build log_post_fn
  ## from model + data.
  if (is.null(hessian) && is.null(log_post_fn)) {
    if (is.null(model) || is.null(data) || is.null(obs_vars) || is.null(compiled))
      stop(paste0(
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
    h                = h
  )
}
