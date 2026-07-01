## R/bayesian-cfcst.R
## --------------------------------------------------------------------------
## bayesian_conditional_forecast() -- conditional forecast distribution from
## posterior draws.  Mirrors the diag_bayesian_irf pattern (R/diag-bayesian-irf.R):
## subsample draws, re-solve perturbation per draw, call conditional_forecast
## per draw, return per-draw paths + quantile summaries.
##
## ROADMAP Tier 8 item 10, Phase 2.
## --------------------------------------------------------------------------


#' Bayesian conditional forecast from posterior draws
#'
#' Re-solves the model at a subsample of posterior parameter draws and
#' collects conditional forecast paths at each draw, returning a distribution
#' over forecast trajectories.  Mirrors the internal
#' \code{diag_bayesian_irf} pattern for impulse responses.
#'
#' @param draws         Posterior draws matrix (\eqn{n\_draws \times n\_params});
#'   column names must match estimated parameter names from the model.
#' @param model         dynhr_mod from \code{parse_mod()}.
#' @param compiled      dynhr_compiled from \code{compile_model()}.
#' @param Y             Observed data matrix (\eqn{T \times n\_obs}) or a
#'   named list of observable variable matrices.  Passed as-is to
#'   \code{\link{conditional_forecast}}.
#' @param plan          A \code{\link{dynhr_plan}} carrying out-of-sample
#'   conditions (required).  The plan's \code{plan_condition} entries specify
#'   which variables and horizons are conditioned.  Tune and shock-scale
#'   entries in the plan are ignored on the forecast side (consistent with
#'   \code{conditional_forecast}'s existing behaviour).
#' @param horizon       Forecast horizon (integer, default \code{8L}).
#' @param obs_names     Character vector of observable variable names
#'   (column order in \code{Y}).  Inferred from \code{colnames(Y)} when
#'   \code{NULL}.
#' @param n_subsample   Maximum number of draws to use (default \code{400L}).
#'   If \code{nrow(draws) > n_subsample} a random subsample is taken.
#' @param ci_bands      Two-element numeric vector of lower/upper quantile
#'   probabilities for the outer credible band (default \code{c(0.10, 0.90)}).
#' @param inner_bands   Inner credible band (default \code{c(0.16, 0.84)}).
#' @param seed          Optional integer seed for the subsample draw.
#' @param ...           Additional arguments forwarded to
#'   \code{\link{conditional_forecast}} (e.g. \code{free_shocks}, \code{Q}).
#'
#' @return An object of class \code{"dynhr_bayesian_cfcst"} with:
#'   \describe{
#'     \item{\code{$paths_array}}{3-D array \eqn{n\_use \times H \times n\_obs}
#'       of per-draw forecast paths (NA for failed draws).}
#'     \item{\code{$quantiles}}{Named list (one entry per observable) of
#'       \eqn{n\_probs \times H} quantile matrices.}
#'     \item{\code{$n_ok}}{Number of draws that produced a valid forecast.}
#'     \item{\code{$n_fail}}{Number of draws that failed (perturbation
#'       did not solve or \code{conditional_forecast} errored).}
#'     \item{\code{$obs_names}}{Observable names.}
#'     \item{\code{$probs}}{Probability levels used for quantiles.}
#'     \item{\code{$plan}}{The plan supplied by the caller.}
#'   }
#'
#' @seealso \code{\link{conditional_forecast}}, \code{\link{dynhr_plan}}
#' @export
bayesian_conditional_forecast <- function(
    draws,
    model,
    compiled,
    Y,
    plan,
    horizon      = 8L,
    obs_names    = NULL,
    n_subsample  = 400L,
    ci_bands     = c(0.10, 0.90),
    inner_bands  = c(0.16, 0.84),
    seed         = NULL,
    ...
) {
  stopifnot(inherits(plan, "dynhr_plan"))

  draws   <- as.matrix(draws)
  n_total <- nrow(draws)
  n_use   <- min(as.integer(n_subsample), n_total)

  if (!is.null(seed)) set.seed(seed)
  draw_idx <- if (n_use < n_total) sample(n_total, n_use) else seq_len(n_total)

  Y <- as.matrix(Y)
  if (is.null(obs_names)) {
    if (!is.null(colnames(Y))) {
      obs_names <- colnames(Y)
    } else {
      stop("bayesian_conditional_forecast: obs_names must be supplied when Y has no colnames.",
           call. = FALSE)
    }
  }

  H         <- as.integer(horizon)
  par_names <- colnames(draws)
  model_pars <- names(model$param_values)

  n_obs <- length(obs_names)

  ## Storage: [draw, horizon, variable]
  paths_array <- array(
    NA_real_,
    dim      = c(n_use, H, n_obs),
    dimnames = list(NULL, paste0("h", seq_len(H)), obs_names)
  )

  n_ok   <- 0L
  n_fail <- 0L
  m_work <- model  ## working copy

  for (ki in seq_len(n_use)) {
    theta_k <- draws[draw_idx[ki], ]

    ## Update model parameters
    pv <- m_work$param_values
    for (nm in par_names) {
      if (nm %in% model_pars) pv[[nm]] <- theta_k[[nm]]
    }
    m_work$param_values <- pv

    ## Re-solve at this draw
    ss_k <- tryCatch(
      suppressWarnings(
        solve_steady(compiled, pv,
                     endo_names = model$var_names,
                     exo_names  = model$varexo_names,
                     verbose    = FALSE)
      ),
      error = function(e) NULL
    )
    if (is.null(ss_k) || !isTRUE(ss_k$converged)) {
      n_fail <- n_fail + 1L
      next
    }

    dr_k <- tryCatch(
      suppressWarnings(
        solve_perturbation(m_work, compiled, ss_k$values, pv, verbose = FALSE)
      ),
      error = function(e) NULL
    )
    if (is.null(dr_k) || !isTRUE(dr_k$bk_satisfied)) {
      n_fail <- n_fail + 1L
      next
    }

    ## Conditional forecast at this draw
    cfcst_k <- tryCatch(
      suppressWarnings(
        conditional_forecast(m_work, dr_k, Y,
                             plan      = plan,
                             horizon   = H,
                             obs_names = obs_names,
                             ...)
      ),
      error   = function(e) NULL,
      warning = function(w) {
        ## Re-try silently: some shock-scale warnings are non-fatal
        tryCatch(
          suppressWarnings(
            conditional_forecast(m_work, dr_k, Y,
                                 plan      = plan,
                                 horizon   = H,
                                 obs_names = obs_names,
                                 ...)
          ),
          error = function(e2) NULL
        )
      }
    )
    if (is.null(cfcst_k) || is.null(cfcst_k$paths_point)) {
      n_fail <- n_fail + 1L
      next
    }

    ## paths_point is H x n_obs (rows = horizons, cols = obs)
    pp <- cfcst_k$paths_point
    if (!is.matrix(pp) || nrow(pp) != H || ncol(pp) != n_obs) {
      n_fail <- n_fail + 1L
      next
    }
    paths_array[ki, , ] <- pp
    n_ok <- n_ok + 1L
  }

  ## Restore model parameters
  m_work$param_values <- model$param_values

  ## ---- Quantile summaries ------------------------------------------------
  probs_all <- sort(unique(c(ci_bands, inner_bands, 0.50)))

  quantiles <- lapply(seq_len(n_obs), function(vi) {
    mat_vi <- paths_array[, , vi, drop = FALSE]
    mat_vi <- matrix(mat_vi, nrow = n_use, ncol = H)  ## n_use x H
    finite_rows <- rowSums(is.finite(mat_vi)) == H
    if (sum(finite_rows) < 2L) return(NULL)
    apply(mat_vi[finite_rows, , drop = FALSE], 2L,
          quantile, probs = probs_all, na.rm = TRUE)
  })
  names(quantiles) <- obs_names

  structure(
    list(
      paths_array = paths_array,
      quantiles   = quantiles,
      n_ok        = n_ok,
      n_fail      = n_fail,
      n_use       = n_use,
      obs_names   = obs_names,
      probs       = probs_all,
      plan        = plan
    ),
    class = c("dynhr_bayesian_cfcst", "list")
  )
}


#' Print method for dynhr_bayesian_cfcst
#' @param x   A \code{dynhr_bayesian_cfcst} object.
#' @param ... Unused.
#' @export
#' @noRd
print.dynhr_bayesian_cfcst <- function(x, ...) {
  cat("<dynhr_bayesian_cfcst>\n")
  cat(sprintf("  Draws used  : %d / %d ok (%d failed)\n",
              x$n_ok, x$n_use, x$n_fail))
  cat(sprintf("  Horizon     : %d\n", dim(x$paths_array)[2L]))
  cat(sprintf("  Observables : %s\n", paste(x$obs_names, collapse = ", ")))
  cat(sprintf("  Quantiles   : %s\n",
              paste(sprintf("%.0f%%", x$probs * 100), collapse = ", ")))
  invisible(x)
}
