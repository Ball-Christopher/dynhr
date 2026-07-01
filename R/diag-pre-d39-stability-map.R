## R/diag-pre-d39-stability-map.R
## --------------------------------------------------------------------------
## D39: BK-feasibility / determinacy stability mapping (Ratto 2008 GSA;
##      Monte-Carlo filtering).
##
## Diagnosis for models where only a fraction of the prior space satisfies
## Blanchard-Kahn conditions, by mapping which parameters drive infeasibility.
## --------------------------------------------------------------------------

#' D39. BK-feasibility stability mapping (Ratto 2008 / Monte-Carlo filtering)
#'
#' Samples \code{n_draws} parameter vectors from the joint prior and
#' classifies each as BK-feasible or infeasible by attempting to solve
#' the model's steady state and first-order perturbation. A two-sample
#' Kolmogorov-Smirnov test (comparing each parameter's marginal distribution
#' across feasible vs infeasible draws) ranks parameters by how strongly they
#' drive infeasibility.
#'
#' This is the Ratto (2008) Monte-Carlo filtering / GSA approach applied to
#' the Blanchard-Kahn determinacy condition. A large KS statistic for a
#' parameter means its value strongly separates the feasible and infeasible
#' regions: the prior for that parameter is poorly positioned relative to the
#' determinacy boundary.
#'
#' @param model      Parsed model object (from \code{parse_mod()}).
#' @param compiled   Compiled model (from \code{compile_model()}). If
#'                   \code{NULL}, compiled internally.
#' @param prior_spec \code{data.frame} with columns \code{name},
#'                   \code{distribution}, \code{p1}, \code{p2},
#'                   \code{lower}, \code{upper} (as returned by
#'                   \code{extract_prior_spec()}). If \code{NULL}, extracted
#'                   from \code{model}.
#' @param n_draws    Number of prior draws (default \code{2000L}).
#' @param seed       Optional integer random seed.
#' @param verbose    Logical; print progress (default \code{FALSE}).
#' @param ...        Passed to \code{compile_model()} if compiling internally.
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{feasible_fraction}}{Scalar in \eqn{[0,1]}: fraction of draws
#'     that are BK-feasible.}
#'   \item{\code{drivers}}{data.frame sorted by \code{ks_stat} (descending):
#'     \code{param}, \code{ks_stat}, \code{p_value}, \code{feasible_lo},
#'     \code{feasible_hi}.}
#'   \item{\code{n_feasible}}{Integer count of feasible draws.}
#'   \item{\code{n_total}}{Integer count of attempted draws (= \code{n_draws}).}
#' }
#'
#' @references
#'   Ratto, M. (2008). Analysing DSGE models with global sensitivity analysis.
#'   \emph{Computational Economics}, 31(2), 115--139.
#'
#'   Saltelli, A., Ratto, M., Andres, T., Campolongo, F., Cariboni, J.,
#'   Gatelli, D., Saisana, M., & Tarantola, S. (2008).
#'   \emph{Global Sensitivity Analysis: The Primer}. John Wiley & Sons.
#'
#' @export
diag_stability_map <- function(model,
                               compiled   = NULL,
                               prior_spec = NULL,
                               n_draws    = 2000L,
                               seed       = NULL,
                               verbose    = FALSE,
                               ...) {

  if (!is.null(seed)) set.seed(seed)

  ## -- Prior spec ----------------------------------------------------------
  if (is.null(prior_spec)) {
    prior_spec <- extract_prior_spec(model, verbose = FALSE)
  }
  if (is.null(prior_spec) || nrow(prior_spec) == 0L) {
    stop("diag_stability_map: no prior_spec available. ",
         "Supply prior_spec= or ensure model has an estimated_params block.")
  }
  param_names <- prior_spec$name
  n_par       <- length(param_names)

  ## -- Compiled model ------------------------------------------------------
  if (is.null(compiled)) {
    if (verbose) message("diag_stability_map: compiling model ...")
    compiled <- compile_model(model, verbose = FALSE, max_order = 1L, ...)
  }

  ## -- Base parameter vector (calibrated values) ---------------------------
  base_params <- model$param_values
  if (is.null(base_params)) base_params <- numeric(0)

  ## -- Build prior sampler -------------------------------------------------
  ## Build one sampler closure per parameter, using p1/p2 directly so that
  ## the distribution-specific semantics are respected regardless of whether
  ## `mean`/`std` helper columns are present (they can shadow p1/p2 in
  ## .smc_make_prior_sampler() for uniform distributions).
  prior_draw_fn <- .d39_make_prior_sampler(prior_spec)

  ## -- Storage for draws and feasibility labels ----------------------------
  draw_mat  <- matrix(NA_real_, nrow = n_draws, ncol = n_par,
                      dimnames = list(NULL, param_names))
  feasible  <- logical(n_draws)

  n_feasible <- 0L

  if (verbose) message(sprintf("diag_stability_map: sampling %d draws ...", n_draws))

  for (i in seq_len(n_draws)) {
    ## Draw from prior: returns a named numeric vector of length n_par
    theta_i <- prior_draw_fn()
    ## Align to param_names order (prior_draw_fn may not match order)
    theta_i <- theta_i[param_names]
    draw_mat[i, ] <- theta_i

    ## Build full parameter vector: start from calibrated values, override
    ## estimated parameters with the draw.
    params_i <- base_params
    for (nm in param_names) {
      if (nm %in% names(params_i)) {
        params_i[nm] <- theta_i[nm]
      } else {
        params_i <- c(params_i, stats::setNames(theta_i[nm], nm))
      }
    }

    ## Try to solve SS + perturbation; classify BK feasibility.
    ## NOTE: solve_perturbation() issues a *warning* (not an error) for BK
    ## violations and returns a DR with bk_satisfied=FALSE.  We must NOT
    ## intercept warnings with tryCatch's warning= handler, because that would
    ## abort before the return value is produced. Use withCallingHandlers to
    ## suppress the BK warning but let execution continue, then read $bk_satisfied.
    ok <- tryCatch({
      ss_i <- solve_steady(compiled, params_i, verbose = FALSE)
      if (is.null(ss_i) || !isTRUE(ss_i$converged)) FALSE
      else {
        dr_i <- withCallingHandlers(
          solve_perturbation(model, compiled, ss_i$values, params_i,
                             verbose = FALSE),
          warning = function(w) invokeRestart("muffleWarning")
        )
        isTRUE(dr_i$bk_satisfied)
      }
    }, error = function(e) FALSE)

    feasible[i] <- ok
    if (ok) n_feasible <- n_feasible + 1L
  }

  feasible_fraction <- n_feasible / n_draws

  if (verbose) {
    message(sprintf("diag_stability_map: feasible fraction = %.3f  (%d / %d)",
                    feasible_fraction, n_feasible, n_draws))
  }

  ## -- KS test: rank parameters by separation power -----------------------
  feas_idx   <- which(feasible)
  infeas_idx <- which(!feasible)

  drivers <- .d39_ks_ranking(draw_mat, feas_idx, infeas_idx, param_names)

  structure(
    list(
      feasible_fraction = feasible_fraction,
      drivers           = drivers,
      n_feasible        = n_feasible,
      n_total           = n_draws
    ),
    class = c("dynhr_stability_map", "list")
  )
}


## --------------------------------------------------------------------------
## .d39_make_prior_sampler: returns function() -> named numeric
##
## Builds samplers from p1/p2 directly to avoid the mean/std column
## shadowing issue in .smc_make_prior_sampler for uniform distributions.
## --------------------------------------------------------------------------

.d39_make_prior_sampler <- function(prior_spec) {
  samplers <- vector("list", nrow(prior_spec))
  names(samplers) <- prior_spec$name

  for (i in seq_len(nrow(prior_spec))) {
    nm   <- prior_spec$name[i]
    dist <- tolower(prior_spec$distribution[i])
    p1   <- prior_spec$p1[i]
    p2   <- prior_spec$p2[i]
    lo   <- prior_spec$lower[i]
    hi   <- prior_spec$upper[i]
    if (is.na(lo))  lo  <- -Inf
    if (is.na(hi))  hi  <- Inf

    samplers[[nm]] <- local({
      d <- dist; a <- p1; b <- p2; lb <- lo; ub <- hi
      switch(d,
        "uniform" =, "unif" = {
          ## p1=lower, p2=upper for uniform
          function() runif(1, a, b)
        },
        "beta" = {
          ## p1=mean, p2=sd (Dynare convention on [lb, ub])
          lo2 <- if (is.finite(lb)) lb else 0
          hi2 <- if (is.finite(ub)) ub else 1
          m01 <- (a - lo2) / (hi2 - lo2)
          s01 <- b / (hi2 - lo2)
          v   <- m01 * (1 - m01) / s01^2 - 1
          v   <- max(v, 2)
          sh1 <- m01 * v; sh2 <- (1 - m01) * v
          function() lo2 + (hi2 - lo2) * rbeta(1, sh1, sh2)
        },
        "gamma" =, "gamm" = {
          shape <- (a / b)^2; rate <- a / b^2
          lo2   <- if (is.finite(lb)) lb else 0
          function() max(lo2, rgamma(1, shape = shape, rate = rate))
        },
        "inv_gamma" =, "invg" =, "inv_gamma1" =, "inv_gamma2" = {
          alpha  <- (a / b)^2 + 2
          beta_p <- a * (alpha - 1)
          lo2    <- if (is.finite(lb)) lb else 0
          function() max(lo2, 1 / rgamma(1, shape = alpha, rate = beta_p))
        },
        ## Default: normal
        {
          function() max(lb, min(ub, rnorm(1, a, b)))
        }
      )
    })
  }

  ## Return a single function that draws ALL params at once.
  function() {
    theta <- vapply(names(samplers), function(nm) samplers[[nm]](), numeric(1))
    stats::setNames(theta, names(samplers))
  }
}


## --------------------------------------------------------------------------
## .d39_ks_ranking: two-sample KS test per parameter
## --------------------------------------------------------------------------

.d39_ks_ranking <- function(draw_mat, feas_idx, infeas_idx, param_names) {

  n_par <- length(param_names)

  ## Guard: need at least 1 obs in each group to run KS.
  if (length(feas_idx) == 0L || length(infeas_idx) == 0L) {
    ## Return a zero-stat table
    return(data.frame(
      param        = param_names,
      ks_stat      = rep(NA_real_, n_par),
      p_value      = rep(NA_real_, n_par),
      feasible_lo  = rep(NA_real_, n_par),
      feasible_hi  = rep(NA_real_, n_par),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(param_names, function(nm) {
    x_feas   <- draw_mat[feas_idx,   nm]
    x_infeas <- draw_mat[infeas_idx, nm]

    ## Remove non-finite values (failed draws may have left NA)
    x_feas   <- x_feas[is.finite(x_feas)]
    x_infeas <- x_infeas[is.finite(x_infeas)]

    if (length(x_feas) < 2L || length(x_infeas) < 2L) {
      return(data.frame(param = nm, ks_stat = NA_real_, p_value = NA_real_,
                        feasible_lo = NA_real_, feasible_hi = NA_real_,
                        stringsAsFactors = FALSE))
    }

    kt <- suppressWarnings(stats::ks.test(x_feas, x_infeas))
    data.frame(
      param        = nm,
      ks_stat      = unname(kt$statistic),
      p_value      = kt$p.value,
      feasible_lo  = min(x_feas),
      feasible_hi  = max(x_feas),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out <- out[order(out$ks_stat, decreasing = TRUE, na.last = TRUE), ]
  rownames(out) <- NULL
  out
}


## --------------------------------------------------------------------------
## print method
## --------------------------------------------------------------------------

#' @export
print.dynhr_stability_map <- function(x, n_top = 10L, ...) {
  cat(sprintf(
    "D39 BK-feasibility stability map\n  feasible: %d / %d (%.1f%%)\n",
    x$n_feasible, x$n_total, 100 * x$feasible_fraction
  ))
  if (!is.null(x$drivers) && nrow(x$drivers) > 0L) {
    cat(sprintf("  top drivers (by KS statistic, showing up to %d):\n", n_top))
    top <- head(x$drivers, n_top)
    for (i in seq_len(nrow(top))) {
      r <- top[i, ]
      cat(sprintf("    %-20s  KS=%.3f  p=%.3g  feasible=[%.4g, %.4g]\n",
                  r$param, r$ks_stat, r$p_value, r$feasible_lo, r$feasible_hi))
    }
  }
  invisible(x)
}
