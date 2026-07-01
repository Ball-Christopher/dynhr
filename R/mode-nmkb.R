## R/mode-nmkb.R
## --------------------------------------------------------------------------
## Phase-2 split from estimate-monolith.R.
##
## nmkb_optimize() -- Nelder-Mead with box constraints (dfoptim::nmkb)
## --------------------------------------------------------------------------

#' Nelder-Mead mode finder with box constraints
#'
#' @param fn         Objective to MINIMISE (return scalar numeric)
#' @param par        Named starting vector
#' @param lower,upper  Box bounds (scalar or vector)
#' @param max_iter   Function-evaluation budget (dfoptim's maxfeval)
#' @param tol        Convergence tolerance
#' @param verbose    Print progress messages
#' @return list(par, value, convergence, iterations, message)
#' @noRd
nmkb_optimize <- function(fn, par, lower = -Inf, upper = Inf,
                          max_iter = 5000, tol = 1e-8,
                          verbose = TRUE, ...) {
  if (!requireNamespace("dfoptim", quietly = TRUE))
    stop("Package 'dfoptim' needed. Install with: install.packages('dfoptim')")

  n         <- length(par)
  par_names <- names(par)

  if (length(lower) == 1) lower <- rep(lower, n)
  if (length(upper) == 1) upper <- rep(upper, n)

  par <- pmax(pmin(par, upper - 1e-8), lower + 1e-8)

  fn_safe <- function(x) {
    names(x) <- par_names
    val <- fn(x)
    if (!is.finite(val)) 1e20 else val
  }

  if (verbose)
    cat(sprintf("  nmkb: n=%d  maxfeval=%d  tol=%.1e\n", n, max_iter, tol))

  res <- dfoptim::nmkb(par, fn_safe, lower = lower, upper = upper,
                       control = list(maxfeval = max_iter, tol = tol))

  best_par <- setNames(res$par, par_names)
  best_val <- res$value
  n_evals  <- res$feval

  if (verbose)
    cat(sprintf("  nmkb done: logpost=%.4f  evals=%d  %s\n",
                -best_val, n_evals, res$message))

  list(par        = best_par,
       value      = best_val,
       convergence = res$convergence,
       iterations  = n_evals,
       message     = res$message %||% "done")
}
