## R/mode-cmaes.R
## --------------------------------------------------------------------------
## Phase-2 split from estimate-monolith.R.
##
## .reflect_bounds()  -- box-constraint reflection (shared with mode-jade.R)
## cmaes_optimize()   -- CMA-ES wrapper over the 'cmaes' package
## --------------------------------------------------------------------------

.reflect_bounds <- function(x, lower, upper) {
  for (i in seq_along(x)) {
    lo <- lower[i]; hi <- upper[i]
    if (!is.finite(lo) && !is.finite(hi)) next
    iter <- 0
    while ((x[i] < lo || x[i] > hi) && iter < 10) {
      if (is.finite(lo) && x[i] < lo) x[i] <- lo + (lo - x[i])
      if (is.finite(hi) && x[i] > hi) x[i] <- hi - (x[i] - hi)
      iter <- iter + 1
    }
    if (is.finite(lo)) x[i] <- max(x[i], lo)
    if (is.finite(hi)) x[i] <- min(x[i], hi)
  }
  x
}


#' CMA-ES mode finder (wraps the 'cmaes' package)
#'
#' @param fn         Objective to MINIMISE (return scalar numeric)
#' @param par        Named starting vector
#' @param lower,upper  Box bounds (scalar or vector)
#' @param max_iter   Generation budget (converted to eval budget internally)
#' @param sigma0     Initial step size (NULL = auto from bounded range)
#' @param lambda,mu  CMA-ES population / parent sizes (NULL = Hansen defaults)
#' @param tol_f,tol_x  Convergence tolerances
#' @param verbose    Print progress messages
#' @param progress   Show cli progress bar
#' @return list(par, value, convergence, iterations, message)
#' @noRd
cmaes_optimize <- function(fn, par, lower = -Inf, upper = Inf,
                           max_iter = 10000, sigma0 = NULL,
                           lambda = NULL, mu = NULL,
                           tol_f = 1e-8, tol_x = 1e-8,
                           verbose = TRUE, progress = TRUE, ...) {
  if (!requireNamespace("cmaes", quietly = TRUE))
    stop("Package 'cmaes' needed. Install with: install.packages('cmaes')")

  n         <- length(par)
  par_names <- names(par)

  if (length(lower) == 1) lower <- rep(lower, n)
  if (length(upper) == 1) upper <- rep(upper, n)

  if (is.null(lambda)) lambda <- 4L + floor(3L * log(n))
  if (is.null(mu))     mu     <- floor(lambda / 2L)

  if (is.null(sigma0)) {
    ranges   <- upper - lower
    finite_r <- is.finite(ranges)
    sigma0   <- if (any(finite_r)) median(ranges[finite_r]) / 6
                else 0.1 * max(abs(par[par != 0]), 1)
  }

  fn_safe <- function(x) {
    names(x) <- par_names
    val <- fn(x)
    if (!is.finite(val)) 1e20 else val
  }

  if (verbose)
    cat(sprintf("  CMA-ES (pkg): n=%d  lambda=%d  mu=%d  sigma0=%.3g  maxit=%d\n",
                n, lambda, mu, sigma0, as.integer(max_iter)))

  ## NB: the 'cmaes' package reads `maxit` (generations) and `stop.tolx`. The
  ## previous control used `stopeval`/`sc.tolx`/`sc.tolf`, which the package
  ## does NOT recognise -- so `max_iter` was silently ignored and cma_es ran to
  ## its default tight convergence. Passing `maxit` makes `max_iter` a real cap
  ## (generations); with the usual large budgets convergence still terminates
  ## first, so default behaviour is unchanged -- the cap only bites for small
  ## explicit budgets. `tol_f` has no package equivalent (only stopfitness).
  ctrl <- list(
    maxit      = as.integer(max_iter),
    sigma      = sigma0,
    lambda     = lambda,
    mu         = mu,
    diag.sigma = FALSE,
    diag.eigen = FALSE,
    diag.pop   = FALSE,
    diag.value = FALSE
  )

  res <- cmaes::cma_es(par, fn_safe, lower = lower, upper = upper,
                       control = ctrl)

  ## cma_es returns par = NULL when it never improves on the start (e.g. a
  ## dispersed start where every nearby point is non-finite). setNames(NULL, .)
  ## would then throw "attempt to set an attribute on NULL"; fall back to the
  ## start point and its objective value instead of crashing the chain.
  if (is.null(res$par)) {
    res$par   <- par
    res$value <- res$value %||% fn_safe(par)
    res$message <- paste0(res$message %||% "", " [no improvement; kept start]")
  }
  best_par <- setNames(res$par, par_names)
  best_val <- res$value
  n_evals  <- res$counts[["function"]] %||% NA_integer_

  if (verbose)
    cat(sprintf("  CMA-ES done: logpost=%.4f  evals=%d\n", -best_val, n_evals))

  list(par        = best_par,
       value      = best_val,
       convergence = res$convergence,
       iterations  = n_evals,
       message     = res$message %||% "done")
}
