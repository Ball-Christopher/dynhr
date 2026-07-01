## R/mode-cmaes.R
## --------------------------------------------------------------------------
## Phase-2 split from estimate-monolith.R.
##
## .reflect_bounds()  -- box-constraint reflection (shared with mode-jade.R)
## cmaes_optimize()   -- CMA-ES wrapper over the 'cmaes' package
## .cmaes_core()      -- legacy custom CMA-ES (kept for reference / fallback)
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


## Legacy custom CMA-ES implementation (Hansen 2016). Kept for reference and
## as a fallback when the 'cmaes' package is unavailable. Not called by
## cmaes_optimize() above; used directly by estimate-monolith.R pipelines
## that pre-date the package wrapper.
.cmaes_core <- function(fn, par, n, sigma, lower, upper,
                        max_iter, tol_f, tol_x,
                        pop_mult = 1,
                        best_val_init = Inf, best_par_init = NULL,
                        verbose = TRUE, progress = TRUE,
                        restart_id = 0) {

  lambda <- max(6, floor((4 + floor(3 * log(n))) * pop_mult))
  mu     <- floor(lambda / 2)

  raw_w   <- log(mu + 0.5) - log(seq_len(mu))
  weights <- raw_w / sum(raw_w)
  mu_eff  <- 1 / sum(weights^2)

  cs    <- (mu_eff + 2) / (n + mu_eff + 5)
  ds    <- 1 + 2 * max(0, sqrt((mu_eff - 1) / (n + 1)) - 1) + cs
  E_chi <- sqrt(n) * (1 - 1/(4*n) + 1/(21*n^2))

  cc    <- (4 + mu_eff / n) / (n + 4 + 2 * mu_eff / n)
  c1    <- 2 / ((n + 1.3)^2 + mu_eff)
  c_mu  <- min(1 - c1, 2 * (mu_eff - 2 + 1/mu_eff) / ((n + 2)^2 + mu_eff))

  m     <- par
  C     <- diag(n)
  ps    <- numeric(n)
  pc    <- numeric(n)
  eigencount <- 0
  B     <- diag(n)
  D     <- rep(1, n)
  invsqrtC <- diag(n)

  best_val <- best_val_init
  best_par <- if (!is.null(best_par_init)) best_par_init else par

  eigen_freq <- max(1, floor(n / (10 * lambda)))

  pb <- NULL
  if (progress) {
    pb <- NULL
    .val <- sprintf("%.4f", -best_val)
    .sigma <- sprintf("%.3g", sigma)
  }
  report_every <- max(1, floor(max_iter / 20))

  converged <- FALSE
  gen <- 0
  f_hist <- numeric(0)

  for (gen in seq_len(max_iter)) {

    samp <- 5*lambda
    arz <- matrix(rnorm(n * samp), n, samp)
    ary <- B %*% diag(D, nrow = n) %*% arz
    arx <- m + sigma * ary

    if (gen == 1) arx[, 1] <- .reflect_bounds(m, lower, upper)

    for (k in 2:samp)
      arx[, k] <- .reflect_bounds(arx[, k], lower, upper)

    f_vals <- numeric(samp)
    for (k in seq_len(samp)) f_vals[k] <- fn(arx[, k])

    idx    <- order(f_vals)
    f_vals <- f_vals[idx]
    arx    <- arx[, idx, drop = FALSE]
    ary    <- ary[, idx, drop = FALSE]
    arz    <- arz[, idx, drop = FALSE]

    if (f_vals[1] < best_val) {
      best_val <- f_vals[1]
      best_par <- arx[, 1]
    }

    m_old <- m
    m     <- drop(arx[, 1:mu, drop = FALSE] %*% weights)
    y_w   <- (m - m_old) / sigma

    ps <- (1 - cs) * ps + sqrt(cs * (2 - cs) * mu_eff) * (invsqrtC %*% y_w)
    h_sig <- as.numeric(sum(ps^2) / n /
                          (1 - (1 - cs)^(2 * gen)) < 2 + 4/(n+1))

    pc <- (1 - cc) * pc + h_sig * sqrt(cc * (2 - cc) * mu_eff) * y_w

    y_sel <- ary[, 1:mu, drop = FALSE]
    C <- (1 - c1 - c_mu) * C +
      c1 * (tcrossprod(pc) + (1 - h_sig) * cc * (2 - cc) * C) +
      c_mu * (y_sel %*% diag(weights, nrow = mu) %*% t(y_sel))

    sigma <- sigma * exp(cs / ds * (sqrt(sum(ps^2)) / E_chi - 1))

    eigencount <- eigencount + 1
    if (eigencount >= eigen_freq) {
      eigencount <- 0
      C <- (C + t(C)) / 2
      eig <- eigen(C, symmetric = TRUE)
      if (all(eig$values > 0)) {
        B <- eig$vectors
        D <- sqrt(eig$values)
        invsqrtC <- B %*% diag(1/D, nrow = n) %*% t(B)
      } else {
        C <- diag(n); B <- diag(n); D <- rep(1, n); invsqrtC <- diag(n)
        sigma <- sigma * 2
      }
    }

    if (!is.null(pb)) {
      .val <- sprintf("%.4f", -best_val)
      .sigma <- sprintf("%.3g", sigma)
      cli::cli_progress_update(id = pb, .envir = environment())
    } else if (verbose && gen %% report_every == 0) {
      cat(sprintf("  gen %5d/%d  logpost=%.4f  sigma=%.4g\n",
                  gen, max_iter, -best_val, sigma))
    }

    f_hist <- c(f_hist, best_val)
    if (length(f_hist) > 30) {
      recent <- tail(f_hist, 30)
      if (max(recent) - min(recent) < tol_f) { converged <- TRUE; break }
    }
    if (sigma * max(D) < tol_x) { converged <- TRUE; break }
    if (max(D) / min(D) > 1e14) {
      C <- diag(n); B <- diag(n); D <- rep(1, n); invsqrtC <- diag(n)
    }
  }

  if (!is.null(pb)) cli::cli_progress_done(id = pb)

  list(par = best_par, value = best_val,
       convergence = as.integer(!converged),
       iterations = gen,
       message = if (converged) "converged" else "max_iter reached")
}
