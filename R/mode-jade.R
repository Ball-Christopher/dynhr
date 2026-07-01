## R/mode-jade.R
## --------------------------------------------------------------------------
## Phase-2 split from estimate-monolith.R.
##
## jade_optimize() -- JADE (adaptive differential evolution) mode finder.
## Uses .reflect_bounds() from mode-cmaes.R (loaded together in the package).
## --------------------------------------------------------------------------

#' JADE (adaptive differential evolution) mode finder
#'
#' @param fn         Objective to MINIMISE (return scalar numeric)
#' @param par        Named starting vector (seeds individual 1 of the population)
#' @param lower,upper  Box bounds (scalar or vector)
#' @param max_iter   Generation budget
#' @param NP         Population size (NULL = max(20, 5*n))
#' @param tol_f      Convergence tolerance on best-value spread
#' @param p_best     Fraction of best individuals for current-to-p-best mutation
#' @param c_adapt    Learning rate for mu_CR / mu_F adaptation
#' @param verbose    Print progress messages
#' @param progress   Show cli progress bar
#' @return list(par, value, convergence, iterations, message)
#' @noRd
jade_optimize <- function(fn, par, lower = -Inf, upper = Inf,
                          max_iter = 10000, NP = NULL,
                          tol_f = 1e-12, p_best = 0.1, c_adapt = 0.1,
                          verbose = TRUE, progress = TRUE) {
  n <- length(par)
  par_names <- names(par)

  if (length(lower) == 1) lower <- rep(lower, n)
  if (length(upper) == 1) upper <- rep(upper, n)
  if (is.null(NP)) NP <- max(20, 5 * n)

  fn_safe <- function(x) {
    names(x) <- par_names
    fn(x)
  }

  pop <- matrix(0, n, NP)
  for (j in seq_len(NP)) {
    for (i in seq_len(n)) {
      lo_i <- if (is.finite(lower[i])) lower[i] else par[i] - 5 * max(abs(par[i]), 1)
      hi_i <- if (is.finite(upper[i])) upper[i] else par[i] + 5 * max(abs(par[i]), 1)
      pop[i, j] <- runif(1, lo_i, hi_i)
    }
  }
  pop[, 1] <- par

  f_pop <- numeric(NP)
  for (j in seq_len(NP)) f_pop[j] <- fn_safe(pop[, j])

  mu_CR <- 0.5
  mu_F  <- 0.5
  archive <- matrix(0, n, 0)

  best_idx <- which.min(f_pop)
  best_val <- f_pop[best_idx]
  best_par <- pop[, best_idx]
  if (verbose) cat(sprintf("  JADE init: f=%.4f (logpost=%.4f) NP=%d\n",
                           best_val, -best_val, NP))

  pb <- NULL
  .val <- sprintf("%.4f", -best_val)
  if (progress) {
    pb <- NULL
  }
  report_every <- max(1, floor(max_iter / 20))

  f_hist    <- numeric(0)
  converged <- FALSE

  for (gen in seq_len(max_iter)) {
    S_CR <- numeric(0)
    S_F  <- numeric(0)

    rank_idx <- order(f_pop)
    p_count  <- max(2, round(p_best * NP))

    new_pop <- pop
    new_f   <- f_pop

    for (j in seq_len(NP)) {
      CR_j <- max(0, min(1, rnorm(1, mu_CR, 0.1)))
      F_j  <- rcauchy(1, mu_F, 0.1)
      iter_f <- 0
      while (F_j <= 0 && iter_f < 20) { F_j <- rcauchy(1, mu_F, 0.1); iter_f <- iter_f + 1 }
      if (F_j <= 0) F_j <- 0.5
      F_j <- min(F_j, 1)

      pbest <- rank_idx[sample.int(p_count, 1)]

      r1 <- j; while (r1 == j) r1 <- sample.int(NP, 1)

      n_arch <- ncol(archive)
      r2 <- j
      attempts <- 0
      while ((r2 == j || r2 == r1) && attempts < 50) {
        r2 <- sample.int(NP + n_arch, 1)
        attempts <- attempts + 1
      }
      if (r2 <= NP) {
        x_r2 <- pop[, r2]
      } else if (n_arch > 0) {
        x_r2 <- archive[, min(r2 - NP, n_arch)]
      } else {
        x_r2 <- pop[, r1]
      }

      v <- pop[, j] + F_j * (pop[, pbest] - pop[, j]) + F_j * (pop[, r1] - x_r2)

      u <- pop[, j]
      j_rand <- sample.int(n, 1)
      for (i in seq_len(n)) {
        if (runif(1) < CR_j || i == j_rand) u[i] <- v[i]
      }

      u <- .reflect_bounds(u, lower, upper)

      f_u <- fn_safe(u)
      if (is.finite(f_u) && f_u <= f_pop[j]) {
        new_pop[, j] <- u
        new_f[j] <- f_u
        if (n_arch < NP) {
          archive <- cbind(archive, pop[, j])
        } else if (NP > 0) {
          archive[, sample.int(NP, 1)] <- pop[, j]
        }
        S_CR <- c(S_CR, CR_j)
        S_F  <- c(S_F, F_j)
      }
    }

    pop   <- new_pop
    f_pop <- new_f

    best_idx <- which.min(f_pop)
    if (f_pop[best_idx] < best_val) {
      best_val <- f_pop[best_idx]
      best_par <- pop[, best_idx]
    }

    if (length(S_CR) > 0) {
      mu_CR <- (1 - c_adapt) * mu_CR + c_adapt * mean(S_CR)
      mu_F  <- (1 - c_adapt) * mu_F + c_adapt * sum(S_F^2) / sum(S_F)
    }

    if (!is.null(pb)) {
      .val <- sprintf("%.4f", -best_val)
      cli::cli_progress_update(id = pb, .envir = environment())
    } else if (verbose && gen %% report_every == 0) {
      cat(sprintf("  gen %5d/%d  logpost=%.4f  mu_F=%.3f  mu_CR=%.3f\n",
                  gen, max_iter, -best_val, mu_F, mu_CR))
    }

    f_hist <- c(f_hist, best_val)
    if (length(f_hist) > 20) {
      recent <- tail(f_hist, 20)
      if (max(recent) - min(recent) < tol_f) { converged <- TRUE; break }
    }
    spread <- max(f_pop) - min(f_pop)
    if (spread < tol_f * 10) { converged <- TRUE; break }
  }

  if (!is.null(pb)) cli::cli_progress_done(id = pb)

  names(best_par) <- par_names
  list(par = best_par, value = best_val,
       convergence = as.integer(!converged),
       iterations = gen,
       message = if (converged) "converged" else "max_iter reached")
}
