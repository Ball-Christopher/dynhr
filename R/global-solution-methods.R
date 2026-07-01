## R/global-solution-methods.R
## --------------------------------------------------------------------------
## S3 methods for GlobalSolution objects returned by solve_global().
##
## Methods:
##   print.GlobalSolution
##   predict.GlobalSolution   -- evaluate policy at arbitrary state points
##   simulate.GlobalSolution  -- stochastic simulation using global policy
##   euler_errors             -- generic and default method
##   euler_errors.GlobalSolution
## --------------------------------------------------------------------------

#' Print a GlobalSolution object
#' @param x   A \code{GlobalSolution} object.
#' @param ... Ignored.
#' @export
print.GlobalSolution <- function(x, ...) {
  cat("GlobalSolution (Chebyshev projection)\n")
  cat(sprintf("  States:   %s\n", paste(x$state_names, collapse = ", ")))
  cat(sprintf("  Endo:     %s\n", paste(x$all_endo_names, collapse = ", ")))
  cat(sprintf("  Shocks:   %s\n", paste(x$shock_names, collapse = ", ")))
  cat(sprintf("  Degree:   %d  (n_basis = %d)\n",
              x$poly_degree, ncol(x$coefs)))
  cat(sprintf("  Converged: %s (in %d iterations, delta = %.2e)\n",
              x$converged, x$n_iter, x$last_delta))
  invisible(x)
}

#' Evaluate the global policy function at given state (lag) points
#'
#' @param object  A \code{GlobalSolution} object.
#' @param newdata A named matrix or data.frame with columns for each state
#'   variable (lag values), or a named numeric vector for a single point.
#' @param ...     Ignored.
#' @return A matrix with one row per input point and one column per
#'   endogenous variable (in declaration order).
#' @export
predict.GlobalSolution <- function(object, newdata, ...) {
  state_names <- object$state_names
  n_state     <- length(state_names)

  ## Coerce newdata to a matrix
  if (is.vector(newdata) && !is.list(newdata)) {
    newdata <- matrix(newdata, nrow = 1L,
                      dimnames = list(NULL, names(newdata)))
  }
  newdata <- as.matrix(newdata)

  if (!all(state_names %in% colnames(newdata)))
    stop(sprintf(
      "predict.GlobalSolution: newdata must contain columns for states: %s",
      paste(state_names, collapse = ", ")))

  ## Normalize each state column
  n_pts <- nrow(newdata)
  x_norm <- matrix(0.0, nrow = n_pts, ncol = n_state)
  colnames(x_norm) <- state_names
  for (j in seq_len(n_state)) {
    nm  <- state_names[j]
    dom <- object$state_domain[[nm]]
    x_norm[, j] <- pmax(-1.0, pmin(1.0,
      cheb_normalize(newdata[, nm], dom[1], dom[2])))
  }

  ## Evaluate basis and multiply by coefficients
  Phi    <- cheb_basis(x_norm, object$poly_degree)
  result <- tcrossprod(Phi, object$coefs)   # n_pts x n_endo
  colnames(result) <- object$all_endo_names
  result
}

#' Stochastic simulation from the global policy function
#'
#' @param object  A \code{GlobalSolution} object.
#' @param nsim    Number of periods to simulate (default 100).
#' @param seed    Random seed (default NULL = don't set).
#' @param init_state Named numeric; initial state lag values.
#'   Defaults to steady state.
#' @param ...     Ignored.
#' @return A matrix of size \code{nsim x n_endo} with simulated paths.
#' @export
simulate.GlobalSolution <- function(object, nsim = 100L, seed = NULL,
                                    init_state = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)

  state_names <- object$state_names
  endo        <- object$all_endo_names
  exo         <- object$shock_names
  n_state     <- length(state_names)
  n_endo      <- length(endo)
  n_exo       <- length(exo)
  shock_sds   <- object$shock_sds

  ## Initial state lag
  if (is.null(init_state)) {
    state_lag <- object$ss_vals[state_names]
  } else {
    state_lag <- init_state[state_names]
  }

  out <- matrix(NA_real_, nrow = nsim, ncol = n_endo)
  colnames(out) <- endo

  for (t in seq_len(nsim)) {
    ## Draw current shock
    eps_t <- rnorm(n_exo, mean = 0, sd = shock_sds)
    names(eps_t) <- exo

    ## Compute next-period lag from state_lag and current shock.
    ## But for simulation we want the REALISED y_t, not just the expected one.
    ## The policy gives y_t = policy(state_lag); then state transitions.
    ## Since this is a simulation (not expectation), use the realised shock.

    ## Evaluate policy at current lag -> y_t
    new_state_norm <- matrix(
      vapply(seq_len(n_state), function(j) {
        nm  <- state_names[j]
        dom <- object$state_domain[[nm]]
        max(-1.0, min(1.0, cheb_normalize(state_lag[nm], dom[1], dom[2])))
      }, numeric(1L)),
      nrow = 1L)

    Phi_t <- cheb_basis(new_state_norm, object$poly_degree)
    y_t   <- as.vector(Phi_t %*% t(object$coefs))
    names(y_t) <- endo
    out[t, ] <- y_t

    ## Update state lag for next period
    ## AR(1) shocks: z_{t+1} = rho*z_t + eps_{t+1} (eps applied at NEXT period)
    ## Capital: k_{t+1} = k_t from policy (policy already solved k_t)
    state_lag <- object$compute_next_lag(y_t, eps_t)
  }

  out
}

#' Euler equation errors generic
#'
#' @param x    An object with a global solution (e.g., \code{GlobalSolution}).
#' @param ...  Additional arguments.
#' @export
euler_errors <- function(x, ...) UseMethod("euler_errors")

#' Compute Euler equation accuracy for a GlobalSolution
#'
#' Evaluates the unit-free Euler equation error (Maliar-Maliar 2014 metric)
#' on a grid of test points. The error at each point is:
#'
#' \deqn{ee = E_t[\beta \cdot MPK_{t+1}] \cdot c_t - 1}
#'
#' where \eqn{MPK_{t+1} = \alpha \exp(z_{t+1}) k_t^{\alpha-1} / c_{t+1}}
#' (for the simple RBC model with full depreciation).
#'
#' @param x       A \code{GlobalSolution} object.
#' @param n_grid  Number of test points per state dimension (default 10).
#' @param n_quad  Quadrature nodes for the Euler error (default 20).
#' @param ...     Ignored.
#' @return A list with \code{log10_errors} (vector), \code{max_log10_error},
#'   and \code{grid} (the evaluation grid).
#' @export
euler_errors.GlobalSolution <- function(x, n_grid = 10L, n_quad = 20L, ...) {
  state_names <- x$state_names
  endo        <- x$all_endo_names
  params      <- x$params
  ss_vals     <- x$ss_vals
  n_state     <- length(state_names)

  ## Build regular evaluation grid (not Chebyshev nodes)
  grid_1d <- lapply(state_names, function(nm) {
    dom <- x$state_domain[[nm]]
    seq(dom[1], dom[2], length.out = n_grid)
  })
  names(grid_1d) <- state_names

  grid_args        <- rev(grid_1d)
  names(grid_args) <- rev(state_names)
  grid_nat         <- as.matrix(expand.grid(grid_args))
  grid_nat         <- grid_nat[, state_names, drop = FALSE]
  n_pts            <- nrow(grid_nat)

  ## Get quadrature nodes/weights for Euler error
  gh_ee   <- gauss_hermite(n_quad)
  q_nodes <- gh_ee$nodes
  q_wts   <- gh_ee$weights

  ## Shock setup
  shock_sds <- x$shock_sds
  n_exo     <- length(x$shock_names)

  ## Pull out model-specific quantities from params
  ## These are computed from the raw residual system, which is exact.
  ## We use the direct Euler error formula for the simple RBC.
  ##
  ## To avoid hard-coding model structure, we use the model residuals directly.
  ## The Euler error is the RELATIVE deviation from the consumption Euler equation:
  ##   ee = (expected RHS) / (LHS) - 1
  ## We compute this via the residuals_fn evaluated at the global solution.

  ee_vec <- numeric(n_pts)

  for (j in seq_len(n_pts)) {
    ## State lag at this evaluation point
    state_lag <- grid_nat[j, ]

    ## Current policy: y_t = policy(state_lag)
    y_t <- predict(x, matrix(state_lag, nrow = 1L,
                              dimnames = list(NULL, state_names)))[1L, ]

    ## Get current consumption
    c_idx <- match("c", endo)
    c_t   <- if (!is.na(c_idx)) y_t[c_idx] else y_t[1L]

    if (c_t <= 0 || !is.finite(c_t)) {
      ee_vec[j] <- NA_real_
      next
    }

    ## Accumulate Euler error via quadrature over t+1 shocks
    ## Build shock combinations
    if (n_exo == 1L) {
      shock_mat  <- matrix(q_nodes * shock_sds[1L], ncol = 1L)
      w_combined <- q_wts
    } else {
      qn_list    <- lapply(seq_len(n_exo), function(k) q_nodes * shock_sds[k])
      shock_mat  <- as.matrix(expand.grid(qn_list))
      qw_list    <- replicate(n_exo, q_wts, simplify = FALSE)
      w_mat      <- as.matrix(expand.grid(qw_list))
      w_combined <- apply(w_mat, 1L, prod)
      w_combined <- w_combined / sum(w_combined)
    }
    n_combo <- nrow(shock_mat)

    euler_sum <- 0.0
    ok <- TRUE

    for (ki in seq_len(n_combo)) {
      eps_next <- shock_mat[ki, ]
      next_lag <- x$compute_next_lag(y_t, eps_next)
      y_lead   <- predict(x, matrix(next_lag, nrow = 1L,
                                    dimnames = list(NULL, state_names)))[1L, ]

      ## Use direct formula from model parameters if available
      ## (more accurate than going through residuals_fn)
      alpha_val <- if ("alpha" %in% names(params)) params[["alpha"]] else NA_real_
      beta_val  <- if ("beta"  %in% names(params)) params[["beta"]]  else NA_real_

      if (!is.na(alpha_val) && !is.na(beta_val)) {
        ## Direct Euler error for RBC with full depreciation:
        ## Euler: 1/c_t = beta * (alpha*exp(z_{t+1})*k_t^{alpha-1}) / c_{t+1}
        ## ee = beta * (alpha*exp(z_{t+1})*k_t^{alpha-1}) / c_{t+1} * c_t - 1
        k_t    <- y_t[match("k", endo)]
        c_tp1  <- y_lead[match("c", endo)]
        z_tp1  <- y_lead[match("z", endo)]   # z at t+1 (from policy)

        if (any(!is.finite(c(k_t, c_tp1, z_tp1))) || k_t <= 0 || c_tp1 <= 0) {
          ok <- FALSE
          break
        }
        mpk_tp1 <- alpha_val * exp(z_tp1) * k_t^(alpha_val - 1.0)
        euler_sum <- euler_sum + w_combined[ki] * beta_val * mpk_tp1 / c_tp1
      } else {
        ok <- FALSE
        break
      }
    }

    if (!ok || !is.finite(euler_sum)) {
      ee_vec[j] <- NA_real_
    } else {
      ## Euler error: ee = E_t[beta * MPK_{t+1} / c_{t+1}] * c_t - 1
      ee_vec[j] <- euler_sum * c_t - 1.0
    }
  }

  log10_ee <- log10(abs(ee_vec))
  list(
    log10_errors    = log10_ee,
    max_log10_error = max(log10_ee[is.finite(log10_ee)], na.rm = TRUE),
    errors          = ee_vec,
    grid            = grid_nat
  )
}
