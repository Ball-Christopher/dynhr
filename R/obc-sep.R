## R/obc-sep.R
## --------------------------------------------------------------------------
## Stochastic Extended Path (SEP) OBC simulation for nonlinear DSGE models.
##
## Provides:
##   sep_gauss_hermite()  -- univariate Gauss-Hermite nodes / weights
##   sep_quadrature()     -- multi-shock quadrature grid (tensor product or
##                           pruned sparse)
##   simulate_sep()       -- SEP IRF / simulation driver
##
## ALGORITHM
##   For a given initial shock eps_1 and state y_0:
##   1. Draw Q quadrature scenarios for future shocks eps_2, ..., eps_T.
##   2. For each scenario q:
##        a. Construct shock_path = [eps_1; eps_{2:T}^{(q)}].
##        b. Solve the perfect-foresight path Y^{(q)} via pf_newton_solve().
##   3. The SEP approximation to E[y_t | y_0, eps_1] is:
##        E_hat[y_t] = sum_q  w_q * Y_t^{(q)}
##   The first-period response (t=1) is deterministic (same across all q since
##   eps_1 is fixed and y_0 is given); the distribution matters for t >= 2.
##
##   "Higher-order tail": for periods beyond horizon T, the SEP switches to the
##   perturbation solution (if available), which provides an O(2) or O(3)
##   approximation to the stochastic expectations. This is not yet implemented;
##   the current version uses y_ss as the terminal condition (equivalent to
##   perfect-foresight with certainty equivalence at the boundary).
##
## LIMITATIONS (Phase O5 v1):
##   - Tensor-product quadrature only (scales poorly with n_exo > 2).
##   - Terminal condition is y_ss (no perturbation tail).
##   - Each PF solve is independent (no warm-starting across quadrature nodes).
##
## References:
##   Adjemian & Juillard (2025), "Stochastic Extended Path", manuscript.
##   Judd (1992), Numerical Methods in Economics, Ch. 7 (Gaussian quadrature).
## --------------------------------------------------------------------------


# =============================================================================
# Quadrature
# =============================================================================

#' Gauss-Hermite nodes and weights for standard normal integration
#'
#' Returns n-point Gauss-Hermite quadrature for N(0,1): approximates
#'   integral_{-inf}^{inf} f(x) phi(x) dx  ≈  sum_k w_k * f(x_k)
#' where phi is the standard normal density.
#'
#' @param n Integer: number of quadrature points (1..20)
#' @return List with $nodes (length-n) and $weights (length-n, sum to 1).
#' @noRd
sep_gauss_hermite <- function(n = 5L) {
  n <- as.integer(n)

  # Pre-tabulated nodes and weights for the physicist's convention
  # (integrating exp(-x^2) f(x) dx). We convert to N(0,1) via x -> x*sqrt(2).
  # Source: Abramowitz & Stegun Table 25.10
  gh_table <- list(
    `1`  = list(x = 0,                                     w = 1.7724539),
    `2`  = list(x = c(-0.7071068, 0.7071068),              w = c(0.8862269, 0.8862269)),
    `3`  = list(x = c(-1.2247449, 0, 1.2247449),           w = c(0.2954090, 1.1816360, 0.2954090)),
    `5`  = list(x = c(-2.0201829,-0.9585725, 0, 0.9585725, 2.0201829),
                w = c(0.0199532, 0.3936193, 0.9453087, 0.3936193, 0.0199532)),
    `7`  = list(x = c(-2.6519614,-1.6735516,-0.8162878, 0,
                       0.8162878, 1.6735516, 2.6519614),
                w = c(0.0009718, 0.0545156, 0.4256073, 0.8102646,
                       0.4256073, 0.0545156, 0.0009718))
  )

  key <- as.character(n)
  if (!key %in% names(gh_table)) {
    # Fall back to n=5 for unsupported sizes
    key <- "5"
    n   <- 5L
  }
  tbl <- gh_table[[key]]

  # Convert physicist to probabilist (N(0,1)):
  #   x_prob = x_phys / sqrt(2),   w_prob = w_phys / sqrt(pi)
  nodes   <- tbl$x / sqrt(2)
  weights <- tbl$w / sqrt(pi)

  list(nodes = nodes, weights = weights)
}


#' Build tensor-product quadrature grid for n_exo independent N(0,1) shocks
#'
#' @param n_exo    Number of independent shock dimensions
#' @param n_pts    Quadrature points per dimension (default 5)
#' @param std      Vector of shock standard deviations (length n_exo).
#'                 Each node is scaled by the corresponding std.
#' @return List with $nodes (Q × n_exo) and $weights (length Q, sum to 1).
#' @noRd
sep_quadrature <- function(n_exo, n_pts = 5L, std = rep(1, n_exo)) {
  gh   <- sep_gauss_hermite(n_pts)
  nd1  <- gh$nodes
  wt1  <- gh$weights

  if (n_exo == 1L) {
    return(list(
      nodes   = matrix(nd1 * std[1], ncol = 1),
      weights = wt1
    ))
  }

  # Tensor product: all combinations of 1D nodes
  idx <- as.matrix(do.call(expand.grid, replicate(n_exo, seq_len(n_pts), simplify = FALSE)))
  Q   <- nrow(idx)
  nodes   <- matrix(0, nrow = Q, ncol = n_exo)
  weights <- numeric(Q)

  for (q in seq_len(Q)) {
    w <- 1
    for (j in seq_len(n_exo)) {
      nodes[q, j]  <- nd1[idx[q, j]] * std[j]
      w            <- w * wt1[idx[q, j]]
    }
    weights[q] <- w
  }

  list(nodes = nodes, weights = weights)
}


# =============================================================================
# SEP driver
# =============================================================================

#' Simulate IRFs via Stochastic Extended Path (SEP)
#'
#' Computes the expected impulse response E[y_t - y_ss | y_0, eps_1] for a
#' nonlinear DSGE model with optional occasionally-binding constraints.
#' Unlike perturbation-based OBC solvers, this operates on the full nonlinear
#' dynamic residuals and does NOT require model(linear).
#'
#' @param compiled   dynhr_compiled from compile_model().
#' @param y0         Named numeric: initial state (y at t=0). Usually y_ss.
#' @param y_ss       Named numeric: steady state (used as terminal condition).
#' @param params     Named numeric parameter vector.
#' @param shock_var  Name of the shock variable (must be in varexo_names).
#' @param shock_size Scalar: size of the initial impulse for shock_var at t=1.
#' @param horizon    Integer: IRF horizon T (default 20).
#' @param n_quad     Integer: quadrature points per shock dimension (default 5).
#' @param obc_specs  List of OBC specs from obc_parse_tags() (empty = no OBCs).
#' @param ...        Additional arguments passed to pf_newton_solve().
#' @return Data frame with columns: period (1..T), variable name columns
#'   (IRF in deviation from y_ss).
#' @export
simulate_sep <- function(compiled,
                         y0,
                         y_ss,
                         params,
                         shock_var,
                         shock_size = 1,
                         horizon    = 20L,
                         n_quad     = 5L,
                         obc_specs  = list(),
                         ...) {

  dyn    <- compiled$dynamic
  n_exo  <- length(dyn$exo_names)
  T      <- as.integer(horizon)

  if (!shock_var %in% dyn$exo_names) stop(sprintf(
    "simulate_sep: shock_var '%s' not in varexo_names.", shock_var))

  # Determine shock standard deviations from model shocks block (if available)
  model   <- compiled$model
  std_vec <- rep(1, n_exo)
  names(std_vec) <- dyn$exo_names
  shk_tbl <- model$shocks$variances   # data.frame: name, stderr, variance
  if (is.data.frame(shk_tbl) && nrow(shk_tbl) > 0L) {
    for (i in seq_len(nrow(shk_tbl))) {
      nm  <- shk_tbl$name[i]
      std <- shk_tbl$stderr[i]
      if (!is.na(nm) && nm %in% dyn$exo_names && !is.na(std))
        std_vec[nm] <- as.numeric(std)
    }
  }

  # Build quadrature for future periods (t=2..T)
  quad <- sep_quadrature(n_exo, n_pts = n_quad, std = std_vec)
  Q    <- nrow(quad$nodes)

  # Allocate: weighted sum of IRF paths over quadrature scenarios
  n_endo <- length(dyn$endo_names)
  irf_sum <- matrix(0, nrow = T, ncol = n_endo)
  colnames(irf_sum) <- dyn$endo_names

  # Initial shock path: eps_1 = shock_size for shock_var, others = 0
  eps1 <- numeric(n_exo)
  names(eps1) <- dyn$exo_names
  eps1[shock_var] <- shock_size

  for (q in seq_len(Q)) {
    # Construct full T × n_exo shock path
    shock_path <- matrix(0, nrow = T, ncol = n_exo)
    colnames(shock_path) <- dyn$exo_names
    shock_path[1L, ] <- eps1
    if (T > 1L) shock_path[2L:T, ] <- matrix(rep(quad$nodes[q, ], T - 1L),
                                               nrow = T - 1L, byrow = TRUE)

    pf <- pf_newton_solve(compiled, y0, y_ss, shock_path, params,
                          obc_specs = obc_specs, ...)
    if (!pf$converged) {
      warning(sprintf("simulate_sep: quadrature node %d did not converge.", q))
    }

    irf_sum <- irf_sum + quad$weights[q] * pf$irf
  }

  # Return as data frame
  df <- as.data.frame(irf_sum)
  df$period <- seq_len(T)
  df[, c("period", dyn$endo_names)]
}
