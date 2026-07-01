## R/dynamic-perturbation-path.R
## --------------------------------------------------------------------------
## Dynamic perturbation around a time-varying path (Mennuni et al. 2025).
##
## Implements the single-regime dynamic perturbation (DP) algorithm
## (Figures 2-4 of the paper) integrated with dynhr's compiled model kernel.
## Supports both the deterministic version (quadrature=NULL) and the stochastic
## version (quadrature != NULL) where the equilibrium expectation is integrated
## over Gauss-Hermite or monomial quadrature nodes of the shock distribution.
##
## Algorithm (deterministic, §2 of scope-C7 brief):
##
##  For each period t = 1 ... T-1:
##    1. Forward auxiliary path (n_aux steps from x_sim[,t] via SS rule).
##    2. Initialize local policy with SS rule (gg_ss, hh_ss from dr$ghx).
##    3. Backward sweep from xaux[:,n] to xaux[:,1]:
##       a. Solve nonlinear system for [xp; y] using nleqslv.
##       b. Build dy adapter vector for jacobian_fn.
##       c. Compute AD Jacobian, split into Eaa/Ebb blocks.
##       d. IFT update: tmp = -solve(Eaa, Ebb); update local hh, gg.
##    4. Advance: x_sim[,t+1] = hh$f0 + ETA %*% shocks[,t+1]
##               y_sim[,t]   = gg$f0
##
## Stochastic extension (quadrature != NULL):
##   Step 3a: the objective is E_t[f(x_t, y_t, x_{t+1}, y_{t+1})] = 0, where
##   x_{t+1} = xp_cand + ETA %*% u, u ~ N(0, Sigma_e).  The expectation is
##   approximated by a weighted sum over quadrature nodes u_k with weights w_k:
##     obj = sum_k w_k * residuals_fn(dy(xp_cand + ETA*u_k, y, xt, yp(xp+ETA*u_k)))
##   The IFT Jacobian at the solution is similarly computed at the quadrature-
##   weighted mean (the deterministic node xp_sol, zero shock), which gives the
##   correct leading-order risk correction to O(sigma^2).
##
## References:
##   Mennuni, Paczos, Walker (2025) "Perturbation around a Time-Varying Path".
##   Standalone R translation in replication/DynamicPerturbation_ReplicationFiles/R/.
## --------------------------------------------------------------------------

#' Dynamic perturbation around a time-varying auxiliary path
#'
#' Implements the single-regime deterministic dynamic perturbation (DP)
#' algorithm from Mennuni, Paczos, Walker (2025) using dynhr's compiled model
#' kernel.  For each simulation period an auxiliary path of length \code{n_aux}
#' is projected forward from the current state via the steady-state linear
#' decision rule; a backward sweep then solves the nonlinear equilibrium at
#' each step and updates a local linear policy via the implicit function theorem
#' (IFT).  The terminal condition is the SS decision rule.
#'
#' @param compiled   A \code{dynhr_compiled} object from \code{compile_model()}.
#' @param ss         Named numeric steady-state vector (from \code{solve_steady()$values}).
#' @param params     Named numeric parameter vector.
#' @param dr         Decision rule list from \code{solve_perturbation(order=1)}.
#'   Must contain \code{dr$ghx} (all-variable linear rule, rows in declaration
#'   order) and \code{dr$state_vars} (names of predetermined state variables).
#' @param n_aux      Length of auxiliary path (integer >= 2; default 15).
#'   Longer auxiliary paths capture more of the nonlinear transition dynamics
#'   but increase computational cost.
#' @param T          Number of simulation periods (integer >= 2; default 2).
#' @param x_0        Initial state vector (length = number of state variables).
#'   Defaults to \code{ss[state_vars]} (start at steady state).
#' @param ETA        Shock loading matrix (\code{n_state x n_shocks}).
#'   Rows correspond to state variables in order.  Pass \code{NULL} or a zero
#'   matrix for a deterministic path.
#' @param shocks     Shock matrix (\code{n_shocks x T}).  Column \code{t}
#'   contains the shocks applied at period \code{t}.  Pass \code{NULL} for
#'   a zero-shock (deterministic) path.
#' @param tol        Convergence tolerance for the nonlinear solver (default 1e-10).
#' @param verbose    Logical; if \code{TRUE}, print progress.
#' @param quadrature Quadrature specification for the stochastic version, or
#'   \code{NULL} (default) for the deterministic path.  Accepted forms:
#'   \describe{
#'     \item{\code{NULL}}{Deterministic path (shocks = 0 in the expectation).
#'       Behaviour is byte-identical to the original single-node version.}
#'     \item{integer scalar \code{n_gh}}{Gauss-Hermite product rule with
#'       \code{n_gh} points per shock dimension (\code{n_gh^n_shock} total nodes).
#'       Use \code{n_gh = 5} for 4th-order accuracy (integrates degree-9 polynomials).}
#'     \item{\code{"monomial"}}{Stroud monomial-2 rule: 2*n_shock nodes, exact for
#'       total degree <= 3.  Captures the O(sigma^2) risk-correction term with
#'       minimal node count.  Preferred for n_shock >= 4.}
#'     \item{list(\code{nodes}, \code{weights})}{Pre-built rule: \code{nodes} is
#'       an (n_nodes x n_shock) matrix of shock vectors and \code{weights} is a
#'       length-n_nodes vector summing to 1.}
#'   }
#'   Requires \code{Sigma_e} when \code{quadrature} is not a pre-built rule.
#'   ACCURACY NOTE: the stochastic version integrates the FIRST-ORDER local
#'   policy over the shock nodes, so it captures the model's f-nonlinearity
#'   \eqn{O(\sigma^2)} risk term but NOT the second-order policy-curvature
#'   (\code{ghxx}/\code{ghuu}) part of the full risk correction.  Its
#'   steady-state risk adjustment has the correct sign but is only a partial
#'   fraction of the order-2 \code{ghss} correction (~1/3 in tests); full
#'   second-order risk accuracy requires a second-order local policy in the
#'   backward sweep (a future increment).
#' @param Sigma_e    Shock covariance matrix (n_shock x n_shock).  Required
#'   when \code{quadrature} is non-NULL.  Typical usage:
#'   \code{Sigma_e = diag(dr$se_estim^2)} or \code{Sigma_e = dr$Sigma_e}.
#' @return A list with:
#'   \describe{
#'     \item{\code{y_sim}}{Matrix (\code{n_jump x T}) of jump variable values.
#'       Rows correspond to \code{jump_vars} (non-state endogenous variables in
#'       declaration order).  Column \code{t} is the period-\code{t} jump.}
#'     \item{\code{x_sim}}{Matrix (\code{n_state x T}) of state variable values.
#'       Rows correspond to \code{dr$state_vars}.  Column 1 is \code{x_0};
#'       column \code{t+1} is the state entering period \code{t+1}.}
#'     \item{\code{state_vars}}{Character vector: names of state variables (rows of x_sim).}
#'     \item{\code{jump_vars}}{Character vector: names of jump variables (rows of y_sim).}
#'   }
#' @export
dynamic_path_perturbation <- function(compiled, ss, params, dr,
                                       n_aux = 15L, T = 2L,
                                       x_0 = NULL, ETA = NULL, shocks = NULL,
                                       tol = 1e-10, verbose = FALSE,
                                       quadrature = NULL, Sigma_e = NULL) {

  # ---- Input validation --------------------------------------------------
  if (!requireNamespace("nleqslv", quietly = TRUE)) {
    stop("dynamic_path_perturbation requires the 'nleqslv' package. ",
         "Install it with: install.packages('nleqslv')")
  }

  n_aux <- as.integer(n_aux)
  T     <- as.integer(T)
  if (n_aux < 2L) stop("n_aux must be >= 2")
  if (T < 2L) stop("T must be >= 2")

  # ---- Quadrature setup -------------------------------------------------
  # Resolve the quadrature specification to a list(nodes, weights) or NULL.
  # When NULL (default), the deterministic path is used (byte-identical to the
  # original implementation).
  q_spec <- NULL
  if (!is.null(quadrature)) {
    if (is.null(Sigma_e)) {
      stop("'Sigma_e' (shock covariance matrix) must be supplied when ",
           "'quadrature' is non-NULL.")
    }
    Sigma_e <- as.matrix(Sigma_e)
    q_spec  <- .dp_resolve_quadrature(quadrature, Sigma_e)
  }

  # Extract model structure
  dyn       <- compiled$dynamic
  dyn_col_map <- dyn$dyn_col_map
  endo      <- compiled$model$var_names
  state_vars <- dr$state_vars

  # Non-state (jump) variables in declaration order
  jump_vars <- setdiff(endo, state_vars)

  nx <- length(state_vars)
  ny <- length(jump_vars)

  # Steady-state subvectors
  x_ss <- ss[state_vars]
  y_ss <- ss[jump_vars]

  # ---- Steady-state linear rule (hh_ss, gg_ss) ---------------------------
  # dr$ghx is (n_endo x n_state): all variable responses to state deviation.
  # hh_ss: state-variable rows (nx x nx)
  # gg_ss: jump-variable rows  (ny x nx)
  ghx <- dr$ghx
  hh_ss <- ghx[state_vars, , drop = FALSE]  # nx x nx
  gg_ss <- ghx[jump_vars,  , drop = FALSE]  # ny x nx

  # ---- Initial state -----------------------------------------------------
  if (is.null(x_0)) {
    x_0 <- x_ss
  } else {
    if (length(x_0) != nx) {
      stop(sprintf("x_0 must have length %d (number of state variables: %s)",
                   nx, paste(state_vars, collapse = ", ")))
    }
    if (is.null(names(x_0))) names(x_0) <- state_vars
  }

  # ---- Shock matrices ----------------------------------------------------
  n_shk <- if (!is.null(ETA)) ncol(ETA) else 0L
  if (is.null(ETA))    ETA    <- matrix(0, nx, 0L)
  if (is.null(shocks)) shocks <- matrix(0, n_shk, T)
  if (n_shk == 0L && nrow(shocks) == 0L) shocks <- matrix(0, 0L, T)

  # ---- Simulation storage ------------------------------------------------
  x_sim <- matrix(0, nx, T)
  y_sim <- matrix(0, ny, T)
  rownames(x_sim) <- state_vars
  rownames(y_sim) <- jump_vars
  x_sim[, 1] <- x_0

  # ---- Auxiliary path storage --------------------------------------------
  xaux <- matrix(0, nx, n_aux + 1L)

  # ---- Residuals closure using dynhr kernel ------------------------------
  # At each backward step: objective = residuals_fn(dy, params, ss)
  # where dy is built from the candidate [xp; y] and the local ghat policy.
  residuals_fn <- dyn$residuals_fn
  jacobian_fn  <- dyn$jacobian_fn

  # ---- Main simulation loop ----------------------------------------------
  t <- 1L
  while (t < T) {

    # 1. Forward auxiliary path (n_aux steps from current state via SS rule)
    xaux[, 1] <- x_sim[, t]
    for (tf in seq_len(n_aux)) {
      dx <- xaux[, tf] - x_ss
      xaux[, tf + 1L] <- x_ss + hh_ss %*% dx
    }

    # 2. Initialize local policy at SS
    gg <- list(df = gg_ss, x0 = x_ss, f0 = y_ss)
    hh <- list(df = hh_ss, x0 = x_ss, f0 = x_ss)

    # 3. Backward sweep from xaux[:,n_aux] back to xaux[:,1]
    tb <- n_aux
    while (tb > 1L) {
      xt <- xaux[, tb - 1L]   # current state (names from rownames(xaux))
      names(xt) <- state_vars

      # Initial guess from current local policy
      xp_guess <- hh$df %*% (xt - hh$x0) + hh$f0
      y_guess  <- gg$df %*% (xt - gg$x0) + gg$f0

      # Apply ghat to get next-period jump approximation
      # (ghat uses the updated local policy gg for yp = gg$df*(xp - gg$x0) + gg$f0)
      # This is evaluated INSIDE the objective with the current gg — so the
      # objective depends on gg, which is fixed for this backward step.

      # Objective function: residuals_fn at candidate z = [xp; y]
      # with yp determined by ghat policy applied to xp.
      #
      # Deterministic (q_spec NULL): evaluate at single node (xp, y, yp=ghat(xp)).
      # Stochastic (q_spec non-NULL): evaluate E_t[f] over quadrature nodes.
      #   x_{t+1} = xp_cand + ETA %*% u_k for shock node u_k.
      #   yp_k = ghat(xp_cand + ETA %*% u_k).
      #   obj = sum_k w_k * f(xt, y_cand, xp_cand + ETA*u_k, yp_k).
      obj_fn <- if (is.null(q_spec)) {
        function(z) {
          xp_cand <- z[seq_len(nx)]
          y_cand  <- z[nx + seq_len(ny)]
          yp_cand <- gg$df %*% (xp_cand - gg$x0) + gg$f0
          dy <- .dp_build_dy(
            xp  = xp_cand,
            y   = y_cand,
            xt  = xt,
            yp  = as.numeric(yp_cand),
            state_vars  = state_vars,
            jump_vars   = jump_vars,
            dyn_col_map = dyn_col_map,
            ss          = ss
          )
          as.numeric(residuals_fn(dy, params, ss))
        }
      } else {
        # Stochastic: integrate residuals over quadrature nodes.
        # q_spec$nodes: (n_nodes x n_shock) shock vectors u_k in R^{n_shock}.
        # ETA: (nx x n_shock) shock-loading matrix mapping shocks to states.
        # xp_{t+1,k} = xp_cand + ETA %*% u_k.
        q_nodes   <- q_spec$nodes    # n_nodes x n_shock
        q_weights <- q_spec$weights  # n_nodes
        n_nodes   <- q_spec$n_nodes

        # Verify ETA dimension matches quadrature shock dimension.
        n_shock_q <- q_spec$n_shock
        if (ncol(ETA) != n_shock_q) {
          stop(sprintf(
            "Dimension mismatch: ETA has %d columns (shocks) but Sigma_e has %d rows. ",
            ncol(ETA), n_shock_q),
            "Ensure ETA and Sigma_e refer to the same shock vector.")
        }

        # ETA_q: precomputed shifted nodes in state space (nx x n_nodes)
        ETA_nodes <- ETA %*% t(q_nodes)   # nx x n_nodes (precomputed outside loop)

        function(z) {
          xp_cand <- z[seq_len(nx)]
          y_cand  <- z[nx + seq_len(ny)]
          # Accumulate weighted residuals over quadrature nodes
          n_eq    <- nx + ny   # number of equilibrium equations
          res_sum <- numeric(n_eq)
          for (k in seq_len(n_nodes)) {
            xp_k  <- xp_cand + ETA_nodes[, k]
            yp_k  <- as.numeric(gg$df %*% (xp_k - gg$x0) + gg$f0)
            dy_k  <- .dp_build_dy(
              xp  = xp_k,
              y   = y_cand,
              xt  = xt,
              yp  = yp_k,
              state_vars  = state_vars,
              jump_vars   = jump_vars,
              dyn_col_map = dyn_col_map,
              ss          = ss
            )
            res_sum <- res_sum + q_weights[k] * as.numeric(residuals_fn(dy_k, params, ss))
          }
          res_sum
        }
      }

      # Solve nonlinear system
      z_guess <- c(xp_guess, y_guess)
      sol <- tryCatch(
        nleqslv::nleqslv(z_guess, obj_fn,
                          control = list(xtol = tol * 1e-2, ftol = tol * 1e-2,
                                         maxit = 10000L, allowSingular = TRUE)),
        error = function(e) NULL
      )

      if (!is.null(sol) && max(abs(sol$fvec)) < tol && all(is.finite(sol$fvec))) {
        z_sol  <- sol$x
        xp_sol <- z_sol[seq_len(nx)]
        y_sol  <- z_sol[nx + seq_len(ny)]

        # 3b. Compute AD Jacobian at solution
        yp_sol <- as.numeric(gg$df %*% (xp_sol - gg$x0) + gg$f0)
        dy_sol <- .dp_build_dy(
          xp  = xp_sol,
          y   = y_sol,
          xt  = xt,
          yp  = yp_sol,
          state_vars  = state_vars,
          jump_vars   = jump_vars,
          dyn_col_map = dyn_col_map,
          ss          = ss
        )
        J_sol <- jacobian_fn(dy_sol, params, ss)

        # 3c. Split into Eaa and Ebb.  gg$df is passed so the xp block of Eaa
        # picks up the chain-rule term J[, jump@+1] %*% gg$df from yp = ghat(xp).
        blocks <- .dp_split_jacobian(J_sol, state_vars, jump_vars, dyn_col_map,
                                     gg_df = gg$df)
        Eaa <- blocks$Eaa
        Ebb <- blocks$Ebb

        # 3d. IFT update: tmp = -solve(Eaa, Ebb)
        # Eaa is (n_eq x (nx+ny)), Ebb is (n_eq x nx)
        # tmp is (nx+ny x nx): [new_hh_df; new_gg_df]
        tmp <- tryCatch(
          -solve(Eaa, Ebb),
          error = function(e) NULL
        )

        if (!is.null(tmp) && all(is.finite(tmp))) {
          hh$df <- tmp[seq_len(nx), , drop = FALSE]
          hh$x0 <- xt
          hh$f0 <- matrix(xp_sol, nrow = nx)
          gg$df <- tmp[nx + seq_len(ny), , drop = FALSE]
          gg$x0 <- xt
          gg$f0 <- matrix(y_sol, nrow = ny)
        } else if (verbose) {
          cat(sprintf("  IFT solve failed at tb=%d, t=%d\n", tb, t))
        }

      } else if (verbose) {
        cat(sprintf("  nleqslv failed at tb=%d, t=%d (max|res|=%g)\n",
                    tb, t,
                    if (!is.null(sol)) max(abs(sol$fvec)) else Inf))
      }

      tb <- tb - 1L
    }

    # 4. Advance one period
    t <- t + 1L
    shock_t <- if (ncol(shocks) >= t) shocks[, t, drop = FALSE] else matrix(0, n_shk, 1L)
    x_sim[, t] <- as.numeric(hh$f0) + ETA %*% shock_t
    y_sim[, t - 1L] <- as.numeric(gg$f0)

    if (verbose) cat(sprintf("dynamic_path_perturbation: period %d done\n", t - 1L))
  }

  list(
    y_sim      = y_sim,
    x_sim      = x_sim,
    state_vars = state_vars,
    jump_vars  = jump_vars
  )
}
