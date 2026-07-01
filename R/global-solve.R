## R/global-solve.R
## --------------------------------------------------------------------------
## Global (projection) solution for DSGE models via Chebyshev polynomial
## collocation and Coleman time iteration.
##
## Algorithm: SIMULTANEOUS-NEWTON approach (brief §3.6 "Alternative").
## At each collocation node (s_lag = state at t-1), we solve for ALL
## current-period endogenous variables y_t = (c_t, k_t, z_t, ...) by
## finding the root of the expected model residuals, where:
##   - y_{t-1} is the collocation node (known lag)
##   - y_t is the unknown (all n_endo vars at current period)
##   - y_{t+1} = policy(state_t) evaluated at the next-period states,
##     averaged over Gauss-Hermite quadrature nodes for eps_{t+1}
##   - eps_t = 0 (expectation is over t+1 shocks; current shock is 0)
##
## dy-assembly rule (from dyn_col_map):
##   var__m1  <- state_lag[var]         (lag vars = collocation grid point)
##   var__0   <- y_cur[var]             (all endo at t = unknowns)
##   var__p1  <- y_lead[var]            (all endo at t+1 from policy + quad)
##   shock__0 <- 0 (or relevant shock)  (shocks enter as exo vars)
##
## The residual at each collocation node is:
##   R(y_cur) = SUM_k w_k * F(y_lag, y_cur, y_lead_k, eps_k)
## where F is the model residual from residuals_fn.
## --------------------------------------------------------------------------

#' Global projection solution for a compiled DSGE model
#'
#' Solves a DSGE model via Chebyshev polynomial collocation with
#' Coleman time iteration and Gauss-Hermite quadrature for expectations.
#' Uses the simultaneous-Newton approach: at each collocation node, all
#' current-period endogenous variables are found jointly via \code{nleqslv}.
#'
#' @param compiled  dynhr_compiled from \code{compile_model()}.
#' @param ss        Steady-state result from \code{solve_steady()} or named
#'   numeric of steady-state values.
#' @param params    Named numeric parameter vector.
#' @param poly_degree Integer; Chebyshev polynomial degree (default 3).
#' @param n_quad    Integer; Gauss-Hermite quadrature nodes per shock dim
#'   (default 5).
#' @param n_nodes   Integer; collocation nodes per state dimension (default 7).
#' @param state_domain Named list; per-state \code{c(lo, hi)} bounds.
#'   If NULL, auto-computed from SS plus coverage based on shock variances.
#' @param tol       Convergence tolerance for max coefficient update (default
#'   1e-7).
#' @param max_iter  Maximum outer iterations (default 500).
#' @param init_from_perturbation Logical; warm-start coefficients from the
#'   order-1 perturbation solution (default TRUE).
#' @param verbose   Logical; print iteration diagnostics (default FALSE).
#' @return A \code{GlobalSolution} S3 object with fields \code{coefs},
#'   \code{state_names}, \code{all_endo_names}, \code{shock_names},
#'   \code{state_domain}, \code{poly_degree}, \code{converged}, etc.
#' @export
solve_global <- function(compiled,
                         ss,
                         params,
                         poly_degree  = 3L,
                         n_quad       = 5L,
                         n_nodes      = 7L,
                         state_domain = NULL,
                         tol          = 1e-7,
                         max_iter     = 500L,
                         init_from_perturbation = TRUE,
                         verbose      = FALSE) {

  ## ------------------------------------------------------------------
  ## 0. Extract model components
  ## ------------------------------------------------------------------
  model   <- compiled$model
  dyn     <- compiled$dynamic
  endo    <- dyn$endo_names      # all endogenous variable names
  exo     <- dyn$exo_names       # shock names
  n_endo  <- length(endo)
  n_exo   <- length(exo)
  n_eq    <- dyn$n_eq
  col_map <- dyn$dyn_col_map     # data.frame: name, lead_lag, col
  tot_cols <- dyn$total_cols

  ss_vals <- if (inherits(ss, "dynhr_steady")) ss$values else ss

  ## Build dy key lookup from col_map
  ## dy_keys[col] = "varname__timing" for each column position
  dy_keys <- character(tot_cols)
  for (k in seq_len(nrow(col_map))) {
    nm  <- col_map$name[k]
    ll  <- col_map$lead_lag[k]
    sfx <- if (ll == 0L) "__0" else if (ll > 0L) paste0("__p", ll) else paste0("__m", abs(ll))
    dy_keys[col_map$col[k]] <- paste0(nm, sfx)
  }

  ## ------------------------------------------------------------------
  ## 1. Classify state (predetermined) vs control (jump) variables
  ## ------------------------------------------------------------------
  lli        <- compiled$lead_lag_incidence %||% model$lead_lag_incidence
  ## LLI rownames are like "t-1", "t", "t+1" — parse to integers
  row_labels <- vapply(rownames(lli), function(rn) {
    rn <- trimws(rn)
    if (rn == "t") return(0L)
    m <- regmatches(rn, regexec("^t([+-]?\\d+)$", rn))[[1]]
    if (length(m) == 2L) return(as.integer(m[2L]))
    0L
  }, integer(1L), USE.NAMES = FALSE)
  sl         <- .structural_lag_lead(lli, row_labels)
  state_names  <- endo[sl$has_lag]    # vars that appear at t-1
  n_state      <- length(state_names)

  if (n_state == 0L)
    stop("solve_global: model has no state variables; projection not applicable.")

  ## ------------------------------------------------------------------
  ## 2. Build index maps for fast dy assembly
  ## ------------------------------------------------------------------
  ## For each (var, timing) we need the column index in the dy vector
  col_of <- function(nm, ll) {
    row <- col_map[col_map$name == nm & col_map$lead_lag == ll, , drop = FALSE]
    if (nrow(row) == 0L) return(NA_integer_)
    row$col[1L]
  }

  ## Pre-build index vectors for all endogenous vars at each timing
  idx_m1 <- vapply(endo, function(nm) col_of(nm, -1L), integer(1))  # lag
  idx_0  <- vapply(endo, function(nm) col_of(nm,  0L), integer(1))  # current
  idx_p1 <- vapply(endo, function(nm) col_of(nm, +1L), integer(1))  # lead
  names(idx_m1) <- endo; names(idx_0) <- endo; names(idx_p1) <- endo

  ## Shock index in dy
  idx_exo <- vapply(exo, function(nm) col_of(nm, 0L), integer(1))
  names(idx_exo) <- exo

  ## ------------------------------------------------------------------
  ## 3. Shock standard deviations
  ## ------------------------------------------------------------------
  shock_sds <- .get_shock_sds(model, params)

  ## ------------------------------------------------------------------
  ## 4. State domain
  ## ------------------------------------------------------------------
  if (is.null(state_domain)) {
    state_domain <- .auto_state_domain(
      state_names = state_names,
      ss_vals     = ss_vals,
      shock_sds   = shock_sds,
      model       = model
    )
  } else {
    for (nm in state_names)
      if (is.null(state_domain[[nm]]))
        stop(sprintf("solve_global: state_domain missing entry for '%s'", nm))
  }

  ## ------------------------------------------------------------------
  ## 5. Gauss-Hermite quadrature
  ## ------------------------------------------------------------------
  gh        <- gauss_hermite(n_quad)
  q_nodes   <- gh$nodes
  q_weights <- gh$weights

  ## ------------------------------------------------------------------
  ## 6. Chebyshev collocation grid
  ## ------------------------------------------------------------------
  ## Tensor-product Chebyshev nodes in each state dimension
  nodes_1d <- lapply(state_names, function(nm) {
    dom <- state_domain[[nm]]
    u   <- cheb_nodes(n_nodes)
    cheb_denormalize(u, dom[1], dom[2])
  })

  ## Grid (n_nodes^n_state rows x n_state cols) — state lag values
  grid_args       <- rev(nodes_1d)
  names(grid_args) <- rev(state_names)
  grid_nat        <- as.matrix(expand.grid(grid_args))
  grid_nat        <- grid_nat[, state_names, drop = FALSE]
  n_coll          <- nrow(grid_nat)

  ## Normalized grid for basis evaluation
  grid_norm <- matrix(0.0, nrow = n_coll, ncol = n_state)
  colnames(grid_norm) <- state_names
  for (j in seq_len(n_state)) {
    nm  <- state_names[j]
    dom <- state_domain[[nm]]
    grid_norm[, j] <- cheb_normalize(grid_nat[, nm], dom[1], dom[2])
  }

  ## ------------------------------------------------------------------
  ## 7. Chebyshev basis at collocation nodes
  ## ------------------------------------------------------------------
  Phi     <- cheb_basis(grid_norm, poly_degree)   # n_coll x n_basis
  n_basis <- ncol(Phi)

  if (verbose)
    cat(sprintf(
      "solve_global: %d states, %d endo, %d coll.nodes, %d basis fns, %d quad pts\n",
      n_state, n_endo, n_coll, n_basis, n_quad))

  ## ------------------------------------------------------------------
  ## 8. Initialize polynomial coefficients (n_endo x n_basis)
  ## ------------------------------------------------------------------
  ## coefs[i, ] gives the Chebyshev coefficients for endogenous var i
  ## as a function of the state LAG values.
  coefs <- matrix(0.0, nrow = n_endo, ncol = n_basis)
  rownames(coefs) <- endo

  if (init_from_perturbation) {
    dr1 <- tryCatch(
      solve_perturbation(model, compiled, ss_vals, params, verbose = FALSE),
      error = function(e) {
        if (verbose) cat("solve_global: perturbation warm-start failed:",
                         conditionMessage(e), "\n")
        NULL
      })
    if (!is.null(dr1)) {
      coefs <- .init_coefs_from_dr1(dr1, endo, state_names, ss_vals,
                                    grid_nat, state_domain, Phi, n_coll, n_basis)
      if (verbose) cat("solve_global: warm-started from order-1 perturbation\n")
    } else {
      ## Constant at SS
      for (i in seq_len(n_endo)) coefs[i, 1L] <- ss_vals[endo[i]]
    }
  } else {
    for (i in seq_len(n_endo)) coefs[i, 1L] <- ss_vals[endo[i]]
  }

  ## ------------------------------------------------------------------
  ## 9. Residuals function reference
  ## ------------------------------------------------------------------
  residuals_fn <- dyn$residuals_fn

  ## ------------------------------------------------------------------
  ## 10. Evaluate policy: (n_pts x n_state) normalized -> (n_pts x n_endo)
  ## ------------------------------------------------------------------
  eval_policy <- function(xnorm_mat) {
    B <- cheb_basis(xnorm_mat, poly_degree)
    tcrossprod(B, coefs)    # n_pts x n_endo
  }

  ## Normalize a state lag vector to [-1,1]^d (clamped)
  norm_state_lag <- function(s) {
    u <- numeric(n_state)
    for (j in seq_len(n_state)) {
      nm  <- state_names[j]
      dom <- state_domain[[nm]]
      u[j] <- max(-1.0, min(1.0, cheb_normalize(s[j], dom[1], dom[2])))
    }
    u
  }

  ## ------------------------------------------------------------------
  ## 11. Time-iteration loop
  ## ------------------------------------------------------------------

  ## Next-period state lag from current endogenous vector.
  ## In Dynare timing, the next-period LAG of states = current state values:
  ##   k(-1) at t+1 = k at t    (i.e., capital chosen at t)
  ##   z(-1) at t+1 = z at t    (since z = rho*z(-1)+eps, z at t IS the state)
  ## So: next_lag[j] = y_cur[state_names[j]]
  next_lag_from_ycur <- function(y_cur) {
    y_cur[state_names]
  }

  ## For the AR(1) shocks (z), next-period z given current z and next eps:
  ## But wait — in the simultaneous Newton approach, y_cur includes z_t,
  ## and the LEAD value z(+1) at time t is z_{t+1}.
  ## z_{t+1} = rho * z_t + eps_{t+1}
  ## Since z_t is in y_cur, we can compute z_{t+1} from the quadrature draw.
  ## But k_{t+1} is given by the policy function evaluated at (k_t, z_{t+1}).
  ##
  ## ACTUALLY: y_cur gives us the current state values AT TIME t.
  ## The LEAD (y_{t+1}) is computed as:
  ##   - first compute next_lag = y_cur[state_names] (the lag for t+1)
  ##   - then for AR(1) shocks: z_{t+1} = rho*z_t + eps_{t+1,k}
  ##     so the z-component of next_lag needs to be UPDATED with the shock
  ##   - the full next_lag for the policy evaluation at t+1 is:
  ##     (k_t, z_{t+1}) = (y_cur["k"], rho*y_cur["z"] + eps_k)
  ## But the policy maps lag -> current, so policy(k_t, z_{t+1}) gives y_{t+1}.
  ## The issue is that z_{t+1} = rho*y_cur["z"] + eps_k is itself the
  ## next-period z, which IS what we feed as the lag to the t+2 step.
  ##
  ## In the model equation z = rho_z*z(-1) + eps_z, the RHS of z at t+1 is
  ## rho_z*z_t + eps_{t+1}. This equals z_{t+1}, so z_{t+1} is the next-period
  ## z value. But z_{t+1}(-1) (lag of z at t+1) = z_t (current z).
  ##
  ## So: to evaluate policy at time t+1, we supply lag = (k_t, z_t),
  ## and the policy function produces (c_{t+1}, k_{t+1}, z_{t+1}).
  ## The model equation z_{t+1} = rho*z_t + eps_{t+1} IS one of the model
  ## residuals, so z_{t+1} is constrained by it.
  ##
  ## However, when computing the LEAD residual at time t, we need y_{t+1}
  ## including z_{t+1}. Since z_{t+1} = rho*z_t + eps_{t+1}, and we're
  ## averaging over eps_{t+1} in quadrature, for each quadrature node k:
  ##   z_{t+1,k} = rho_z * z_t + eps_{t+1,k}   (AR(1) propagation)
  ##   lag for policy at t+1: (k_t, z_{t+1,k})
  ##   y_{t+1,k} = policy(k_t, z_{t+1,k})
  ## But y_{t+1,k} includes k_{t+1}, c_{t+1}, z_{t+1} — all from policy.
  ## z from policy at (k_t, z_{t+1,k}) should equal z_{t+1,k} itself
  ## (the z equation is in the model), so it's consistent.
  ##
  ## PRACTICAL APPROACH: compute lag for t+1 from current y and each eps draw.
  ## Lag for t+1 = (y_cur["k"], z_{t+1,k}) where z_{t+1,k} = rho*y_cur["z"]+eps_k.
  ##
  ## This requires knowing which state is an AR(1) shock process and which is
  ## a capital-like accumulation. We detect this from the model structure:
  ## - A state nm with a corresponding shock "eps_<nm>" is AR(1) propagated
  ## - Other states have their lag carried forward as-is (policy determines them)
  ##
  ## For the rbc_simple model:
  ##   state_names = c("k", "z")
  ##   k: no eps_k; lag for t+1 = current k from policy = y_cur["k"]
  ##   z: eps_z exists; lag for t+1 = rho_z * y_cur["z"] + eps_{t+1,k}
  ##
  ## So the lag vector fed to the t+1 policy evaluation is:
  ##   next_lag_k[j] = y_cur[state_names[j]]     if no corresponding shock
  ##   next_lag_k[j] = rho_j*y_cur[j] + eps_k    if AR(1) shock state

  ## Get AR(1) parameters for each state (NULL if not AR(1) shock state)
  ar1_rho       <- setNames(vector("list", n_state), state_names)
  ar1_shock_idx <- setNames(rep(NA_integer_, n_state), state_names)
  for (j in seq_len(n_state)) {
    nm       <- state_names[j]
    shock_nm <- paste0("eps_", nm)
    si       <- match(shock_nm, exo)
    if (!is.na(si)) {
      rho_nm <- paste0("rho_", nm)
      rho_val <- if (rho_nm %in% names(params)) params[[rho_nm]] else
                 if ("rho" %in% names(params)) params[["rho"]] else 0.9
      ar1_rho[[nm]]       <- rho_val
      ar1_shock_idx[nm]   <- si
    }
  }

  ## Compute the "next lag" vector given y_cur and shock_eps_next (next-period shocks).
  ##
  ## For the Euler expectation at time t, we need y_{t+1} for each draw eps_{t+1}.
  ## y_{t+1} = policy(next_lag_{t+1}) where next_lag_{t+1} is what we feed into
  ## the policy function to get time-t+1 values.
  ##
  ## The policy maps (k_lag, z_lag) -> y_cur where y_cur["z"] = rho*z_lag (shock=0).
  ## To get y_{t+1}["z"] = rho*z_t + eps_{t+1}, we need to feed:
  ##   next_lag["z"] = (rho*z_t + eps_{t+1}) / rho = z_t + eps_{t+1}/rho
  ## because then policy(..., next_lag["z"])["z"] = rho*(z_t + eps_{t+1}/rho) = rho*z_t + eps_{t+1}.
  ## Similarly, c_{t+1} from the policy at (k_t, z_t + eps_{t+1}/rho) will depend on
  ## the correct z_{t+1} via the resource constraint.
  ##
  ## For capital-like states (k): next_lag["k"] = y_cur["k"] (current k is lag at t+1).
  compute_next_lag <- function(y_cur, shock_eps_next) {
    next_lag <- y_cur[state_names]
    for (j in seq_len(n_state)) {
      nm <- state_names[j]
      si <- ar1_shock_idx[nm]
      if (!is.na(si)) {
        rho_j <- ar1_rho[[nm]]
        ## Feed z_lag that, when processed by the policy, gives
        ## z_{t+1} = rho*z_t + eps_{t+1}:
        ## z_lag_feed = (rho*y_cur[nm] + eps) / rho = y_cur[nm] + eps/rho
        next_lag[nm] <- y_cur[nm] + shock_eps_next[si] / rho_j
      }
      ## For non-AR(1) states (like k): next_lag[nm] = y_cur[nm] (carry forward)
    }
    next_lag
  }

  ## ------------------------------------------------------------------
  ## Build the dy template once (named zeros, length tot_cols)
  ## ------------------------------------------------------------------
  dy_template <- setNames(numeric(tot_cols), dy_keys)

  ## Assemble dy given y_lag, y_cur, y_lead, shock_cur.
  ## y_lag, y_cur, y_lead must be NAMED numeric vectors (over endo names).
  ## y_lag only needs to contain state-variable entries (non-state lag
  ## positions are not in the dyn_col_map and idx_m1 will be NA for them).
  assemble_dy <- function(y_lag, y_cur, y_lead, shock_cur) {
    dy <- dy_template
    ## Lag values (var__m1): only non-NA idx_m1 entries exist in col_map
    for (j in seq_len(n_endo)) {
      ci <- idx_m1[j]
      if (!is.na(ci)) {
        v <- y_lag[endo[j]]   # named lookup; NA if not in y_lag
        if (!is.na(v)) dy[ci] <- v
      }
    }
    ## Current values (var__0): all endo vars have a t=0 entry
    for (j in seq_len(n_endo)) {
      ci <- idx_0[j]
      if (!is.na(ci)) dy[ci] <- y_cur[endo[j]]
    }
    ## Lead values (var__p1): only vars with a future appearance
    for (j in seq_len(n_endo)) {
      ci <- idx_p1[j]
      if (!is.na(ci)) {
        v <- y_lead[endo[j]]
        if (!is.na(v)) dy[ci] <- v
      }
    }
    ## Shock values (shock__0)
    for (s in seq_len(n_exo)) {
      ci <- idx_exo[s]
      if (!is.na(ci)) dy[ci] <- shock_cur[exo[s]]
    }
    dy
  }

  ## ------------------------------------------------------------------
  ## Compute the expected residual at a collocation node
  ## given y_lag (lag values) and y_cur_vec (candidate current values)
  ## Returns a numeric vector of length n_eq.
  ## ------------------------------------------------------------------
  expected_residual <- function(y_lag_vec, y_cur_vec) {
    ## Build shock combinations (tensor product for multi-shock models)
    ## Each draw = sigma_k * standard_node_k
    if (n_exo == 1L) {
      shock_mat  <- matrix(q_nodes * shock_sds[1L], ncol = 1L)
      w_combined <- q_weights
    } else {
      qn_list <- lapply(seq_len(n_exo), function(k) q_nodes * shock_sds[k])
      shock_mat  <- as.matrix(expand.grid(qn_list))
      qw_list    <- replicate(n_exo, q_weights, simplify = FALSE)
      w_mat      <- as.matrix(expand.grid(qw_list))
      w_combined <- apply(w_mat, 1L, prod)
      w_combined <- w_combined / sum(w_combined)
    }
    n_combo <- nrow(shock_mat)

    acc <- numeric(n_eq)

    for (ki in seq_len(n_combo)) {
      ## Next-period shock draw (for t+1)
      eps_next <- shock_mat[ki, ]

      ## Next-period state lag = f(y_cur, eps_next)
      next_lag_vec <- compute_next_lag(y_cur_vec, eps_next)

      ## Evaluate policy at next-period lag -> y_{t+1}
      next_lag_norm <- matrix(norm_state_lag(next_lag_vec), nrow = 1L)
      y_lead_mat    <- eval_policy(next_lag_norm)
      y_lead_vec    <- as.vector(y_lead_mat)
      names(y_lead_vec) <- endo

      ## Current-period shock is 0 (expectations at t; shock at t is integrated out
      ## by setting it to 0 for the Euler equation approach).
      ## For the contemporaneous equations (e.g., resource constraint, AR(1) shock),
      ## the current shock IS the realised shock. But we are computing E_t[F(...)],
      ## where the CURRENT shock eps_t is already realised.
      ## Convention: set current shock = 0 (expected residual approach).
      shock_cur <- setNames(rep(0.0, n_exo), exo)

      ## Assemble dy
      dy <- assemble_dy(y_lag_vec, y_cur_vec, y_lead_vec, shock_cur)

      ## Guard: non-finite values in dy
      bad <- which(!is.finite(dy))
      if (length(bad) > 0L) {
        stop(sprintf(
          "solve_global: non-finite dy entries [%s] at ki=%d.\n  y_lag=%s\n  y_cur=%s",
          paste(dy_keys[bad], collapse = ","), ki,
          paste(sprintf("%s=%.4g", endo, y_lag_vec), collapse = " "),
          paste(sprintf("%s=%.4g", endo, y_cur_vec), collapse = " ")))
      }

      ## Evaluate residual
      res_k <- residuals_fn(dy, params, ss_vals)

      if (any(!is.finite(res_k))) {
        bad_eq <- which(!is.finite(res_k))
        stop(sprintf(
          "solve_global: residuals_fn returned non-finite value(s) at eq(s) [%s].\n  ki=%d, dy=%s",
          paste(bad_eq, collapse = ","), ki,
          paste(sprintf("%s=%.4g", dy_keys, dy), collapse = " ")))
      }

      acc <- acc + w_combined[ki] * res_k
    }
    acc
  }

  ## ------------------------------------------------------------------
  ## Inner solve: find y_cur at a collocation node by rooting residual
  ## ------------------------------------------------------------------
  solve_at_node <- function(y_lag_vec, y_init) {
    nleqslv::nleqslv(
      x   = y_init,
      fn  = function(y) expected_residual(y_lag_vec, y),
      method  = "Broyden",
      control = list(maxit = 300L, ftol = 1e-9, xtol = 1e-9,
                     allowSingular = TRUE, trace = 0L)
    )
  }

  ## ------------------------------------------------------------------
  ## Main iteration
  ## ------------------------------------------------------------------
  last_delta <- Inf
  converged  <- FALSE

  for (iter in seq_len(max_iter)) {

    new_ycur_mat <- matrix(NA_real_, nrow = n_coll, ncol = n_endo)
    colnames(new_ycur_mat) <- endo

    n_fail <- 0L
    for (j in seq_len(n_coll)) {
      ## State lag at this collocation node
      y_lag_j <- grid_nat[j, ]   # named: state_names

      ## Initial guess for y_cur from current policy
      phi_j   <- Phi[j, , drop = FALSE]            # 1 x n_basis
      y_init  <- as.vector(phi_j %*% t(coefs))     # n_endo
      names(y_init) <- endo

      ## Solve
      sol <- tryCatch(
        solve_at_node(y_lag_j, y_init),
        error = function(e) list(termcd = 99L, x = y_init)
      )

      if (sol$termcd <= 3L) {
        new_ycur_mat[j, ] <- sol$x
      } else {
        n_fail <- n_fail + 1L
        new_ycur_mat[j, ] <- y_init   # keep old guess
      }
    }

    if (verbose && n_fail > 0L)
      cat(sprintf("  iter %3d: %d/%d nodes failed nleqslv\n",
                  iter, n_fail, n_coll))

    ## Update coefficients via least squares
    new_coefs <- matrix(NA_real_, nrow = n_endo, ncol = n_basis)
    rownames(new_coefs) <- endo
    for (i in seq_len(n_endo)) {
      fit <- lm.fit(Phi, new_ycur_mat[, i])
      new_coefs[i, ] <- fit$coefficients
    }

    last_delta <- max(abs(new_coefs - coefs))
    if (verbose)
      cat(sprintf("  iter %3d: max_delta_coef = %.3e  failures = %d/%d\n",
                  iter, last_delta, n_fail, n_coll))

    coefs <- new_coefs

    if (last_delta < tol && n_fail == 0L) {
      converged <- TRUE
      if (verbose)
        cat(sprintf("solve_global: converged in %d iterations (delta=%.2e)\n",
                    iter, last_delta))
      break
    }
  }

  if (!converged && verbose)
    cat(sprintf(
      "solve_global: did NOT converge in %d iter (delta=%.2e, failures=%d)\n",
      max_iter, last_delta, n_fail))

  ## ------------------------------------------------------------------
  ## Return GlobalSolution
  ## ------------------------------------------------------------------
  structure(
    list(
      coefs          = coefs,          # n_endo x n_basis
      poly_degree    = poly_degree,
      state_names    = state_names,
      all_endo_names = endo,           # all endogenous vars
      shock_names    = exo,
      shock_sds      = shock_sds,
      state_domain   = state_domain,
      ss_vals        = ss_vals,
      params         = params,
      compiled       = compiled,
      n_quad         = n_quad,
      n_nodes        = n_nodes,
      converged      = converged,
      n_iter         = if (converged) iter else max_iter,
      last_delta     = last_delta,
      tol            = tol,
      assemble_dy    = assemble_dy,     # exposed for Euler error computation
      residuals_fn   = dyn$residuals_fn,
      compute_next_lag = compute_next_lag,
      norm_state_lag = norm_state_lag,
      eval_policy    = eval_policy
    ),
    class = "GlobalSolution"
  )
}

## =========================================================================
## Internal helpers
## =========================================================================

## Get shock standard deviations from the model's shocks block.
## Returns a named numeric vector of length n_exo.
.get_shock_sds <- function(model, params) {
  exo <- model$varexo_names
  sds <- setNames(rep(0.01, length(exo)), exo)
  if (!is.null(model$shocks)) {
    for (sh in model$shocks) {
      nm <- sh$name
      if (!is.null(nm) && nm %in% names(sds)) {
        val <- sh$stderr
        if (is.character(val)) {
          val <- if (val %in% names(params)) params[[val]] else 0.01
        }
        if (!is.null(val) && is.numeric(val) && is.finite(val) && val > 0)
          sds[nm] <- val
      }
    }
  }
  sds
}

## Compute default state domain from SS ± coverage based on shock variances.
.auto_state_domain <- function(state_names, ss_vals, shock_sds, model) {
  domain <- vector("list", length(state_names))
  names(domain) <- state_names
  params <- model$param_values

  for (nm in state_names) {
    ss_val <- ss_vals[nm]
    shock_nm <- paste0("eps_", nm)
    if (shock_nm %in% names(shock_sds)) {
      ## AR(1) shock state: stationary variance = sigma^2 / (1 - rho^2)
      rho_nm  <- paste0("rho_", nm)
      rho_val <- if (!is.null(params) && rho_nm %in% names(params))
        params[[rho_nm]] else 0.9
      sigma_stat <- shock_sds[shock_nm] / sqrt(max(1 - rho_val^2, 0.01))
      lo <- ss_val - 3.5 * sigma_stat
      hi <- ss_val + 3.5 * sigma_stat
    } else {
      ## Capital-like: ±40% of SS value
      half_width <- max(0.4 * abs(ss_val), 0.1)
      lo <- ss_val - half_width
      hi <- ss_val + half_width
    }
    domain[[nm]] <- c(lo, hi)
  }
  domain
}

## Initialize polynomial coefficients from order-1 perturbation.
## At each collocation node (state_lag), the linear policy gives
## y_cur ≈ ss + ghx * (state_lag - ss_state).
## We regress these predicted values on the Chebyshev basis.
.init_coefs_from_dr1 <- function(dr1, endo, state_names, ss_vals,
                                 grid_nat, state_domain, Phi, n_coll, n_basis) {
  n_endo  <- length(endo)
  n_state <- length(state_names)
  ghx     <- dr1$ghx   # n_endo x n_state_dr; rownames = endo, colnames = dr1$state_vars

  coefs <- matrix(0.0, nrow = n_endo, ncol = n_basis)
  rownames(coefs) <- endo

  ## At each grid node, compute the order-1 approximation for each endo var
  ycur_approx <- matrix(0.0, nrow = n_coll, ncol = n_endo)
  for (j in seq_len(n_coll)) {
    for (i in seq_len(n_endo)) {
      nm  <- endo[i]
      val <- ss_vals[nm]
      ## Add linear correction from state deviations
      for (s in seq_len(n_state)) {
        snm <- state_names[s]
        if (snm %in% colnames(ghx) && nm %in% rownames(ghx)) {
          dev_s <- grid_nat[j, snm] - ss_vals[snm]
          val   <- val + ghx[nm, snm] * dev_s
        }
      }
      ycur_approx[j, i] <- val
    }
  }

  ## Regress each endo var's linear approximation onto the basis
  for (i in seq_len(n_endo)) {
    fit <- lm.fit(Phi, ycur_approx[, i])
    coefs[i, ] <- fit$coefficients
  }
  coefs
}
