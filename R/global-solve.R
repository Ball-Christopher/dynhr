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
#'   \code{state_domain}, \code{poly_degree}, \code{converged}, the order-1
#'   stationary state sds \code{state_sd}, and two domain-clipping
#'   diagnostics.  \code{predict.GlobalSolution()} CLIPS a state outside
#'   \code{state_domain} to the boundary, silently, so the Euler quadrature's
#'   own \eqn{t+1} points leaving the box would corrupt the solution
#'   invisibly; both diagnostics report the FRACTION of those \eqn{t+1} state
#'   coordinates that do so, at the converged policy.
#'   \code{domain_clip_solve} measures it over the collocation grid (whose
#'   corners are the extremes of the box, so a small positive value is
#'   normal); \code{domain_clip_frac} measures it over the ergodic box
#'   (steady state \eqn{\pm 3} order-1 stationary sds of each state), where
#'   it is expected to be ZERO and a positive value raises a warning.  Both
#'   are \code{NA} when no order-1 solution is available.
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

  ## ---- Argument sanity, BEFORE anything can fail obscurely ---------------
  ## The signature is solve_global(compiled, ss, params, ...). A two-argument
  ## call -- solve_global(cm, model$param_values) -- silently binds the
  ## PARAMETER vector to `ss` and leaves `params` a missing promise, which
  ## used to surface hundreds of Newton solves later as "NA/NaN/Inf in 'x'"
  ## out of lm.fit. Say what is actually wrong instead.
  if (missing(params) || is.null(params) || !is.numeric(params) ||
      is.null(names(params)))
    stop("solve_global: `params` must be a NAMED numeric parameter vector. ",
         "Note the argument order -- solve_global(compiled, ss, params, ...) ",
         "-- so a two-argument call passes the parameters as `ss`.",
         call. = FALSE)

  ss_vals <- if (inherits(ss, "dynhr_steady")) ss$values else ss
  if (!is.numeric(ss_vals) || is.null(names(ss_vals)) ||
      !all(endo %in% names(ss_vals)))
    stop("solve_global: `ss` must be a solve_steady() result or a NAMED ",
         "numeric steady state covering every endogenous variable (missing: ",
         paste(setdiff(endo, names(ss_vals)), collapse = ", "), ").",
         call. = FALSE)

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
  state_names  <- .global_state_names(compiled, endo)
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
  ## 2b. VALIDATE the model class, before anything can fail obscurely.
  ## ------------------------------------------------------------------
  ## The projection scheme represents only a narrow class of models (see
  ## .global_shock_pairing() below for the exact contract). Anything outside
  ## it used to produce either a converged-but-WRONG deterministic policy
  ## (when the name-based `eps_<state>` pairing silently found nothing) or an
  ## opaque "NA/NaN/Inf in 'x'" from nleqslv/lm.fit several hundred Newton
  ## solves later. Fail here instead, naming the offending shock/equation.
  pairing <- .global_shock_pairing(compiled, params, ss_vals, state_names)

  ## ------------------------------------------------------------------
  ## 3. Shock standard deviations
  ## ------------------------------------------------------------------
  shock_sds <- .get_shock_sds(model, params)

  ## ------------------------------------------------------------------
  ## 4. State domain
  ## ------------------------------------------------------------------
  ## The order-1 solution is needed TWICE -- for the shock-aware default
  ## domain (the stationary sd of the endogenous, capital-like states) and,
  ## further down, for the warm start -- so it is solved once, here, and
  ## carried forward.  A failure is non-fatal in both places, and a caller
  ## who supplies `state_domain` AND opts out of the warm start pays nothing.
  dr1 <- if (is.null(state_domain) || isTRUE(init_from_perturbation)) {
    tryCatch(
      solve_perturbation(model, compiled, ss_vals, params, verbose = FALSE),
      error = function(e) {
        if (verbose) cat("solve_global: order-1 perturbation failed:",
                         conditionMessage(e), "\n")
        NULL
      })
  } else NULL

  ## The order-1 stationary state sds size BOTH the default box and the
  ## post-solve clipping diagnostic (§ "domain clipping" below), so they are
  ## computed once here rather than in the `is.null(state_domain)` branch --
  ## a caller who supplies an explicit box still gets the diagnostic.
  sd_state <- .dr_state_sd(dr1, model, params)

  if (is.null(state_domain)) {
    state_domain <- .auto_state_domain(
      state_names = state_names,
      ss_vals     = ss_vals,
      shock_sds   = shock_sds,
      model       = model,
      pairing     = pairing,
      exo         = exo,
      state_sd    = sd_state,
      n_quad      = n_quad
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

  ## AR(1) parameters per state, from .global_shock_pairing() -- read off the
  ## compiled Jacobian, NOT off the parameter names. The old code required a
  ## shock literally called `eps_<state>` and a parameter called `rho_<state>`
  ## or `rho`, silently DEFAULTING rho to 0.9 when neither existed, and
  ## silently leaving the state unshocked when the name did not match.
  ## `ar1_psi` is the shock loading (1 for `z = rho*z(-1) + eps_z`, sigma for
  ## `nu = rho*nu(-1) + sigma*eps_nu`), which the old code assumed was 1.
  ar1_rho       <- pairing$rho
  ar1_psi       <- pairing$psi
  ar1_shock_idx <- pairing$shock_idx

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
        ## z_{t+1} = rho*z_t + psi*eps_{t+1}:
        ## z_lag_feed = (rho*y_cur[nm] + psi*eps) / rho = y_cur[nm] + psi*eps/rho
        ## psi == 1 for the plain `z = rho*z(-1) + eps_z` form, in which case
        ## `(1*eps)/rho` is bit-identical to the previous `eps/rho`.
        next_lag[nm] <- y_cur[nm] + ar1_psi[[nm]] * shock_eps_next[si] / rho_j
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
  ## 12. Domain-clipping diagnostic (F3-B)
  ## ------------------------------------------------------------------
  ## `norm_state_lag()` CLAMPS its argument to `state_domain`, so whenever
  ## the Euler quadrature's own t+1 feed point s_t + psi*eps_{t+1}/rho leaves
  ## the box, `eval_policy()` returns the value AT the boundary and the
  ## residual driven to zero is not the model's.  The clamp stays (an
  ## unclamped Chebyshev extrapolation diverges), but it is no longer
  ## invisible: the fraction of t+1 state COORDINATES that leave the box is
  ## measured at the converged policy and reported.
  ##
  ## Two boxes are measured, because they answer different questions:
  ##   `domain_clip_solve` -- over the COLLOCATION grid, i.e. what the solver
  ##      itself did.  Its corners are the extremes of the box, so a small
  ##      positive value is normal and is not a defect.
  ##   `domain_clip_frac`  -- over the ERGODIC box (steady state +- `n_sd`
  ##      order-1 stationary sds of each state), i.e. where the data live.
  ##      This one is expected to be ZERO; anything positive means the box is
  ##      too narrow for its own quadrature and a warning is issued.
  clip_at <- .gs_clip_frac(state_names = state_names,
                           state_domain = state_domain,
                           eval_policy = eval_policy,
                           norm_state_lag = norm_state_lag,
                           compute_next_lag = compute_next_lag,
                           endo = endo, exo = exo,
                           shock_sds = shock_sds, n_quad = n_quad)
  domain_clip_solve <- clip_at(grid_nat)

  erg_grid <- .gs_ergodic_grid(state_names, ss_vals, sd_state,
                               n_sd = .gs_clip_n_sd, n_grid = 5L)
  domain_clip_frac <- if (is.null(erg_grid)) NA_real_ else clip_at(erg_grid)

  if (isTRUE(is.finite(domain_clip_frac) && domain_clip_frac > 0))
    .gs_warn_once("domain_clip", sprintf(
      paste0("solve_global: %.2f%% of the Euler-quadrature t+1 state ",
             "coordinates leave `state_domain` over the ergodic +-%g ",
             "stationary-sd box, where predict() CLIPS them silently. The ",
             "collocation box is too narrow for its own quadrature; widen ",
             "it (`state_domain=`, or a larger `endo_cover`/`domain_cover` ",
             "through .auto_state_domain()). See `$domain_clip_frac`. ",
             "This warning fires at most once per session."),
      100 * domain_clip_frac, .gs_clip_n_sd))

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
      ## Domain-clipping diagnostics; see § 12 above.
      domain_clip_frac  = domain_clip_frac,
      domain_clip_solve = domain_clip_solve,
      state_sd          = sd_state,
      assemble_dy    = assemble_dy,     # exposed for Euler error computation
      residuals_fn   = dyn$residuals_fn,
      ## AR(1) "feed" parameters, per state (NA for capital-like states).
      ## compute_next_lag() applies them one state vector at a time; the
      ## global particle filter (R/global-likelihood.R) needs the SAME map
      ## applied to an N x n_state particle matrix at once, so the pieces are
      ## exposed rather than re-derived (a second derivation is a second
      ## chance to disagree with the solver about the timing convention).
      ar1_rho        = ar1_rho,
      ar1_psi        = ar1_psi,
      ar1_shock_idx  = ar1_shock_idx,
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
##
## Delegates to `.get_shock_stderr()` -- the package's canonical, estimation-
## aware resolver (parameter-named stderr > shocks-block *_expr re-evaluated
## against the CURRENT params > frozen parse-time numerics).  Two bugs are
## fixed by that delegation:
##   * The old body looped `for (sh in model$shocks)`, but `model$shocks` is
##     `list(variances = <data.frame>, correlations = <data.frame>)`, so `sh`
##     was a WHOLE data.frame and `sh$name` its whole name COLUMN.  With one
##     shock the length-1 column happened to work; with two or more,
##     `nm %in% names(sds)` is length > 1 and the condition ERRORED -- i.e.
##     solve_global() could not run on any multi-shock model.
##   * Even on the single-shock path it read the frozen numeric `stderr`, so
##     an ESTIMATED shock standard deviation (`stderr sig_a;`) had no effect
##     on the projection solution -- the P0 bug, in the global solver.
## The 0.01 fallback is kept for shocks whose std resolves to 0 (no shocks
## block at all): a zero-variance shock would collapse the quadrature.
.get_shock_sds <- function(model, params) {
  exo <- model$varexo_names
  sds <- setNames(rep(0.01, length(exo)), exo)
  se  <- tryCatch(.get_shock_stderr(model, exo, params),
                  error = function(e) NULL)
  if (!is.null(se)) {
    ok <- is.finite(se) & se > 0
    sds[ok] <- se[ok]
  }
  sds
}

## Unconditional (stationary) standard deviation of each state under the
## ORDER-1 solution, as a named vector over `dr$state_vars`.
##
## Used to size the default collocation box: the Chebyshev nodes should sit
## where the data live, and for an endogenous, capital-like state the only
## model-derived scale available before the projection is solved is the
## linear model's own stationary sd.
##
## `shock_sds` optionally overrides the model's shock standard deviations
## (the SBC harness sizes its box at a prior-UPPER sd).  The override
## rescales the shock covariance, so any declared correlations survive.
## Returns NULL when the rule cannot be applied (no order-1 solution, or a
## non-stationary Lyapunov solve).
#' @noRd
.dr_state_sd <- function(dr, model, params, shock_sds = NULL) {
  if (is.null(dr) || is.null(dr$state_vars) || !length(dr$state_vars))
    return(NULL)
  TT <- dr$ghx[dr$state_idx, , drop = FALSE]
  RR <- dr$ghu[dr$state_idx, , drop = FALSE]
  Se <- tryCatch(.get_shock_cov(model, dr$exo_names, params),
                 error = function(e) NULL)
  if (is.null(Se)) return(NULL)
  if (!is.null(shock_sds)) {
    s0 <- sqrt(pmax(diag(as.matrix(Se)), 0))
    sn <- shock_sds[dr$exo_names]
    ok <- is.finite(s0) & s0 > 0 & is.finite(sn) & sn > 0
    if (any(ok)) {
      f      <- rep(1, length(s0))
      f[ok]  <- sn[ok] / s0[ok]
      D      <- diag(f, nrow = length(f))
      Se     <- D %*% as.matrix(Se) %*% D
    }
  }
  P0 <- tryCatch(kf_stationary_init(TT, RR, Se), error = function(e) NULL)
  if (is.null(P0)) return(NULL)
  stats::setNames(sqrt(pmax(diag(P0), 0)), dr$state_vars)
}


## =========================================================================
## Domain-clipping diagnostic (F3-B)
## =========================================================================

## Half-width, in order-1 stationary state sds, of the box over which
## `domain_clip_frac` is measured (and over which the warning fires).  Three
## sds of every state jointly is already a rare corner of the ergodic set;
## the point of the measure is that even THERE the quadrature's t+1 feed must
## stay inside the collocation box.
.gs_clip_n_sd <- 3

## Package-private store for one-time warnings (the projection solve is
## re-run at EVERY theta inside make_log_posterior_global_pf(), so a
## per-solve warning would spam an entire MCMC run).  Same pattern as
## `.hank_egm_warn_once()` / `.cumulant_warn_once()`.
.gs_warn_env <- new.env(parent = emptyenv())

#' Emit a warning at most once per session, keyed by \code{key}.
#' @keywords internal
#' @noRd
.gs_warn_once <- function(key, ...) {
  if (isTRUE(.gs_warn_env[[key]])) return(invisible(NULL))
  .gs_warn_env[[key]] <- TRUE
  warning(paste0(...), call. = FALSE)
}

## Largest Gauss-Hermite abscissa at `n_quad` nodes, in shock-sd units: how
## far into the tail of eps_{t+1} the Euler quadrature actually reaches.  The
## default AR(1) cover is sized so that even THAT node stays inside the box.
#' @noRd
.gs_quad_reach <- function(n_quad) {
  q <- tryCatch(gauss_hermite(n_quad)$nodes, error = function(e) NULL)
  if (is.null(q) || !length(q) || !all(is.finite(q))) return(3)
  max(abs(q))
}

## Tensor grid of steady state +- `n_sd` order-1 stationary sds per state.
## NULL when the stationary sds are unavailable (no order-1 solution, or a
## non-stationary Lyapunov solve), in which case the diagnostic is NA rather
## than silently computed on some other box.
#' @noRd
.gs_ergodic_grid <- function(state_names, ss_vals, state_sd, n_sd = 3,
                             n_grid = 5L) {
  if (is.null(state_sd) || !all(state_names %in% names(state_sd)))
    return(NULL)
  sdv <- state_sd[state_names]
  if (!all(is.finite(sdv)) || any(sdv <= 0)) return(NULL)
  g1 <- lapply(state_names, function(nm)
    ss_vals[[nm]] + seq(-n_sd, n_sd, length.out = n_grid) * state_sd[[nm]])
  names(g1) <- state_names
  args <- rev(g1); names(args) <- rev(state_names)
  out <- as.matrix(expand.grid(args))
  out[, state_names, drop = FALSE]
}

## Build the clip-fraction measurement for a solved policy.  Returns a
## function of a matrix of state LAGS; for each row it evaluates the policy,
## forms the Gauss-Hermite t+1 feed points the Euler equation integrates over
## (the solver's OWN `compute_next_lag()`, so the two cannot disagree about
## the timing) and returns the fraction of feed COORDINATES outside the box.
##
#' @noRd
.gs_clip_frac <- function(state_names, state_domain, eval_policy,
                          norm_state_lag, compute_next_lag, endo, exo,
                          shock_sds, n_quad) {
  n_state <- length(state_names)
  n_exo   <- length(exo)
  lo <- vapply(state_names, function(nm) state_domain[[nm]][1L], numeric(1))
  hi <- vapply(state_names, function(nm) state_domain[[nm]][2L], numeric(1))

  ## The SAME tensor-product Gauss-Hermite node set expected_residual() uses.
  gh <- gauss_hermite(n_quad)
  shock_mat <- if (n_exo == 1L) matrix(gh$nodes * shock_sds[1L], ncol = 1L) else
    as.matrix(expand.grid(lapply(seq_len(n_exo),
                                 function(k) gh$nodes * shock_sds[k])))
  n_combo <- nrow(shock_mat)

  function(grid_nat) {
    if (is.null(grid_nat) || !nrow(grid_nat)) return(NA_real_)
    n_out <- 0L
    n_tot <- 0L
    for (j in seq_len(nrow(grid_nat))) {
      s_lag <- grid_nat[j, state_names]
      y <- as.vector(eval_policy(matrix(norm_state_lag(s_lag), nrow = 1L)))
      names(y) <- endo
      if (!all(is.finite(y))) next
      for (ki in seq_len(n_combo)) {
        nl <- compute_next_lag(y, shock_mat[ki, ])[state_names]
        n_out <- n_out + sum(nl < lo | nl > hi, na.rm = TRUE)
        n_tot <- n_tot + n_state
      }
    }
    if (n_tot == 0L) NA_real_ else n_out / n_tot
  }
}


## Compute default state domain from SS ± coverage based on shock variances.
##
## `pairing` is .global_shock_pairing()'s structural rho/psi/shock_idx map, so
## the persistence and the effective innovation std |psi|*sigma come from the
## MODEL rather than from a `rho_<state>` name lookup with a 0.9 default.
##
## TWO shock-awareness fixes (F2-C):
##
##  * AR(1) states are covered at the FEED point s_{t-1} + psi*eps_t/rho --
##    the argument the filter and simulator actually evaluate the policy at
##    -- whose stationary sd is `sigma_stat / rho`, STRICTLY wider than the
##    state's own.  The old rule covered `domain_cover` sds of the STATE, so
##    predict()'s silent clip bit in the tails.
##
##  * Endogenous, capital-like states (no paired shock) used a fixed
##    ±40% of the steady-state LEVEL, regardless of shock size.  Measured on
##    the F1-C two-shock RBC, the resulting box is 59 stationary sds wide at
##    `sd_scale = 0.1` and only 5.9 at `sd_scale = 1` -- i.e. the SAME rule
##    lands on opposite sides of the accuracy optimum depending on the
##    calibration, which is exactly what a shock-blind rule does.  Both ends
##    hurt, for different reasons:
##      - too WIDE and a degree-3 Chebyshev fit is spread over curvature that
##        the ergodic set never visits, so the fit near the data is loose;
##      - too NARROW and the Euler quadrature's own t+1 states leave the box,
##        where `predict()` CLIPS silently, so the residual being driven to
##        zero is not the model's.
##    `state_sd` (the order-1 stationary sd, `.dr_state_sd()`) puts the box
##    at a FIXED number of stationary sds -- `endo_cover` -- so it tracks the
##    shocks instead of the level.  The old ±40% rule survives as a CAP
##    (`level_cap_frac`) and as the fallback when no order-1 solution is
##    available.  `floor_frac` keeps a degenerate (zero-sd) state from
##    collapsing the grid to a point.
##
## WHAT F3-B CHANGED, and why the F2-C numbers above no longer read the same.
## F2-C set `endo_cover = 12` from ONE fixture at ONE shock scale, and the
## U-shape it fitted had a NARROW arm that was not a property of the fit at
## all: with `domain_cover = 3.5` the AR(1) box did not reach the largest
## Gauss-Hermite node, so 2.7-4.0% of the Euler quadrature's own t+1 state
## coordinates over the ergodic ±3 sd box were being CLIPPED -- on all four
## F3-B fixtures, at every `endo_cover`, and invisibly.  Once the AR(1)
## branch is floored at the quadrature's reach (see the `quad_reach` term
## below) the narrow arm collapses and the optimum moves to ~6 stationary
## sds and stays there.  Measured max |k_proj - k_order3| over a ±3
## stationary-sd box in units of the stationary sd of k, on the F1-C
## two-shock RBC (n_quad 5, n_nodes 5), cover in stationary sds of k:
##
##             4 sd     6 sd     8 sd    10 sd    12 sd    20 sd    30 sd
##   deg 3, sd_scale 1
##            2.8e-2   9.4e-3   2.0e-2   4.0e-2   8.1e-2   1.4e-1   2.8e-1
##   deg 3, sd_scale 0.1
##            5.9e-3   6.7e-5   7.3e-5   8.1e-5   9.1e-5   1.6e-4   3.2e-4
##   deg 4, sd_scale 1
##            1.8e-2   3.3e-3   4.4e-3   1.0e-2   3.0e-2      --       --
##
## So `endo_cover = 6` is at or within 1.1x of the per-fixture optimum at
## both shock scales and both degrees, whereas 12 is 8.7x (degree 3) and
## 9.3x (degree 4) off it at `sd_scale = 1`.  It is a CONSTANT, not a rule:
## the capital-like optimum did not move with the calibration once the AR(1)
## clipping was gone, which is exactly what a "fixed number of stationary
## sds" is supposed to deliver and what the shock-blind ±40% rule could not.
.auto_state_domain <- function(state_names, ss_vals, shock_sds, model,
                               pairing, exo, domain_cover = 3.5,
                               endo_cover = 6, state_sd = NULL,
                               floor_frac = 1e-3, level_cap_frac = 0.4,
                               level_cap_min = 0.1, pos_floor_frac = 0.05,
                               boundary_frac = 0.4,
                               n_quad = 5L, erg_cover = .gs_clip_n_sd,
                               reach_margin = 1.05) {
  domain <- vector("list", length(state_names))
  names(domain) <- state_names

  for (nm in state_names) {
    ss_val <- ss_vals[nm]
    si     <- pairing$shock_idx[[nm]]
    ## Hoisted: the boundary bound below needs it as well as the cover rule.
    sd_s   <- if (!is.null(state_sd) && nm %in% names(state_sd))
                state_sd[[nm]] else NA_real_
    if (!is.na(si)) {
      ## AR(1) shock state: stationary variance = (psi*sigma)^2 / (1 - rho^2);
      ## the FEED point adds (psi*sigma/rho)^2 on top.
      rho_val    <- pairing$rho[[nm]]
      sd_eff     <- abs(pairing$psi[[nm]]) * shock_sds[[exo[si]]]
      sigma_stat <- sd_eff / sqrt(max(1 - rho_val^2, 0.01))
      sigma_feed <- sqrt(sigma_stat^2 + (sd_eff / rho_val)^2)
      ## F3-B: `domain_cover` FEED sds is a distributional cover -- it says
      ## how much of the stationary feed distribution the box holds -- and it
      ## is NOT the quantity the Euler quadrature needs.  At a state lag
      ## s = erg_cover*sigma_stat the solver evaluates the policy at
      ##   rho*s + psi*q*sigma/rho,   q = the LARGEST Gauss-Hermite abscissa,
      ## a DETERMINISTIC point, and `predict()` clips it silently if the box
      ## does not reach it.  Measured on all four F3-B fixtures, the shipped
      ## `domain_cover = 3.5` did not: 2.7-4.0% of the t+1 coordinates over
      ## the ergodic +-3 sd box were clipped on EVERY fixture and at EVERY
      ## `endo_cover`, i.e. a defect of the AR(1) branch that no capital-like
      ## cover constant could have repaired.  The half-width is therefore the
      ## MAX of the distributional cover and the quadrature's own reach, so
      ## the default box is clip-free BY CONSTRUCTION and the
      ## `domain_clip_frac` warning below means what it says.  The floor is
      ## `n_quad`-aware because the reach is (q grows like sqrt(2*n_quad)).
      ## `reach_margin` is not decoration.  The t+1 feed is built from the
      ## POLICY value at the ergodic edge, not from rho*s exactly, so the
      ## interpolation error decides a bare equality: measured without a
      ## margin, `kpr` and `persist` still clipped 10 of 125 coordinates by
      ## ~1e-9 while `rbc2s` (the same rule, the error of the other sign)
      ## clipped none.  A knife-edge that a numerical noise term resolves is
      ## not a guarantee.
      quad_reach <- reach_margin *
        (abs(rho_val) * erg_cover * sigma_stat +
           .gs_quad_reach(n_quad) * sd_eff / abs(rho_val))
      half_width <- max(domain_cover * sigma_feed, quad_reach)
    } else {
      ## The old, shock-blind ±40%-of-level rule.  Since F3-B (follow-up) it
      ## is ONLY the no-order-1-solution fallback; it is no longer a CAP on
      ## the sd-based rule.
      ##
      ## WHY THE CAP HAD TO GO.  `max(level_cap_frac*|ss|, level_cap_min)` is
      ## a LEVEL bound: it does not scale with the shocks, so for a large
      ## enough sigma it truncates the box BELOW the ergodic set, which is
      ## precisely the silent clipping this wave exists to remove.  Measured
      ## on test-global-solve.R's `.rbc_mod_string(0.05)` fixture (the same
      ## full-depreciation RBC as `smallk`, alpha = 0.33, sigma_z = 0.05, so
      ## sd_k = 0.0311 against k_ss = 0.1883): `endo_cover * sd_k` asks for
      ## 0.1865 but the cap allowed only 0.1, i.e. 3.22 stationary sds -- a
      ## box barely wider than the +-3 sd ergodic set it has to hold.  The
      ## policy at the (k, z) = (+3 sd, +3 sd) corner returns k_t = 0.3008
      ## against a box top of 0.2883, and 2.0% of the Euler quadrature's t+1
      ## coordinates were clipped.  At sigma_z = 0.01 the same rule asks for
      ## 0.0373, the cap never binds, and nothing clips -- the defect is
      ## exactly the cap's shock-blindness.
      ##
      ## Note the reach cannot be floored the way the AR(1) branch is: the
      ## order-1 estimate of the policy's image of the ergodic box,
      ## erg_cover * sum_j |ghx[k,j]| * sigma_j = 0.089 here, UNDERSTATES the
      ## true image (0.1125) by 26%, because k = alpha*beta*exp(z)*k^alpha is
      ## strongly convex in z at sigma_z = 0.05. There is no cheap linear
      ## bound to floor against, so the cap is dropped rather than repaired.
      ##
      ## The cap's two stated jobs are both covered elsewhere now:
      ##   * "the new rule can only ever NARROW the box" was an F2-C
      ##     migration-safety argument, not a correctness one, and it is what
      ##     reintroduced the clipping;
      ##   * "keeps a positive-valued state off zero" is now the EXPLICIT
      ##     `pos_floor_frac` bound below, which is what actually does that
      ##     job here (k_ss - endo_cover*sd_k = 0.0018 is floored to 0.0094).
      ## On all four F3-B fixtures the cap sat at 7.8-21.8 stationary sds,
      ## i.e. above `endo_cover = 6`, so dropping it is a NO-OP for them and
      ## changes only the large-shock regime where it was doing harm.
      cap  <- max(level_cap_frac * abs(ss_val), level_cap_min)
      half_width <- if (is.finite(sd_s) && sd_s > 0) {
        max(endo_cover * sd_s, floor_frac * max(abs(ss_val), 1))
      } else cap
    }
    lo <- ss_val - half_width
    hi <- ss_val + half_width

    ## F3-B: an EXPLICIT lower bound for a positive-valued capital-like state.
    ## Before F3-B nothing stopped the box from crossing k = 0 -- the ±40%
    ## cap kept it off zero only by accident, and only while
    ## `level_cap_min` (0.1) happened to be smaller than the steady-state
    ## level. On a calibration whose k_ss is BELOW that floor (the F3-B
    ## `smallk` fixture, k_ss = 0.0888, alpha = 0.12) the capped half-width
    ## 0.1 puts collocation nodes at k < 0, where k^(alpha-1) is not real and
    ## the Newton solve returns NaN.  A capital-like state (no paired shock)
    ## with a strictly positive steady state is therefore floored at
    ## `pos_floor_frac` of its own level.  The box becomes ASYMMETRIC, which
    ## is correct: the natural boundary is one-sided, and so is the ergodic
    ## distribution near it.  AR(1) states are exempt -- a log-deviation
    ## process is centred at zero and has no such boundary.
    ##
    ## F3-B FOLLOW-UP: `pos_floor_frac` alone was set at 5%, which turned out
    ## to be far too PERMISSIVE, and that -- not the cover, and not the
    ## quadrature -- is what a large shock exposes.  On
    ## test-global-solve.R's `.rbc_mod_string(0.05)` fixture (sd_k/k_ss =
    ## 0.165) `endo_cover = 6` puts the lower edge 5.76 sd below the steady
    ## state, at 5% of k_ss: the box then spans a factor 40 in k, and a
    ## degree-3 Chebyshev fit is spread over the k -> 0 curvature (k^(alpha-1)
    ## diverges) that the ergodic set never visits.  Measured against this
    ## model's CLOSED-FORM policy, max |k_proj - k_exact| over the +-3 sd box
    ## in sd_k units, holding the +6 sd upper edge fixed and moving only the
    ## lower edge:
    ##     lower edge   5.76 sd   5.45 sd   4.85 sd   4.00 sd   3.03 sd
    ##     |proj-exact|  0.2326    0.1594    0.0920    0.0432    0.0244
    ## -- a 5.4x accuracy swing that is ENTIRELY the low edge (sweeping the
    ## UPPER edge at a fixed -4 sd lower edge gives 0.0725 / 0.0415 / 0.0432 /
    ## 0.0572 at +4 / +5 / +6 / +8 sd, i.e. flat, with the optimum at the
    ## shipped `endo_cover = 6`).
    ##
    ## So the boundary bound is TWO-sided, and the second half is what was
    ## missing: the floor may lift the lower edge, but it must never lift it
    ## into the ergodic set, and equally the box must never run closer to the
    ## boundary than the fit can afford.  `erg_cover + 1` (one stationary sd
    ## of margin below the +-erg_cover box the clipping diagnostic measures
    ## on) is that limit.
    ##
    ## WHEN IT BINDS.  At the shipped `endo_cover = 6`, `erg_cover = 3` and
    ## `boundary_frac = 0.4` the bound is active exactly when
    ## sd > ss/10 -- via the 0.4*ss term for 0.10 < sd/ss <= 0.15, and via
    ## the (erg_cover + 1) sd term above that.  So it is a NO-OP on all four
    ## original F3-B fixtures (sd/ss = 0.045-0.067) and on the same RBC at
    ## sigma_z = 0.01, and moves only the sigma_z = 0.05 case (sd/ss =
    ## 0.165), where it pulls the lower edge from 5.76 sd back to 4 sd and
    ## the accuracy from 0.233 to 0.043 sd.  (The threshold scales with the
    ## cover: at a hand-passed `endo_cover = 20` it would bind from
    ## sd > 0.03*ss.  That is intended -- the wider the requested box, the
    ## sooner it would otherwise reach the boundary.)
    if (is.na(si) && is.finite(ss_val) && ss_val > 0) {
      lo <- max(lo, pos_floor_frac * ss_val)
      if (is.finite(sd_s) && sd_s > 0)
        lo <- max(lo, min(boundary_frac * ss_val,
                          ss_val - (erg_cover + 1) * sd_s))
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


## =========================================================================
## Model-class validation (E1-B follow-up)
## -------------------------------------------------------------------------
## solve_global() can only represent a NARROW class of models, and until this
## function existed it never said so.  Its time-iteration convention sets the
## CURRENT shock to zero in the residual (see expected_residual()) and injects
## a realised shock only through the AR(1) "feed" lag inside
## compute_next_lag().  A shock that does not sit in an AR(1) state process
## therefore never enters the solution AT ALL -- silently.
##
## Worse, the pairing used to be by NAME: a state `nm` was treated as an
## AR(1) shock state only when a shock literally called `eps_<nm>` existed,
## and its persistence was read from a parameter called `rho_<nm>`, `rho`, or
## -- failing both -- DEFAULTED TO 0.9.  On
##
##     x = rho*x(-1) + e;   y = beta*x + 0.2*x(+1) + u;
##
## (shocks `e` and `u`; no `eps_x`) nothing paired, both shocks vanished from
## the policy, and the solver returned a converged, deterministic, WRONG
## solution.  Downstream, make_log_posterior_global_pf() then produced a
## perfectly finite log-likelihood of -425 where the Kalman filter says 151.
##
## The pairing is now STRUCTURAL, read from the compiled Jacobian rather than
## from names, and anything outside the representable class is a hard error
## naming the offending shock and equation.  For each shock s the equation it
## appears in must be exactly
##
##     a0 * nm(t) + a1 * nm(t-1) + b * s(t) = 0     (linear, no other terms)
##
## for a single state `nm`, giving rho = -a1/a0 and psi = -b/a0, i.e.
## nm_t = rho*nm_{t-1} + psi*s_t.  The feed that reproduces a realised shock
## through a policy solved at s_t = 0 is then
##
##     feed = nm_lag + psi * s / rho
##
## which generalises the old (psi = 1) form bit-identically.
##
## @return list(rho, psi, shock_idx) -- each a vector over `state_names`, with
##   NA for a state that carries no shock (a capital-like state, which IS
##   representable: its lag is simply carried forward).
#' @noRd
.global_shock_pairing <- function(compiled, params, ss_vals, state_names,
                                  context = "solve_global") {
  map  <- .acc_dy_map(compiled)
  dyn  <- compiled$dynamic
  exo  <- map$exo

  if (length(exo) == 0L)
    stop(context, ": the model declares no exogenous shocks, so there is ",
         "nothing for the projection solver to integrate over.", call. = FALSE)

  col0  <- map$col0_of
  colm1 <- stats::setNames(map$m1$col, map$m1$name)
  colex <- stats::setNames(map$ux$col, map$ux$name)

  ## Human-readable name for a set of dy columns, for the error messages.
  lbl <- function(cols) {
    out <- map$keys[cols]
    out <- sub("__p1$", "(t+1)", sub("__m1$", "(t-1)", sub("__0$", "(t)", out)))
    paste(out, collapse = ", ")
  }

  ## Two probe points, so that a coefficient which merely LOOKS constant at
  ## one point (a nonlinear AR process such as log(a) = rho*log(a(-1)) + eps,
  ## for which the feed formula is simply wrong) is detected.
  probe <- function(h) {
    dy <- map$template
    for (k in seq_along(map$keys)) {
      nm    <- sub("__(0|p1|m1)$", "", map$keys[k])
      v     <- if (nm %in% names(ss_vals)) ss_vals[[nm]] else 0
      dy[k] <- v + h * (1 + abs(v))
    }
    dy
  }
  J1 <- tryCatch(dyn$jacobian_fn(probe(1e-3), params, ss_vals),
                 error = function(e) NULL)
  J2 <- tryCatch(dyn$jacobian_fn(probe(-7e-3), params, ss_vals),
                 error = function(e) NULL)
  if (is.null(J1) || is.null(J2) || !all(is.finite(J1)) || !all(is.finite(J2)))
    stop(context, ": could not evaluate the model Jacobian near the steady ",
         "state, so the shock/state structure cannot be validated.",
         call. = FALSE)

  A   <- pmax(abs(J1), abs(J2))
  rsc <- pmax(apply(A, 1L, max), .Machine$double.eps)
  nzr <- function(q) which(A[q, ] > 1e-10 * rsc[q])

  rho       <- stats::setNames(rep(NA_real_, length(state_names)), state_names)
  psi       <- stats::setNames(rep(NA_real_, length(state_names)), state_names)
  shock_idx <- stats::setNames(rep(NA_integer_, length(state_names)),
                               state_names)

  for (si in seq_along(exo)) {
    s  <- exo[si]
    cs <- colex[[s]]
    if (is.null(cs) || is.na(cs))
      stop(context, ": shock '", s, "' has no column in the compiled dynamic ",
           "system, so it cannot enter the projection solution.",
           call. = FALSE)

    rows <- which(A[, cs] > 1e-10 * rsc)
    if (length(rows) != 1L)
      stop(context, ": shock '", s, "' enters ", length(rows), " equations (",
           if (length(rows)) paste0("equation(s) ",
                                    paste(rows, collapse = ", ")) else "none",
           "). The projection solver sets the CURRENT shock to zero in the ",
           "residual and re-injects it through an AR(1) state's feed lag, so ",
           "every shock must appear in exactly ONE equation, and that ",
           "equation must be an AR(1) state process `state = rho*state(-1) + ",
           "psi*", s, "`.", call. = FALSE)
    q <- rows

    nz  <- nzr(q)
    ## The single endogenous variable dated t in this equation.
    at0 <- names(col0)[match(intersect(nz, col0), col0)]
    at0 <- at0[!is.na(at0)]
    if (length(at0) != 1L)
      stop(context, ": shock '", s, "' enters equation ", q,
           ", which is not an AR(1) state process -- it involves ",
           if (length(at0)) paste0("the current-dated variables ",
                                   paste(at0, collapse = ", "))
           else "no current-dated endogenous variable",
           ". Full equation-", q, " terms: ", lbl(nz),
           ". solve_global() can only represent a shock that enters a single ",
           "equation of the form `state = rho*state(-1) + psi*", s, "`.",
           call. = FALSE)
    nm <- at0

    if (!nm %in% state_names)
      stop(context, ": shock '", s, "' drives '", nm, "' (equation ", q,
           "), which is not a state variable (it never appears at t-1). ",
           "The projection policy is a function of the state lags only, so ",
           "such a shock cannot be represented.", call. = FALSE)
    if (!is.na(shock_idx[[nm]]))
      stop(context, ": state '", nm, "' is driven by more than one shock ('",
           exo[shock_idx[[nm]]], "' and '", s, "'). solve_global() feeds one ",
           "shock per AR(1) state.", call. = FALSE)

    c0 <- col0[[nm]]
    c1 <- colm1[[nm]]
    if (is.null(c1) || is.na(c1))
      stop(context, ": shock '", s, "' drives '", nm, "' (equation ", q,
           ") but '", nm, "' has no t-1 term there, so it is not an AR(1) ",
           "process. Equation-", q, " terms: ", lbl(nz), ".", call. = FALSE)

    extra <- setdiff(nz, c(c0, c1, cs))
    if (length(extra))
      stop(context, ": equation ", q, " (the process for '", nm,
           "' driven by shock '", s, "') also involves ", lbl(extra),
           ". solve_global() requires the shock-carrying equation to be ",
           "exactly `", nm, " = rho*", nm, "(-1) + psi*", s, "`.",
           call. = FALSE)

    ## Linearity: the three coefficients must not move between probe points.
    co1 <- J1[q, c(c0, c1, cs)]
    co2 <- J2[q, c(c0, c1, cs)]
    if (max(abs(co1 - co2)) > 1e-8 * max(1, max(abs(co1))))
      stop(context, ": the process for '", nm, "' driven by shock '", s,
           "' (equation ", q, ") is NONLINEAR in its own terms (for example ",
           "`log(", nm, ") = rho*log(", nm, "(-1)) + ", s,
           "`). solve_global() injects the shock as ", nm,
           "(-1) + psi*", s, "/rho, which is only correct for a LINEAR AR(1) ",
           "process.", call. = FALSE)

    a0 <- J1[q, c0]; a1 <- J1[q, c1]; b <- J1[q, cs]
    r  <- -a1 / a0
    ps <- -b  / a0
    if (!is.finite(r) || abs(r) < 1e-8)
      stop(context, ": the process for '", nm, "' driven by shock '", s,
           "' has persistence rho = ", format(r), ". The shock is injected ",
           "as ", nm, "(-1) + psi*", s, "/rho, so a zero (or non-finite) rho ",
           "cannot be represented; an i.i.d. state needs a different solver.",
           call. = FALSE)
    if (!is.finite(ps) || abs(ps) < 1e-12)
      stop(context, ": shock '", s, "' has a zero loading on '", nm,
           "' (equation ", q, ").", call. = FALSE)
    if (abs(r) >= 1)
      stop(context, ": the process for '", nm, "' driven by shock '", s,
           "' has |rho| = ", format(abs(r)), " >= 1. The projection grid is ",
           "built from a STATIONARY distribution, so a unit-root state has no ",
           "bounded state domain.", call. = FALSE)

    rho[[nm]]       <- r
    psi[[nm]]       <- ps
    shock_idx[[nm]] <- si
  }

  list(rho = rho, psi = psi, shock_idx = shock_idx)
}


## State (predetermined) variables of a compiled model: the endogenous
## variables that appear at t-1.  Factored out of solve_global() so that
## make_log_posterior_global_pf() can run .global_shock_pairing() at FACTORY
## time -- one derivation, not two that can drift apart.
#' @noRd
.global_state_names <- function(compiled, endo = NULL) {
  endo <- endo %||% compiled$dynamic$endo_names
  lli  <- compiled$lead_lag_incidence %||% compiled$model$lead_lag_incidence
  if (is.null(lli))
    stop("solve_global: the compiled model carries no lead_lag_incidence.",
         call. = FALSE)
  ## LLI rownames are like "t-1", "t", "t+1" -- parse to integers.
  row_labels <- vapply(rownames(lli), function(rn) {
    rn <- trimws(rn)
    if (rn == "t") return(0L)
    m <- regmatches(rn, regexec("^t([+-]?\\d+)$", rn))[[1]]
    if (length(m) == 2L) return(as.integer(m[2L]))
    0L
  }, integer(1L), USE.NAMES = FALSE)
  endo[.structural_lag_lead(lli, row_labels)$has_lag]
}
