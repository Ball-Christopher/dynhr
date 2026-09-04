## R/cumulant-gradient.R
## --------------------------------------------------------------------------
## Analytic/semi-analytic gradient of the cumulant log-likelihood.
##
## Two routes available via the `deriv` argument:
##
##   "fd"  (default): For each parameter theta_k, re-solve steady-state +
##     order-2 perturbation at theta_k ± h (central FD of the full re-solve).
##
##   "implicit": Implicit differentiation of the order-2 Kronecker/linear
##     systems (Foundation A one order up). The model-Hessian FD is the only
##     per-parameter primitive FD; K_xx and A_L are factorized once and reused
##     for all parameters.  Validated for orders 1-2 (mean + variance).
##     Orders 3-4 fall back to FD-of-resolve (the cumulant-Lyapunov sensitivity
##     was not validated; see agent-E-implicit-order2.md).
##
## Reference: Mutschler (2015), Sections 3.1-3.2.
## --------------------------------------------------------------------------


# ============================================================================
# Internal: re-solve order-2 DR at perturbed params
# ============================================================================

#' Re-solve steady state + order-2 perturbation at perturbed params.
#' The returned DecisionRules2 carries dr$order = 2L intrinsically, which
#' .dr_perturbation_order() / .build_moment_vector() use to activate the
#' order-3/4 terms. We additionally propagate any explicit "order" attribute
#' from `dr_template` so a manual override on the base DR is honoured.
#' Returns NULL on failure (BK not satisfied, non-convergence, etc.).
#' @noRd
.resolve_order2 <- function(model, compiled, params_perturbed, dr_template) {
  ss <- tryCatch(
    solve_steady_state(model, compiled, params_perturbed, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(ss) || !isTRUE(ss$converged)) return(NULL)

  dr2 <- tryCatch(
    solve_perturbation(model, compiled, ss$values, params_perturbed,
                       order = 2L, verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(dr2) || !isTRUE(dr2$bk_satisfied)) return(NULL)

  # Propagate the "order" attribute from the base DR so moment assembly matches
  ord_attr <- attr(dr_template, "order")
  if (!is.null(ord_attr)) attr(dr2, "order") <- ord_attr

  dr2
}


# ============================================================================
# Internal: extract the flat model-moment vector from a dr2 object
# ============================================================================

#' Build the model moment vector m(theta) from a DR object, for given
#' obs_vars and orders.  Replicates the EXACT assembly logic of
#' .cumulant_loglik() — including the attr(dr, "order") guard for orders 3-4
#' — so that the gradient is consistent with the function being differentiated.
#'
#' `lags` (E2-A) appends the model-implied autocovariances Gamma(h) =
#' Cov(y_t, y_{t-h}) for the strictly-positive lags, AFTER the contemporaneous
#' cumulant blocks -- the single definition of the moment ordering shared by
#' `.cumulant_loglik()`, `estimate_gmm_weight_matrix()` and
#' `method_of_moments()`. `lags = integer(0)` (the default) is byte-identical
#' to the pre-E2-A function.
#'
#' Returns NULL if any quantity is non-finite.
#' @noRd
.build_moment_vector <- function(dr2, model, params, obs_vars, orders,
                                 me_variance = 0, lags = integer(0)) {
  n_obs   <- length(obs_vars)
  obs_idx <- match(obs_vars, dr2$endo_names)
  if (any(is.na(obs_idx))) return(NULL)

  n_endo <- length(dr2$endo_names)
  lags   <- .mom_lags(lags)

  ## compute_moments() is the ONE source of Sigma_y / Sigma_state; it solves a
  ## Lyapunov equation, so when BOTH the order-2 block and the autocovariance
  ## block are requested it is evaluated once here and shared.
  mom_cache <- if (2L %in% orders || length(lags))
    compute_moments(dr2, model, params = params) else NULL

  m_model <- numeric(0)

  # ---- Order 1: mean — mirrors .cumulant_loglik lines 706-711 ----
  if (1L %in% orders) {
    mean_model <- dr2$ys[obs_vars]
    if (!is.null(dr2$ghss)) {
      mean_model <- mean_model + 0.5 * dr2$ghss[obs_vars]
    }
    m_model <- c(m_model, mean_model)
  }

  # ---- Order 2: variance — mirrors .cumulant_loglik lines 704-715 ----
  if (2L %in% orders) {
    moments <- mom_cache
    Sigma_y <- moments$var_cov[obs_vars, obs_vars, drop = FALSE]
    if (me_variance > 0) diag(Sigma_y) <- diag(Sigma_y) + me_variance
    m_model <- c(m_model, as.numeric(Sigma_y))
  }

  # ---- Orders 3-4: mirrors .cumulant_loglik ----
  # Use the same robust order detector as .cumulant_loglik so the gradient
  # activates (or skips) orders 3-4 under exactly the same condition. A
  # standard DecisionRules2 from solve_perturbation(order = 2) carries
  # dr$order = 2L, so the higher cumulants are now activated automatically.
  dr_order <- .dr_perturbation_order(dr2)

  c3_model <- NULL
  c4_model <- NULL

  if (any(orders >= 3L) && dr_order >= 2L) {
    c3_result <- tryCatch(
      compute_third_cumulant(dr2, model, params),
      error = function(e) NULL
    )
    if (is.null(c3_result)) return(NULL)
    # Mirror .cumulant_loglik EXACTLY: project rows AND (j,k) columns onto the
    # observables, giving n_obs × n_obs^2 in the sample_cumulants()$c3 layout.
    # (Was a two-level for(a)/for(b) loop keeping only the (i,i,k) slice, so
    # the model side was structurally zero off that slice -- see
    # .project_c3_obs.)
    c3_model <- .project_c3_obs(c3_result$c3_obs, obs_idx, n_endo)
  }

  if (any(orders >= 4L) && dr_order >= 2L) {
    c4_result <- tryCatch(
      compute_fourth_cumulant(dr2, model, params),
      error = function(e) NULL
    )
    if (is.null(c4_result)) return(NULL)
    # Mirror .cumulant_loglik EXACTLY: project rows AND (j,k,l) columns onto
    # the observables, giving n_obs × n_obs^3 in the sample_cumulants()$c4
    # layout.  (Was a row-only subset -> n_obs × n_endo^3, which made the
    # loglik's `m_emp - m_model` recycle when n_obs < n_endo.)
    c4_model <- .project_c4_obs(c4_result$c4_obs, obs_idx, n_endo)
  }

  if (3L %in% orders && !is.null(c3_model)) {
    m_model <- c(m_model, as.numeric(c3_model))
  }

  if (4L %in% orders && !is.null(c4_model)) {
    m_model <- c(m_model, as.numeric(c4_model))
  }

  ## ---- Autocovariances Gamma(h), h >= 1 (E2-A) ----
  ## Measurement error is i.i.d., so it enters Gamma(0) (handled above) and
  ## NOT Gamma(h) for h >= 1 -- do not propagate me_variance here.
  if (length(lags)) {
    ac <- .mom_model_autocov(dr2, model, params, obs_vars, lags,
                             moments = mom_cache)
    if (is.null(ac)) return(NULL)
    for (A in ac) m_model <- c(m_model, as.numeric(A))
  }

  if (!all(is.finite(m_model))) return(NULL)
  m_model
}


# ============================================================================
# Public gradient function
# ============================================================================

#' Analytic/semi-analytic gradient of the cumulant log-likelihood
#'
#' Computes the gradient of \code{.cumulant_loglik} with respect to a set of
#' structural parameters using one of two routes controlled by \code{deriv}:
#'
#' \describe{
#'   \item{\code{"fd"} (default)}{For each parameter theta_k, re-solves steady
#'     state + order-2 perturbation at theta_k ± h and central-differences the
#'     full model-moment vector.  Robust; works for all cumulant orders.}
#'   \item{\code{"implicit"}}{Implicit differentiation of the order-2 Kronecker
#'     and linear systems.  K_xx and A_L are factorized once and reused for all
#'     parameters.  Orders 1-2 (mean + variance) differentiated analytically via
#'     \code{solution_derivatives_order2()}.  Orders 3-4 (third/fourth cumulant)
#'     are differentiated by CENTRAL FINITE DIFFERENCES of the forward cumulant
#'     functions along the analytic order-1/2 solution-derivative direction
#'     (a deliberate choice — a closed-form tensor-Lyapunov sensitivity of the
#'     chain term is not implemented; see the history note in
#'     R/cumulant-cumulant-deriv.R).}
#' }
#'
#' @param model     dynhr_mod from \code{parse_mod()}.
#' @param compiled  dynhr_compiled (max_order >= 2) from
#'   \code{compile_model()}.
#' @param dr        DecisionRules2 (order 2) at the base parameter point.
#' @param params    Named numeric: full parameter vector at the base point.
#' @param param_names Character: names of parameters to differentiate.
#' @param obs_vars  Character: observed variable names.
#' @param data      T × n_obs matrix (columns = obs_vars).
#' @param orders    Integer vector: cumulant orders to match (default 1:4).
#' @param me_variance Scalar measurement-error variance (default 0).
#' @param h_rel     Relative step size for central differences (default 1e-4).
#' @param deriv     Character: \code{"fd"} (default) for FD-of-order-2-solve,
#'   or \code{"implicit"} for implicit-differentiation of the order-2 systems.
#'
#' @return Named numeric vector of length \code{length(param_names)}.
#'   Entries are \code{NA_real_} for parameters where a required solve failed.
#'
#' @seealso \code{.cumulant_loglik}, \code{compute_third_cumulant},
#'   \code{compute_fourth_cumulant}, \code{solution_derivatives_order2}
#' @noRd
cumulant_loglik_grad <- function(model, compiled, dr, params, param_names,
                                 obs_vars, data, orders = 1:4,
                                 me_variance = 0, h_rel = 1e-4,
                                 deriv = c("fd", "implicit")) {
  deriv <- match.arg(deriv)

  ## Dispatch to implicit path for orders 1-2 only; fall back to FD for 3-4.
  if (deriv == "implicit") {
    return(.cumulant_loglik_grad_implicit(model, compiled, dr, params,
                                          param_names, obs_vars, data,
                                          orders = orders,
                                          me_variance = me_variance,
                                          h_rel = h_rel))
  }

  ## ---- FD path (original implementation) ----

  if (!all(param_names %in% names(params)))
    stop("cumulant_loglik_grad: params missing: ",
         paste(param_names[!param_names %in% names(params)], collapse = ", "))

  n_obs   <- length(obs_vars)
  T_obs   <- nrow(data)
  max_ord <- max(orders)

  # ---- 1. Sample cumulants (same as .cumulant_loglik) ----
  data <- .cumulant_subset_obs(data, obs_vars)
  sc <- sample_cumulants(data, max_order = max_ord)

  # ---- 2. Base model-moment vector ----
  m_base <- suppressWarnings(
    .build_moment_vector(dr, model, params, obs_vars, orders, me_variance))
  if (is.null(m_base)) {
    warning("cumulant_loglik_grad: base moment vector is NULL / non-finite")
    return(setNames(rep(NA_real_, length(param_names)), param_names))
  }

  # ---- 3. Build empirical moment vector (mirrors .cumulant_loglik m_emp exactly) ----
  # Replicate .cumulant_loglik lines 752-780
  dr_order_base <- .dr_perturbation_order(dr)
  m_emp <- numeric(0)
  if (1L %in% orders) m_emp <- c(m_emp, sc$mean[obs_vars])
  if (2L %in% orders) {
    m_emp <- c(m_emp, as.numeric(sc$var_cov[obs_vars, obs_vars, drop = FALSE]))
  }
  # For orders 3-4 in m_emp: only append if the model side also appended them
  # (i.e., dr_order_base >= 2). Mirrors .cumulant_loglik exactly.
  c3_active <- 3L %in% orders && dr_order_base >= 2L
  c4_active <- 4L %in% orders && dr_order_base >= 2L

  if (c3_active) {
    # m_base for order 3 includes c3_model (n_obs × n_obs^2 = n_obs*n_obs^2 entries)
    # .cumulant_loglik m_emp: sc_c3_obs = sc$c3 (n_obs × n_obs^2)
    if (!is.null(sc$c3)) m_emp <- c(m_emp, as.numeric(sc$c3))
  }
  if (c4_active) {
    # .cumulant_loglik m_emp: sc$c4 (n_obs × n_obs^3).  Since .build_moment_vector
    # now projects the model c4 onto the observables too, the two blocks have
    # matching length for any n_obs <= n_endo (they used to only line up when
    # every endogenous variable was observed).
    if (!is.null(sc$c4)) m_emp <- c(m_emp, as.numeric(sc$c4))
  }

  n_moments <- length(m_base)
  # Defensive: with the corrected order-4 layout length(m_emp) == n_moments,
  # so rep_len is the identity.  Kept so a future block-length regression
  # degrades the same way .cumulant_loglik does rather than erroring here.
  m_emp_matched <- rep_len(m_emp, n_moments)
  delta_base <- m_emp_matched - m_base     # m_hat - m(theta), length = n_moments

  # ---- 4. For each param: central-difference m(theta) ----
  grad <- setNames(rep(NA_real_, length(param_names)), param_names)

  for (k in seq_along(param_names)) {
    pnm <- param_names[k]
    pval <- params[pnm]

    h <- h_rel * (abs(pval) + 1e-8)   # absolute step

    params_p <- params; params_p[pnm] <- pval + h
    params_m <- params; params_m[pnm] <- pval - h

    dr_p <- .resolve_order2(model, compiled, params_p, dr)
    dr_m <- .resolve_order2(model, compiled, params_m, dr)

    if (is.null(dr_p) || is.null(dr_m)) {
      # Fall back to one-sided if one direction failed
      if (!is.null(dr_p)) {
        m_p <- suppressWarnings(
          .build_moment_vector(dr_p, model, params_p, obs_vars, orders, me_variance))
        m_c <- suppressWarnings(
          .build_moment_vector(dr,   model, params,   obs_vars, orders, me_variance))
        if (!is.null(m_p) && !is.null(m_c)) {
          dm <- (m_p - m_c) / h
          grad[k] <- sum(delta_base * dm) * T_obs / n_moments
        }
      } else if (!is.null(dr_m)) {
        m_m <- suppressWarnings(
          .build_moment_vector(dr_m, model, params_m, obs_vars, orders, me_variance))
        m_c <- suppressWarnings(
          .build_moment_vector(dr,   model, params,   obs_vars, orders, me_variance))
        if (!is.null(m_m) && !is.null(m_c)) {
          dm <- (m_c - m_m) / h
          grad[k] <- sum(delta_base * dm) * T_obs / n_moments
        }
      }
      next
    }

    m_p <- suppressWarnings(
      .build_moment_vector(dr_p, model, params_p, obs_vars, orders, me_variance))
    m_m <- suppressWarnings(
      .build_moment_vector(dr_m, model, params_m, obs_vars, orders, me_variance))

    if (is.null(m_p) || is.null(m_m)) next

    dm <- (m_p - m_m) / (2 * h)   # central difference of moment vector

    # dL/dtheta_k = (m_hat - m)' * dm * T / n_moments
    grad[k] <- sum(delta_base * dm) * T_obs / n_moments
  }

  grad
}


# ============================================================================
# Reverse-mode (adjoint_solution) gradient — Consumer 1, O(1) in P
# ============================================================================

#' Reverse-mode ("adjoint_solution") gradient of the cumulant log-likelihood.
#'
#' Orders 1-3 are differentiated by a single reverse pass: the moment-vector
#' cotangent is back-propagated to cotangents on the solution blocks (ghx, ghu,
#' ghxx, ghxu, ghuu, ghss, Sigma_e) — through the observable projection and the
#' third-cumulant tensor-Lyapunov (via .solve_third_cross_cumulant_adjoint) —
#' then ONE .solution_adjoint (first-order ghx/ghu/ys channel) + ONE
#' .solution_adjoint_order2 (the four order-2 blocks) + a per-parameter
#' d(Sigma_e) contraction turn them into the structural-parameter gradient.
#' This is O(1) in P in the number of factorizations/solution-solves (both
#' adjoint kernels share their factorizations across parameters).
#'
#' Order 4 (kurtosis) is NOT reversed in closed form (its symbolic reverse is a
#' large, error-prone build; see R/cumulant-cumulant-deriv.R STATUS). When
#' \code{4 \%in\% orders}, the order-4 block's gradient contribution is returned
#' as \code{NA} for every parameter, so the caller falls back to exact
#' FD-of-forward for the WHOLE parameter (keeping the gradient consistent with
#' the forward loglik). Orders 1-3 always go through the reverse path.
#'
#' @return Named numeric vector (length param_names). \code{NA} for a parameter
#'   whose reverse contribution could not be formed (e.g. order 4 requested, or
#'   an adjoint kernel returned not-ok).
#' @noRd
.cumulant_loglik_grad_adjoint <- function(model, compiled, dr, params,
                                          param_names, obs_vars, data,
                                          orders = 1:4, me_variance = 0,
                                          h_rel = 1e-4) {

  np <- length(param_names)
  na_out <- setNames(rep(NA_real_, np), param_names)

  if (!inherits(dr, "DecisionRules2")) return(na_out)

  n_obs   <- length(obs_vars)
  T_obs   <- nrow(data)
  endo    <- dr$endo_names
  exo     <- dr$exo_names
  n_endo  <- length(endo)
  n_exo   <- length(exo)
  state_idx <- dr$state_idx
  n_s     <- length(state_idx)
  obs_idx <- match(obs_vars, endo)
  if (any(is.na(obs_idx))) return(na_out)

  ## Order 4: the reverse chain does not cover kurtosis. If it is requested and
  ## actually active (order-2 DR), decline entirely so the caller FD-fallbacks.
  dr_order <- .dr_perturbation_order(dr)
  if (4L %in% orders && dr_order >= 2L) return(na_out)

  ## ---- 1. Sample cumulants + base moment vector + delta -------------------
  max_ord <- max(orders)
  data <- .cumulant_subset_obs(data, obs_vars)
  sc <- sample_cumulants(data, max_order = max_ord)
  m_base <- suppressWarnings(
    .build_moment_vector(dr, model, params, obs_vars, orders, me_variance))
  if (is.null(m_base)) return(na_out)

  m_emp <- numeric(0)
  if (1L %in% orders) m_emp <- c(m_emp, sc$mean[obs_vars])
  if (2L %in% orders)
    m_emp <- c(m_emp, as.numeric(sc$var_cov[obs_vars, obs_vars, drop = FALSE]))
  c3_active <- 3L %in% orders && dr_order >= 2L
  if (c3_active && !is.null(sc$c3)) m_emp <- c(m_emp, as.numeric(sc$c3))
  n_moments <- length(m_base)
  m_emp_matched <- rep_len(m_emp, n_moments)
  delta_base <- m_emp_matched - m_base

  ## dL/dm_model = delta_base * T / n_moments  (loglik = -0.5 sum(delta^2)/n * T)
  bar_m <- delta_base * T_obs / n_moments

  ## ---- 2. Split bar_m into per-order blocks -------------------------------
  ghx <- dr$ghx; ghu <- dr$ghu; ghss <- dr$ghss
  Sigma_e <- diag(.get_shock_stderr(model, exo, params)^2, n_exo)
  hx <- ghx[state_idx, , drop = FALSE]
  hu <- ghu[state_idx, , drop = FALSE]
  Sigma_state <- .state_covariance(hx, hu, Sigma_e)   # n_s x n_s

  ## Block cotangents on the solution matrices (full endo-row layout).
  bar_ghx  <- matrix(0, n_endo, ncol(ghx))
  bar_ghu  <- matrix(0, n_endo, n_exo)
  bar_ghss <- numeric(n_endo)
  bar_ys   <- setNames(numeric(n_endo), endo)
  bar_ghxx <- matrix(0, n_endo, ncol(dr$ghxx))
  bar_ghxu <- if (!is.null(dr$ghxu)) matrix(0, n_endo, ncol(dr$ghxu)) else NULL
  bar_ghuu <- if (!is.null(dr$ghuu)) matrix(0, n_endo, ncol(dr$ghuu)) else NULL
  bar_Sigma_state <- matrix(0, n_s, n_s)
  bar_Sigma_e     <- matrix(0, n_exo, n_exo)

  idx <- 0L

  ## ---- Order 1 (mean = ys[obs] + 0.5 ghss[obs]) ---------------------------
  if (1L %in% orders) {
    bm <- bar_m[seq_len(n_obs) + idx]
    bar_ys[obs_vars] <- bar_ys[obs_vars] + bm
    if (!is.null(ghss)) bar_ghss[obs_idx] <- bar_ghss[obs_idx] + 0.5 * bm
    idx <- idx + n_obs
  }

  ## ---- Order 2 (variance Sigma_y[obs,obs]) --------------------------------
  ## Sigma_y = ghx Sigma_state ghx' + ghu Sigma_e ghu'  (me_variance adds to diag,
  ## which is param-free so contributes nothing to the block cotangents).
  if (2L %in% orders) {
    bS_obs <- matrix(bar_m[seq_len(n_obs^2) + idx], n_obs, n_obs)
    idx <- idx + n_obs^2
    ## scatter obs-block cotangent into full n_endo x n_endo
    bS <- matrix(0, n_endo, n_endo)
    bS[obs_idx, obs_idx] <- bS_obs
    ## reverse ghx Sigma_state ghx'
    bar_ghx <- bar_ghx +
      bS %*% ghx %*% t(Sigma_state) + t(bS) %*% ghx %*% Sigma_state
    bar_Sigma_state <- bar_Sigma_state + t(ghx) %*% bS %*% ghx
    ## reverse ghu Sigma_e ghu'
    bar_ghu <- bar_ghu +
      bS %*% ghu %*% t(Sigma_e) + t(bS) %*% ghu %*% Sigma_e
    bar_Sigma_e <- bar_Sigma_e + t(ghu) %*% bS %*% ghu
  }

  ## ---- Order 3 (third cumulant) -------------------------------------------
  if (c3_active) {
    ## bar_m order-3 block is n_obs x n_obs^2 (dst-col layout of .cumulant_loglik).
    n34 <- n_obs * n_obs * n_obs
    b3_obs <- bar_m[seq_len(n34) + idx]
    idx <- idx + n34
    ## Scatter into full c3_obs (n_endo x n_endo^2).  The forward gather is
    ##   c3_model <- c3_obs[obs_idx, .c3_obs_col_index(obs_idx, n_endo)]
    ## (a pure sub-selection with no repeated source entry), so its reverse is
    ## the transposed scatter onto exactly those rows/columns.  It used to be a
    ## two-level for(a)/for(b) loop hitting only the (i,i,k) slice, which left
    ## the (i,j,k), j != i, cotangents at zero -- the adjoint mirror of the
    ## forward E4-B bug.
    bar_c3_obs <- matrix(0, n_endo, n_endo * n_endo)
    b3_full <- matrix(b3_obs, n_obs, n_obs * n_obs)
    bar_c3_obs[obs_idx, .c3_obs_col_index(obs_idx, n_endo)] <- b3_full
    rev3 <- tryCatch(
      .compute_third_cumulant_adjoint(dr, model, params, bar_c3_obs),
      error = function(e) NULL)
    if (is.null(rev3)) return(na_out)
    bar_ghx  <- bar_ghx  + rev3$bar_ghx
    bar_ghu  <- bar_ghu  + rev3$bar_ghu
    bar_ghxx <- bar_ghxx + rev3$bar_ghxx
    if (!is.null(bar_ghxu) && !is.null(rev3$bar_ghxu))
      bar_ghxu <- bar_ghxu + rev3$bar_ghxu
    if (!is.null(bar_ghuu) && !is.null(rev3$bar_ghuu))
      bar_ghuu <- bar_ghuu + rev3$bar_ghuu
    bar_Sigma_state <- bar_Sigma_state + rev3$bar_Sigma_x
    bar_Sigma_e     <- bar_Sigma_e     + rev3$bar_Sigma_e
  }

  ## ---- 3. Fold bar_Sigma_state through the Lyapunov solve -----------------
  ## Sigma_state = A Sigma_state A' + Q, A = hx, Q = hu Sigma_e hu'.
  if (n_s > 0 && any(bar_Sigma_state != 0)) {
    lyr <- .lyap_solve_adjoint(hx, Sigma_state, bar_Sigma_state)
    ## NULL => the transposed Lyapunov had no stationary solution (explosive or
    ## unit-root hx). `.solve_lyapunov` signals that with NaN, not a condition,
    ## so without this test a non-stationary draw silently yields a NaN
    ## gradient. Mirror the adjoint-KF siblings: return the all-NA gradient.
    if (is.null(lyr)) return(na_out)
    bar_hx_ly <- lyr$bar_A                       # n_s x n_s
    bar_Q     <- lyr$bar_Q                       # n_s x n_s
    ## Q = hu Sigma_e hu'
    bar_hu_ly <- bar_Q %*% hu %*% t(Sigma_e) + t(bar_Q) %*% hu %*% Sigma_e
    bar_Sigma_e <- bar_Sigma_e + t(hu) %*% bar_Q %*% hu
    ## fold hx/hu cotangents into full ghx/ghu state rows/cols
    bar_ghx[state_idx, seq_len(n_s)] <-
      bar_ghx[state_idx, seq_len(n_s)] + bar_hx_ly
    bar_ghu[state_idx, ] <- bar_ghu[state_idx, , drop = FALSE] + bar_hu_ly
  }

  ## ---- 4. First-order channel: ghx/ghu/ys -> structural params ------------
  ## .solution_adjoint reverses the whole first-order fixed point + ys. Its
  ## block-extraction reverse is ADDITIVE: V_G[state_idx,] += G_TT and
  ## V_G[obs_idx,] += G_ZZ. Using obs_vars = ALL endo makes G_ZZ/G_DD cover
  ## every row, so we route the COMPLETE bar_ghx/bar_ghu through G_ZZ/G_DD and
  ## set G_TT/G_RR = 0 to avoid double-counting the state rows.
  fo <- tryCatch(
    .solution_adjoint(
      model, compiled, dr, params, param_names,
      obs_vars = endo,
      bars = list(
        G_TT = matrix(0, n_s, n_s),
        G_RR = matrix(0, n_s, n_exo),
        G_ZZ = bar_ghx[, seq_len(n_s), drop = FALSE],
        G_DD = bar_ghu,
        g_d  = as.numeric(bar_ys))),
    error = function(e) NULL)

  ## ---- 5. Order-2 channel: ghxx/ghxu/ghuu/ghss -> structural params -------
  o2 <- tryCatch(
    .solution_adjoint_order2(
      model, compiled, dr, params, param_names,
      bars = list(bar_ghxx = bar_ghxx, bar_ghxu = bar_ghxu,
                  bar_ghuu = bar_ghuu, bar_ghss = bar_ghss)),
    error = function(e) NULL)

  ## ---- 6. Sigma_e direct channel (per-parameter d(Sigma_e)) ---------------
  ## The cumulant moments depend on Sigma_e directly (not only through the
  ## solution). The order-2 adjoint accounts for the Sigma_e that flows THROUGH
  ## the solution (its own dvSe); the DIRECT dependence (bar_Sigma_e above) is
  ## contracted here with each parameter's central-FD d(Sigma_e) — one cheap 2x
  ## FD of .get_shock_cov per parameter, no solution solve.
  grad <- setNames(rep(NA_real_, np), param_names)
  bar_vSe <- as.numeric(bar_Sigma_e)
  for (pnm in param_names) {
    g_fo <- if (!is.null(fo) && isTRUE(fo$ok[[pnm]])) fo$grad[[pnm]] else NA_real_
    g_o2 <- if (!is.null(o2) && isTRUE(o2$ok[[pnm]])) o2$grad[[pnm]] else NA_real_
    if (is.na(g_fo) || is.na(g_o2)) { grad[pnm] <- NA_real_; next }
    gj <- g_fo + g_o2
    if (any(bar_vSe != 0)) {
      pval <- params[[pnm]]
      h <- max(h_rel * abs(pval), 1e-7)
      dSe <- .o2sd_dSigma_e(model, params, pnm, h, n_exo)
      gj <- gj + sum(bar_vSe * as.numeric(dSe))
    }
    grad[pnm] <- gj
  }
  grad
}


# ============================================================================
# Implicit-differentiation gradient (orders 1-2 analytic; 3-4 FD fallback)
# ============================================================================

#' Implicit-differentiation gradient of the cumulant log-likelihood (orders 1-4).
#'
#' Uses \code{solution_derivatives_order2()} to obtain d(ghxx)/dθ, d(ghss)/dθ,
#' d(Sigma_x)/dθ analytically, then chains through the moment formulas for
#' orders 1 (mean) and 2 (variance) analytically; orders 3 (third cumulant) and
#' 4 (fourth cumulant) are obtained by CENTRAL FINITE DIFFERENCES of the forward
#' cumulant functions along that analytic solution-derivative direction (a
#' closed-form tensor-Lyapunov sensitivity of the chain term is deliberately not
#' implemented). See R/cumulant-cumulant-deriv.R.
#' @noRd
.cumulant_loglik_grad_implicit <- function(model, compiled, dr, params,
                                            param_names, obs_vars, data,
                                            orders = 1:4,
                                            me_variance = 0,
                                            h_rel = 1e-4) {

  if (!inherits(dr, "DecisionRules2"))
    stop(".cumulant_loglik_grad_implicit: dr must be a DecisionRules2 object")

  n_obs   <- length(obs_vars)
  T_obs   <- nrow(data)
  np      <- length(param_names)

  ## ---- 1. Sample cumulants ----
  max_ord <- max(orders)
  data <- .cumulant_subset_obs(data, obs_vars)
  sc <- sample_cumulants(data, max_order = max_ord)

  ## ---- 2. Base model-moment vector ----
  m_base <- suppressWarnings(
    .build_moment_vector(dr, model, params, obs_vars, orders, me_variance))
  if (is.null(m_base)) {
    warning(".cumulant_loglik_grad_implicit: base moment vector NULL/non-finite")
    return(setNames(rep(NA_real_, np), param_names))
  }

  ## ---- 3. Empirical moment vector (same assembly as FD path) ----
  dr_order_base <- .dr_perturbation_order(dr)
  m_emp <- numeric(0)
  if (1L %in% orders) m_emp <- c(m_emp, sc$mean[obs_vars])
  if (2L %in% orders) {
    m_emp <- c(m_emp, as.numeric(sc$var_cov[obs_vars, obs_vars, drop = FALSE]))
  }
  c3_active <- 3L %in% orders && dr_order_base >= 2L
  c4_active <- 4L %in% orders && dr_order_base >= 2L
  if (c3_active && !is.null(sc$c3)) m_emp <- c(m_emp, as.numeric(sc$c3))
  if (c4_active && !is.null(sc$c4)) m_emp <- c(m_emp, as.numeric(sc$c4))

  n_moments <- length(m_base)
  m_emp_matched <- rep_len(m_emp, n_moments)
  delta_base <- m_emp_matched - m_base

  ## ---- 4. Determine which orders are present ----
  orders_12 <- intersect(orders, 1:2)      # orders 1-2: analytic via solution_deriv_order2
  orders_34 <- intersect(orders, 3:4)      # orders 3-4: fully analytic via cumulant-Lyapunov

  has_12 <- length(orders_12) > 0L
  has_34 <- length(orders_34) > 0L

  ## ---- 5. Implicit derivatives of the order-2 solution ----
  ## Always call solution_derivatives_order2 — orders 3-4 need d_ghxx / d_Sigma_x
  ## even if orders 1-2 are not requested.
  o2d <- tryCatch(
    solution_derivatives_order2(model, compiled, dr, params, param_names,
                                 h_rel = h_rel, h_hess = max(h_rel * 100, 1e-4)),
    error = function(e) {
      warning(sprintf(".cumulant_loglik_grad_implicit: solution_derivatives_order2 failed: %s",
                      conditionMessage(e)))
      NULL
    }
  )
  if (is.null(o2d)) {
    warning(".cumulant_loglik_grad_implicit: falling back fully to FD route")
    return(cumulant_loglik_grad(model, compiled, dr, params, param_names,
                                 obs_vars, data, orders = orders,
                                 me_variance = me_variance, h_rel = h_rel,
                                 deriv = "fd"))
  }

  ## ---- 6. Dimensions ----
  endo_names <- dr$endo_names
  state_idx  <- dr$state_idx
  exo_names  <- dr$exo_names
  n          <- length(endo_names)
  n_s        <- length(state_idx)
  n_u        <- length(exo_names)
  obs_idx    <- match(obs_vars, endo_names)

  ghx     <- dr$ghx
  ghu     <- dr$ghu
  Sigma_e <- dr$Sigma_e
  hx      <- ghx[state_idx, , drop = FALSE]
  hu      <- ghu[state_idx, , drop = FALSE]

  ## Base Sigma_x and Sigma_y (for order-2 derivative)
  Sigma_x <- .state_covariance(hx, hu, Sigma_e)   # n_s x n_s
  Sigma_y_full <- ghx %*% Sigma_x %*% t(ghx) + ghu %*% Sigma_e %*% t(ghu)  # n x n
  if (me_variance > 0) {
    Sigma_y_full_obs <- Sigma_y_full[obs_idx, obs_idx, drop = FALSE]
    diag(Sigma_y_full_obs) <- diag(Sigma_y_full_obs) + me_variance
  } else {
    Sigma_y_full_obs <- Sigma_y_full[obs_idx, obs_idx, drop = FALSE]
  }

  ## ---- 6b. Analytic order-3/4 cumulant derivatives ----
  ## Pre-compute base cumulants and call cumulant_moment_derivs_3_4 once for
  ## all parameters (all parameter loops share the same base quantities).
  cd34 <- NULL   # will hold list(pnm = list(d_c3_obs_sub, d_c4_obs_sub, ok))
  if (has_34 && c3_active) {
    c3_base <- tryCatch(
      compute_third_cumulant(dr, model, params),
      error = function(e) NULL
    )
    c4_base <- tryCatch(
      compute_fourth_cumulant(dr, model, params),
      error = function(e) NULL
    )
    if (!is.null(c3_base) && !is.null(c4_base)) {
      cd34 <- tryCatch(
        cumulant_moment_derivs_3_4(
          dr2          = dr,
          model        = model,
          params       = params,
          o2d          = o2d,
          obs_idx      = obs_idx,
          obs_vars     = obs_vars,
          Sigma_x      = Sigma_x,
          Sigma_y_full = Sigma_y_full,
          c3_base      = c3_base,
          c4_base      = c4_base,
          compiled     = compiled,
          h_rel        = h_rel
        ),
        error = function(e) {
          warning(sprintf(".cumulant_loglik_grad_implicit: cumulant_moment_derivs_3_4 failed: %s",
                          conditionMessage(e)))
          NULL
        }
      )
    }
  }

  ## ---- 7. Per-parameter gradient ----
  grad <- setNames(rep(NA_real_, np), param_names)

  for (k in seq_along(param_names)) {
    pnm  <- param_names[k]
    pval <- params[[pnm]]
    h    <- h_rel * (abs(pval) + 1e-8)

    dm <- numeric(n_moments)
    idx <- 0L

    ## ---- Orders 1-2: analytic ----
    if (has_12 && !is.null(o2d)) {
      d2 <- o2d$derivs[[pnm]]
      d1 <- o2d$first$derivs[[pnm]]

      if (!isTRUE(d2$ok) || !isTRUE(d1$ok)) {
        ## Fail gracefully: NA for this param
        next
      }

      dG   <- d1$dG     # n x n_s
      dH   <- d1$dH     # n x n_u
      dys  <- d1$dys    # n
      dhx  <- dG[state_idx, , drop = FALSE]   # n_s x n_s
      dhu  <- dH[state_idx, , drop = FALSE]   # n_s x n_u

      d_Sigma_x <- d2$d_Sigma_x   # n_s x n_s
      d_ghss    <- d2$d_ghss      # n
      ## Sigma_e's OWN theta-dependence (estimated shock stds wired via
      ## stderr_expr): without the ghu*dSigma_e*ghu' term below,
      ## d(Sigma_y)/dtheta silently drops the shock-covariance channel (C5
      ## sibling of the .o2sd_dSigma_x gap). Zero for params that do not
      ## enter Sigma_e.
      d_Sigma_e <- d2$d_Sigma_e
      if (is.null(d_Sigma_e)) d_Sigma_e <- matrix(0, nrow(Sigma_e), ncol(Sigma_e))

      ## d(Sigma_y)[obs,obs]/dθ = dG[obs,] Σ_x G[obs,]' + G[obs,] dΣ_x G[obs,]' + ... + sym
      dSigma_y_full <- dG %*% Sigma_x %*% t(ghx) + ghx %*% d_Sigma_x %*% t(ghx) +
                       ghx %*% Sigma_x %*% t(dG) +
                       dH %*% Sigma_e %*% t(ghu)  + ghu %*% Sigma_e %*% t(dH) +
                       ghu %*% d_Sigma_e %*% t(ghu)
      dSigma_y_obs <- dSigma_y_full[obs_idx, obs_idx, drop = FALSE]

      if (1L %in% orders_12) {
        ## d(mean)/dθ = dys[obs] + 0.5 * d_ghss[obs]
        d_mean <- dys[obs_vars]
        if (!is.null(dr$ghss)) d_mean <- d_mean + 0.5 * d_ghss[obs_vars]
        dm[seq_len(n_obs) + idx] <- d_mean
        idx <- idx + n_obs
      }
      if (2L %in% orders_12) {
        dm[seq_len(n_obs^2) + idx] <- as.numeric(dSigma_y_obs)
        idx <- idx + n_obs^2
      }
    }

    ## ---- Orders 3-4: fully analytic via tensor-Lyapunov sensitivity ----
    if (has_34 && c3_active) {
      if (!is.null(cd34) && !is.null(cd34[[pnm]]) && isTRUE(cd34[[pnm]]$ok)) {
        cd_k <- cd34[[pnm]]

        if (3L %in% orders_34) {
          dm34_c3 <- as.numeric(cd_k$d_c3_obs_sub)
          n34_c3  <- length(dm34_c3)
          dm[seq_len(n34_c3) + idx] <- dm34_c3
          idx <- idx + n34_c3
        }

        if (4L %in% orders_34) {
          ## cumulant_moment_derivs_3_4() returns d_c4_obs_sub row-subset only
          ## (n_obs x n_endo^3).  The moment block is n_obs x n_obs^3, so
          ## project the (j,k,l) columns here with the SAME index map
          ## .build_moment_vector / .cumulant_loglik use for the forward c4.
          dm34_c4 <- as.numeric(
            cd_k$d_c4_obs_sub[, .c4_obs_col_index(obs_idx, n), drop = FALSE])
          n34_c4  <- length(dm34_c4)
          dm[seq_len(n34_c4) + idx] <- dm34_c4
          idx <- idx + n34_c4
        }
      }
      ## If cd34 failed for this param, dm stays 0 for orders 3-4 (safe)
    }

    grad[k] <- sum(delta_base * dm) * T_obs / n_moments
  }

  grad
}
