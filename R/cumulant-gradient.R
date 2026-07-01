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
#' Returns NULL if any quantity is non-finite.
#' @noRd
.build_moment_vector <- function(dr2, model, params, obs_vars, orders,
                                 me_variance = 0) {
  n_obs   <- length(obs_vars)
  obs_idx <- match(obs_vars, dr2$endo_names)
  if (any(is.na(obs_idx))) return(NULL)

  n_endo <- length(dr2$endo_names)

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
    moments <- compute_moments(dr2, model, params = params)
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
    # Mirror .cumulant_loglik lines 727-738 EXACTLY (including row-then-col subset)
    c3_model_raw <- c3_result$c3_obs[obs_idx, , drop = FALSE]  # n_obs × n_endo^2
    c3_obs_only  <- matrix(0, n_obs, n_obs * n_obs)
    for (a in seq_len(n_obs)) {
      for (b in seq_len(n_obs)) {
        src_col <- (obs_idx[a] - 1L) * n_endo + obs_idx[b]
        dst_col <- (a - 1L) * n_obs + b
        c3_obs_only[a, dst_col] <- c3_model_raw[a, src_col]
      }
    }
    c3_model <- c3_obs_only
  }

  if (any(orders >= 4L) && dr_order >= 2L) {
    c4_result <- tryCatch(
      compute_fourth_cumulant(dr2, model, params),
      error = function(e) NULL
    )
    if (is.null(c4_result)) return(NULL)
    # Mirror .cumulant_loglik line 744 EXACTLY (row subset only, NOT col subset)
    # NB: this means c4_model is n_obs × n_endo^3 (intentionally matches loglik)
    c4_model <- c4_result$c4_obs[obs_idx, , drop = FALSE]
  }

  if (3L %in% orders && !is.null(c3_model)) {
    m_model <- c(m_model, as.numeric(c3_model))
  }

  if (4L %in% orders && !is.null(c4_model)) {
    m_model <- c(m_model, as.numeric(c4_model))
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
#'     parameters.  Orders 1-2 (mean + variance) differentiated via
#'     \code{solution_derivatives_order2()}.  Orders 3-4 (third/fourth cumulant)
#'     differentiated analytically via tensor-Lyapunov sensitivity in
#'     \code{cumulant_moment_derivs_3_4()} — see R/cumulant-cumulant-deriv.R.}
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
    # .cumulant_loglik m_emp: sc$c4 (n_obs × n_obs^3)
    # Note: length mismatch with m_model for order 4 when n_endo > n_obs;
    # we replicate the exact same assembly so that delta = m_emp - m_base
    # uses R's recycling just as .cumulant_loglik does.
    if (!is.null(sc$c4)) m_emp <- c(m_emp, as.numeric(sc$c4))
  }

  n_moments <- length(m_base)
  # Replicate .cumulant_loglik's delta = m_emp - m_model:
  # when lengths differ, R recycles m_emp. Mirror that exactly.
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
# Implicit-differentiation gradient (orders 1-2 analytic; 3-4 FD fallback)
# ============================================================================

#' Implicit-differentiation gradient of the cumulant log-likelihood (orders 1-4).
#'
#' Uses \code{solution_derivatives_order2()} to obtain d(ghxx)/dθ, d(ghss)/dθ,
#' d(Sigma_x)/dθ analytically, then chains through the moment formulas for
#' orders 1 (mean), 2 (variance), 3 (third cumulant, fully analytic via
#' tensor-Lyapunov sensitivity), and 4 (fourth cumulant, fully analytic via
#' the approximate Lyapunov construction in compute_fourth_cumulant).
#' See R/cumulant-cumulant-deriv.R and agent-F-tensor-lyapunov.md.
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

      ## d(Sigma_y)[obs,obs]/dθ = dG[obs,] Σ_x G[obs,]' + G[obs,] dΣ_x G[obs,]' + ... + sym
      dSigma_y_full <- dG %*% Sigma_x %*% t(ghx) + ghx %*% d_Sigma_x %*% t(ghx) +
                       ghx %*% Sigma_x %*% t(dG) +
                       dH %*% Sigma_e %*% t(ghu)  + ghu %*% Sigma_e %*% t(dH)
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
          dm34_c4 <- as.numeric(cd_k$d_c4_obs_sub)
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
