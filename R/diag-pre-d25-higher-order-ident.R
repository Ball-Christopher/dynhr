## R/diag-pre-d25-higher-order-ident.R
## --------------------------------------------------------------------------
## D25: Higher-order (pruned second-order) local identification.
##
## Compares the rank of the Jacobian d m(theta) / d theta' of the observables'
## unconditional moments (mean, covariance, autocovariances) computed from
##   (i)  the first-order solution, and
##   (ii) the pruned second-order solution (AFVRR 2018 augmented state space),
## both RE-SOLVED at each perturbed theta. Parameters that only move the
## second-order terms (risk / curvature corrections) are invisible to (i) but
## identified by (ii); parameters entering only through a product are
## unidentified in both.
##
## Key references:
##   Mutschler, W. (2015). Identification of DSGE models -- the effect of
##     higher-order approximation and pruning. Journal of Economic Dynamics
##     and Control, 56, 34-54.
##   Iskrev, N. (2010). Local identification in DSGE models.
##     Journal of Monetary Economics, 57(2), 189-202.
##   Andreasen, M. M., Fernandez-Villaverde, J., & Rubio-Ramirez, J. F.
##     (2018). The pruned state-space system for non-linear DSGE models.
##     Review of Economic Studies, 85(1), 1-49.
## --------------------------------------------------------------------------

#' D25. Higher-order (pruned order-2) identification diagnostic
#'
#' Local identification check in the spirit of Mutschler (2015): the
#' Jacobian of the observables' unconditional moments with respect to the
#' parameters is computed twice from the same re-solved model -- once from
#' the first-order solution and once from the pruned second-order solution
#' (exact AFVRR 2018 moments) -- and the two ranks are compared with the
#' shared equilibrated, finite-difference-noise-aware rank rule
#' (\code{.ident_equilibrated_rank}). The moment vector at each order is the
#' observables' mean (optional), \code{vech} of their covariance and the full
#' autocovariance matrices at lags \code{1..n_lags}.
#'
#' The diagnostic needs a re-solve closure: with a fixed decision rule the
#' moments do not move with theta and every Jacobian is zero. Supply
#' \code{dr_solve_fn}, or \code{model} + \code{compiled} (the closure is then
#' built internally; a \code{compiled} object of order < 2 is recompiled at
#' \code{max_order = 2}).
#'
#' @param dr_solve_fn Function \code{theta -> DecisionRules2} (or
#'   \code{NULL} where the model does not solve). The returned object should
#'   carry \code{$Sigma_e} evaluated at theta so estimated shock standard
#'   errors are differentiated; otherwise \code{Sigma_e} is used (fixed).
#' @param params     Named numeric vector: the parameters to check (the
#'   point at which the Jacobian is evaluated).
#' @param model,compiled Parsed and compiled model, used to build
#'   \code{dr_solve_fn} when it is \code{NULL}.
#' @param obs_names  Observed endogenous variable names (default: all
#'   endogenous variables). Ignored when \code{obs_mat} is given.
#' @param obs_mat    Optional \code{n_obs x n_endo} observation matrix.
#' @param Sigma_e    Shock covariance used when a solved rule carries no
#'   \code{$Sigma_e}. There is no identity default.
#' @param param_names Optional parameter labels (default \code{names(params)}).
#' @param n_lags     Number of autocovariance lags in the moment vector.
#' @param include_mean Include the observables' unconditional mean (set
#'   \code{FALSE} for demeaned data).
#' @param ramsey_result Optional \code{dynhr_ramsey_result2}. A Ramsey
#'   solution cannot be re-solved at perturbed theta here, so without an
#'   explicit \code{dr_solve_fn} the diagnostic returns INFO.
#' @param eps        Central finite-difference step (the rank tolerance also
#'   uses the Jacobian at \code{2 * eps}).
#' @param verbose    Print progress messages.
#' @param meta       Plot provenance descriptor.
#'
#' @return A \code{dynhr_diagnostic}; \code{pass} is \code{TRUE} when the
#'   parameters are locally identified by the pruned second-order moments,
#'   \code{FALSE} when they are not, \code{NA} when the check could not run.
#'   \code{result} holds \code{order1_jacobian}, \code{order2_jacobian},
#'   \code{order1_rank}, \code{order2_rank}, \code{identification_gain},
#'   \code{order1_singular_values}, \code{order2_singular_values},
#'   \code{order1_tol}, \code{order2_tol}, \code{param_identification}
#'   (per-parameter status and distinctness at each order),
#'   \code{higher_order_only} (parameters identified only at order 2),
#'   \code{unidentified_order2}, \code{moments}, \code{moment_vector_info}.
#'
#' @references
#'   Mutschler, W. (2015). Identification of DSGE models -- the effect of
#'     higher-order approximation and pruning. \emph{Journal of Economic
#'     Dynamics and Control}, 56, 34-54.
#'
#' @noRd
d25_higher_order_identification <- function(dr_solve_fn = NULL,
                                            params = NULL,
                                            model = NULL,
                                            compiled = NULL,
                                            obs_names = NULL,
                                            obs_mat = NULL,
                                            Sigma_e = NULL,
                                            param_names = NULL,
                                            n_lags = 4L,
                                            include_mean = TRUE,
                                            ramsey_result = NULL,
                                            eps = 1e-5,
                                            verbose = FALSE,
                                            meta = NULL) {
  skip <- function(reason, text) {
    .make_result(
      pass = NA,
      summary = paste0("D25 Higher-order identification: ", text),
      llm_summary = sprintf(
        "D25 | Higher-Order Identification | INFO\n  status=skipped reason=%s\n  note: %s",
        reason, text))
  }

  if (is.null(params) || length(params) == 0L)
    return(skip("no_params", "no parameters supplied."))
  if (is.null(dr_solve_fn) && !is.null(ramsey_result))
    return(skip("ramsey_no_resolve", paste0(
      "a Ramsey solution cannot be re-solved at perturbed parameters here ",
      "(supply dr_solve_fn); the rank was not checked.")))
  if (is.null(dr_solve_fn)) {
    if (is.null(model) || is.null(compiled))
      return(skip("no_resolve_closure", paste0(
        "needs a re-solve closure (dr_solve_fn, or model + compiled); with a ",
        "fixed decision rule the moment Jacobian is identically zero.")))
    dr_solve_fn <- .d25_dr2_solve_fn(model, compiled)
  }
  if (!is.function(dr_solve_fn))
    .dynhr_abort("`dr_solve_fn` must be a function theta -> DecisionRules2.")

  if (is.null(param_names)) param_names <- names(params)
  if (is.null(param_names)) param_names <- paste0("theta_", seq_along(params))
  if (length(param_names) != length(params))
    .dynhr_abort(sprintf("`param_names` has %d entries but `params` has %d.",
                         length(param_names), length(params)))
  n_par <- length(params)
  n_lags <- as.integer(n_lags)

  # ---- 1. Solve at the evaluation point and fix the observation map ----
  dr0 <- dr_solve_fn(params)
  if (is.null(dr0))
    return(skip("no_solution", "the model does not solve at params."))
  if (is.null(dr0$ghxx) || is.null(dr0$ghuu) || is.null(dr0$ghxu) ||
      is.null(dr0$ghss))
    .dynhr_abort("`dr_solve_fn` must return a second-order (DecisionRules2) solution.")
  endo_names <- dr0$endo_names %||% rownames(dr0$ghx)
  n_endo <- nrow(dr0$ghx)
  if (is.null(obs_mat)) {
    if (is.null(obs_names)) obs_names <- endo_names
    idx <- match(obs_names, endo_names)
    if (anyNA(idx))
      .dynhr_abort(sprintf("observables not among the endogenous variables: %s",
                           paste(obs_names[is.na(idx)], collapse = ", ")))
    obs_mat <- diag(n_endo)[idx, , drop = FALSE]
  } else {
    obs_mat <- as.matrix(obs_mat)
    if (ncol(obs_mat) != n_endo)
      .dynhr_abort(sprintf("`obs_mat` has %d columns; the model has %d endogenous variables.",
                           ncol(obs_mat), n_endo))
    obs_names <- rownames(obs_mat) %||% paste0("obs_", seq_len(nrow(obs_mat)))
  }
  if (is.null(dr0$Sigma_e) && is.null(Sigma_e))
    .dynhr_abort(paste0("no shock covariance: the solved rule carries no ",
                        "$Sigma_e and `Sigma_e` is NULL (no identity default)."))

  # Both moment vectors come from ONE solve per theta (a DecisionRules2 holds
  # the first-order solution in ghx/ghu).
  moment_fn <- function(th) {
    names(th) <- names(params)
    d <- if (identical(unname(th), unname(params))) dr0 else dr_solve_fn(th)
    if (is.null(d)) return(NULL)
    Se <- as.matrix(d$Sigma_e %||% Sigma_e)
    m1 <- .d25_moments_order1(d, Se, obs_mat, n_lags, include_mean)
    m2 <- .d25_moments_order2(d, Se, obs_mat, n_lags, include_mean)
    c(m1, m2)
  }
  m0 <- moment_fn(params)
  n_m <- .d25_moment_count(nrow(obs_mat), n_lags, include_mean)
  if (length(m0) != 2L * n_m || !all(is.finite(m0)))
    return(skip("nonfinite_moments", paste0(
      "the unconditional moments are not finite at params (non-stationary ",
      "or explosive pruned system); the rank was not checked.")))
  mnames <- .d25_moment_names(obs_names, n_lags, include_mean)

  # ---- 2. Jacobians at h and 2h ----
  if (verbose) .dynhr_inform(sprintf("[d25] %d parameters, %d re-solves.",
                                     n_par, 4L * n_par))
  J  <- .numerical_jacobian(moment_fn, params, eps = eps)
  J2 <- .numerical_jacobian(moment_fn, params, eps = 2 * eps)
  if (!all(is.finite(J)) || !all(is.finite(J2))) {
    bad <- param_names[colSums(!is.finite(J)) > 0 | colSums(!is.finite(J2)) > 0]
    return(skip("nonfinite_jacobian", sprintf(
      "the model does not solve (or is non-stationary) near params for: %s; the rank was not checked.",
      paste(bad, collapse = ", "))))
  }
  dimnames(J) <- dimnames(J2) <- list(c(paste0("o1:", mnames), paste0("o2:", mnames)),
                                      param_names)
  r1 <- seq_len(n_m); r2 <- n_m + seq_len(n_m)
  J_o1 <- J[r1, , drop = FALSE]; J_o2 <- J[r2, , drop = FALSE]
  rk1 <- .ident_equilibrated_rank(J_o1, J2[r1, , drop = FALSE])
  rk2 <- .ident_equilibrated_rank(J_o2, J2[r2, , drop = FALSE])
  rank_o1 <- as.integer(rk1$rank); rank_o2 <- as.integer(rk2$rank)
  gain <- rank_o2 - rank_o1

  # ---- 3. Per-parameter status ----
  status <- function(rk) ifelse(param_names %in% rk$unidentified_params,
                                "Unidentified", "Identified")
  param_id <- data.frame(
    parameter        = param_names,
    order1_status    = status(rk1),
    order2_status    = status(rk2),
    order1_distinct  = .d25_distinctness(rk1$Je, n_par),
    order2_distinct  = .d25_distinctness(rk2$Je, n_par),
    stringsAsFactors = FALSE
  )
  if (rank_o1 == n_par) param_id$order1_status <- "Identified"
  if (rank_o2 == n_par) param_id$order2_status <- "Identified"
  ho_only <- param_id$parameter[param_id$order1_status == "Unidentified" &
                                param_id$order2_status == "Identified"]
  unid2 <- param_id$parameter[param_id$order2_status == "Unidentified"]

  pass <- rank_o2 == n_par

  # ---- 4. Plots ----
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plots$sv_comparison <- .apply_meta(
      .d25_plot_sv(rk1, rk2, rank_o1, rank_o2, n_par), meta)
    plots$param_identification <- .apply_meta(
      .d25_plot_params(param_id), meta)
  }

  # ---- 5. Text ----
  fmt_set <- function(x) if (length(x)) paste(x, collapse = ", ") else "none"
  verdict <- if (pass && gain > 0) {
    sprintf("Identified at pruned order 2, but NOT at order 1 (higher-order-only: %s).",
            fmt_set(ho_only))
  } else if (pass) {
    "Identified at both orders."
  } else if (gain > 0) {
    sprintf("Order 2 adds %d direction(s) (higher-order-only: %s) but %d direction(s) remain unidentified (%s).",
            gain, fmt_set(ho_only), n_par - rank_o2, fmt_set(unid2))
  } else {
    sprintf("Not identified even at pruned order 2 (%d direction(s) missing: %s).",
            n_par - rank_o2, fmt_set(unid2))
  }
  summary_text <- sprintf(
    "D25 Higher-order identification: rank order-1 = %d/%d, pruned order-2 = %d/%d (gain %+d). %s",
    rank_o1, n_par, rank_o2, n_par, gain, verdict)

  llm_summary <- paste(c(
    sprintf("D25 | Higher-Order Identification | %s", if (pass) "PASS" else "FAIL"),
    sprintf("  rank_o1=%d rank_o2=%d gain=%+d n_par=%d n_moments_per_order=%d n_obs=%d n_lags=%d mean=%s",
            rank_o1, rank_o2, gain, n_par, n_m, nrow(obs_mat), n_lags, include_mean),
    sprintf("  tol_o1=%.3g (%s) tol_o2=%.3g (%s)", rk1$tol, rk1$tol_source,
            rk2$tol, rk2$tol_source),
    sprintf("  higher_order_only: %s", fmt_set(ho_only)),
    sprintf("  unidentified_order2: %s", fmt_set(unid2)),
    sprintf("  action: %s",
            if (pass && gain > 0)
              "Estimate at order >= 2 (pruned); a first-order likelihood cannot pin down the higher-order-only parameters."
            else if (pass)
              "None."
            else
              "Fix or reparameterise the unidentified combination, or add observables; higher-order moments do not resolve it.")
  ), collapse = "\n")

  .make_result(
    result = list(
      order1_jacobian = J_o1,
      order2_jacobian = J_o2,
      order1_rank = rank_o1,
      order2_rank = rank_o2,
      identification_gain = gain,
      order1_singular_values = rk1$singular_values,
      order2_singular_values = rk2$singular_values,
      order1_tol = rk1$tol,
      order2_tol = rk2$tol,
      param_identification = param_id,
      higher_order_only = ho_only,
      unidentified_order2 = unid2,
      order2_null_space = rk2$null_space,
      moments = list(order1 = stats::setNames(m0[r1], mnames),
                     order2 = stats::setNames(m0[r2], mnames)),
      moment_vector_info = list(n_obs = nrow(obs_mat), obs_names = obs_names,
                                n_lags = n_lags, include_mean = include_mean,
                                n_moments_per_order = n_m)
    ),
    pass = pass,
    plots = plots,
    summary = summary_text,
    llm_summary = llm_summary
  )
}


## theta -> DecisionRules2 closure for a parsed + compiled model. Re-solves the
## steady state, the first-order solution (BK + stationarity guard) and the
## second-order solution, and attaches Sigma_e evaluated at theta (estimated
## stderr entries are therefore differentiated). NULL where it does not solve.
.d25_dr2_solve_fn <- function(model, compiled) {
  if ((compiled$dynamic$max_order %||% 1L) < 2L)
    compiled <- compile_model(model, verbose = FALSE, max_order = 2L)
  if (is.null(compiled$lead_lag_incidence) &&
      !is.null(compiled$model$lead_lag_incidence))
    compiled$lead_lag_incidence <- compiled$model$lead_lag_incidence
  sys_cache <- cache_system_structure(compiled)
  state <- new.env(parent = emptyenv())
  state$ss_warm <- NULL
  function(theta) {
    sol <- .solve_dr_for_theta(model, compiled, sys_cache, theta, state,
                               lik_init = "stationary")
    if (is.null(sol)) return(NULL)
    Se <- .get_shock_cov(model, sol$dr$exo_names, sol$params)
    dr2 <- solve_perturbation_order2(model, compiled, state$ss_warm,
                                     sol$params, sol$dr, Sigma_e = Se)
    dr2$Sigma_e <- Se
    dr2
  }
}

.d25_moment_count <- function(n_obs, n_lags, include_mean) {
  (if (include_mean) n_obs else 0L) + n_obs * (n_obs + 1L) / 2L +
    n_lags * n_obs * n_obs
}

.d25_moment_names <- function(obs, n_lags, include_mean) {
  lt <- lower.tri(diag(length(obs)), diag = TRUE)
  rn <- matrix(obs, length(obs), length(obs))
  cn <- t(rn)
  acv <- unlist(lapply(seq_len(n_lags), function(k)
    sprintf("acov%d(%s,%s)", k, as.vector(rn), as.vector(cn))))
  c(if (include_mean) paste0("mean(", obs, ")"),
    sprintf("cov(%s,%s)", rn[lt], cn[lt]), acv)
}

## [mean, vech(Gamma_0), vec(Gamma_1), ..., vec(Gamma_K)], where
## Gamma_k = Cov(y_{t+k}, y_t) and, for y_t = Dx s_{t-1} + Dv r_t,
## s_{t+1} = A s_t + B r_t (r_t white, Cov Cr):
##   Gamma_0 = Dx S Dx' + Dv Cr Dv',  Gamma_k = Dx A^{k-1} (A S Dx' + B Cr Dv').
.d25_stack_moments <- function(mu, Dx, Dv, A, B, S, Cr, Z, n_lags, include_mean) {
  G0 <- Z %*% (Dx %*% S %*% t(Dx) + Dv %*% Cr %*% t(Dv)) %*% t(Z)
  C  <- (A %*% S %*% t(Dx) + B %*% Cr %*% t(Dv)) %*% t(Z)
  ZD <- Z %*% Dx
  acv <- vector("list", n_lags)
  for (k in seq_len(n_lags)) {
    acv[[k]] <- as.vector(ZD %*% C)
    C <- A %*% C
  }
  c(if (include_mean) as.vector(Z %*% mu),
    G0[lower.tri(G0, diag = TRUE)], unlist(acv))
}

## First-order moments: y_t = ys + ghx x_{t-1} + ghu e_t, x_t = hx x_{t-1} + hu e_t.
.d25_moments_order1 <- function(dr, Sigma_e, Z, n_lags, include_mean) {
  si <- dr$state_idx
  hx <- dr$ghx[si, , drop = FALSE]; hu <- dr$ghu[si, , drop = FALSE]
  Sx <- solve_lyapunov(hx, hu %*% Sigma_e %*% t(hu))
  .d25_stack_moments(dr$ys, dr$ghx, dr$ghu, hx, hu, Sx, Sigma_e, Z,
                     n_lags, include_mean)
}

## Pruned second-order moments from the AFVRR augmented system (exact).
.d25_moments_order2 <- function(dr, Sigma_e, Z, n_lags, include_mean) {
  sys <- .order2_aug_system(dr, Sigma_e)
  Sx  <- solve_lyapunov(sys$hx, sys$hu %*% Sigma_e %*% t(sys$hu))
  Cr  <- .order2_cov_r(numeric(sys$n_s), Sx, Sigma_e)
  Sxi <- solve_lyapunov(sys$Tlin, sys$G %*% Cr %*% t(sys$G))
  mu_xi <- as.numeric(solve(diag(sys$d) - sys$Tlin, sys$cc + sys$c_u))
  mu <- sys$ys + as.numeric(sys$Dxi %*% mu_xi) + 0.5 * sys$ghss + sys$c_v
  .d25_stack_moments(mu, sys$Dxi, sys$Gv, sys$Tlin, sys$G, Sxi, Cr, Z,
                     n_lags, include_mean)
}

## Per-parameter distinctness: norm of the residual of equilibrated column j
## projected off the other columns (sine of the angle to their span); 0 when
## the parameter's moment effect is a combination of the others' (or zero).
.d25_distinctness <- function(Je, n_par) {
  if (nrow(Je) == 0L) return(rep(0, n_par))
  vapply(seq_len(n_par), function(j) {
    x <- Je[, j]
    if (n_par == 1L) return(sqrt(sum(x^2)))
    sqrt(sum(qr.resid(qr(Je[, -j, drop = FALSE]), x)^2))
  }, numeric(1)) |> pmin(1)
}

.d25_plot_sv <- function(rk1, rk2, rank_o1, rank_o2, n_par) {
  labs_ord <- c(sprintf("First order (rank %d/%d)", rank_o1, n_par),
                sprintf("Pruned order 2 (rank %d/%d)", rank_o2, n_par))
  sv <- c(rk1$singular_values, rk2$singular_values)
  tol <- c(rk1$tol, rk2$tol)
  pos <- sv[sv > 0]
  floor_val <- if (length(pos)) min(c(pos, tol[tol > 0])) / 100 else 1e-16
  df <- data.frame(
    index = rep(seq_len(n_par), 2L) + rep(c(-0.08, 0.08), each = n_par),
    order = factor(rep(labs_ord, each = n_par), levels = labs_ord),
    value = pmax(sv, floor_val),
    status = ifelse(sv > rep(tol, each = n_par), "Above tolerance", "At/below tolerance"),
    stringsAsFactors = FALSE
  )
  tol_df <- data.frame(order = factor(labs_ord, levels = labs_ord),
                       tol = pmax(tol, floor_val))
  cols <- stats::setNames(dynhr_palette_vibrant[c(1L, 3L)], labs_ord)
  ggplot2::ggplot(df, ggplot2::aes(x = index, y = value, colour = order)) +
    ggplot2::geom_hline(data = tol_df,
                        ggplot2::aes(yintercept = tol, colour = order),
                        linetype = "dashed", linewidth = 0.5) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(shape = status), size = 3) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(breaks = seq_len(n_par)) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_shape_manual(values = c("Above tolerance" = 16,
                                           "At/below tolerance" = 1),
                                name = NULL, drop = FALSE) +
    theme_dynhr() +
    ggplot2::theme(legend.position = "bottom", legend.box = "vertical") +
    ggplot2::labs(
      title = "D25: Moment-Jacobian singular values, first order vs pruned order 2",
      subtitle = paste0("Equilibrated Jacobians; dashed = rank tolerance (from finite-difference noise).\n",
                        "Open points = unidentified directions",
                        if (any(sv <= floor_val)) "; exact zeros are drawn at the plot floor." else "."),
      x = "Singular value index (largest first)",
      y = "Singular value (log scale)")
}

.d25_plot_params <- function(param_id) {
  n <- nrow(param_id)
  lvl <- param_id$parameter[order(param_id$order2_distinct, param_id$order1_distinct)]
  labs_ord <- c("First order", "Pruned order 2")
  df <- data.frame(
    parameter = factor(rep(param_id$parameter, 2L), levels = lvl),
    order = factor(rep(labs_ord, each = n), levels = labs_ord),
    distinct = c(param_id$order1_distinct, param_id$order2_distinct),
    status = c(param_id$order1_status, param_id$order2_status),
    stringsAsFactors = FALSE
  )
  df$y <- as.integer(df$parameter) + ifelse(df$order == labs_ord[1L], 0.12, -0.12)
  seg <- data.frame(y0 = match(param_id$parameter, lvl) + 0.12,
                    y1 = match(param_id$parameter, lvl) - 0.12,
                    x0 = param_id$order1_distinct, x1 = param_id$order2_distinct)
  cols <- stats::setNames(dynhr_palette_vibrant[c(1L, 3L)], labs_ord)
  ggplot2::ggplot(df, ggplot2::aes(x = distinct, y = y)) +
    ggplot2::geom_segment(data = seg,
                          ggplot2::aes(x = x0, xend = x1, y = y0, yend = y1),
                          colour = dynhr_na_colour, linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(colour = order, shape = status), size = 3.2,
                        stroke = 1.1) +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::scale_shape_manual(values = c(Identified = 16, Unidentified = 1),
                                name = NULL, drop = FALSE) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(breaks = seq_along(lvl), labels = lvl,
                                expand = ggplot2::expansion(add = 0.5)) +
    theme_dynhr() +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(
      title = "D25: Per-parameter identification, first order vs pruned order 2",
      subtitle = paste0("Distinctness = sine of the angle between a parameter's moment effect\n",
                        "and the span of the other parameters' effects (0 = not separable)."),
      x = "Distinctness (0 = unidentified, 1 = orthogonal to other parameters)",
      y = NULL)
}
