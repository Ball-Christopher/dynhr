## R/diag-pre-d26-calibration-sensitivity.R
## --------------------------------------------------------------------------
## D26 Calibration sensitivity: how the estimated parameters depend on the
## calibrated (fixed) ones (Iskrev 2019b).
##
## Two statistics, both on the moment map m(theta_c, theta_e) and a weight W:
##   1. The implicit-function derivative of the minimum-distance estimator
##        theta_e_hat = argmin (m* - m(theta_c, theta_e))' W (m* - m(.))
##      at the truth, d theta_e_hat / d theta_c = -(J_e' W J_e)^-1 J_e' W J_c:
##      the first-order shift of the estimates caused by a miscalibration.
##      Reported as a derivative and as an elasticity.
##   2. The stability of the estimated parameters' asymptotic precision
##      (SE_i = sqrt(diag((J_e' W J_e)^-1)), the D20 Fisher bound) and of the
##      rank of J_e when each calibrated parameter is moved by +/- 10%.
## --------------------------------------------------------------------------

#' D26. Calibration sensitivity diagnostic (Iskrev 2019b)
#'
#' Measures how the estimates of the free parameters \eqn{\theta_e} depend on
#' the values at which the other parameters \eqn{\theta_c} were calibrated.
#'
#' \strong{Miscalibration derivative.} For the minimum-distance estimator
#' \eqn{\hat\theta_e = \arg\min (m^* - m(\theta_c, \theta_e))' W
#' (m^* - m(\theta_c, \theta_e))}, evaluated where the model fits
#' (\eqn{m^* = m(\theta_c, \theta_e)}), the implicit-function theorem gives
#' \deqn{D = \partial\hat\theta_e / \partial\theta_c' =
#'   -(J_e' W J_e)^{-1} J_e' W J_c,}
#' with \eqn{J_e}, \eqn{J_c} the moment Jacobians (central differences). A
#' calibration error \eqn{\Delta} shifts the estimates by \eqn{D\Delta} to
#' first order (exactly, for a linear moment map). \code{result$elasticity}
#' is \eqn{E_{ij} = D_{ij}\,\theta_{c,j} / \theta_{e,i}} (the percent shift in
#' \eqn{\hat\theta_{e,i}} per percent error in \eqn{\theta_{c,j}}); it is
#' \code{NA} where either value is zero, and such cells are not gated. With an
#' exactly identified model \eqn{D} does not depend on \eqn{W}.
#'
#' \strong{Precision stability.} Each calibrated parameter is moved by the
#' relative amounts in \code{perturbation_grid} (a zero-valued parameter by
#' the same amounts in absolute terms), and \eqn{J_e} is recomputed. The
#' standard errors \eqn{SE_i = \sqrt{[(J_e' W J_e)^{-1}]_{ii}}} and the rank of
#' \eqn{W^{1/2} J_e} (shared D1/D20 helper) are recorded. \code{se_ratio} is
#' \eqn{\max SE_i / \min SE_i} over the baseline and all perturbations; it is
#' invariant to the units of the moments and to the value of
#' \eqn{\theta_{e,i}}. With the default \code{W = I}, SE is not an asymptotic
#' standard error (no sample size); only its ratio is used.
#'
#' \strong{Gate (0.9.4).} An elasticity breach
#' (\eqn{\max_j |E_{ij}| >} \code{elasticity_tol}) is a \strong{WARN} and never
#' a FAIL: the derivative is a descriptive sensitivity measure (Andrews,
#' Gentzkow & Shapiro 2017) -- it says the estimate inherits a number you chose,
#' which is worth reporting, not evidence that the estimation is wrong. FAIL is
#' reserved for a loss of precision (\code{se_ratio >} \code{fragility_tol}) or
#' of identification (rank drop) at some perturbation, which do invalidate the
#' inference. PASS iff neither fires. The result is
#' informational (\code{pass = NA}) when \eqn{J_e} is rank deficient at the
#' baseline (the derivative does not exist; see D1) or when no calibrated
#' parameter moves any moment -- the signature of a \code{model_solve_fn}
#' that does not re-solve the model (e.g. a fixed decision rule), which
#' would otherwise report every estimate as perfectly robust.
#'
#' @param model_solve_fn Function: full named parameter vector
#'   (calibrated + estimated) -> numeric vector of moments. The moment names
#'   (if any) must be the same at every call; a call may signal failure by
#'   returning non-finite values.
#' @param theta_c,theta_e Numeric vectors of calibrated / estimated parameter
#'   values. When named they are matched to \code{param_names_c} /
#'   \code{param_names_e} by NAME (not position).
#' @param param_names_c,param_names_e Character vectors of calibrated /
#'   estimated parameter names (default: \code{names(theta_c)} /
#'   \code{names(theta_e)}). They must not overlap.
#' @param weight Optional symmetric positive semi-definite
#'   n_moment x n_moment weight matrix W (matched to the moment names when both
#'   are named). \code{NULL} (default) = identity.
#' @param perturbation_grid Relative perturbations applied to each calibrated
#'   parameter for the precision-stability check (default \code{c(-0.1, 0.1)});
#'   the baseline (0) is always included.
#' @param eps Finite-difference step.
#' @param elasticity_tol WARN when \eqn{|E_{ij}|} exceeds this (default 1: a
#'   1\% calibration error moves the estimate by more than 1\%). \strong{The
#'   value 1 is a package choice}: Andrews, Gentzkow & Shapiro (2017) define
#'   the sensitivity measure but set no threshold on it. Unit elasticity is
#'   simply the natural break-even point -- below it the calibration error is
#'   attenuated, above it amplified.
#' @param fragility_tol Flag \code{se_ratio} above this (default 10).
#' @param tol_rank \code{NULL} (finite-difference-aware rank tolerance, as D1)
#'   or a relative rank tolerance.
#' @param calibration_choice Which parameters to calibrate at all
#'   (Alegre Canton 2026, arXiv:2606.25688): a named list of options, or
#'   \code{NULL}/\code{FALSE} to skip. From the Jacobians above, every split of
#'   \eqn{(\theta_c, \theta_e)} into estimated and fixed sets is ranked by the
#'   sensitivity statistic \eqn{K_S} (worst-case first-order bias of an object
#'   of interest per normalised calibration error); see
#'   \code{d26_calibration_choice} for the method. Options: \code{target}
#'   (parameter names -- default the estimated \eqn{\theta_e}, each measured
#'   relative to its value and always estimated -- or a
#'   \code{function(eta)} of the full named parameter vector), \code{ranges}
#'   (normalisation widths \eqn{\Delta_j}: named widths or a (min, max) matrix
#'   with parameter row names; default \eqn{|\eta_j|}, 1 at zero, so errors
#'   are relative), \code{epsilon} (default 0.05), \code{n_obs} (enables the
#'   \eqn{\sqrt{\log n/n}} weak-identification cutoff), \code{must_estimate},
#'   \code{must_calibrate}, \code{min_estimated}, \code{max_estimated},
#'   \code{max_partitions} (default 4096). Informational: it never changes the
#'   D26 badge.
#' @param meta Plot provenance (see \code{diag_meta}).
#' @return dynhr_diagnostic; \code{result} holds \code{derivative} (D),
#'   \code{elasticity} (E), \code{jacobian_e}, \code{jacobian_c},
#'   \code{rank_e}, \code{inert_calibrated}, \code{se_array}
#'   (n_est x n_cal x n_pert), \code{strength_array} (\eqn{|\theta_e|/SE}),
#'   \code{rank_array} (n_cal x n_pert), \code{perturbation_values},
#'   \code{fragility_table}, \code{fragile_params}, \code{status},
#'   \code{calibration_choice} (see \code{d26_calibration_choice}; \code{NULL}
#'   when switched off).
#' @references
#'   Iskrev, N. (2019). What to expect when you're calibrating: measuring the
#'   effect of calibration on the estimation of macroeconomic models.
#'   \emph{Journal of Economic Dynamics and Control}, 99, 54--81.
#'
#'   Andrews, I., Gentzkow, M., & Shapiro, J. M. (2017). Measuring the
#'   sensitivity of parameter estimates to estimation moments.
#'   \emph{Quarterly Journal of Economics}, 132(4), 1553--1592.
#'
#'   Alegre Canton, J. (2026). Choosing what to calibrate and what to estimate
#'   in structural models. arXiv:2606.25688 [econ.EM], June 2026.
#' @noRd
d26_calibration_sensitivity <- function(model_solve_fn,
                                        theta_c,
                                        theta_e,
                                        param_names_c = NULL,
                                        param_names_e = NULL,
                                        weight = NULL,
                                        perturbation_grid = c(-0.10, 0.10),
                                        eps = 1e-5,
                                        elasticity_tol = 1,
                                        fragility_tol = 10,
                                        tol_rank = NULL,
                                        calibration_choice = list(),
                                        meta = NULL) {
  if (!is.function(model_solve_fn))
    .dynhr_abort("D26: `model_solve_fn` must be a function.")
  theta_c <- .d26_align_theta(theta_c, param_names_c, "theta_c", "theta_c_")
  theta_e <- .d26_align_theta(theta_e, param_names_e, "theta_e", "theta_e_")
  names_c <- names(theta_c)
  names_e <- names(theta_e)
  both <- intersect(names_c, names_e)
  if (length(both))
    .dynhr_abort("D26: parameter(s) both calibrated and estimated: ",
                 paste(both, collapse = ", "), ".")
  n_cal <- length(theta_c)
  n_est <- length(theta_e)
  perturbation_grid <- sort(unique(c(0, as.numeric(perturbation_grid))))
  if (!all(is.finite(perturbation_grid)) || any(perturbation_grid <= -1))
    .dynhr_abort("D26: `perturbation_grid` must be finite relative changes > -1.")
  n_pert <- length(perturbation_grid)
  i0 <- which(perturbation_grid == 0)

  ## ---- moment map with a fixed moment order ------------------------------
  f0 <- model_solve_fn(c(theta_c, theta_e))
  if (!is.numeric(f0) || length(f0) == 0L || !all(is.finite(f0)))
    .dynhr_abort("D26: `model_solve_fn` returned no / non-finite moments at ",
                 "the baseline parameters.")
  n_mom <- length(f0)
  mom_names <- names(f0) %||% paste0("m_", seq_len(n_mom))
  if (anyDuplicated(mom_names))
    .dynhr_abort("D26: `model_solve_fn` returned duplicated moment names.")
  named_out <- !is.null(names(f0))
  moments <- function(th_c, th_e) {
    out <- model_solve_fn(c(th_c, th_e))
    if (named_out) {
      if (!is.numeric(out) || !all(mom_names %in% names(out)))
        return(rep(NA_real_, n_mom))
      out <- out[mom_names]
    }
    if (!is.numeric(out) || length(out) != n_mom) return(rep(NA_real_, n_mom))
    as.numeric(out)
  }

  ## ---- weight: Jw = W^{1/2} J so that J'WJ = crossprod(Jw) ----------------
  Wh <- .d26_weight_sqrt(weight, mom_names)
  weigh <- function(J) {
    J <- if (is.null(Wh)) J else Wh %*% J
    dimnames(J) <- list(NULL, colnames(J))
    J
  }
  jac <- function(fn, at, h, cn) {
    J <- .numerical_jacobian(fn, at, eps = h)
    dimnames(J) <- list(mom_names, cn)
    J
  }

  ## ---- baseline Jacobians (step h and 2h for the FD error) -----------------
  fe <- function(th_e) moments(theta_c, th_e)
  fc <- function(th_c) moments(th_c, theta_e)
  J_e  <- jac(fe, theta_e, eps, names_e)
  J_e2 <- jac(fe, theta_e, 2 * eps, names_e)
  J_c  <- jac(fc, theta_c, eps, names_c)
  J_c2 <- jac(fc, theta_c, 2 * eps, names_c)
  if (!all(is.finite(J_e)) || !all(is.finite(J_c)))
    .dynhr_abort("D26: the moment Jacobian is non-finite at the baseline ",
                 "(the model does not solve within +/- eps of it).")

  ## A calibrated parameter is inert when its moment column is within 10x of
  ## its own finite-difference error (exactly zero for a moment function that
  ## ignores it).
  fd_c <- if (all(is.finite(J_c2))) apply(abs(J_c - J_c2), 2, max) else rep(0, n_cal)
  inert <- names_c[apply(abs(J_c), 2, max) <= 10 * fd_c]

  rk0 <- .ident_equilibrated_rank(weigh(J_e),
                                  if (all(is.finite(J_e2))) weigh(J_e2),
                                  tol_rank = tol_rank)
  rel_tol <- rk0$tol / max(rk0$sv_max, 1e-300)
  identified0 <- rk0$rank == n_est

  ## ---- miscalibration derivative -----------------------------------------
  D <- matrix(NA_real_, n_est, n_cal, dimnames = list(estimated = names_e,
                                                      calibrated = names_c))
  if (identified0) {
    Jew <- weigh(J_e)
    cs <- sqrt(colSums(Jew^2))
    D[] <- -qr.coef(qr(sweep(Jew, 2, cs, "/")), weigh(J_c)) / cs
  }
  E <- D * outer(1 / theta_e, theta_c)
  E[outer(theta_e == 0, theta_c == 0, "|")] <- NA_real_
  dimnames(E) <- dimnames(D)

  ## ---- precision stability under +/- perturbations -------------------------
  ## SE from the inverse on the identified subspace (as D20); a parameter that
  ## loads on a null-space direction has SE = Inf.
  se_of <- function(Jw) {
    rk <- .ident_equilibrated_rank(Jw, tol_rank = rel_tol)
    unid <- union(rk$unidentified_params,
                  names_e[rowSums(rk$null_space^2) >= 0.01])
    cs <- sqrt(colSums(Jw^2))
    cs[cs <= 0] <- 1
    k <- rk$rank
    sv <- svd(sweep(Jw, 2, cs, "/"), nu = 0, nv = n_est)
    Vk <- sv$v[, seq_len(k), drop = FALSE]
    se <- sqrt(rowSums(sweep(Vk, 2, sv$d[seq_len(k)], "/")^2)) / cs
    se[names_e %in% unid] <- Inf
    list(se = se, rank = k, unid = unid)
  }
  base <- se_of(weigh(J_e))
  pert_lab <- sprintf("%+g%%", 100 * perturbation_grid)
  se_array <- array(NA_real_, dim = c(n_est, n_cal, n_pert),
                    dimnames = list(estimated = names_e, calibrated = names_c,
                                    perturbation = pert_lab))
  rank_array <- matrix(NA_integer_, n_cal, n_pert,
                       dimnames = list(calibrated = names_c, perturbation = pert_lab))
  pert_values <- matrix(NA_real_, n_cal, n_pert,
                        dimnames = list(names_c, pert_lab))
  lost_id <- stats::setNames(rep(FALSE, n_est), names_e)
  n_failed <- 0L
  for (j in seq_len(n_cal)) {
    step_j <- if (theta_c[[j]] == 0) 1 else abs(theta_c[[j]])
    for (k in seq_len(n_pert)) {
      th_c <- theta_c
      th_c[j] <- theta_c[[j]] + perturbation_grid[k] * step_j
      pert_values[j, k] <- th_c[[j]]
      if (k == i0) {
        s_k <- base
      } else {
        Jk <- jac(function(th_e) moments(th_c, th_e), theta_e, eps, names_e)
        if (!all(is.finite(Jk))) { n_failed <- n_failed + 1L; next }
        s_k <- se_of(weigh(Jk))
        lost_id[names_e %in% s_k$unid] <- TRUE
      }
      se_array[, j, k] <- s_k$se
      rank_array[j, k] <- s_k$rank
    }
  }
  strength_array <- abs(theta_e) / se_array

  ## ---- per-parameter table and gate -----------------------------------------
  fragility_table <- do.call(rbind, lapply(seq_len(n_est), function(i) {
    se_i <- as.numeric(se_array[i, , ])
    se_i <- se_i[!is.na(se_i)]
    fin <- se_i[is.finite(se_i)]
    ratio <- if (length(fin) < length(se_i)) Inf
             else if (length(fin)) max(fin) / min(fin) else NA_real_
    e_i <- abs(E[i, ])
    max_el <- if (all(is.na(e_i))) NA_real_ else max(e_i, na.rm = TRUE)
    worst <- if (is.na(max_el)) NA_character_ else names_c[which.max(e_i)]
    fr_el <- isTRUE(max_el > elasticity_tol)
    fr_se <- isTRUE(ratio > fragility_tol) || lost_id[[i]]
    data.frame(
      parameter = names_e[i], theta = theta_e[[i]],
      se_baseline = base$se[i],
      s_baseline = abs(theta_e[[i]]) / base$se[i],
      se_min = if (length(fin)) min(fin) else NA_real_,
      se_max = if (length(fin)) max(fin) else NA_real_,
      se_ratio = ratio, lost_identification = lost_id[[i]],
      max_abs_elasticity = max_el, worst_calibrated = worst,
      fragile_elasticity = fr_el, fragile_precision = fr_se,
      fragile = identified0 && (fr_el || fr_se),
      stringsAsFactors = FALSE, row.names = NULL)
  }))
  fragile_params <- fragility_table$parameter[fragility_table$fragile]

  ## An elasticity breach says "your estimate moves with a number you chose",
  ## which is information about the calibration, not evidence the estimation is
  ## broken -- so it WARNs and never FAILs (Andrews, Gentzkow & Shapiro 2017
  ## present exactly this derivative as a descriptive sensitivity measure).
  ## Losing precision or rank across the perturbation grid IS a failure.
  ft <- fragility_table
  fail_params <- ft$parameter[identified0 & ft$fragile_precision]
  sensitive_params <- ft$parameter[identified0 & ft$fragile_elasticity &
                                     !ft$fragile_precision]

  all_inert <- n_cal > 0L && length(inert) == n_cal
  status <- if (all_inert) "inert"
            else if (!identified0) "unidentified"
            else if (length(fail_params)) "fragile"
            else if (length(sensitive_params)) "sensitive" else "robust"
  pass <- switch(status, robust = TRUE, sensitive = TRUE, fragile = FALSE, NA)
  warn <- identical(status, "sensitive")

  ## ---- which parameters to calibrate (Alegre Canton 2026) -------------------
  ## Informational: ranks every estimated/fixed split from the Jacobians above;
  ## the badge stays D26's own (the paper gives a ranking, not a pass/fail).
  choice <- .d26_choice_from_d26(calibration_choice, theta_c, theta_e,
                                 J_e, J_e2, J_c, J_c2, weight, tol_rank, eps,
                                 all_inert)
  choice_txt <- if (is.null(choice)) "" else
    paste0(" Calibration choice (INFO): ", choice$message)

  ## ---- text ----------------------------------------------------------------
  w_lab <- if (is.null(Wh)) "W = I (unit moment weights)" else "user W"
  inert_txt <- if (length(inert) && !all_inert)
    sprintf(" Calibrated parameters that move no moment (INFO): %s.",
            paste(inert, collapse = ", ")) else ""
  fail_txt <- if (n_failed)
    sprintf(" %d perturbed calibration(s) did not solve and were skipped.", n_failed) else ""
  unid_txt <- paste(rk0$unidentified_params, collapse = ", ")
  body <- switch(status,
    unidentified = sprintf(
      "INFO -- estimated parameters are not identified at the baseline (rank %d/%d%s); the miscalibration derivative does not exist (see D1).",
      rk0$rank, n_est, if (nzchar(unid_txt)) paste0("; ", unid_txt) else ""),
    inert = paste0(
      "INFO -- no calibrated parameter moves any moment, so the calibration ",
      "cannot matter for these moments. If the calibrated parameters are ",
      "structural, `model_solve_fn` is probably not re-solving the model ",
      "(e.g. it reuses a fixed decision rule); D26 cannot be assessed."),
    robust = sprintf(
      "PASS -- max |elasticity| = %.2f (tol %.2g), max SE ratio = %.2f (tol %.2g).",
      .d26_max(fragility_table$max_abs_elasticity), elasticity_tol,
      .d26_max(fragility_table$se_ratio), fragility_tol),
    sensitive = sprintf(
      "WARN -- %d estimated parameter(s) move more than the calibration does: %s (max |elasticity| = %.2f, tol %.2g). Precision is stable (max SE ratio = %.2f, tol %.2g).",
      length(sensitive_params), paste(sensitive_params, collapse = ", "),
      .d26_max(fragility_table$max_abs_elasticity), elasticity_tol,
      .d26_max(fragility_table$se_ratio), fragility_tol),
    fragile = sprintf(
      "FAIL -- %d estimated parameter(s) lose precision or identification across the perturbation grid: %s (max SE ratio = %.2f, tol %.2g; max |elasticity| = %.2f, tol %.2g).",
      length(fail_params), paste(fail_params, collapse = ", "),
      .d26_max(fragility_table$se_ratio), fragility_tol,
      .d26_max(fragility_table$max_abs_elasticity), elasticity_tol))
  summary_text <- sprintf(
    "D26 Calibration sensitivity: %d estimated, %d calibrated, %d moments, %s. %s%s%s%s",
    n_est, n_cal, n_mom, w_lab, body, inert_txt, fail_txt, choice_txt)

  top <- character(0)
  if (identified0) {
    ord <- order(-abs(E), na.last = NA)
    ord <- utils::head(ord, 5L)
    top <- sprintf("%s<-%s: d=%.3g elast=%.2f",
                   names_e[row(E)[ord]], names_c[col(E)[ord]], D[ord], E[ord])
  }
  llm <- paste(c(
    sprintf("D26 | Calibration Sensitivity (Iskrev 2019b) | %s",
            .badge_str(list(pass = pass, errored = FALSE, warn = warn))),
    sprintf("  status=%s n_est=%d n_cal=%d n_mom=%d rank_e=%d weight=%s",
            status, n_est, n_cal, n_mom, rk0$rank, if (is.null(Wh)) "I" else "user"),
    "  stat: D = -(Je'WJe)^-1 Je'WJc (miscalibration shift); elasticity = D*theta_c/theta_e",
    sprintf("  gate: WARN if max|elasticity|>%.2g; FAIL if SE ratio over +/-perturbations>%.2g OR rank loss",
            elasticity_tol, fragility_tol),
    if (length(top)) sprintf("  largest: %s", paste(top, collapse = "; ")),
    if (length(inert)) sprintf("  inert_calibrated: %s", paste(inert, collapse = ", ")),
    if (length(fragile_params)) sprintf("  fragile: %s", paste(fragile_params, collapse = ", ")),
    if (!is.null(choice)) sprintf("  calibration_choice: %s", choice$message),
    sprintf("  action: %s", switch(status,
      unidentified = "Estimated parameters not identified (D1); fix identification first.",
      inert = "Supply a model_solve_fn that re-solves the model at the calibrated values.",
      robust = "Estimates are robust to the calibration.",
      sensitive = "Estimates move more than 1-for-1 with the listed calibrated parameters; report the sensitivity (Andrews-Gentzkow-Shapiro 2017) or estimate them.",
      fragile = "Precision or identification of the listed estimates collapses under plausible recalibration; estimate the calibrated parameters or report results across calibrations."))
  ), collapse = "\n")

  plots <- .d26_plots(D, E, fragility_table, status, elasticity_tol,
                      fragility_tol, perturbation_grid, w_lab, meta)

  .make_result(
    result = list(
      status = status,
      derivative = D,
      elasticity = E,
      jacobian_e = J_e,
      jacobian_c = J_c,
      rank_e = rk0$rank,
      rank_tolerance = rel_tol,
      unidentified_params = rk0$unidentified_params,
      inert_calibrated = inert,
      se_baseline = stats::setNames(base$se, names_e),
      se_array = se_array,
      strength_array = strength_array,
      rank_array = rank_array,
      perturbation_values = pert_values,
      n_failed = n_failed,
      fragility_table = fragility_table,
      fragile_params = fragile_params,
      fail_params = fail_params,
      sensitive_params = sensitive_params,
      weighted = !is.null(Wh),
      calibration_choice = choice
    ),
    pass = pass,
    warn = warn,
    plots = plots,
    summary = summary_text,
    llm_summary = llm
  )
}


## Max of the finite values (NA when there are none), for text output.
.d26_max <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) max(x) else NA_real_
}


## Resolve a parameter vector against its names BY NAME.
.d26_align_theta <- function(theta, nms, what, prefix) {
  if (!is.numeric(theta) || length(theta) == 0L || !all(is.finite(theta)))
    .dynhr_abort("D26: `", what, "` must be a non-empty finite numeric vector.")
  if (is.null(nms)) nms <- names(theta) %||% paste0(prefix, seq_along(theta))
  nms <- as.character(nms)
  if (anyDuplicated(nms) || any(!nzchar(nms)))
    .dynhr_abort("D26: names for `", what, "` must be unique and non-empty.")
  if (!is.null(names(theta))) {
    miss <- setdiff(nms, names(theta))
    extra <- setdiff(names(theta), nms)
    if (length(miss) || length(extra) || length(theta) != length(nms))
      .dynhr_abort("D26: names(", what, ") do not match its parameter names",
                   if (length(miss)) paste0("; missing: ", paste(miss, collapse = ", ")),
                   if (length(extra)) paste0("; unexpected: ", paste(extra, collapse = ", ")),
                   ".")
    theta <- theta[nms]
  } else if (length(theta) != length(nms)) {
    .dynhr_abort("D26: `", what, "` has ", length(theta), " values but ",
                 length(nms), " names.")
  }
  stats::setNames(as.numeric(theta), nms)
}


## Symmetric square root of a PSD weight matrix aligned to the moments
## (NULL for the identity).
.d26_weight_sqrt <- function(weight, mom_names) {
  if (is.null(weight)) return(NULL)
  W <- as.matrix(weight)
  n <- length(mom_names)
  if (!is.numeric(W) || !all(dim(W) == n) || !all(is.finite(W)))
    .dynhr_abort("D26: `weight` must be a finite ", n, " x ", n, " matrix.")
  if (!is.null(rownames(W)) && !is.null(colnames(W))) {
    if (!all(mom_names %in% rownames(W)) || !all(mom_names %in% colnames(W)))
      .dynhr_abort("D26: dimnames(weight) do not match the moment names.")
    W <- W[mom_names, mom_names, drop = FALSE]
  }
  if (max(abs(W - t(W))) > 1e-10 * max(abs(W), 1e-300))
    .dynhr_abort("D26: `weight` must be symmetric.")
  ev <- eigen((W + t(W)) / 2, symmetric = TRUE)
  if (min(ev$values) < -1e-10 * max(abs(ev$values)))
    .dynhr_abort("D26: `weight` must be positive semi-definite.")
  ev$vectors %*% (sqrt(pmax(ev$values, 0)) * t(ev$vectors))
}


## D26 plots: elasticity heatmap and SE-ratio bars.
.d26_plots <- function(D, E, ft, status, elasticity_tol, fragility_tol,
                       grid, w_lab, meta) {
  plots <- list()
  if (!requireNamespace("ggplot2", quietly = TRUE) || status == "unidentified")
    return(plots)
  est <- rownames(D)
  cal <- colnames(D)
  df <- expand.grid(Estimated = est, Calibrated = cal, stringsAsFactors = FALSE)
  df$elasticity <- as.numeric(E)
  df$derivative <- as.numeric(D)
  df$Estimated <- factor(df$Estimated, levels = rev(est))
  df$Calibrated <- factor(df$Calibrated, levels = cal)
  if (status == "inert") df$elasticity <- NA_real_
  df$over <- !is.na(df$elasticity) & abs(df$elasticity) > elasticity_tol
  df$label <- if (status == "inert") "not assessable"
              else ifelse(is.na(df$elasticity),
                          sprintf("n/a\n(D = %.2g)", df$derivative),
                          sprintf("%.2f", df$elasticity))
  lim <- max(c(abs(df$elasticity), elasticity_tol), na.rm = TRUE)
  df$dark <- !is.na(df$elasticity) & abs(df$elasticity) > 0.6 * lim
  small <- length(est) * length(cal) <= 80L
  p_el <- ggplot2::ggplot(df, ggplot2::aes(x = Calibrated, y = Estimated,
                                           fill = elasticity)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_tile(data = df[df$over, , drop = FALSE], fill = NA,
                       colour = "#1A1A1A", linewidth = 1) +
    scale_fill_dynhr_sunset(limits = c(-lim, lim), name = "Elasticity",
                            guide = ggplot2::guide_colourbar(barwidth = 12,
                                                             barheight = 0.5)) +
    theme_dynhr_diagnostic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = "D26: Estimate shift per calibration error (elasticity)",
      subtitle = sprintf(
        "%% change in estimate per 1%% error in calibrated value; outlined: abs > %.2g. %s.",
        elasticity_tol, w_lab),
      x = "Calibrated parameter", y = "Estimated parameter")
  if (small)
    p_el <- p_el +
      ggplot2::geom_text(ggplot2::aes(label = label, colour = dark), size = 3.2,
                         lineheight = 0.85, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = c(`FALSE` = "#1A1A1A", `TRUE` = "#FFFFFF"),
                                   guide = "none")
  if (status == "inert")
    p_el <- p_el + ggplot2::labs(
      subtitle = "No calibrated parameter moves any moment: sensitivity cannot be assessed.")
  plots$elasticity_heatmap <- .apply_meta(p_el, meta)

  if (status != "inert") {
    bd <- ft
    bd$ratio_plot <- ifelse(is.finite(bd$se_ratio), bd$se_ratio, NA_real_)
    top <- max(c(bd$ratio_plot, fragility_tol), na.rm = TRUE) * 2
    bd$ratio_plot[bd$lost_identification | is.infinite(bd$se_ratio)] <- top
    bd$state <- ifelse(bd$lost_identification, "Loses identification (drawn at top)",
                ifelse(bd$fragile_precision, "Fragile precision", "Stable precision"))
    bd$parameter <- factor(bd$parameter, levels = bd$parameter[order(bd$ratio_plot)])
    p_se <- ggplot2::ggplot(bd, ggplot2::aes(x = parameter, y = ratio_plot,
                                             colour = state)) +
      ggplot2::geom_segment(ggplot2::aes(xend = parameter, y = 1, yend = ratio_plot),
                            linewidth = 0.8) +
      ggplot2::geom_point(size = 2.6) +
      ggplot2::geom_hline(yintercept = fragility_tol, linetype = "dashed",
                          colour = dynhr_colours$grey, linewidth = 0.5) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_log10() +
      ggplot2::scale_colour_manual(
        values = c("Stable precision" = dynhr_colours$mid_blue,
                   "Fragile precision" = dynhr_colours$orange,
                   "Loses identification (drawn at top)" = dynhr_colours$red),
        name = NULL) +
      theme_dynhr_diagnostic() +
      ggplot2::labs(
        title = "D26: Precision stability across calibrations",
        subtitle = sprintf(
          "Max/min SE over calibrated values moved by %s; dashed: fragility tol %.2g.",
          paste(sprintf("%+g%%", 100 * grid[grid != 0]), collapse = ", "),
          fragility_tol),
        x = NULL, y = "SE ratio, max / min (log scale; 1 = unaffected)")
    plots$se_ratio <- .apply_meta(p_se, meta)
  }
  plots
}
