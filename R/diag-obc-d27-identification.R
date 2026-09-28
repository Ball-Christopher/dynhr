## R/diag-obc-d27-identification.R
## --------------------------------------------------------------------------
## Phase H: D27 — OBC/Piecewise-Linear Identification Diagnostic
##
## Assesses how occasionally binding constraints (OBCs) affect local parameter
## identification. Each OccBin regime (slack, or a subset of constraints
## binding) has its own linear policy, hence its own model-implied moments
## and identification Jacobian.
##
## Algorithm:
##   1. At every parameter point the model is re-solved: steady state, slack
##      first-order policy, system matrices, and (when a tag bound names a
##      parameter) the OBC specs.  Each binding regime's policy is the OccBin
##      policy obc_solve_binding() at THAT point, so the regime Jacobian moves
##      with theta.  (Pre-0.9.4 the binding policy was cached at the
##      calibration and the binding-regime Jacobian was identically zero.)
##   2. Regime-conditional moments: mean, variance, contemporaneous covariance
##      and lag-1 autocovariance of the selected variables, treating the
##      regime policy as a stationary linear process.  The regime is HELD FIXED
##      while differentiating, so the moments are smooth in theta and central
##      differences never straddle the constraint's kink.
##   3. Regime weights = shadow probability of the regime under the slack
##      regime's stationary Gaussian distribution (probability that exactly
##      the regime's constraints would be violated by the unconstrained
##      solution), or user-supplied frequencies.
##   4. Rank per regime and of the probability-weighted pooled Fisher
##      information, with the FD-noise-aware tolerance shared with D1/D20
##      (.ident_equilibrated_rank).  Parameters identified only through
##      regimes with probability < rare_prob are flagged.
##
## References:
##   Guerrieri & Iacoviello (2015). OccBin: A toolkit for solving dynamic
##     models with occasionally binding constraints easily. JME 70.
##   Iskrev, N. (2010). Local identification in DSGE models. JME 57.
## --------------------------------------------------------------------------

#' D27. OBC/Piecewise-Linear Identification Diagnostic
#'
#' Assesses how occasionally binding constraints affect local parameter
#' identification. Each OBC regime (all slack, or a subset of constraints
#' binding) induces a different linear policy (OccBin, Guerrieri-Iacoviello
#' 2015). For each regime the diagnostic differentiates the regime-conditional
#' moments (mean, variance, contemporaneous covariance and lag-1
#' autocovariance of the selected variables) with respect to the parameters,
#' re-solving the model and the binding-regime policy at every
#' finite-difference point. The regime itself is held fixed while
#' differentiating, so the moments are smooth and the difference quotient
#' never crosses the constraint's kink.
#'
#' \strong{Regime weights.} By default each regime is weighted by its shadow
#' probability: the probability, under the slack regime's stationary Gaussian
#' distribution, that exactly that regime's constraints are violated by the
#' unconstrained solution. This approximates the OccBin binding frequency (it
#' is exact when the constrained variable's shadow value does not feed back
#' from the constraint). Pass \code{regime_probs} (e.g. simulated or smoothed
#' binding frequencies) to override it.
#'
#' \strong{Verdict.} Rank is decided per regime and for the pooled Fisher
#' information \eqn{\sum_r w_r J_r'J_r}, with the finite-difference-aware
#' tolerance shared with D1/D20. FAIL when a parameter is unidentified even
#' after pooling all regimes with positive weight. WARN when every parameter is
#' identified but some are identified only through regimes whose probability is
#' below \code{rare_prob} (e.g. a bound that matters only in a rarely binding
#' regime): that is a statement about how informative the data are, the OBC
#' analogue of weak instruments, and the OccBin literature treats it
#' qualitatively rather than as a non-identification result.
#' A parameter that is unidentified inside one regime but identified
#' elsewhere (e.g. Taylor-rule coefficients at the ZLB) is expected under
#' OccBin and is reported, not failed. Strength
#' \eqn{|\theta_i| / SE_i} uses unit moment weights (D20 \code{"none"}), so it
#' is unit-dependent and informational.
#'
#' @param model           A dynhr_mod object with OBC (MCP) tags.
#' @param params          Named full parameter vector at the evaluation point.
#' @param compiled        Compiled model (compiled here when NULL).
#' @param obc_specs       List of OBC specs. When NULL they are parsed from
#'   \code{model} at every parameter point, so a bound given by a parameter
#'   name moves with that parameter; user-supplied specs are held fixed.
#' @param regime_subset   Integer regime indices (bit j set = constraint j
#'   binds). NULL (default): all slack plus every single-binding regime.
#' @param regime_probs    Optional regime weights, one per element of
#'   \code{regime_subset} (in that order). NULL: shadow probabilities.
#' @param param_names     Parameters to analyse (default \code{names(params)}).
#' @param obs_names       Variables whose moments are used (default: all
#'   endogenous variables).
#' @param eps             Finite-difference step (a second Jacobian at
#'   \code{2 * eps} calibrates the rank tolerance).
#' @param strength_threshold \eqn{|\theta|/SE} below which a parameter is
#'   listed as weak in a regime (informational).
#' @param strength_ratio_threshold Max/min cross-regime strength ratio above
#'   which a parameter is listed as regime-dependent (informational).
#' @param rare_prob       Regimes with weight below this are "rare"; a
#'   parameter identified only through rare regimes WARNs. \strong{The default
#'   0.05 is a package choice with no literature source}: neither Iskrev (2010)
#'   nor the OccBin literature (Guerrieri & Iacoviello 2015; Cuba-Borda et al.
#'   2019) proposes any numeric rare-regime cutoff. It borrows the familiar
#'   5\% tail convention and is meant to be adjusted per application.
#' @param verbose         Print progress messages.
#' @param meta            Optional plot metadata.
#'
#' @return A \code{dynhr_diagnostic} whose \code{result} holds
#'   \code{regime_jacobians}, \code{regime_strength} (n_param x (regimes +
#'   pooled)), \code{regime_ranks}, \code{regime_probs},
#'   \code{regime_unidentified}, \code{pooled_rank},
#'   \code{pooled_unidentified}, \code{rare_only_params},
#'   \code{regime_dependent_params}, \code{weak_params_by_regime},
#'   \code{regime_labels}, \code{regime_indices}, \code{n_obc_specs}.
#' @noRd
d27_obc_identification <- function(model,
                                    params,
                                    compiled = NULL,
                                    obc_specs = NULL,
                                    regime_subset = NULL,
                                    regime_probs = NULL,
                                    param_names = NULL,
                                    obs_names = NULL,
                                    eps = 1e-5,
                                    strength_threshold = 1.0,
                                    strength_ratio_threshold = 3.0,
                                    rare_prob = 0.05,
                                    verbose = FALSE,
                                    meta = NULL) {
  # ---- 1. Validate inputs ----
  if (!inherits(model, "dynhr_mod")) {
    return(.make_result(
      pass    = NA,
      summary = "D27 OBC Identification: model must be a dynhr_mod object."
    ))
  }
  if (is.null(params) || is.null(names(params))) {
    return(.make_result(
      pass    = NA,
      summary = "D27 OBC Identification: a named params vector is required."
    ))
  }
  params <- unlist(params)
  param_names <- param_names %||% names(params)
  if (!all(param_names %in% names(params))) {
    .dynhr_abort("d27: param_names not in names(params): ",
                 paste(setdiff(param_names, names(params)), collapse = ", "))
  }
  n_par <- length(param_names)

  fixed_specs <- obc_specs
  specs0 <- fixed_specs %||% .d27_parse_specs(model, params)
  k <- length(specs0)
  if (k == 0L) {
    return(.make_result(
      pass    = NA,
      summary = paste(
        "D27 OBC Identification: No OBC constraints detected or provided.",
        "Add [mcp = 'var > bound'] equation tags or pass obc_specs."
      ),
      llm_summary = paste(
        "[INFO] D27 OBC Identification status=skipped reason=no_obc_specs",
        "action: define OBC constraints via mcp tags"
      )
    ))
  }
  n_regimes_total <- 2L ^ k
  regime_subset <- regime_subset %||% c(0L, 2L ^ (seq_len(k) - 1L))
  regime_subset <- unique(as.integer(regime_subset))
  if (anyNA(regime_subset) || any(regime_subset < 0L | regime_subset >= n_regimes_total)) {
    .dynhr_abort(sprintf("d27: regime_subset must lie in [0, %d] for %d constraint(s).",
                         n_regimes_total - 1L, k))
  }
  if (!is.null(regime_probs) &&
      (length(regime_probs) != length(regime_subset) ||
       any(!is.finite(regime_probs)) || any(regime_probs < 0))) {
    .dynhr_abort("d27: regime_probs must be finite, >= 0, one per regime in regime_subset.")
  }
  n_regimes <- length(regime_subset)

  spec_names <- vapply(specs0, function(s) s$name %||% s$var_name %||% "obc",
                       character(1))
  regime_labels <- vapply(regime_subset, function(idx) {
    if (idx == 0L) return("all slack")
    paste(sprintf("%s binds", spec_names[obc_regime_flags(idx, k)]), collapse = " + ")
  }, character(1))

  compiled <- compiled %||% compile_model(model, verbose = FALSE)
  sc <- cache_system_structure(compiled)
  var_names <- model$var_names
  obs_names <- obs_names %||% var_names
  if (!all(obs_names %in% var_names)) {
    .dynhr_abort("d27: obs_names not among the endogenous variables: ",
                 paste(setdiff(obs_names, var_names), collapse = ", "))
  }
  sel <- match(obs_names, var_names)

  if (verbose) .dynhr_cat(sprintf("[d27] %d OBC spec(s); regimes: %s\n", k,
                                  paste(regime_labels, collapse = ", ")))

  # ---- 2. Moments for every regime at one parameter point ----
  # One model solve per point serves all regimes; a failed regime gives NA
  # rows (fixed length), a failed solve gives an all-NA vector.
  n_sel <- length(sel)
  n_mom <- 3L * n_sel + n_sel * (n_sel - 1L) / 2L
  eval_point <- function(th) {
    p_full <- params
    p_full[names(th)] <- th
    sol <- .d27_solve_point(model, compiled, sc, p_full, fixed_specs)
    out <- rep(NA_real_, n_mom * n_regimes)
    if (is.null(sol)) return(out)
    for (r in seq_len(n_regimes)) {
      m <- .d27_regime_moments(sol, regime_subset[r], sel)
      if (!is.null(m)) out[(r - 1L) * n_mom + seq_len(n_mom)] <- m$moments
    }
    out
  }

  theta0 <- params[param_names]
  sol0 <- .d27_solve_point(model, compiled, sc, params, fixed_specs)
  if (is.null(sol0)) {
    return(.make_result(
      pass = NA,
      summary = paste("D27 OBC Identification: the slack model could not be solved",
                      "(steady state or Blanchard-Kahn failure) at params.")
    ))
  }
  slack0 <- .d27_regime_moments(sol0, 0L, sel)
  if (is.null(slack0)) {
    return(.make_result(
      result = list(regime_labels = regime_labels, regime_indices = regime_subset),
      pass = NA,
      summary = paste("D27 OBC Identification: the slack regime has no finite",
                      "stationary moments at params (unit/explosive root).")
    ))
  }
  mom0 <- lapply(regime_subset, function(idx) .d27_regime_moments(sol0, idx, sel))
  mom_names <- .d27_moment_names(obs_names)

  # ---- 3. Regime weights ----
  prob_source <- if (is.null(regime_probs)) "shadow" else "user"
  if (is.null(regime_probs)) {
    regime_probs <- .d27_shadow_probs(slack0, sol0, regime_subset, k)
  }
  regime_probs <- as.numeric(regime_probs)
  names(regime_probs) <- regime_labels

  # ---- 4. Jacobians (step h and 2h, regime held fixed) ----
  J_all  <- .numerical_jacobian(eval_point, theta0, eps = eps)
  J2_all <- .numerical_jacobian(eval_point, theta0, eps = 2 * eps)

  regime_jacobians <- vector("list", n_regimes)
  regime_J2 <- vector("list", n_regimes)
  names(regime_jacobians) <- regime_labels
  regime_ok <- logical(n_regimes)
  for (r in seq_len(n_regimes)) {
    rows <- (r - 1L) * n_mom + seq_len(n_mom)
    J  <- matrix(J_all[rows, ],  n_mom, n_par, dimnames = list(mom_names, param_names))
    J2 <- matrix(J2_all[rows, ], n_mom, n_par, dimnames = list(mom_names, param_names))
    regime_ok[r] <- !is.null(mom0[[r]]) && all(is.finite(J))
    if (regime_ok[r]) {
      regime_jacobians[[r]] <- J
      regime_J2[[r]] <- if (all(is.finite(J2))) J2
    }
  }
  if (!any(regime_ok)) {
    return(.make_result(
      result = list(regime_labels = regime_labels, regime_indices = regime_subset),
      pass = NA,
      summary = paste("D27 OBC Identification: no regime has a finite Jacobian",
                      "(failed solve or non-stationary policy within eps of params).")
    ))
  }

  # ---- 5. Rank + strength per regime and pooled ----
  cols <- c(regime_labels, "pooled")
  regime_strength <- matrix(NA_real_, n_par, n_regimes + 1L,
                            dimnames = list(param_names, cols))
  regime_ranks <- rep(NA_integer_, n_regimes)
  names(regime_ranks) <- regime_labels
  regime_unid <- stats::setNames(vector("list", n_regimes), regime_labels)
  for (r in which(regime_ok)) {
    rk <- .ident_equilibrated_rank(regime_jacobians[[r]], regime_J2[[r]])
    regime_ranks[r] <- rk$rank
    regime_unid[[r]] <- rk$unidentified_params
    regime_strength[, r] <- .d27_strength(regime_jacobians[[r]], rk, theta0)
  }

  pool <- function(use) {
    use <- use[regime_ok[use]]
    if (length(use) == 0L) return(NULL)
    Jw  <- do.call(rbind, lapply(use, function(r) sqrt(regime_probs[r]) * regime_jacobians[[r]]))
    J2_ok <- all(vapply(regime_J2[use], Negate(is.null), logical(1)))
    J2w <- if (J2_ok) do.call(rbind, lapply(use, function(r) sqrt(regime_probs[r]) * regime_J2[[r]]))
    list(Jw = Jw, rk = .ident_equilibrated_rank(Jw, J2w))
  }
  probs_ok <- all(is.finite(regime_probs))
  pooled <- if (probs_ok) pool(which(regime_probs > 0)) else NULL
  common <- if (probs_ok) pool(which(regime_probs >= rare_prob)) else NULL
  pooled_unid <- if (!is.null(pooled)) pooled$rk$unidentified_params else NA_character_
  pooled_rank <- if (!is.null(pooled)) pooled$rk$rank else NA_integer_
  rare_only <- character(0)
  if (!is.null(pooled)) {
    regime_strength[, "pooled"] <- .d27_strength(pooled$Jw, pooled$rk, theta0)
    common_unid <- if (is.null(common)) param_names else common$rk$unidentified_params
    rare_only <- setdiff(common_unid, pooled_unid)
  }

  # ---- 6. Cross-regime (informational) ----
  regime_dependent <- character(0)
  weak_by_regime <- list()
  for (pname in param_names) {
    s_i <- regime_strength[pname, seq_len(n_regimes)]
    weak <- which(is.finite(s_i) & s_i < strength_threshold)
    if (length(weak) > 0L) weak_by_regime[[pname]] <- regime_labels[weak]
    s_fin <- s_i[is.finite(s_i)]
    if (length(s_fin) >= 2L &&
        (min(s_fin) <= 0 && max(s_fin) > 0 ||
         min(s_fin) > 0 && max(s_fin) / min(s_fin) > strength_ratio_threshold)) {
      regime_dependent <- c(regime_dependent, pname)
    }
  }

  ## A parameter no regime identifies is a hard identification failure. A
  ## parameter identified only through a rare regime is a weak-identification /
  ## data-informativeness statement -- the OBC literature treats it the way the
  ## IV literature treats weak instruments, qualitatively -- so it WARNs.
  pass <- if (is.null(pooled)) NA else length(pooled_unid) == 0L
  warn <- isTRUE(pass) && length(rare_only) > 0L

  # ---- 7. Plot ----
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plots$regime_strength <- .apply_meta(
      .d27_plot_strength(regime_strength, regime_probs, regime_unid,
                         pooled_unid, rare_only, prob_source, k),
      meta)
  }

  result <- list(
    regime_jacobians        = regime_jacobians,
    regime_strength         = regime_strength,
    regime_ranks            = regime_ranks,
    regime_probs            = regime_probs,
    prob_source             = prob_source,
    regime_unidentified     = regime_unid,
    pooled_rank             = pooled_rank,
    pooled_unidentified     = pooled_unid,
    rare_only_params        = rare_only,
    regime_dependent_params = regime_dependent,
    weak_params_by_regime   = weak_by_regime,
    regime_labels           = regime_labels,
    regime_indices          = stats::setNames(regime_subset, regime_labels),
    n_obc_specs             = k,
    n_params                = n_par
  )

  fmt_list <- function(x) if (length(x)) paste(x, collapse = ", ") else "none"
  regime_str <- paste(sprintf("%s: p=%.3g rank=%s", regime_labels, regime_probs,
                              ifelse(is.na(regime_ranks), "failed",
                                     sprintf("%d/%d", regime_ranks, n_par))),
                      collapse = "; ")
  verdict <- if (is.na(pass)) {
    "Regime weights unavailable; pooled identification not assessed."
  } else if (warn) {
    sprintf(paste0("Pooled (probability-weighted) rank %d/%d: all parameters ",
                   "identified, but %s identified only via regimes with ",
                   "p < %.3g -- weakly informed by the data."),
            pooled_rank, n_par, fmt_list(rare_only), rare_prob)
  } else if (pass) {
    sprintf("Pooled (probability-weighted) rank %d/%d: all parameters identified.",
            pooled_rank, n_par)
  } else {
    paste0(sprintf("Pooled rank %d/%d.", pooled_rank, n_par),
           if (length(pooled_unid)) sprintf(" Unidentified in every regime: %s.", fmt_list(pooled_unid)),
           if (length(rare_only)) sprintf(" Identified only via regimes with p < %.3g: %s.",
                                          rare_prob, fmt_list(rare_only)))
  }
  summary_str <- sprintf(
    "D27 OBC Identification: %d spec(s), %d regime(s) (%s weights) [%s]. %s Regime-dependent: %s.",
    k, n_regimes, prob_source, regime_str, verdict, fmt_list(regime_dependent))

  .make_result(
    result  = result,
    pass    = pass,
    warn    = warn,
    plots   = plots,
    summary = summary_str,
    llm_summary = sprintf(
      "[%s] D27 OBC Identification n_obc_specs=%d n_regimes=%d pooled_rank=%s/%d ranks=[%s] probs=[%s] pooled_unidentified=%s rare_only=%s regime_dependent=%s",
      .badge_str(list(pass = pass, errored = FALSE, warn = warn)),
      k, n_regimes, pooled_rank, n_par,
      paste(regime_ranks, collapse = ","),
      paste(sprintf("%.3g", regime_probs), collapse = ","),
      fmt_list(pooled_unid), fmt_list(rare_only), fmt_list(regime_dependent)
    )
  )
}


# ==========================================================================
# Internal helpers for D27
# ==========================================================================

#' Parse MCP tags with the model's parameter values set to `params`
#'
#' A bound given by a parameter name then follows that parameter. Returns an
#' empty list when the model has no MCP tags.
#' @noRd
.d27_parse_specs <- function(model, params) {
  tags <- vapply(model$equations, function(e) {
    t <- e$tag_raw %||% e$tag
    if (is.null(t) || is.na(t)) "" else as.character(t)
  }, character(1))
  if (!any(grepl("mcp", tags, fixed = TRUE))) return(list())
  pv <- model$param_values
  common <- intersect(names(params), names(pv))
  pv[common] <- params[common]
  model$param_values <- pv
  obc_parse_tags(model)
}

#' Solve the model (steady state, slack policy, system matrices) at a point
#'
#' @return list(dr, sys, specs, Sigma_e, ys) or NULL when the steady state
#'   fails or Blanchard-Kahn does not hold.
#' @noRd
.d27_solve_point <- function(model, compiled, sc, p_full, fixed_specs) {
  ss <- solve_steady(compiled, params = p_full, verbose = FALSE)
  if (is.null(ss) || !isTRUE(ss$converged)) return(NULL)
  dr <- solve_perturbation(model, compiled, ss$values, p_full, order = 1L,
                           verbose = FALSE)
  if (is.null(dr) || !isTRUE(dr$bk_satisfied)) return(NULL)
  ys <- if (inherits(dr$ys, "dynhr_steady")) dr$ys$values else dr$ys
  list(
    dr      = dr,
    sys     = extract_system_matrices_fast(sc, ss$values, p_full),
    specs   = fixed_specs %||% .d27_parse_specs(model, p_full),
    # Sigma_e from the model's shock block at THIS point (shock-std parameters
    # move it); never an identity default.
    Sigma_e = .get_shock_cov(model, dr$exo_names, p_full),
    ys      = as.numeric(ys)
  )
}

#' Stationary moments of one regime's linear policy
#'
#' Policy y_t = ys + c + G s_{t-1} + H e_t with s = y[state_idx] (deviation).
#' @param sel Indices of the variables whose moments are returned.
#' @return list(moments, mean, V) or NULL (singular binding system or
#'   non-stationary regime policy).
#' @noRd
.d27_regime_moments <- function(sol, regime_idx, sel) {
  dr <- sol$dr
  si <- dr$state_idx
  n  <- nrow(dr$ghx)
  if (regime_idx == 0L) {
    G <- dr$ghx; H <- dr$ghu; cc <- numeric(n)
  } else {
    flags <- obc_regime_flags(regime_idx, length(sol$specs))
    b <- obc_solve_binding(sol$sys, dr, sol$specs[flags], obs_idx = sel)
    if (is.null(b)) return(NULL)
    G <- b$dr$ghx; H <- b$dr$ghu; cc <- as.numeric(b$c_full)
  }
  G <- matrix(G, n); H <- matrix(H, n)
  ns <- length(si)
  if (ns > 0L) {
    Ts <- G[si, , drop = FALSE]
    if (max(Mod(eigen(Ts, only.values = TRUE)$values)) >= 1 - 1e-8) return(NULL)
    Rs <- H[si, , drop = FALSE]
    Vs <- solve_lyapunov(Ts, Rs %*% sol$Sigma_e %*% t(Rs))
    ms <- solve(diag(ns) - Ts, cc[si])
    mu <- sol$ys + cc + as.numeric(G %*% ms)
    V  <- G %*% Vs %*% t(G) + H %*% sol$Sigma_e %*% t(H)
    G1 <- G %*% V[si, , drop = FALSE]          # Cov(y_t, y_{t-1})
  } else {
    mu <- sol$ys + cc
    V  <- H %*% sol$Sigma_e %*% t(H)
    G1 <- matrix(0, n, n)
  }
  Vsel <- V[sel, sel, drop = FALSE]
  m <- c(mu[sel], diag(Vsel), Vsel[upper.tri(Vsel)],
         diag(G1[sel, sel, drop = FALSE]))
  if (!all(is.finite(m))) return(NULL)
  list(moments = m, mean = mu, V = V)
}

#' Moment labels matching .d27_regime_moments() order
#' @noRd
.d27_moment_names <- function(nm) {
  pairs <- outer(nm, nm, paste, sep = ",")
  c(paste0("mean(", nm, ")"), paste0("var(", nm, ")"),
    paste0("cov(", pairs[upper.tri(pairs)], ")"), paste0("acov1(", nm, ")"))
}

#' Shadow probability of each regime
#'
#' Under the slack regime's stationary Gaussian distribution, the probability
#' that exactly the regime's constraints are violated. Deterministic
#' (mvtnorm Miwa algorithm; no RNG draw).
#' @param slack Output of .d27_regime_moments() for regime 0 (full mean / V).
#' @noRd
.d27_shadow_probs <- function(slack, sol, regime_subset, k) {
  specs <- sol$specs
  vi  <- vapply(specs, function(s) s$var_idx, integer(1))
  sgn <- vapply(specs, function(s) if (identical(s$op, ">")) -1 else 1, numeric(1))
  bnd <- vapply(specs, function(s) as.numeric(s$bound), numeric(1))
  # u_j = sgn_j * (x_j - bound_j) > 0  <=>  constraint j violated
  mu_u <- sgn * (slack$mean[vi] - bnd)
  V_u  <- (sgn %o% sgn) * slack$V[vi, vi, drop = FALSE]
  sd_u <- sqrt(pmax(diag(V_u), 0))
  vapply(regime_subset, function(idx) {
    viol <- obc_regime_flags(idx, k)
    det <- sd_u <= 1e-12 * max(1, abs(mu_u))
    # Constraints with no variance are violated deterministically or never.
    if (any(det & (viol != (mu_u > 0)))) return(0)
    st <- !det
    if (!any(st)) return(1)
    lo <- ifelse(viol[st], 0, -Inf)
    hi <- ifelse(viol[st], Inf, 0)
    if (sum(st) == 1L) {
      return(stats::pnorm(hi, mu_u[st], sd_u[st]) - stats::pnorm(lo, mu_u[st], sd_u[st]))
    }
    as.numeric(mvtnorm::pmvnorm(lower = lo, upper = hi, mean = mu_u[st],
                                sigma = V_u[st, st, drop = FALSE],
                                algorithm = mvtnorm::Miwa()))
  }, numeric(1))
}

#' |theta| / SE from I = J'J on the identified subspace (D20 "none" weighting)
#'
#' Parameters loading on an unidentified direction get strength 0.
#' @noRd
.d27_strength <- function(J, rk, theta) {
  n_par <- ncol(J)
  out <- rep(0, n_par)
  names(out) <- colnames(J)
  if (rk$rank == 0L) return(out)
  cs <- sqrt(colSums(J^2))
  cs[cs <= 0] <- 1
  sj <- svd(sweep(J, 2, cs, "/"), nu = 0, nv = n_par)
  q <- rk$rank
  Vq <- sj$v[, seq_len(q), drop = FALSE]
  se <- sqrt(pmax(rowSums(sweep(Vq, 2, sj$d[seq_len(q)], "/")^2), 0)) / cs
  out <- abs(as.numeric(theta)) / se
  names(out) <- colnames(J)
  unid <- union(rk$unidentified_params,
                colnames(J)[rowSums(rk$null_space^2) >= 0.01])
  out[unid] <- 0
  out
}

#' Heatmap of identification strength by regime (+ pooled column)
#' @noRd
.d27_plot_strength <- function(S, probs, regime_unid, pooled_unid, rare_only,
                               prob_source, k) {
  col_lab <- c(sprintf("%s\np = %.3g", colnames(S)[-ncol(S)], probs),
               "pooled\n(p-weighted)")
  unid <- c(regime_unid, list(pooled = pooled_unid))
  df <- data.frame(
    Parameter = factor(rep(rownames(S), ncol(S)), levels = rev(rownames(S))),
    Regime    = factor(rep(col_lab, each = nrow(S)), levels = col_lab),
    Strength  = as.vector(S),
    stringsAsFactors = FALSE
  )
  unid_flag <- unlist(lapply(seq_len(ncol(S)), function(j)
    rownames(S) %in% (unid[[j]] %||% character(0))))
  df$label <- ifelse(is.na(df$Strength), "failed",
              ifelse(unid_flag, "unid.", formatC(df$Strength, digits = 2, format = "g")))
  df$fill <- ifelse(is.finite(df$Strength) & df$Strength > 0, log10(df$Strength), NA_real_)
  # Light text on the dark (low) half of cividis, dark text on the light half.
  rng <- range(df$fill, na.rm = TRUE)
  pos <- if (all(is.finite(rng)) && diff(rng) > 0) (df$fill - rng[1]) / diff(rng) else 0
  df$txt <- ifelse(is.finite(df$fill) & pos < 0.55, "white", "grey10")
  df$rare <- df$Regime == col_lab[length(col_lab)] & rownames(S)[
    (seq_len(nrow(df)) - 1L) %% nrow(S) + 1L] %in% rare_only
  subtitle <- paste0(
    sprintf("%d OBC constraint(s); p = %s regime probability.\n", k,
            if (prob_source == "shadow") "shadow (slack-Gaussian)" else "user-supplied"),
    "'unid.' = not identified in that regime (FD-aware rank test).",
    if (length(rare_only)) paste0("\nRed outline: identified only via rare regimes (",
                                  paste(rare_only, collapse = ", "), ")") else "")
  ggplot2::ggplot(df, ggplot2::aes(x = .data$Regime, y = .data$Parameter)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$fill), colour = "white", linewidth = 0.6) +
    ggplot2::geom_tile(data = df[df$rare, , drop = FALSE], fill = NA,
                       colour = tol_vibrant[["red"]], linewidth = 1.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, colour = .data$txt), size = 3.4) +
    ggplot2::scale_colour_identity() +
    ggplot2::geom_vline(xintercept = ncol(S) - 0.5, colour = "grey20", linewidth = 0.8) +
    scale_fill_dynhr_cividis(name = "log10 |t|\n(unit weights)") +
    theme_dynhr_diagnostic() +
    ggplot2::labs(
      title = "D27: identification strength by OBC regime",
      subtitle = subtitle,
      x = NULL, y = NULL
    )
}
