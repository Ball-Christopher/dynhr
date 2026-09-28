## R/diag-obc-d31-scenario.R
## --------------------------------------------------------------------------
## Phase H: D31 — OBC Scenario Comparison (Constrained vs Unconstrained)
##
## Compares a constrained (occasionally-binding) path with the unconstrained
## (always-slack) path of the same scenario.  Either the two paths are
## supplied, or D31 generates them itself from the slack decision rule: one
## impulse, simulated once under the all-slack policy and once through the
## OccBin/Boehl complementarity solver (boehl_solve_regime_path).
##
## Everything is compared in LEVEL units (steady state + deviation), which is
## the unit an OBC bound is written in (obc_solve_binding() uses
## bound - ys).  The bound is drawn at that level in the constrained
## variable's panel, and binding periods are shaded there.
##
## Badge: FAIL when the constrained path violates a bound (or the solver did
## not converge); PASS when every bound is respected; NA when no bound is
## known (supplied paths without constraint information).
## --------------------------------------------------------------------------

#' D31. OBC Scenario Comparison (Constrained vs Unconstrained)
#'
#' Compares a constrained path with the unconstrained path of the same
#' scenario.  Rows are periods, columns are variables.  Supplied paths are
#' aligned by column name (and must have identical row names / length when
#' both carry them); they are never recycled.
#'
#' Metrics are per variable (max |constrained - unconstrained|, the period
#' where it occurs, and that max relative to the variable's own unconstrained
#' range), because a single max across variables with different units is not
#' meaningful.  Per constraint: periods where the unconstrained path violates
#' the bound, periods where the constrained path sits at the bound, and the
#' largest violation of the bound by the constrained path.
#'
#' @param pf_constrained    T x n matrix / data frame — constrained path
#'   (rows = periods).  A single \code{Date} column is used as \code{dates}.
#' @param pf_unconstrained  T x n matrix / data frame — unconstrained path.
#' @param var_names         Optional variable names (when paths are unnamed).
#' @param constraint_name   Label for the constraint (e.g. "ZLB").
#' @param constraint_var    Character vector: constrained variable(s) in the
#'   supplied paths.  Defaults to the \code{obc_specs} variables.
#' @param bound             Numeric vector (same length as
#'   \code{constraint_var}), in the SAME units as the supplied paths.
#' @param op                ">" (lower bound) or "<" (upper bound), recycled
#'   to \code{constraint_var}.
#' @param dates             Optional x-axis values (length T).
#' @param sys,dr_slack      System matrices and the slack-regime rule solved
#'   at the SAME parameters.  If either is NULL both are rebuilt from
#'   \code{compiled} at \code{params} (overlaid on the calibration).
#' @param obc_specs         OBC spec list (auto-generation; also supplies the
#'   bound for supplied paths, which are then taken to be LEVELS).
#' @param obs_idx           Observable indices (passed to the OBC solver).
#' @param model             dynhr_mod (shock standard deviations).
#' @param compiled          Compiled model (to rebuild \code{sys}/\code{dr_slack}).
#' @param params            Parameter values (may be a subset, e.g. the
#'   estimated ones); also used for the shock standard deviations.
#' @param shock             Shock to impulse (name or index).  NULL: the
#'   shock that moves the first constrained variable furthest toward its
#'   bound per standard deviation.
#' @param n_periods         Horizon of the generated scenario (default 20).
#' @param shock_scale       Impulse size in standard deviations.  NULL: the
#'   larger of 1 and the size at which the unconstrained path overshoots the
#'   bound by the steady-state distance to it (so the constraint binds).
#' @param tol               Relative tolerance for "at the bound" / "no
#'   effect" (scaled by the path magnitude).
#' @param meta              Optional dynhr_diag_meta.
#' @return dynhr_diagnostic list.
#' @noRd
d31_obc_scenario_comparison <- function(pf_constrained    = NULL,
                                         pf_unconstrained  = NULL,
                                         var_names         = NULL,
                                         constraint_name   = "OBC",
                                         constraint_var    = NULL,
                                         bound             = NULL,
                                         op                = ">",
                                         dates             = NULL,
                                         sys               = NULL,
                                         dr_slack          = NULL,
                                         obc_specs         = NULL,
                                         obs_idx           = NULL,
                                         model             = NULL,
                                         compiled          = NULL,
                                         params            = NULL,
                                         shock             = NULL,
                                         n_periods         = 20L,
                                         shock_scale       = NULL,
                                         tol               = 1e-8,
                                         meta              = NULL) {

  have_c <- !is.null(pf_constrained)
  have_u <- !is.null(pf_unconstrained)
  if (xor(have_c, have_u)) {
    .dynhr_abort("D31: supply both `pf_constrained` and `pf_unconstrained` ",
                 "(or neither, to generate the scenario).")
  }

  scenario <- NULL
  y_label  <- "Value"
  if (!have_c) {
    gen <- .d31_generate(dr_slack = dr_slack, sys = sys, obc_specs = obc_specs,
                         obs_idx = obs_idx, model = model, compiled = compiled,
                         params = params, shock = shock,
                         n_periods = n_periods, shock_scale = shock_scale)
    if (is.null(gen$paths_c)) {
      return(.make_result(
        result  = NULL,
        pass    = NA,
        plots   = list(),
        summary = sprintf("D31 OBC scenario comparison (%s): not run -- %s",
                          constraint_name, gen$reason),
        llm_summary = sprintf(
          "[INFO] D31 OBC scenario status=not_run reason=%s", gen$reason_code)
      ))
    }
    pf_c <- gen$paths_c
    pf_u <- gen$paths_u
    constraints <- gen$constraints
    scenario <- gen$scenario
    y_label  <- "Level (steady state + deviation)"
  } else {
    al <- .d31_align_paths(pf_constrained, pf_unconstrained, var_names, dates)
    pf_c  <- al$c
    pf_u  <- al$u
    dates <- al$dates
    constraints <- .d31_constraints_supplied(colnames(pf_c), constraint_var,
                                             bound, op, obc_specs)
  }

  n_T  <- nrow(pf_c)
  vars <- colnames(pf_c)
  diff_mat <- pf_c - pf_u
  scale <- max(1, abs(pf_c), abs(pf_u))
  eps_t <- tol * scale

  # ---- per-variable effect (each in its own units) ----
  max_abs  <- apply(abs(diff_mat), 2, max)
  t_max    <- apply(abs(diff_mat), 2, which.max)
  u_range  <- apply(abs(pf_u), 2, max)
  rel      <- ifelse(u_range > eps_t, max_abs / u_range, NA_real_)
  affected <- max_abs > eps_t
  effects <- data.frame(variable = vars, max_abs_effect = unname(max_abs),
                        period_of_max = unname(t_max),
                        rel_effect = unname(rel), affected = unname(affected),
                        stringsAsFactors = FALSE)

  # ---- per-constraint binding / violation ----
  cons_tab <- NULL
  binding  <- list()
  if (length(constraints$var) > 0L) {
    rows <- lapply(seq_along(constraints$var), function(j) {
      v <- constraints$var[j]; b <- constraints$bound[j]; o <- constraints$op[j]
      sgn <- if (o == ">") 1 else -1          # feasible side: sgn*(x - b) >= 0
      gap_c <- sgn * (pf_c[, v] - b)
      gap_u <- sgn * (pf_u[, v] - b)
      at_b  <- abs(gap_c) <= eps_t
      list(tab = data.frame(constraint = constraints$label[j], variable = v,
                            op = o, bound = b,
                            n_unconstrained_violations = sum(gap_u < -eps_t),
                            n_at_bound = sum(at_b),
                            first_at_bound = if (any(at_b)) which(at_b)[1] else NA_integer_,
                            last_at_bound = if (any(at_b)) max(which(at_b)) else NA_integer_,
                            max_violation = max(0, -gap_c),
                            stringsAsFactors = FALSE),
           at_b = at_b)
    })
    cons_tab <- do.call(rbind, lapply(rows, `[[`, "tab"))
    binding  <- stats::setNames(lapply(rows, `[[`, "at_b"), constraints$label)
  }

  binds_any <- !is.null(cons_tab) &&
    any(cons_tab$n_unconstrained_violations > 0L | cons_tab$n_at_bound > 0L)
  violated  <- !is.null(cons_tab) && any(cons_tab$max_violation > eps_t)
  converged <- if (is.null(scenario)) NA else isTRUE(scenario$converged)
  # With no violation of the unconstrained path, the all-slack solution is a
  # valid constrained solution, so the paths must coincide.
  spurious  <- !is.null(cons_tab) && !binds_any && any(affected)

  pass <- if (is.null(cons_tab)) NA
          else !violated && !isFALSE(converged) && !spurious

  # ---- text ----
  top <- effects[order(-effects$max_abs_effect), , drop = FALSE]
  top <- top[top$affected, , drop = FALSE]
  eff_txt <- if (nrow(top) == 0L) "no variable path is altered"
    else paste(sprintf("%s %.4g (t=%d)", utils::head(top$variable, 4L),
                       utils::head(top$max_abs_effect, 4L),
                       utils::head(top$period_of_max, 4L)), collapse = ", ")
  cons_txt <- if (is.null(cons_tab)) "no bound supplied, binding not assessed"
    else paste(sprintf("%s%s %s %.4g: unconstrained violates in %d/%d periods, constrained at bound in %d%s, max violation %.2g",
                       ifelse(cons_tab$constraint == constraints$label_default, "",
                              paste0(cons_tab$constraint, " ")),
                       cons_tab$variable, cons_tab$op,
                       cons_tab$bound, cons_tab$n_unconstrained_violations, n_T,
                       cons_tab$n_at_bound,
                       ifelse(is.na(cons_tab$first_at_bound), "",
                              sprintf(" [t=%d..%d]", cons_tab$first_at_bound,
                                      cons_tab$last_at_bound)),
                       cons_tab$max_violation), collapse = "; ")
  scen_txt <- if (is.null(scenario)) ""
    else sprintf(" Scenario: %+.3g sd %s impulse, %d periods%s.",
                 scenario$shock_sd, scenario$shock, n_T,
                 if (isTRUE(scenario$converged)) "" else ", SOLVER NOT CONVERGED")
  verdict <- if (is.na(pass)) ""
    else if (violated) " FAIL: constrained path violates the bound."
    else if (isFALSE(converged)) " FAIL: OBC solver did not converge."
    else if (spurious) " FAIL: constraint never binds yet paths differ."
    else if (!binds_any) " Constraint never binds in this scenario."
    else ""
  summary <- sprintf("D31 OBC scenario comparison (%s):%s Largest effects: %s. %s.%s",
                     constraint_name, scen_txt, eff_txt, cons_txt, verdict)
  llm <- sprintf(
    "[%s] D31 OBC scenario | constraint=%s binds=%s n_affected=%d/%d max_violation=%s converged=%s",
    .badge_str(list(pass = pass, errored = FALSE)), constraint_name,
    if (is.null(cons_tab)) "NA" else binds_any, sum(affected), length(vars),
    if (is.null(cons_tab)) "NA" else sprintf("%.3g", max(cons_tab$max_violation)),
    converged)

  plots <- .d31_plots(pf_c, pf_u, diff_mat, effects, constraints, binding,
                      dates, constraint_name, y_label, scenario, meta)

  .make_result(
    result  = list(paths_constrained = pf_c, paths_unconstrained = pf_u,
                   diff_mat = diff_mat, effects = effects,
                   constraints = cons_tab, binds = binds_any,
                   n_vars_affected = sum(affected), scenario = scenario),
    pass    = pass,
    plots   = plots,
    summary = summary,
    llm_summary = llm
  )
}


# ---------------------------------------------------------------------------
# Supplied paths: validate and align (never recycle)
# ---------------------------------------------------------------------------
.d31_align_paths <- function(pc, pu, var_names, dates) {
  prep <- function(p, nm) {
    if (is.data.frame(p)) {
      is_date <- vapply(p, function(col) inherits(col, c("Date", "POSIXt")), logical(1))
      d <- NULL
      if (any(is_date)) {
        if (sum(is_date) > 1L)
          .dynhr_abort("D31: `", nm, "` has more than one date column.")
        d <- p[[which(is_date)]]
        p <- p[, !is_date, drop = FALSE]
      }
      num <- vapply(p, is.numeric, logical(1))
      if (!all(num))
        .dynhr_abort("D31: `", nm, "` has non-numeric column(s): ",
                     paste(names(p)[!num], collapse = ", "), ".")
      m <- as.matrix(p)
      if (is.null(d)) rn <- if (is.character(attr(p, "row.names"))) rownames(p) else NULL
      else rn <- as.character(d)
      rownames(m) <- rn
      return(list(m = m, dates = d))
    }
    m <- as.matrix(p)
    if (!is.numeric(m)) .dynhr_abort("D31: `", nm, "` must be numeric.")
    list(m = m, dates = NULL)
  }
  a <- prep(pc, "pf_constrained")
  b <- prep(pu, "pf_unconstrained")
  mc <- a$m; mu <- b$m
  if (nrow(mc) != nrow(mu)) {
    .dynhr_abort(sprintf(paste0(
      "D31: the paths have different lengths (constrained %d periods, ",
      "unconstrained %d). Align them to a common period range first."),
      nrow(mc), nrow(mu)))
  }
  if (!is.null(rownames(mc)) && !is.null(rownames(mu)) &&
      !identical(rownames(mc), rownames(mu))) {
    .dynhr_abort("D31: the paths cover different periods (row names / dates ",
                 "differ). Align them to a common period range first.")
  }
  cn_c <- colnames(mc); cn_u <- colnames(mu)
  if (!is.null(cn_c) && !is.null(cn_u)) {
    if (!setequal(cn_c, cn_u) || anyDuplicated(cn_c) || anyDuplicated(cn_u)) {
      .dynhr_abort("D31: the paths have different variables: ",
                   paste(union(setdiff(cn_c, cn_u), setdiff(cn_u, cn_c)),
                         collapse = ", "), ".")
    }
    mu <- mu[, cn_c, drop = FALSE]
  } else {
    if (ncol(mc) != ncol(mu)) {
      .dynhr_abort(sprintf("D31: the paths have %d and %d columns.",
                           ncol(mc), ncol(mu)))
    }
    nm <- if (!is.null(var_names)) var_names else if (!is.null(cn_c)) cn_c
          else if (!is.null(cn_u)) cn_u else paste0("var_", seq_len(ncol(mc)))
    if (length(nm) != ncol(mc))
      .dynhr_abort("D31: `var_names` has length ", length(nm), " but the paths have ",
                   ncol(mc), " columns.")
    colnames(mc) <- nm; colnames(mu) <- nm
  }
  if (nrow(mc) == 0L || ncol(mc) == 0L)
    .dynhr_abort("D31: the supplied paths are empty.")
  if (!all(is.finite(mc)) || !all(is.finite(mu)))
    .dynhr_abort("D31: the supplied paths contain non-finite values.")
  if (is.null(dates)) {
    dates <- if (!is.null(a$dates)) a$dates else b$dates
    if (!is.null(a$dates) && !is.null(b$dates) && !identical(a$dates, b$dates))
      .dynhr_abort("D31: the paths have different date columns.")
  }
  if (!is.null(dates) && length(dates) != nrow(mc))
    .dynhr_abort(sprintf("D31: `dates` has length %d, paths have %d periods.",
                         length(dates), nrow(mc)))
  rownames(mc) <- NULL; rownames(mu) <- NULL
  list(c = mc, u = mu, dates = dates)
}

.d31_label <- function(var, op, bound) {
  sprintf("%s %s %s", var, op, vapply(bound, format, character(1), digits = 4))
}

.d31_constraints_supplied <- function(vars, constraint_var, bound, op, obc_specs) {
  if (is.null(constraint_var) && is.null(bound) && length(obc_specs) > 0L) {
    constraint_var <- vapply(obc_specs, function(s) s$var_name, character(1))
    bound <- vapply(obc_specs, function(s) as.numeric(s$bound), numeric(1))
    op    <- vapply(obc_specs, function(s) s$op, character(1))
  }
  if (is.null(constraint_var) && is.null(bound))
    return(list(var = character(0), bound = numeric(0), op = character(0),
                label = character(0)))
  if (is.null(constraint_var) || is.null(bound) ||
      length(bound) != length(constraint_var)) {
    .dynhr_abort("D31: `constraint_var` and `bound` must be supplied together ",
                 "with equal lengths.")
  }
  if (!all(is.finite(bound))) .dynhr_abort("D31: `bound` must be finite.")
  miss <- setdiff(constraint_var, vars)
  if (length(miss) > 0L)
    .dynhr_abort("D31: constraint variable(s) not in the paths: ",
                 paste(miss, collapse = ", "), ".")
  op <- rep_len(as.character(op), length(constraint_var))
  op[op %in% c(">", ">=")] <- ">"
  op[op %in% c("<", "<=")] <- "<"
  if (!all(op %in% c(">", "<")))
    .dynhr_abort("D31: `op` must be '>' or '<'.")
  lab <- .d31_label(constraint_var, op, bound)
  list(var = constraint_var, bound = as.numeric(bound), op = op,
       label = lab, label_default = lab)
}


# ---------------------------------------------------------------------------
# Auto-generation: one impulse, all-slack vs OccBin/Boehl complementarity
# ---------------------------------------------------------------------------
.d31_generate <- function(dr_slack, sys, obc_specs, obs_idx, model, compiled,
                          params, shock, n_periods, shock_scale) {
  none <- function(reason, code) list(paths_c = NULL, reason = reason,
                                      reason_code = code)
  if (length(obc_specs) == 0L)
    return(none("needs obc_specs (or supply both paths)", "no_specs"))
  if (is.null(model))
    return(none("needs `model` for the shock standard deviations", "no_model"))
  # `params` may be a subset (e.g. estimated parameters); overlay it on the
  # calibration so the rule, the system and the shock stds share ONE vector.
  full <- model$param_values
  if (!is.null(params)) {
    full[intersect(names(params), names(full))] <-
      params[intersect(names(params), names(full))]
    params <- c(full, params[setdiff(names(params), names(full))])
  } else {
    params <- full
  }
  if (is.null(sys) || is.null(dr_slack)) {
    # Re-solve at `full` rather than pairing a supplied rule with matrices
    # built at other parameter values.
    if (is.null(compiled))
      return(none("needs `sys` + `dr_slack`, or `compiled` to build them",
                  "no_sys"))
    ss_res <- solve_steady_state(model, compiled, full, verbose = FALSE)
    if (!isTRUE(ss_res$converged))
      return(none("steady state did not converge", "ss_failed"))
    sys <- extract_system_matrices_fast(cache_system_structure(compiled),
                                        ss_res$ss, full)
    dr_slack <- .solve_from_system(sys, model, compiled, ss_res$ss, full, FALSE)
    if (is.null(dr_slack) || !isTRUE(dr_slack$bk_satisfied))
      return(none("slack-regime solution violates Blanchard-Kahn", "bk_failed"))
  }
  ys <- if (inherits(dr_slack$ys, "dynhr_steady")) dr_slack$ys$values else dr_slack$ys
  endo <- dr_slack$endo_names
  exo  <- dr_slack$exo_names
  n_T  <- as.integer(n_periods)
  if (length(n_T) != 1L || is.na(n_T) || n_T < 2L)
    .dynhr_abort("D31: `n_periods` must be an integer >= 2.")
  ys_v <- if (is.null(ys)) rep(0, length(endo)) else as.numeric(ys)[seq_along(endo)]
  names(ys_v) <- endo

  cvar  <- vapply(obc_specs, function(s) endo[s$var_idx], character(1))
  cbnd  <- vapply(obc_specs, function(s) as.numeric(s$bound), numeric(1))
  cop   <- vapply(obc_specs, function(s) s$op, character(1))
  lab0 <- .d31_label(cvar, cop, cbnd)
  constraints <- list(
    var = cvar, bound = cbnd, op = cop, label_default = lab0,
    label = vapply(seq_along(obc_specs), function(j)
      obc_specs[[j]]$name %||% lab0[j], character(1)))
  # steady-state distance to each bound on its feasible side
  dist <- ifelse(cop == ">", ys_v[cvar] - cbnd, cbnd - ys_v[cvar])
  if (any(dist < 0)) {
    .dynhr_abort("D31: the steady state violates constraint(s) ",
                 paste(constraints$label[dist < 0], collapse = ", "), ".")
  }

  sd_e <- .get_shock_stderr(model, exo, params)
  cache <- new.env(parent = emptyenv(), hash = TRUE)
  if (is.null(obs_idx)) obs_idx <- seq_along(endo)
  obc_ensure_policy(0L, cache, sys, dr_slack, obc_specs, obs_idx)

  slack_path <- function(k, size) {
    e <- matrix(0, length(exo), n_T)
    e[k, 1L] <- size
    boehl_simulate(e, dr_slack, cache, integer(n_T))$paths
  }
  # movement of the first constrained variable toward its bound per +1 sd
  v1 <- cvar[1]; s1 <- if (cop[1] == ">") -1 else 1
  depth <- vapply(seq_along(exo), function(k) {
    x <- s1 * slack_path(k, sd_e[k])[v1, ]
    c(max(x), max(-x))
  }, numeric(2))
  if (is.null(shock)) {
    k <- which.max(apply(depth, 2, max))
  } else {
    k <- if (is.character(shock)) match(shock, exo) else as.integer(shock)
    if (length(k) != 1L || is.na(k) || k < 1L || k > length(exo))
      .dynhr_abort("D31: `shock` must name or index one of: ",
                   paste(exo, collapse = ", "), ".")
  }
  if (sd_e[k] <= 0)
    return(none(sprintf("shock %s has zero standard deviation", exo[k]),
                "zero_shock_sd"))
  sgn <- if (depth[1, k] >= depth[2, k]) 1 else -1
  dep <- max(depth[, k])
  if (is.null(shock_scale)) {
    shock_scale <- if (dep > 0) max(1, 2 * dist[1] / dep) else 1
  }
  if (!is.numeric(shock_scale) || length(shock_scale) != 1L ||
      !is.finite(shock_scale))
    .dynhr_abort("D31: `shock_scale` must be a finite number.")

  size <- unname(sgn * shock_scale * sd_e[k])
  e <- matrix(0, length(exo), n_T)
  e[k, 1L] <- size
  u <- boehl_simulate(e, dr_slack, cache, integer(n_T))$paths
  sol <- boehl_solve_regime_path(e, dr_slack, sys, obc_specs, obs_idx = obs_idx,
                                 single_spell = FALSE, max_iter = 50L)
  to_level <- function(p) {
    m <- t(p) + matrix(ys_v, n_T, length(endo), byrow = TRUE)
    colnames(m) <- endo
    m
  }
  list(paths_c = to_level(sol$paths), paths_u = to_level(u),
       constraints = constraints,
       scenario = list(shock = exo[k], shock_sd = sgn * shock_scale,
                       shock_size = size, converged = isTRUE(sol$converged),
                       regime_path = sol$regime_path, fb_norm = sol$fb_norm))
}


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

# The bound LEVEL is the point of D31, so keep y-axis values (the compact
# theme hides them); only tighten facet spacing.
.d31_theme <- function(.gg) {
  theme_dynhr() +
    ggplot2::theme(panel.spacing.x = ggplot2::unit(0.8, "lines"),
                   panel.spacing.y = ggplot2::unit(0.5, "lines"),
                   strip.text = ggplot2::element_text(size = ggplot2::rel(0.85),
                                                      face = "bold", family = "sans"))
}

.d31_plots <- function(pf_c, pf_u, diff_mat, effects, constraints, binding,
                       dates, constraint_name, y_label, scenario, meta) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(list())
  .gg <- getNamespace("ggplot2")
  n_T  <- nrow(pf_c)
  xval <- if (!is.null(dates)) dates else seq_len(n_T)
  x_lab <- if (!is.null(dates)) "Date" else "Period"
  # shade whole periods: half the typical spacing either side of a point
  half <- if (n_T > 1L) stats::median(diff(as.numeric(xval))) / 2 else 0.5
  ttl <- if (identical(constraint_name, "OBC")) "OBC"
         else sprintf("OBC (%s)", constraint_name)

  # constrained variables first, then by effect size; at most 12 panels
  ord  <- order(!(effects$variable %in% constraints$var), -effects$max_abs_effect)
  pick <- effects$variable[ord][seq_len(min(12L, nrow(effects)))]
  n_dropped <- nrow(effects) - length(pick)
  lvl <- pick

  plot_df <- data.frame(
    x = rep(xval, 2L * length(pick)),
    Value = c(as.vector(pf_c[, pick, drop = FALSE]),
              as.vector(pf_u[, pick, drop = FALSE])),
    Variable = factor(rep(rep(pick, each = n_T), 2L), levels = lvl),
    Scenario = rep(c("Constrained", "Unconstrained"), each = n_T * length(pick)))

  sub <- if (!is.null(scenario))
    sprintf("%+.3g sd impulse to %s at period 1%s", scenario$shock_sd,
            scenario$shock,
            if (scenario$converged) "" else " -- OBC solver NOT converged")
    else NULL
  notes <- c(if (length(constraints$var) > 0L)
               "dotted line = bound; shading = constrained path at bound",
             if (n_dropped > 0L) sprintf("%d more variable(s) not shown", n_dropped))
  sub <- paste(c(sub, notes), collapse = "; ")
  if (!nzchar(sub)) sub <- NULL

  p <- .gg$ggplot(plot_df, .gg$aes(x = x, y = Value))
  shade_df <- NULL
  bound_df <- NULL
  if (length(constraints$var) > 0L) {
    keep <- constraints$var %in% pick
    bound_df <- data.frame(Variable = factor(constraints$var[keep], levels = lvl),
                           bound = constraints$bound[keep])
    shade <- lapply(which(keep), function(j) {
      at <- binding[[j]]
      if (!any(at)) return(NULL)
      r <- rle(at); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
      s <- starts[r$values]; e <- ends[r$values]
      data.frame(Variable = factor(constraints$var[j], levels = lvl),
                 xmin = xval[s] - half, xmax = xval[e] + half)
    })
    shade_df <- do.call(rbind, shade)
    if (!is.null(shade_df) && nrow(shade_df) > 0L) {
      p <- p + .gg$geom_rect(data = shade_df, inherit.aes = FALSE,
                             .gg$aes(xmin = xmin, xmax = xmax),
                             ymin = -Inf, ymax = Inf,
                             fill = dynhr_colours$light_grey, alpha = 0.6)
    }
    if (nrow(bound_df) > 0L) {
      p <- p + .gg$geom_hline(data = bound_df, .gg$aes(yintercept = bound),
                              colour = dynhr_colours$red, linetype = "dotted",
                              linewidth = 0.6)
    }
  }
  p <- p +
    .gg$geom_line(.gg$aes(colour = Scenario, linetype = Scenario), linewidth = 0.7) +
    .gg$facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
    .gg$scale_colour_manual(values = c("Constrained" = dynhr_colours$mid_blue,
                                       "Unconstrained" = dynhr_colours$orange)) +
    .gg$scale_linetype_manual(values = c("Constrained" = "solid",
                                         "Unconstrained" = "dashed")) +
    .gg$scale_y_continuous(n.breaks = 4L) +
    .d31_theme(.gg) +
    .gg$labs(title = sprintf("D31: %s -- constrained vs unconstrained path", ttl),
             subtitle = sub, x = x_lab, y = y_label, colour = NULL, linetype = NULL)
  plots <- list(comparison = .apply_meta(p, meta))

  # ---- difference plot: only variables the constraint actually moves ----
  act <- effects$variable[effects$affected]
  act <- act[order(match(act, pick))]
  act <- act[!is.na(match(act, pick))]
  if (length(act) == 0L) {
    msg <- if (length(constraints$var) > 0L &&
               !any(vapply(binding, any, logical(1))))
      "The constraint never binds in this scenario:\nconstrained and unconstrained paths coincide."
    else "Constrained and unconstrained paths coincide\n(no variable is altered)."
    pd <- .gg$ggplot() +
      .gg$annotate("text", x = 0, y = 0, label = msg, size = 4.2,
                   colour = dynhr_colours$dark_blue) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 14),
                     plot.margin = ggplot2::margin(8, 8, 8, 8)) +
      .gg$labs(title = sprintf("D31: %s -- constraint effect", ttl))
  } else {
    n_flat <- length(pick) - length(act)
    diff_df <- data.frame(
      x = rep(xval, length(act)),
      Diff = as.vector(diff_mat[, act, drop = FALSE]),
      Variable = factor(rep(act, each = n_T), levels = act))
    pd <- .gg$ggplot(diff_df, .gg$aes(x = x, y = Diff))
    if (!is.null(shade_df) && nrow(shade_df) > 0L) {
      sd2 <- shade_df[as.character(shade_df$Variable) %in% act, , drop = FALSE]
      if (nrow(sd2) > 0L) {
        sd2$Variable <- factor(as.character(sd2$Variable), levels = act)
        pd <- pd + .gg$geom_rect(data = sd2, inherit.aes = FALSE,
                                 .gg$aes(xmin = xmin, xmax = xmax),
                                 ymin = -Inf, ymax = Inf,
                                 fill = dynhr_colours$light_grey, alpha = 0.6)
      }
    }
    pd <- pd +
      .gg$geom_hline(yintercept = 0, colour = dynhr_colours$grey, linewidth = 0.3) +
      .gg$geom_line(colour = dynhr_colours$mid_blue, linewidth = 0.7) +
      .gg$facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
      .gg$scale_y_continuous(n.breaks = 4L) +
      .d31_theme(.gg) +
      .gg$labs(title = sprintf("D31: %s -- constraint effect", ttl),
               subtitle = paste0("constrained minus unconstrained, in each variable's own units",
                                 if (n_flat > 0L)
                                   sprintf("; %d unaffected variable(s) omitted", n_flat)
                                 else ""),
               x = x_lab, y = "Constrained - unconstrained")
  }
  plots$difference <- .apply_meta(pd, meta)
  plots
}
