## R/diag-obc-d32-binding.R
## --------------------------------------------------------------------------
## Phase H: D32 -- OBC Binding Periods Summary
##
## How often, when, and for how long each occasionally binding constraint
## binds along a regime path.
##
## Regime encoding: every dynhr OBC solver reports the regime as an integer
## BITFIELD (bit b-1 set <=> constraint b binds; see R/obc-regime-path.R), so
## regime value 3 means "constraints 1 AND 2 bind", not "constraint 3". This
## diagnostic decodes the bitfield into a T x n_constr logical matrix and
## computes every statistic per constraint from that matrix.
##
## (This file replaces the former duplicate d_obc_binding_summary() in
## R/diag-obc.R, which treated regime values as constraint indices.)
## --------------------------------------------------------------------------

## Maximal runs of TRUE in a logical vector: data.frame(start, end, duration).
.d32_spells <- function(x) {
  x <- as.logical(x)
  if (!length(x) || !any(x))
    return(data.frame(start = integer(0), end = integer(0),
                      duration = integer(0)))
  r      <- rle(x)
  ends   <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  k      <- which(r$values)
  data.frame(start = as.integer(starts[k]), end = as.integer(ends[k]),
             duration = as.integer(r$lengths[k]))
}

## Normalise the regime_path argument to a T x n_constr logical matrix.
.d32_binding_matrix <- function(regime_path, n_constr) {
  if (inherits(regime_path, "obc_regime_path"))
    return(regime_path$binding)
  if (is.list(regime_path) && !is.data.frame(regime_path)) {
    ## OBC solver result: logical $regime matrix, or a bitfield
    regime_path <- regime_path$regime %||% regime_path$regime_path %||%
      regime_path$active_set
    if (is.null(regime_path))
      .dynhr_abort(paste0("d32_obc_binding_summary(): the solver result has ",
                          "no 'regime'/'regime_path'/'active_set' field."))
  }
  if (is.matrix(regime_path) || is.data.frame(regime_path)) {
    B <- as.matrix(regime_path)
    if (!is.logical(B) && !all(B %in% c(0, 1, NA)))
      .dynhr_abort(paste0("d32_obc_binding_summary(): a matrix regime_path ",
                          "must be a T x n_constr logical (or 0/1) matrix."))
    if (anyNA(B))
      .dynhr_abort("d32_obc_binding_summary(): regime_path contains NA.")
    return(matrix(as.logical(B), nrow = nrow(B),
                  dimnames = list(NULL, colnames(B))))
  }
  if (anyNA(regime_path))
    .dynhr_abort("d32_obc_binding_summary(): regime_path contains NA.")
  codes <- as.integer(regime_path)
  if (any(codes != regime_path) || any(codes < 0L))
    .dynhr_abort(paste0("d32_obc_binding_summary(): regime_path must be a ",
                        "non-negative integer bitfield (0 = slack)."))
  n_bits <- if (length(codes) && max(codes) > 0L)
    as.integer(floor(log2(max(codes)))) + 1L else 1L
  if (is.null(n_constr)) {
    n_constr <- n_bits
  } else if (n_bits > n_constr) {
    .dynhr_abort(sprintf(paste0(
      "d32_obc_binding_summary(): regime_path contains bitfield value %d, ",
      "which needs %d constraints, but only %d are declared."),
      max(codes), n_bits, n_constr))
  }
  .obc_bitfield_to_binding(codes, n_constr)
}


#' D32. OBC Binding Periods Summary
#'
#' How often, when, and for how long each occasionally binding constraint
#' binds along a regime path. The regime path is decoded as a bitfield (bit
#' \code{b-1} set means constraint \code{b} binds), so simultaneous binding
#' of several constraints is attributed to each of them.
#'
#' Reported per constraint: number of binding periods, binding frequency,
#' number of spells (maximal runs of consecutive binding periods), mean and
#' maximum spell length, and a spell table flagging spells that are censored
#' by the sample (start at \code{t = 1} or run to \code{t = T}). An
#' "episode" is a maximal run in which AT LEAST ONE constraint binds.
#'
#' The diagnostic is informational (badge INFO): binding frequency has no
#' pass/fail criterion.
#'
#' When \code{regime_path} is \code{NULL} and \code{sys}, \code{dr_slack} and
#' \code{obc_specs} are supplied, the path is generated with
#' \code{boehl_solve_regime_path()} for a one-off impulse of \code{shock_size}
#' standard deviations to shock \code{shock} in period 1 (an IRF path, not a
#' historical frequency).
#'
#' @param regime_path Regime path: an integer bitfield vector (length T,
#'   0 = slack), a \code{T x n_constr} logical matrix, an
#'   \code{obc_regime_path} object, or an OBC solver result list.
#' @param constraint_names Character vector, one name per constraint (bit).
#'   Defaults to the spec \code{name} / \code{var_name}, else
#'   \code{"Constraint k"}.
#' @param dates Optional Date or numeric vector (length T) for the time axis.
#' @param var_paths Optional \code{n_endo x T} matrix of endogenous paths.
#' @param var_names Optional character vector naming the rows of
#'   \code{var_paths}.
#' @param constraint_var Optional character vector (one per constraint): the
#'   constrained variables. Defaults to the specs' \code{var_name}.
#' @param sys,dr_slack,obc_specs,obs_idx,model Inputs for auto-generation.
#'   \code{obc_specs} also supplies names, variables and bounds.
#' @param n_periods Periods for the auto-generated path (default 20).
#' @param shock Shock (name or index) for the auto-generated path.
#' @param shock_size Impulse size in standard deviations (default 1).
#' @param Sigma_e Shock covariance used to scale the auto-generated impulse.
#'   \code{NULL} means a unit innovation (\code{eps = shock_size}).
#' @param meta Optional dynhr_diag_meta.
#' @return dynhr_diagnostic list. \code{$result} holds \code{binding}
#'   (T x n_constr logical), \code{per_constraint}, \code{spells},
#'   \code{episodes}, \code{n_episodes}, \code{total_binding},
#'   \code{binding_frac}, \code{avg_duration}, \code{max_duration},
#'   \code{n_simultaneous}, \code{always_binding}.
#' @noRd
d32_obc_binding_summary <- function(regime_path      = NULL,
                                    constraint_names = NULL,
                                    dates            = NULL,
                                    var_paths        = NULL,
                                    var_names        = NULL,
                                    constraint_var   = NULL,
                                    sys              = NULL,
                                    dr_slack         = NULL,
                                    obc_specs        = NULL,
                                    obs_idx          = NULL,
                                    model            = NULL,
                                    n_periods        = 20L,
                                    shock            = 1L,
                                    shock_size       = 1,
                                    Sigma_e          = NULL,
                                    meta             = NULL) {

  source_note <- ""
  n_specs <- length(obc_specs)

  # ---- Auto-generate regime path if data provided ----
  if (is.null(regime_path) && !is.null(sys) && !is.null(dr_slack) &&
      n_specs > 0L) {
    n_shock <- NCOL(dr_slack$ghu)
    exo_nm  <- dr_slack$exo_names %||% paste0("shock_", seq_len(n_shock))
    k <- if (is.character(shock)) match(shock, exo_nm) else as.integer(shock)
    if (length(k) != 1L || is.na(k) || k < 1L || k > n_shock)
      .dynhr_abort(sprintf("d32_obc_binding_summary(): unknown shock '%s'.",
                           paste(shock, collapse = ",")))
    sd_k <- 1
    if (!is.null(Sigma_e)) {
      S <- as.matrix(Sigma_e)
      if (!all(dim(S) == n_shock))
        .dynhr_abort(sprintf(paste0("d32_obc_binding_summary(): Sigma_e is ",
                                    "%d x %d but the model has %d shocks."),
                             nrow(S), ncol(S), n_shock))
      sd_k <- sqrt(S[k, k])
    }
    shock_seq <- matrix(0, nrow = n_shock, ncol = n_periods)
    shock_seq[k, 1L] <- shock_size * sd_k
    if (is.null(obs_idx)) obs_idx <- seq_len(nrow(dr_slack$ghx))
    occ <- boehl_solve_regime_path(shock_seq, dr_slack, sys, obc_specs,
                                   obs_idx = obs_idx, max_iter = 30L)
    regime_path <- occ$regime_path
    if (is.null(var_paths)) var_paths <- occ$paths
    if (is.null(var_names))
      var_names <- dr_slack$endo_names %||% model$var_names
    source_note <- sprintf(" [auto IRF path: %+g %s impulse to %s at t=1%s]",
                           shock_size,
                           if (is.null(Sigma_e)) "unit" else "s.d.",
                           exo_nm[k],
                           if (isTRUE(occ$converged)) "" else ", NOT converged")
  }

  # ---- Placeholder ----
  if (is.null(regime_path)) {
    return(.make_result(
      result  = NULL,
      pass    = NA,
      plots   = list(),
      summary = paste(
        "D32 OBC binding summary: NOT YET AVAILABLE.",
        "Provide regime_path (bitfield vector, logical matrix or",
        "obc_regime_path) or supply sys + dr_slack + obc_specs."
      ),
      llm_summary = "[INFO] D32 OBC binding status=placeholder reason=no_data"
    ))
  }

  # ---- Decode ----
  n_constr_decl <- if (n_specs > 0L) n_specs
                   else if (length(constraint_names)) length(constraint_names)
                   else NULL
  B <- .d32_binding_matrix(regime_path, n_constr_decl)
  T_total  <- nrow(B)
  n_constr <- ncol(B)

  if (is.null(constraint_names)) {
    constraint_names <- if (n_specs == n_constr && n_specs > 0L) {
      vapply(obc_specs, function(s) s$name %||% s$var_name %||% NA_character_,
             character(1))
    } else if (inherits(regime_path, "obc_regime_path") &&
               length(regime_path$constraints) == n_constr) {
      regime_path$constraints
    } else if (!is.null(colnames(B))) {
      colnames(B)
    } else {
      rep(NA_character_, n_constr)
    }
    miss <- is.na(constraint_names) | !nzchar(constraint_names)
    constraint_names[miss] <- sprintf("Constraint %d", which(miss))
  }
  constraint_names <- as.character(constraint_names)
  if (length(constraint_names) != n_constr)
    .dynhr_abort(sprintf(paste0("d32_obc_binding_summary(): %d constraint ",
                                "name(s) for %d constraint(s)."),
                         length(constraint_names), n_constr))
  colnames(B) <- constraint_names

  if (!is.null(dates) && length(dates) != T_total)
    .dynhr_abort(sprintf(paste0("d32_obc_binding_summary(): dates has length ",
                                "%d but the regime path has %d periods."),
                         length(dates), T_total))
  date_at <- function(i) if (is.null(dates)) rep(NA, length(i)) else dates[i]

  # ---- Per-constraint spells ----
  spell_list <- lapply(seq_len(n_constr), function(j) {
    sp <- .d32_spells(B[, j])
    if (!nrow(sp)) return(NULL)
    data.frame(constraint = constraint_names[j], spell = seq_len(nrow(sp)),
               start_idx = sp$start, end_idx = sp$end, duration = sp$duration,
               start_date = date_at(sp$start), end_date = date_at(sp$end),
               left_censored = sp$start == 1L,
               right_censored = sp$end == T_total,
               stringsAsFactors = FALSE)
  })
  spells <- do.call(rbind, spell_list)
  if (is.null(spells))
    spells <- data.frame(constraint = character(0), spell = integer(0),
                         start_idx = integer(0), end_idx = integer(0),
                         duration = integer(0), start_date = date_at(integer(0)),
                         end_date = date_at(integer(0)),
                         left_censored = logical(0), right_censored = logical(0))

  per_constraint <- do.call(rbind, lapply(seq_len(n_constr), function(j) {
    d <- spells$duration[spells$constraint == constraint_names[j]]
    w <- which(B[, j])
    data.frame(constraint = constraint_names[j],
               n_binding = length(w),
               binding_frac = length(w) / T_total,
               n_spells = length(d),
               mean_spell = if (length(d)) mean(d) else NA_real_,
               max_spell = if (length(d)) max(d) else NA_integer_,
               first_binding = if (length(w)) min(w) else NA_integer_,
               last_binding = if (length(w)) max(w) else NA_integer_,
               stringsAsFactors = FALSE)
  }))

  # ---- Episodes: runs where at least one constraint binds ----
  any_bind <- rowSums(B) > 0L
  ep <- .d32_spells(any_bind)
  episodes <- data.frame(
    episode    = seq_len(nrow(ep)),
    start_idx  = ep$start, end_idx = ep$end, duration = ep$duration,
    start_date = date_at(ep$start), end_date = date_at(ep$end),
    constraints = vapply(seq_len(nrow(ep)), function(i) {
      rows <- ep$start[i]:ep$end[i]
      paste(constraint_names[colSums(B[rows, , drop = FALSE]) > 0L],
            collapse = "+")
    }, character(1)),
    stringsAsFactors = FALSE)

  n_episodes     <- nrow(episodes)
  total_binding  <- sum(any_bind)
  binding_frac   <- total_binding / T_total
  avg_duration   <- if (n_episodes) mean(episodes$duration) else NA_real_
  max_duration   <- if (n_episodes) max(episodes$duration) else NA_integer_
  n_simultaneous <- sum(rowSums(B) >= 2L)
  always_binding <- T_total > 0L && all(any_bind)

  # ---- Plots ----
  plots <- list()
  if (T_total > 0L && requireNamespace("ggplot2", quietly = TRUE)) {
    .gg <- getNamespace("ggplot2")
    x_all <- if (is.null(dates)) seq_len(T_total) else dates
    x_lab <- if (is.null(dates)) "Period" else "Date"
    step  <- if (T_total >= 2L) stats::median(diff(as.numeric(x_all))) else 1
    ## One bar per constraint (first constraint on top): a slack background
    ## over the whole sample plus one rectangle per binding spell, each
    ## period occupying [x_t - step/2, x_t + step/2] (step in x-axis units,
    ## i.e. days on a Date axis).
    x_num <- as.numeric(x_all)
    ypos  <- rev(seq_len(n_constr))
    as_x  <- function(v) if (inherits(dates, "Date"))
      as.Date(v, origin = "1970-01-01") else v
    rect_of <- function(j, s, e, state)
      data.frame(xmin = as_x(x_num[s] - step / 2),
                 xmax = as_x(x_num[e] + step / 2),
                 ymin = ypos[j] - 0.4, ymax = ypos[j] + 0.4,
                 state = factor(state, levels = c("Slack", "Binding")))
    tl_df <- do.call(rbind, c(
      lapply(seq_len(n_constr), function(j) rect_of(j, 1L, T_total, "Slack")),
      lapply(seq_len(n_constr), function(j) {
        sp <- .d32_spells(B[, j])
        if (nrow(sp)) rect_of(j, sp$start, sp$end, "Binding")
      })))
    subtitle <- sprintf(
      "Any constraint binding in %d of %d periods (%.1f%%); %d episode(s)%s",
      total_binding, T_total, 100 * binding_frac, n_episodes,
      if (n_constr > 1L) sprintf("; %d period(s) with >= 2 binding",
                                 n_simultaneous) else "")
    plots$regime_timeline <- .apply_meta(
      .gg$ggplot(tl_df) +
        .gg$geom_rect(.gg$aes(xmin = .data$xmin, xmax = .data$xmax,
                              ymin = .data$ymin, ymax = .data$ymax,
                              fill = .data$state)) +
        .gg$scale_fill_manual(values = c(Slack = dynhr_na_fill,
                                         Binding = dynhr_colours$red),
                              drop = FALSE) +
        .gg$scale_y_continuous(breaks = ypos, labels = constraint_names) +
        ## a rect-only layer does not auto-select a date scale
        (if (inherits(dates, "Date")) .gg$scale_x_date() else NULL) +
        .gg$labs(title = "D32: OBC binding periods by constraint",
                 subtitle = subtitle, x = x_lab, y = NULL, fill = NULL) +
        theme_dynhr_diagnostic() +
        .gg$theme(panel.grid.major.y = .gg$element_blank(),
                  panel.grid.minor = .gg$element_blank(),
                  legend.position = "top"),
      meta)

    ## Constrained variables with binding periods shaded
    cvar <- constraint_var %||%
      (if (n_specs == n_constr)
         vapply(obc_specs, function(s) s$var_name %||% NA_character_,
                character(1)))
    bounds <- if (n_specs == n_constr)
      vapply(obc_specs, function(s) as.numeric(s$bound %||% NA_real_),
             numeric(1)) else rep(NA_real_, n_constr)
    if (!is.null(var_paths) && !is.null(cvar) && !is.null(var_names)) {
      vp <- as.matrix(var_paths)
      if (length(cvar) == 1L && n_constr > 1L) cvar <- rep(cvar, n_constr)
      idx <- match(cvar, var_names)
      ok  <- which(!is.na(idx) & idx <= nrow(vp))
      if (ncol(vp) == T_total && length(ok)) {
        lab <- ifelse(constraint_names[ok] == cvar[ok], cvar[ok],
                      sprintf("%s (%s)", constraint_names[ok], cvar[ok]))
        lab_f <- function(v) factor(v, levels = lab)
        line_df <- data.frame(
          x = rep(x_all, length(ok)),
          value = as.vector(t(vp[idx[ok], , drop = FALSE])),
          panel = lab_f(rep(lab, each = T_total)))
        rect_df <- do.call(rbind, lapply(seq_along(ok), function(m) {
          sp <- .d32_spells(B[, ok[m]])
          if (!nrow(sp)) return(NULL)
          data.frame(xmin = as_x(x_num[sp$start] - step / 2),
                     xmax = as_x(x_num[sp$end] + step / 2),
                     panel = lab_f(lab[m]))
        }))
        hl_df <- data.frame(bound = bounds[ok], panel = lab_f(lab))
        hl_df <- hl_df[is.finite(hl_df$bound), , drop = FALSE]
        p <- .gg$ggplot(line_df, .gg$aes(x = .data$x, y = .data$value))
        if (!is.null(rect_df))
          p <- p + .gg$geom_rect(
            data = rect_df, inherit.aes = FALSE,
            .gg$aes(xmin = .data$xmin, xmax = .data$xmax,
                    ymin = -Inf, ymax = Inf, fill = "Binding"),
            alpha = 0.35)
        if (nrow(hl_df))
          p <- p + .gg$geom_hline(data = hl_df,
                                  .gg$aes(yintercept = .data$bound,
                                          linetype = "Bound"),
                                  colour = "grey30")
        p <- p +
          .gg$geom_line(colour = dynhr_primary_colour, linewidth = 0.7) +
          .gg$facet_wrap(~ panel, ncol = 1L, scales = "free_y") +
          .gg$scale_fill_manual(values = c(Binding = dynhr_colours$red),
                                name = NULL) +
          .gg$scale_linetype_manual(values = c(Bound = "dashed"), name = NULL) +
          .gg$labs(title = "D32: Constrained variables (binding periods shaded)",
                   x = x_lab, y = "Level (model units)") +
          theme_dynhr_diagnostic() +
          .gg$theme(legend.position = "top")
        plots$constraint_var <- .apply_meta(p, meta)
      }
    }
  }

  # ---- Summary ----
  degenerate_note <- if (always_binding)
    " [DEGENERATE: a constraint binds in every period -- check the bound or shock size.]"
    else ""
  detail <- paste(vapply(seq_len(n_constr), function(j) {
    r <- per_constraint[j, ]
    if (r$n_binding == 0L) return(sprintf("%s: never binds", r$constraint))
    sprintf("%s: %d/%d periods (%.1f%%), %d spell(s), mean %.1f, max %d",
            r$constraint, r$n_binding, T_total, 100 * r$binding_frac,
            r$n_spells, r$mean_spell, r$max_spell)
  }, character(1)), collapse = "; ")

  summary_str <- if (n_episodes == 0L) {
    sprintf("D32 OBC binding summary: no constraint binds in %d period(s).%s",
            T_total, source_note)
  } else {
    sprintf(paste0("D32 OBC binding summary: %d episode(s), %d/%d periods ",
                   "with any constraint binding (%.1f%%). %s.%s%s"),
            n_episodes, total_binding, T_total, 100 * binding_frac, detail,
            degenerate_note, source_note)
  }

  .make_result(
    result  = list(binding = B, per_constraint = per_constraint,
                   spells = spells, episodes = episodes,
                   n_episodes = n_episodes, total_binding = total_binding,
                   binding_frac = binding_frac, avg_duration = avg_duration,
                   max_duration = max_duration,
                   n_simultaneous = n_simultaneous,
                   always_binding = always_binding),
    pass    = NA,
    plots   = plots,
    summary = summary_str,
    llm_summary = sprintf(
      "[INFO] D32 OBC binding summary | n_episodes=%d binding_frac=%.3f %s%s",
      n_episodes, binding_frac,
      paste(sprintf("%s=%d/%d", per_constraint$constraint,
                    per_constraint$n_binding, T_total), collapse = " "),
      if (always_binding) " degenerate=always_binding" else "")
  )
}
