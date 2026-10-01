## R/diag-deep-d33-structural.R
## --------------------------------------------------------------------------
## D33. Structural-vs-reduced-form gap ("borrowed identification").
##
## A deep primitive theta is only truly "deep" if the data speaks to *it*, not
## merely to a reduced-form coefficient kappa = g(theta) that a calibration map
## then back-solves. The canonical case is the New Keynesian Phillips Curve:
## the data identifies the slope kappa, while the Calvo probability theta behind
## it is recovered only through
##     kappa = (1 - theta) (1 - beta theta) / theta .
## (Mavroeidis, Plagborg-Moller & Stock 2014.)
##
## This diagnostic reads the deep->reduced-form map declared in @dynhr:deep
## (`reduced_form=` on each primitive), reconstructs g() from the .mod
## calibration expression (or a supplied map function), and measures how much
## identification the map *transfers* from kappa to theta via the local
## elasticity  E = (theta/kappa) dkappa/dtheta.  |E| < 1 means the map dilutes
## identification (a tight band on kappa becomes a loose band on theta): the
## primitive's apparent precision is "borrowed".  By the delta method the
## relative variance inflation is 1/E^2.
##
## It needs only the model + taxonomy (no estimation), so it is a Group-B
## (pre-estimation) diagnostic. It exposes `$result$passport_axis` so the
## Deep-Parameter Passport picks up the "structural" verdict automatically.
##
## References:
##   Mavroeidis, S., Plagborg-Moller, M., & Stock, J. H. (2014). Empirical
##     evidence on inflation expectations in the New Keynesian Phillips Curve.
##     Journal of Economic Literature, 52(1), 124-188.
##   Canova, F., & Sala, L. (2009). Back to square one: identification issues
##     in DSGE models. Journal of Monetary Economics, 56(4), 431-449.
## --------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# Calibration-program reader + restricted evaluator.
#
# parse_mod() evaluates the top-level `name = expr;` assignments IN SOURCE
# ORDER (after macro expansion, comment stripping and block removal) and keeps
# only the resulting numbers. To differentiate a reduced-form coefficient kappa
# with respect to a primitive theta, D33 must re-run the part of that program
# kappa depends on, with theta perturbed. The reader below mirrors the parser's
# pre-processing exactly, so that
#   * commented-out assignments are ignored,
#   * the LAST assignment before use wins (as in the parser), and
#   * chained derivations (lam = f(theta); kappa = lam * c) propagate the
#     perturbation instead of silently giving dkappa/dtheta = 0.
#
# Expressions are evaluated in a sandbox whose only functions are arithmetic
# and elementary maths (no baseenv, so a .mod cannot run system() through
# D33), after an AST check that every call is whitelisted and every symbol is
# defined; unsupported expressions are reported as unevaluable, not as NA
# numbers that masquerade as a flat map.
# ---------------------------------------------------------------------------

# 0.9.4: D33's private sandbox was folded into the package-wide .mod
# expression sandbox (`.dynhr_safe_fn_names` / `.dynhr_sandbox_env()` in
# R/parse-blocks.R), which the PARSER now uses too. One allowlist, one
# constructor -- D33 and parse_mod() can no longer drift apart.

# Cheap lexical screen run BEFORE parse(): only arithmetic characters, balanced
# parentheses, no dangling/doubled binary operators, no juxtaposed operands
# (MATLAB `[1 2]`-style text). Whatever passes is valid R arithmetic.
.rf_text_ok <- function(s) {
  if (!nzchar(s) || !grepl("^[-A-Za-z0-9_.+*/^(), ]+$", s)) return(FALSE)
  ch <- strsplit(s, "", fixed = TRUE)[[1]]
  depth <- cumsum((ch == "(") - (ch == ")"))
  if (any(depth < 0) || depth[length(depth)] != 0) return(FALSE)
  if (grepl("[-+*/^(,]\\s*$", s) || grepl("^\\s*[*/^),]", s)) return(FALSE)
  if (grepl("[*/^]\\s*[*/^),]", s) || grepl("[-+]\\s*[*/^),]", s)) return(FALSE)
  if (grepl("\\(\\s*[),]", s) || grepl(",\\s*[),]", s)) return(FALSE)
  if (grepl("[A-Za-z0-9_.)]\\s+[A-Za-z0-9_.]", s)) return(FALSE)
  if (grepl("\\)\\s*\\(", s)) return(FALSE)
  if (grepl("\\.\\.|[0-9]\\.[0-9]*\\.", s)) return(FALSE)
  # an identifier/number directly followed by "(" is a call; a number is not
  if (grepl("(^|[^A-Za-z0-9_.])[0-9.][0-9.eE]*\\s*\\(", s)) return(FALSE)
  TRUE
}

# Ordered calibration program of a model parsed from a .mod file:
# data.frame(name, text, ok) + list columns `expr` and `syms`.
# Returns NULL when the model has no readable source file.
.rf_calibration_program <- function(model) {
  sf <- model$source_file
  if (is.null(sf) || length(sf) != 1L || is.na(sf) || !file.exists(sf)) return(NULL)
  txt <- paste(readLines(sf, warn = FALSE), collapse = "\n")
  txt <- iconv(txt, from = "LATIN1", to = "UTF-8", sub = "?")
  txt <- sub(paste0("^", intToUtf8(0xFEFF)), "", txt)   # BOM, as parse_mod
  txt <- expand_macros(txt, mod_dir = dirname(sf))
  txt <- remove_blocks(strip_comments_and_macros(txt))
  pat <- "([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*([^;\\n]+)\\s*;"
  hits <- regmatches(txt, gregexpr(pat, txt, perl = TRUE))[[1]]
  n <- length(hits)
  nm <- character(n); tx <- character(n); ok <- logical(n)
  exprs <- vector("list", n); syms <- vector("list", n)
  for (i in seq_len(n)) {
    parts <- regmatches(hits[i], regexec(pat, hits[i], perl = TRUE))[[1]]
    nm[i] <- parts[2]
    tx[i] <- gsub("\\s+", " ", trimws(parts[3]))
    ok[i] <- .rf_text_ok(tx[i])
    if (ok[i]) {
      e <- parse(text = tx[i], keep.source = FALSE)[[1]]
      calls <- setdiff(all.names(e, functions = TRUE), all.vars(e))
      ok[i] <- all(calls %in% .dynhr_safe_fn_names)
      exprs[[i]] <- e
      syms[[i]] <- setdiff(all.vars(e), "pi")
    }
  }
  prog <- data.frame(name = nm, text = tx, ok = ok, stringsAsFactors = FALSE)
  prog$expr <- exprs
  prog$syms <- syms
  prog
}

# Indices (source order) of the assignments `target` transitively depends on,
# resolving each symbol to its LAST valid assignment before the point of use.
# Returns list(idx, deps); idx is empty if `target` is never (validly) assigned.
.rf_closure <- function(prog, target) {
  last_def <- function(sym, before) {
    w <- which(prog$name == sym & prog$ok & seq_len(nrow(prog)) < before)
    if (length(w)) max(w) else NA_integer_
  }
  start <- last_def(target, nrow(prog) + 1L)
  if (is.na(start)) return(list(idx = integer(0), deps = character(0)))
  idx <- integer(0); todo <- start; deps <- character(0)
  while (length(todo)) {
    i <- todo[1]; todo <- todo[-1]
    if (i %in% idx) next
    idx <- c(idx, i)
    for (s in prog$syms[[i]]) {
      deps <- union(deps, s)
      j <- last_def(s, i)
      if (!is.na(j)) todo <- c(todo, j)
    }
  }
  list(idx = sort(idx), deps = deps)
}

# Evaluator for one reduced-form coefficient from the .mod program.
# `pinned` (the primitives) are inputs: never recomputed from their own .mod
# expressions, so a perturbation of theta is not overwritten. Constant
# assignments of names already in `pv` keep the `pv` value (so user-supplied
# param_values are respected); derived assignments are recomputed.
.rf_program_evaluator <- function(prog, idx, target, pinned) {
  force(prog); force(idx); force(target); force(pinned)
  sandbox <- .dynhr_sandbox_env()
  function(pv) {
    env <- list2env(as.list(pv), parent = sandbox)
    for (i in idx) {
      nm <- prog$name[i]
      if (nm %in% pinned) next
      sy <- prog$syms[[i]]
      if (length(sy) == 0L && exists(nm, envir = env, inherits = FALSE)) next
      if (!all(vapply(sy, exists, logical(1), envir = env, inherits = FALSE)))
        return(NA_real_)
      v <- suppressWarnings(eval(prog$expr[[i]], envir = env))
      if (!is.numeric(v) || length(v) != 1L) return(NA_real_)
      assign(nm, v, envir = env)
    }
    if (!exists(target, envir = env, inherits = FALSE)) return(NA_real_)
    v <- get(target, envir = env, inherits = FALSE)
    if (is.numeric(v) && length(v) == 1L) as.numeric(v) else NA_real_
  }
}

# Evaluator from a user-supplied reduced_form_fn (named vector -> named vector).
# A failing map must signal it by returning non-finite values or omitting the
# name; errors propagate.
.rf_fn_evaluator <- function(reduced_form_fn, target) {
  force(reduced_form_fn); force(target)
  function(pv) {
    out <- reduced_form_fn(unlist(pv))
    if (!is.numeric(out) || is.null(names(out)) || !(target %in% names(out)))
      return(NA_real_)
    as.numeric(out[[target]])
  }
}

# Central-difference d kappa / d theta at `param_vec`.
# Step h = max(|theta|, 1e-3) * eps (eps = 1e-5 ~ cube root of machine eps,
# the truncation/round-off optimum for a central difference).
.rf_jacobian <- function(eval_fn, param_vec, theta_name, eps = 1e-5) {
  th <- suppressWarnings(as.numeric(param_vec[[theta_name]]))
  if (length(th) != 1L || !is.finite(th)) return(NA_real_)
  h  <- max(abs(th), 1e-3) * eps
  vp <- param_vec; vp[[theta_name]] <- th + h
  vm <- param_vec; vm[[theta_name]] <- th - h
  (eval_fn(vp) - eval_fn(vm)) / (2 * h)
}


# ---------------------------------------------------------------------------
#' D33. Structural-vs-reduced-form gap ("borrowed identification")
#'
#' For each deep primitive that declares a reduced-form image
#' (\code{reduced_form=} in the \code{@dynhr:deep} block), measures how much
#' identification the calibration map \eqn{\kappa = g(\theta)} transfers from
#' the data-facing coefficient \eqn{\kappa} to the primitive \eqn{\theta}, via
#' the local log-elasticity
#' \eqn{E = \partial\log|\kappa|/\partial\log|\theta| = (\theta/\kappa)\,\partial\kappa/\partial\theta}
#' (central finite difference, all other parameters held fixed).
#' \eqn{|E| < } \code{elasticity_threshold} flags a primitive whose apparent
#' precision is "borrowed" from its reduced-form image (delta-method relative
#' variance inflation \eqn{1/E^2}); two or more primitives declared on the same
#' \eqn{\kappa} are flagged "confounded".
#'
#' The map is rebuilt from the .mod calibration program (the assignments
#' \eqn{\kappa} transitively depends on, last assignment before use wins, as in
#' \code{parse_mod}); every declared primitive is held as an input. Only
#' arithmetic and elementary functions are evaluated. Rows whose map cannot be
#' evaluated (model parsed from text, unsupported expression) are verdict
#' \code{"unevaluable"}; a declared \eqn{\kappa} whose calibration does not
#' involve \eqn{\theta} is \code{"unmapped"}; and \eqn{\theta = 0} or
#' \eqn{\kappa = 0} (log-elasticity undefined) is \code{"undefined (zero level)"}.
#' These rows give \code{NA} on the Passport axis.
#'
#' @param model     Parsed model (from \code{parse_mod}); supplies the
#'   \code{@dynhr:deep} map, \code{source_file} and \code{param_values}.
#' @param deep_spec Optional \code{\link{build_deep_spec}} (built from
#'   \code{model} if omitted).
#' @param param_values Optional named numeric parameter vector (defaults to
#'   \code{model$param_values}); the point at which the map is differentiated.
#' @param reduced_form_fn Optional override: a function mapping a named
#'   parameter vector to a named vector of reduced-form coefficients.  Used
#'   instead of the .mod calibration expressions when supplied (required for a
#'   model parsed from text). Signal failure by returning non-finite values.
#' @param info_kappa Optional named numeric vector: the data's precision
#'   (Fisher information, \eqn{1/\mathrm{var}}) for each reduced-form
#'   coefficient, matched by name. Reported only: the table gains
#'   \code{kappa_info} and the implied delta-method precision of the primitive,
#'   \code{theta_info} \eqn{= (\partial\kappa/\partial\theta)^2} \code{kappa_info}.
#' @param elasticity_threshold \eqn{|E|} below which a primitive is flagged
#'   diluted (default 1).
#' @param meta      Optional \code{\link{diag_meta}} provenance descriptor.
#' @return A \code{dynhr_diagnostic}; \code{$result$passport_axis} is a named
#'   logical (TRUE = direct, FALSE = confounded/diluted, NA = not assessable)
#'   for the Passport. \code{pass} is FALSE if any primitive is borrowed, NA if
#'   none is borrowed but some could not be assessed, TRUE otherwise.
#' @noRd
d33_structural_vs_reduced_form <- function(model = NULL,
                                           deep_spec = NULL,
                                           param_values = NULL,
                                           reduced_form_fn = NULL,
                                           info_kappa = NULL,
                                           elasticity_threshold = 1.0,
                                           meta = NULL) {
  if (!is.null(reduced_form_fn) && !is.function(reduced_form_fn))
    .dynhr_abort("d33_structural_vs_reduced_form(): `reduced_form_fn` must be a function.")
  if (!is.numeric(elasticity_threshold) || length(elasticity_threshold) != 1L ||
      !is.finite(elasticity_threshold) || elasticity_threshold < 0)
    .dynhr_abort("d33_structural_vs_reduced_form(): `elasticity_threshold` must be a single finite number >= 0.")
  if (is.null(deep_spec)) deep_spec <- build_deep_spec(model = model)
  if (is.null(param_values)) param_values <- model$param_values
  if (is.null(param_values) || length(param_values) == 0)
    return(.make_result(pass = NA,
      summary = "D33 structural-vs-reduced-form: no parameter values available."))
  param_values <- as.list(param_values)

  # primitives that declare a reduced-form image
  prim <- deep_spec[deep_spec$role == "primitive" &
                    !is.na(deep_spec$reduced_form), , drop = FALSE]
  if (nrow(prim) == 0)
    return(.make_result(pass = NA,
      summary = paste("D33 structural-vs-reduced-form: no reduced-form maps",
                      "declared (add reduced_form= to a @dynhr:deep primitive).")))

  kappa_names <- unique(prim$reduced_form)

  # Per-kappa evaluator (+ the .mod symbols kappa depends on, NULL if unknown).
  eval_fns <- list(); kappa_deps <- list(); why_na <- list()
  if (!is.null(reduced_form_fn)) {
    for (kn in kappa_names) eval_fns[[kn]] <- .rf_fn_evaluator(reduced_form_fn, kn)
  } else {
    prog <- if (is.null(model)) NULL else .rf_calibration_program(model)
    for (kn in kappa_names) {
      if (is.null(prog)) {
        why_na[[kn]] <- "no readable .mod source (supply reduced_form_fn)"
        eval_fns[[kn]] <- function(pv) NA_real_
        next
      }
      cl <- .rf_closure(prog, kn)
      if (length(cl$idx) == 0L) {
        why_na[[kn]] <- "no evaluable calibration expression in the .mod"
        eval_fns[[kn]] <- function(pv) NA_real_
        next
      }
      kappa_deps[[kn]] <- cl$deps
      # the primitives declared on this kappa are the inputs being perturbed
      eval_fns[[kn]] <- .rf_program_evaluator(
        prog, cl$idx, kn, pinned = setdiff(prim$param[prim$reduced_form == kn], kn))
    }
  }

  # Per-primitive geometry of the map.
  rows <- lapply(seq_len(nrow(prim)), function(i) {
    th_name <- prim$param[i]
    kn      <- prim$reduced_form[i]
    ef      <- eval_fns[[kn]]
    th_cal  <- suppressWarnings(as.numeric(param_values[[th_name]]))
    th_cal  <- if (length(th_cal) == 1L) th_cal else NA_real_
    k_cal   <- ef(param_values)
    J       <- .rf_jacobian(ef, param_values, th_name)
    ok      <- is.finite(J) && is.finite(k_cal) && is.finite(th_cal)
    elas    <- if (ok && th_cal != 0 && k_cal != 0) (th_cal / k_cal) * J else NA_real_
    data.frame(
      param         = th_name,
      reduced_form  = kn,
      theta_cal     = th_cal,
      kappa_cal     = k_cal,
      dkappa_dtheta = J,
      elasticity    = elas,
      rel_var_infl  = if (is.finite(elas) && elas != 0) 1 / elas^2 else NA_real_,
      mapped        = is.null(kappa_deps[[kn]]) || th_name %in% kappa_deps[[kn]],
      evaluable     = ok,
      stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)

  # Optional: data precision on kappa -> implied precision on theta (reported).
  tab$kappa_info <- NA_real_
  if (!is.null(info_kappa) && !is.null(names(info_kappa))) {
    hit <- match(tab$reduced_form, names(info_kappa))
    tab$kappa_info <- as.numeric(info_kappa)[hit]
  }
  tab$theta_info <- tab$dkappa_dtheta^2 * tab$kappa_info

  # Verdict. Identification through a reduced-form coefficient fails in two
  # distinct ways:
  #   (1) CONFOUNDING -- more than one deep primitive is declared on the same
  #       kappa, so the single data constraint on kappa cannot separate them.
  #   (2) DILUTION    -- the map compresses identification, |E| < threshold,
  #       so a tight band on kappa still leaves a loose band on theta (a flat
  #       map, E = 0, is the extreme case).
  # Rows that cannot be assessed are reported as such, never as borrowed.
  n_co <- stats::ave(seq_len(nrow(tab)), tab$reduced_form, FUN = length) - 1L
  tab$coprimitives <- as.integer(n_co)
  zero_lvl <- tab$evaluable & !is.finite(tab$elasticity)
  tab$verdict <- ifelse(tab$coprimitives >= 1L, "confounded",
                 ifelse(!tab$mapped,            "unmapped",
                 ifelse(!tab$evaluable,         "unevaluable",
                 ifelse(zero_lvl,               "undefined (zero level)",
                 ifelse(abs(tab$elasticity) < elasticity_threshold,
                                                "diluted", "direct")))))

  axis <- ifelse(tab$verdict == "direct", TRUE,
          ifelse(tab$verdict %in% c("confounded", "diluted"), FALSE, NA))
  passport_axis <- stats::setNames(axis, tab$param)
  borrowed <- tab$param[axis %in% FALSE]
  unassessed <- tab$param[is.na(axis)]
  pass <- if (length(borrowed)) FALSE else if (length(unassessed)) NA else TRUE

  # --- plots ---
  plots <- list()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plots$elasticity <- .plot_d33_elasticity(tab, elasticity_threshold, meta)
    mc <- .plot_d33_map_curves(tab, eval_fns, param_values, meta,
                               threshold = elasticity_threshold)
    if (!is.null(mc)) plots$map_curves <- mc
  }

  una_txt <- if (length(unassessed))
    sprintf(" Not assessable: %s.", paste(sprintf("%s (%s%s)", unassessed,
      tab$verdict[is.na(axis)],
      vapply(tab$reduced_form[is.na(axis)], function(k)
        if (is.null(why_na[[k]])) "" else paste0(": ", why_na[[k]]), character(1))),
      collapse = ", ")) else ""
  summary_txt <- sprintf(
    "D33 Structural-vs-reduced-form: %d primitive(s) with a reduced-form map. %s%s",
    nrow(tab),
    if (length(borrowed))
      sprintf("BORROWED: %s (data identifies the reduced-form coefficient, not the primitive).",
              paste(borrowed, collapse = ", "))
    else if (length(unassessed)) "None borrowed among those assessed."
    else sprintf("All identified directly (|elasticity| >= %.2f).", elasticity_threshold),
    una_txt)

  conf_set <- tab$param[tab$verdict == "confounded"]
  dil_set  <- tab$param[tab$verdict == "diluted"]
  badge <- if (isTRUE(pass)) "PASS" else if (isFALSE(pass)) "FAIL" else "INFO"
  fin <- tab[is.finite(tab$elasticity), , drop = FALSE]
  fin <- utils::head(fin[order(abs(fin$elasticity)), , drop = FALSE], 4)
  llm <- paste(c(
    sprintf("D33 | Structural-vs-Reduced-Form | %s", badge),
    sprintf("  primitives_with_map=%d borrowed=%d (confounded=%d diluted=%d) unassessed=%d threshold=|E|<%.2f",
            nrow(tab), length(borrowed), length(conf_set), length(dil_set),
            length(unassessed), elasticity_threshold),
    if (nrow(fin))
      sprintf("  elasticities (smallest |E| first): %s",
              paste(sprintf("%s->%s E=%.3g", fin$param, fin$reduced_form,
                            fin$elasticity), collapse = ", ")),
    if (length(conf_set))
      sprintf("  confounded: %s (multiple primitives share one reduced-form coeff)",
              paste(conf_set, collapse = ", ")),
    if (length(dil_set))
      sprintf("  diluted: %s (map compresses identification, |E|<%.2f)",
              paste(dil_set, collapse = ", "), elasticity_threshold),
    if (length(unassessed))
      sprintf(" %s", una_txt),
    sprintf("  action: %s",
            if (isTRUE(pass)) "Each primitive maps 1:1 to a distinct reduced-form coeff with |E|>=threshold; identification is not borrowed."
            else if (length(conf_set))
              sprintf("%s share a reduced-form coefficient: add an equation/observable that separates them, or calibrate all but one.",
                      paste(utils::head(conf_set, 3), collapse = ", "))
            else if (length(dil_set))
              sprintf("%s: the data identifies the reduced-form coefficient, not the primitive. Estimate the slope directly, calibrate the primitive, or report identification-robust bands.",
                      paste(utils::head(dil_set, 3), collapse = ", "))
            else "Fix the reduced_form= declarations or supply reduced_form_fn so every map can be evaluated.")
  ), collapse = "\n")

  .make_result(
    result  = list(table = tab, borrowed = borrowed, unassessed = unassessed,
                   passport_axis = passport_axis),
    pass    = pass,
    plots   = plots,
    summary = summary_txt,
    llm_summary = llm)
}


.d33_verdict_colours <- c(direct = unname(tol_vibrant["teal"]),
                          diluted = unname(tol_vibrant["orange"]),
                          confounded = unname(tol_vibrant["red"]),
                          `not assessed` = unname(tol_vibrant["grey"]))

.d33_verdict_group <- function(v) {
  factor(ifelse(v %in% c("direct", "diluted", "confounded"), v, "not assessed"),
         levels = names(.d33_verdict_colours))
}

# Bar chart of |elasticity| with the dilution threshold. Rows without a finite
# elasticity are listed in the subtitle, not drawn as zero-length "borrowed" bars.
.plot_d33_elasticity <- function(tab, threshold, meta) {
  df <- tab[is.finite(tab$elasticity), , drop = FALSE]
  missing <- tab[!is.finite(tab$elasticity), , drop = FALSE]
  note <- if (nrow(missing))
    paste0("Not drawn (no finite elasticity): ",
           paste(sprintf("%s -> %s [%s]", missing$param, missing$reduced_form,
                         missing$verdict), collapse = "; "))
  else NULL
  if (nrow(df) == 0L) {
    p <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No finite elasticity to plot",
                        size = 5) +
      ggplot2::labs(title = "D33: structural-vs-reduced-form identification transfer",
                    subtitle = note, x = NULL, y = NULL) +
      theme_dynhr() +
      ggplot2::theme(axis.text = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank())
    return(.apply_meta(p, meta))
  }
  df$abs_E <- abs(df$elasticity)
  df$lab   <- sprintf("%s -> %s", df$param, df$reduced_form)
  df$lab   <- factor(df$lab, levels = unique(df$lab[order(df$abs_E)]))
  df$Verdict <- .d33_verdict_group(df$verdict)
  df$txt   <- sprintf("E = %.3g", df$elasticity)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = abs_E, y = lab, fill = Verdict)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = txt), hjust = -0.1, size = 3.5) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                        colour = tol_vibrant[["red"]], linewidth = 0.5) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::scale_fill_manual(values = .d33_verdict_colours, drop = TRUE,
                               name = NULL) +
    ggplot2::labs(
      title    = "D33: structural-vs-reduced-form identification transfer",
      subtitle = paste(c(sprintf("Elasticity of reduced-form coeff w.r.t. primitive; dashed: |E| = %.2f (below = diluted)",
                                 threshold), note), collapse = "\n"),
      x = expression("|" * E * "|  =  |(" * theta * "/" * kappa * ")  d" * kappa * "/d" * theta * "|"),
      y = NULL)
  p <- p + theme_dynhr()
  .apply_meta(p, meta)
}

# Theta grid for a map curve: +/-50% around the calibrated value (same sign),
# capped inside (0, 1) for probability-like params; always contains theta.
.d33_theta_grid <- function(th, n = 81) {
  lo <- th - abs(th) / 2
  hi <- th + abs(th) / 2
  if (th > 0 && th < 1) hi <- min(hi, (1 + th) / 2)
  sort(unique(c(seq(lo, hi, length.out = n), th)))
}

# Small multiples of the map on log-log axes RELATIVE to the calibration:
# x = theta/theta_cal, y = kappa/kappa_cal. The slope at (1, 1) is exactly the
# elasticity E, so the dashed guides y = x^(+/-threshold) are the dilution
# boundary: a curve flatter than the guides is diluted. Rows without a finite
# elasticity (theta or kappa = 0, unevaluable) have no panel.
.plot_d33_map_curves <- function(tab, eval_fns, param_values, meta,
                                 threshold = 1) {
  curves <- list(); guides <- list()
  for (i in seq_len(nrow(tab))) {
    if (!is.finite(tab$elasticity[i])) next
    th_name <- tab$param[i]; kn <- tab$reduced_form[i]
    th_cal  <- tab$theta_cal[i]; k_cal <- tab$kappa_cal[i]
    facet <- sprintf("%s -> %s", th_name, kn)
    grid <- .d33_theta_grid(th_cal)
    ef <- eval_fns[[kn]]
    kv <- vapply(grid, function(g) {
      pv <- param_values; pv[[th_name]] <- g; ef(pv)
    }, numeric(1))
    xr <- grid / th_cal; yr <- kv / k_cal
    keep <- is.finite(yr) & yr > 0
    if (!any(keep)) next
    curves[[length(curves) + 1]] <- data.frame(
      facet = facet, x = xr[keep], y = yr[keep], stringsAsFactors = FALSE)
    guides[[length(guides) + 1]] <- data.frame(
      facet = facet, x = rep(range(xr), 2),
      y = c(range(xr)^threshold, range(xr)^(-threshold)),
      grp = rep(c("up", "down"), each = 2), stringsAsFactors = FALSE)
  }
  if (length(curves) == 0) return(NULL)
  cdf <- do.call(rbind, curves)
  gdf <- do.call(rbind, guides)
  pts <- data.frame(
    facet = sprintf("%s -> %s", tab$param, tab$reduced_form),
    x = 1, y = 1,
    Verdict = .d33_verdict_group(tab$verdict),
    lab = sprintf("E = %.3g", tab$elasticity),
    stringsAsFactors = FALSE)
  pts <- pts[pts$facet %in% cdf$facet, , drop = FALSE]

  p <- ggplot2::ggplot(cdf, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(data = gdf, ggplot2::aes(x = x, y = y, group = grp),
                       linetype = "dashed", colour = "grey45", linewidth = 0.4,
                       inherit.aes = FALSE) +
    ggplot2::geom_line(colour = dynhr_primary_colour, linewidth = 0.8) +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = x, y = y, colour = Verdict),
                        size = 2.8, inherit.aes = FALSE) +
    ggplot2::geom_label(data = pts, ggplot2::aes(x = x, y = y, label = lab),
                        vjust = -0.5, size = 3.2, linewidth = 0,
                        fill = "white", alpha = 0.8, inherit.aes = FALSE) +
    ggplot2::scale_x_log10(breaks = c(0.5, 0.75, 1, 1.25, 1.5)) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_colour_manual(values = .d33_verdict_colours, drop = TRUE,
                                 name = NULL) +
    ggplot2::facet_wrap(~ facet, scales = "free_y") +
    ggplot2::labs(
      title    = "D33: reduced-form map relative to the calibration (log-log)",
      subtitle = sprintf("Slope at the dot = elasticity E; dashed: |E| = %.2f (flatter = diluted). Other parameters held at calibration",
                         threshold),
      x = expression(theta / theta[cal] ~ "(structural primitive, log scale)"),
      y = expression(kappa / kappa[cal] ~ "(reduced form, log scale)"))
  p <- p + theme_dynhr()
  .apply_meta(p, meta)
}
