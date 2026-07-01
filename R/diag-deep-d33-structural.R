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


# Extract the RHS calibration expression for each reduced-form parameter from
# the .mod source (the parser evaluates and discards these, so we re-read them).
# Returns a named list: kappa_name -> expression text.
.extract_reduced_form_exprs <- function(model, kappa_names) {
  if (length(kappa_names) == 0) return(list())
  sf <- model$source_file
  if (is.null(sf) || length(sf) != 1 || is.na(sf) || !file.exists(sf)) return(list())
  txt <- paste(readLines(sf, warn = FALSE), collapse = "\n")
  pat <- "([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*([^;\\n]+)\\s*;"
  hits <- regmatches(txt, gregexpr(pat, txt, perl = TRUE))[[1]]
  out <- list()
  for (m in hits) {
    parts <- regmatches(m, regexec(pat, m, perl = TRUE))[[1]]
    nm <- parts[2]
    if (!(nm %in% kappa_names)) next
    if (!is.null(out[[nm]])) next                 # first assignment wins
    out[[nm]] <- gsub("\\s+", " ", trimws(parts[3]))
  }
  out
}

# Evaluate a calibration expression at a named parameter vector.
.rf_eval <- function(expr_text, param_vec) {
  env <- list2env(as.list(param_vec), parent = baseenv())
  val <- tryCatch(eval(parse(text = expr_text), envir = env),
                  error = function(e) NA_real_)
  if (is.numeric(val) && length(val) == 1) val else NA_real_
}

# Central-difference d kappa / d theta_j at the calibrated point.
.rf_jacobian <- function(eval_fn, param_vec, theta_name, eps = 1e-5) {
  th <- param_vec[[theta_name]]
  if (is.null(th) || !is.finite(th)) return(NA_real_)
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
#' the local elasticity \eqn{E = (\theta/\kappa)\,\partial\kappa/\partial\theta}.
#' \eqn{|E| < 1} flags a primitive whose apparent precision is "borrowed" from
#' its reduced-form image (relative variance inflation \eqn{1/E^2}).
#'
#' @param model     Parsed model (from \code{parse_mod}); supplies the
#'   \code{@dynhr:deep} map, \code{source_file} and \code{param_values}.
#' @param deep_spec Optional \code{\link{build_deep_spec}} (built from
#'   \code{model} if omitted).
#' @param param_values Optional named numeric parameter vector (defaults to
#'   \code{model$param_values}); the point at which the map is differentiated.
#' @param reduced_form_fn Optional override: a function mapping a named
#'   parameter vector to a named vector of reduced-form coefficients.  Used
#'   instead of the .mod calibration expressions when supplied.
#' @param info_kappa Optional named numeric vector of the data's identification
#'   strength (or precision) for each reduced-form coefficient, used to scale
#'   the borrowed verdict by how informative the data actually is about
#'   \eqn{\kappa}.
#' @param elasticity_threshold \eqn{|E|} below which a primitive is flagged
#'   borrowed (default 1).
#' @param meta      Optional \code{\link{diag_meta}} provenance descriptor.
#' @return A \code{dynhr_diagnostic}; \code{$result$passport_axis} is a named
#'   logical (TRUE = direct / structural, FALSE = borrowed) for the Passport.
#' @noRd
d33_structural_vs_reduced_form <- function(model = NULL,
                                           deep_spec = NULL,
                                           param_values = NULL,
                                           reduced_form_fn = NULL,
                                           info_kappa = NULL,
                                           elasticity_threshold = 1.0,
                                           meta = NULL) {
  tryCatch({
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

    # Build the per-kappa evaluator: a supplied fn, else .mod expressions.
    eval_fns <- list()
    if (!is.null(reduced_form_fn)) {
      for (kn in kappa_names)
        eval_fns[[kn]] <- (function(k) { force(k); function(pv) {
          out <- tryCatch(reduced_form_fn(unlist(pv)), error = function(e) NULL)
          if (!is.null(out) && k %in% names(out)) as.numeric(out[[k]]) else NA_real_
        } })(kn)
    } else {
      exprs <- .extract_reduced_form_exprs(model, kappa_names)
      for (kn in kappa_names) {
        et <- exprs[[kn]]
        # force(e): without it the promise resolves -- after the loop -- to the
        # LAST expression, so every kappa would share one map (a real bug that
        # only bites when there are >=2 distinct reduced-form coefficients).
        eval_fns[[kn]] <- if (is.null(et)) function(pv) NA_real_
                          else (function(e) { force(e); function(pv) .rf_eval(e, pv) })(et)
      }
    }

    # Per-primitive geometry of the map.
    rows <- lapply(seq_len(nrow(prim)), function(i) {
      th_name <- prim$param[i]
      kn      <- prim$reduced_form[i]
      ef      <- eval_fns[[kn]]
      th_cal  <- suppressWarnings(as.numeric(param_values[[th_name]]))
      k_cal   <- ef(param_values)
      J       <- .rf_jacobian(ef, param_values, th_name)
      elas    <- if (isTRUE(is.finite(J)) && isTRUE(is.finite(k_cal)) &&
                     isTRUE(is.finite(th_cal)) && abs(k_cal) > 0)
        (th_cal / k_cal) * J else NA_real_
      data.frame(
        param        = th_name,
        reduced_form = kn,
        theta_cal    = th_cal,
        kappa_cal    = k_cal,
        dkappa_dtheta = J,
        elasticity   = elas,
        rel_var_infl = if (isTRUE(is.finite(elas)) && abs(elas) > 0) 1 / elas^2 else NA_real_,
        stringsAsFactors = FALSE)
    })
    tab <- do.call(rbind, rows)

    # Optional refinement: scale by how informative the data is about kappa.
    tab$kappa_info <- NA_real_
    if (!is.null(info_kappa) && !is.null(names(info_kappa))) {
      hit <- match(tab$reduced_form, names(info_kappa))
      tab$kappa_info <- as.numeric(info_kappa)[hit]
    }

    # Verdict. Identification through a reduced-form coefficient fails in two
    # distinct ways:
    #   (1) CONFOUNDING -- more than one deep primitive is routed through the
    #       same kappa, so the single data constraint on kappa cannot separate
    #       them (the canonical NKPC case when beta is not pinned down too).
    #   (2) DILUTION    -- the map compresses identification, |E| < threshold,
    #       so a tight band on kappa still leaves a loose band on theta. A
    #       near-flat / sign-changing map (|E| not finite) is the extreme case.
    # Confounding is the primary, structural signal; dilution is secondary.
    n_co <- stats::ave(seq_len(nrow(tab)), tab$reduced_form, FUN = length) - 1L
    tab$coprimitives <- as.integer(n_co)
    flat       <- !is.finite(tab$elasticity)
    diluted    <- !flat & abs(tab$elasticity) < elasticity_threshold
    confounded <- tab$coprimitives >= 1L
    tab$verdict <- ifelse(confounded, "confounded",
                   ifelse(flat,       "borrowed (flat map)",
                   ifelse(diluted,    "diluted", "direct")))

    is_direct <- tab$verdict == "direct"
    borrowed  <- tab$param[!is_direct]
    passport_axis <- stats::setNames(is_direct, tab$param)

    pass <- length(borrowed) == 0

    # --- plots ---
    plots <- list()
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      plots$elasticity <- .plot_d33_elasticity(tab, elasticity_threshold, meta)
      mc <- tryCatch(.plot_d33_map_curves(tab, eval_fns, param_values, meta),
                     error = function(e) NULL)
      if (!is.null(mc)) plots$map_curves <- mc
    }

    summary_txt <- sprintf(
      "D33 Structural-vs-reduced-form: %d primitive(s) with a reduced-form map. %s",
      nrow(tab),
      if (pass) "All identified directly (|elasticity| >= threshold)."
      else sprintf("BORROWED: %s (data identifies the slope, not the primitive).",
                   paste(borrowed, collapse = ", ")))

    conf_set <- tab$param[tab$verdict == "confounded"]
    dil_set  <- tab$param[grepl("dilut|flat", tab$verdict)]
    llm <- paste(c(
      sprintf("D33 | Structural-vs-Reduced-Form | %s", if (pass) "PASS" else "FAIL"),
      sprintf("  primitives_with_map=%d borrowed=%d (confounded=%d diluted=%d) threshold=|E|<%.2f",
              nrow(tab), length(borrowed), length(conf_set), length(dil_set),
              elasticity_threshold),
      {
        worst <- tab[order(abs(tab$elasticity)), , drop = FALSE]
        worst <- utils::head(worst, 4)
        sprintf("  elasticities: %s",
                paste(sprintf("%s->%s E=%.2f", worst$param, worst$reduced_form,
                              worst$elasticity), collapse = ", "))
      },
      if (length(conf_set))
        sprintf("  confounded: %s (multiple primitives share one reduced-form coeff)",
                paste(conf_set, collapse = ", ")),
      if (length(dil_set))
        sprintf("  diluted: %s (map compresses identification, |E|<%.2f)",
                paste(dil_set, collapse = ", "), elasticity_threshold),
      sprintf("  action: %s",
              if (pass) "Each primitive maps 1:1 to a distinct reduced-form coeff with |E|>=threshold; identification is not borrowed."
              else if (length(conf_set))
                sprintf("%s share a reduced-form coefficient: add an equation/observable that separates them, or calibrate all but one.",
                        paste(utils::head(conf_set, 3), collapse = ", "))
              else sprintf("%s: the data identifies the reduced-form coefficient, not the primitive. Estimate the slope directly, calibrate the primitive, or report identification-robust bands.",
                           paste(utils::head(borrowed, 3), collapse = ", ")))
    ), collapse = "\n")

    .make_result(
      result  = list(table = tab, borrowed = borrowed,
                     passport_axis = passport_axis),
      pass    = pass,
      plots   = plots,
      summary = summary_txt,
      llm_summary = llm)
  }, error = function(e) {
    .make_result(pass = NA, errored = TRUE,
                 summary = paste("D33 structural-vs-reduced-form: ERROR --",
                                 conditionMessage(e)))
  })
}


# Bar chart of |elasticity| with the borrowed threshold.
.plot_d33_elasticity <- function(tab, threshold, meta) {
  df <- tab
  df$abs_E <- abs(df$elasticity)
  df$abs_E[!is.finite(df$abs_E)] <- 0
  df$lab   <- sprintf("%s -> %s", df$param, df$reduced_form)
  df$lab   <- factor(df$lab, levels = df$lab[order(df$abs_E)])
  df$Verdict <- ifelse(df$verdict == "direct", "direct", "borrowed")

  p <- ggplot2::ggplot(df, ggplot2::aes(x = abs_E, y = lab, fill = Verdict)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                        colour = dynhr_colours$red, linewidth = 0.5) +
    ggplot2::scale_fill_manual(
      values = c(direct = dynhr_colours$teal, borrowed = dynhr_colours$red),
      drop = FALSE, name = NULL) +
    ggplot2::labs(
      title    = "D33: Structural-vs-reduced-form identification transfer",
      subtitle = sprintf("|elasticity| of the deep->reduced-form map; |E| < %.2f (dashed) = borrowed",
                         threshold),
      x = expression("|" * E * "|  =  |(" * theta * "/" * kappa * ")  d" * kappa * "/d" * theta * "|"),
      y = NULL)
  p <- p + theme_dynhr()
  .apply_meta(p, meta)
}

# Small-multiples of the map curve kappa(theta) with the calibrated point.
.plot_d33_map_curves <- function(tab, eval_fns, param_values, meta) {
  curves <- list()
  for (i in seq_len(nrow(tab))) {
    th_name <- tab$param[i]; kn <- tab$reduced_form[i]
    th_cal  <- tab$theta_cal[i]
    if (!is.finite(th_cal)) next
    # plausible range: probability-like params (in (0,1)) span (0.05, 0.97),
    # otherwise +/-50% around the calibrated value.
    if (th_cal > 0 && th_cal < 1) {
      grid <- seq(0.05, 0.97, length.out = 60)
    } else {
      grid <- seq(max(1e-6, th_cal * 0.5), th_cal * 1.5, length.out = 60)
    }
    ef <- eval_fns[[kn]]
    kv <- vapply(grid, function(g) {
      pv <- param_values; pv[[th_name]] <- g; ef(pv)
    }, numeric(1))
    curves[[length(curves) + 1]] <- data.frame(
      facet = sprintf("%s -> %s", th_name, kn),
      theta = grid, kappa = kv, stringsAsFactors = FALSE)
  }
  if (length(curves) == 0) return(NULL)
  cdf <- do.call(rbind, curves)
  pts <- data.frame(
    facet = sprintf("%s -> %s", tab$param, tab$reduced_form),
    theta = tab$theta_cal, kappa = tab$kappa_cal,
    verdict = ifelse(tab$verdict == "direct", "direct", "borrowed"),
    stringsAsFactors = FALSE)
  pts <- pts[is.finite(pts$theta) & is.finite(pts$kappa), , drop = FALSE]

  p <- ggplot2::ggplot(cdf, ggplot2::aes(x = theta, y = kappa)) +
    ggplot2::geom_line(colour = dynhr_colours$mid_blue, linewidth = 0.7) +
    ggplot2::geom_point(data = pts,
                        ggplot2::aes(x = theta, y = kappa, colour = verdict),
                        size = 2.6, inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(
      values = c(direct = dynhr_colours$teal, borrowed = dynhr_colours$red),
      drop = FALSE, name = NULL) +
    ggplot2::facet_wrap(~ facet, scales = "free") +
    ggplot2::labs(
      title    = "D33: reduced-form map kappa(theta) with the calibrated point",
      subtitle = "A flat / compressing map turns a tight kappa band into a loose theta band",
      x = expression(theta ~ "(structural primitive)"),
      y = expression(kappa ~ "(reduced-form coefficient)"))
  p <- p + theme_dynhr()
  .apply_meta(p, meta)
}
