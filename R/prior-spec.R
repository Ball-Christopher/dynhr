## R/prior-spec.R
## --------------------------------------------------------------------------
## Phase-2 split from estimate-monolith.R.
##
## extract_prior_spec(): build the prior_spec data.frame from a parsed model.
## Helper functions for shock-to-sigma renaming and distribution defaults.
## --------------------------------------------------------------------------

#' Extract prior specification from a parsed model
#'
#' Reads the \code{estimated_params} block from the parsed model and returns a
#' data.frame suitable for \code{make_log_posterior()}, following Dynare's
#' \code{estimated_params} semantics (Dynare 7.1 reference manual, "Estimation
#' based on likelihood"; \code{estimation/set_prior.m}, \code{estimation/prior_bounds.m},
#' \code{estimation/dynare_estimation_init.m}):
#'
#' * Long form \code{NAME, INITVAL, LB, UB, PRIOR_SHAPE, P1, P2 [, P3, P4, JSCALE]}
#'   and short form \code{NAME, PRIOR_SHAPE, P1, P2 [, P3, P4, JSCALE]} (no
#'   INITVAL/LB/UB).
#' * \code{lower}/\code{upper} are the prior's support (after any P3/P4
#'   generalisation) INTERSECTED with \code{[LB, UB]} -- \code{dynare_estimation_init.m}:
#'   \code{bounds = prior_bounds(...); bounds.lb = max(bounds.lb, lb);
#'   bounds.ub = min(bounds.ub, ub)}. A \code{corr} row's \code{[LB, UB]} is clamped to
#'   \code{[-1, 1]} as in \code{set_prior.m}. An \code{estimated_params_bounds} block
#'   (\code{model$estimated_params_bounds}) replaces LB/UB before the intersection,
#'   as in Dynare.
#' * P3/P4 by shape: BETA -> generalised beta on \code{[P3, P4]} (P1/P2 are the
#'   mean/sd on that interval); GAMMA, INV_GAMMA (1 and 2) -> P3 is a shift
#'   (lower end of the support) and P4 is ignored with a warning (as
#'   in Dynare); NORMAL -> P3/P4 are truncation bounds.
#' * UNIFORM follows Dynare's \code{uniform_specification()}: **P1 is the MEAN and
#'   P2 the STANDARD DEVIATION**, so the support is
#'   \code{[P1 - sqrt(3) P2, P1 + sqrt(3) P2]}. To give the BOUNDS instead, leave
#'   P1/P2 empty and put them in P3/P4:
#'   \code{NAME, INITVAL, LB, UB, uniform_pdf, , , lb, ub;} (or the short form
#'   \code{NAME, uniform_pdf, , , lb, ub;}). P1/P2 together with P3/P4 is an error,
#'   as in Dynare. \code{uniform_pdf, 0, 1} is therefore NOT the unit interval but
#'   \code{[-1.732, 1.732]}; rows that look like a bounds-style spec read as
#'   mean/sd raise a \code{dynhr_warning_uniform_parameterisation} warning (the
#'   calibrated / INITVAL value lies in \code{[P1, P2]} but outside the implied
#'   support; P2 > P1 >= 0; or the support crosses the hard boundary of a
#'   \code{stderr} (below 0) or \code{corr} (outside -1 to 1) parameter). Rows written
#'   with P3/P4 never warn.
#' * The returned \code{p1}/\code{p2} are the MEAN/SD for every shape (for a uniform
#'   given by P3/P4: \code{(P3 + P4)/2} and \code{(P4 - P3)/sqrt(12)}), and the \code{p3}/\code{p4}
#'   columns carry the beta support / shift (NA = the default 0/1 resp. 0) or
#'   the uniform's support \code{[a, b]} -- Dynare's \code{bayestopt_} layout; see
#'   R/prior-density.R.
#' * \code{init} is INITVAL, overridden by an \code{estimated_params_init} block
#'   (\code{model$estimated_params_init}); NA when neither gives one.
#'   \code{run_mode_finding()} starts from it (see \code{.rmf_default_theta_init()}).
#'
#' A \code{stderr eps_X} prior is renamed to the model parameter \code{sig_X} only when
#' the shocks block's stderr / variance expression for \code{eps_X} references
#' \code{sig_X} (e.g. \code{var eps_X; stderr sig_X;}), so the draw reaches the
#' likelihood through that expression. Otherwise the prior keeps the shock name
#' and the draw is injected as the shock's own standard deviation.
#'
#' @param model   dynhr_mod (from parse_mod())
#' @param verbose Logical. Print extracted prior info to console (default TRUE).
#' @return data.frame with columns: name, distribution, p1, p2, lower, upper,
#'         mean, std, p3, p4, init
#' @export
extract_prior_spec <- function(model, verbose = TRUE) {
  ep <- model$estimated_params
  if (is.null(ep) || nrow(ep) == 0) {
    stop("No estimated_params found in model (add an estimated_params block ",
         "to the .mod file, or supply model$estimated_params directly)")
  }

  cnames <- tolower(names(ep))
  ## `corr a, b` rows carry the SECOND shock in `name2` and the row kind in
  ## `type`. Reading only `name` collapsed such a row to the FIRST shock's
  ## name, so the likelihood read the correlation draw as that shock's STDERR
  ## -- or, if a `stderr` prior for the same shock was
  ## also present, tripped the duplicate-name guard.
  col_of <- function(nm) if (nm %in% cnames) which(cnames == nm)[1] else NA_integer_
  col_type  <- col_of("type")
  col_name2 <- col_of("name2")
  col_init  <- col_of("init")
  col_lb    <- col_of("lb")
  col_ub    <- col_of("ub")
  if ("name" %in% cnames && "prior" %in% cnames) {
    col_name <- which(cnames == "name")[1]
    col_dist <- which(cnames == "prior")[1]
    col_p1   <- col_of("p1")
    col_p2   <- col_of("p2")
    col_p3   <- col_of("p3")
    col_p4   <- col_of("p4")
  } else {
    col_name <- 1; col_dist <- 2; col_p1 <- 3; col_p2 <- 4
    col_p3 <- if (ncol(ep) >= 5) 5 else NA
    col_p4 <- if (ncol(ep) >= 6) 6 else NA
    col_init <- NA_integer_; col_lb <- NA_integer_; col_ub <- NA_integer_
  }
  cell_num <- function(i, col) {
    if (is.na(col)) return(NA_real_)
    suppressWarnings(as.numeric(as.character(ep[i, col])))
  }
  abort_row <- function(nm, ...)
    .dynhr_abort("extract_prior_spec: prior '", nm, "': ", ...,
                 class = "dynhr_error_prior_spec")

  n <- nrow(ep)
  key <- character(n); dist_v <- character(n); is_corr <- logical(n)
  p1v <- p2v <- p3v <- p4v <- rep(NA_real_, n)
  nat_lo <- rep(-Inf, n); nat_hi <- rep(Inf, n)
  lbv <- rep(-Inf, n); ubv <- rep(Inf, n); initv <- rep(NA_real_, n)
  uni_meansd <- logical(n); row_kind <- rep("parameter", n)

  for (i in seq_len(n)) {
    nm        <- trimws(as.character(ep[i, col_name]))
    row_type  <- if (!is.na(col_type))
      tolower(trimws(as.character(ep[i, col_type]))) else NA_character_

    if (identical(row_type, "corr")) {
      nm2 <- if (!is.na(col_name2))
        trimws(as.character(ep[i, col_name2])) else NA_character_
      if (is.na(nm2) || !nzchar(nm2) || identical(nm2, "NA")) {
        .dynhr_abort(sprintf(
          paste0("extract_prior_spec: the `corr` prior for '%s' has no second ",
                 "shock. Write it as `corr <shock_a>, <shock_b>, <prior>, ...;` ",
                 "in the estimated_params block."), nm),
          class = "dynhr_error_prior_spec")
      }
      nm <- .corr_key(nm, nm2)
      is_corr[i] <- TRUE
    } else if (identical(row_type, "skew")) {
      ## `skew <shock>` rows also collapsed to the bare shock name, making an
      ## estimated skewness indistinguishable from an estimated stderr; the
      ## draw was then consumed as the STDERR and alpha never moved.
      nm <- .skew_key(nm)
    }
    key[i] <- nm
    if (!is.na(row_type)) row_kind[i] <- row_type

    dist_base <- .normalize_dist(as.character(ep[i, col_dist]))
    dist_v[i] <- dist_base
    p1 <- cell_num(i, col_p1); p2 <- cell_num(i, col_p2)
    p3 <- cell_num(i, col_p3); p4 <- cell_num(i, col_p4)

    if (identical(dist_base, "uniform")) {
      ## Dynare uniform_specification() (estimation/uniform_specification.m):
      ## with P3 AND P4 given the support is [P3, P4] and P1/P2 must be empty
      ## (error_indicator 1); otherwise P1/P2 are the MEAN/SD, must be finite
      ## (error_indicator 2), and the support is P1 -+ sqrt(3) P2. A lone P3
      ## or P4 is ignored, as in Dynare. The row is stored in Dynare's
      ## bayestopt_ layout: p1/p2 = mean/sd, p3/p4 = the support.
      if (!is.na(p3) && !is.na(p4)) {
        if (!is.na(p1) || !is.na(p2))
          abort_row(nm, "a UNIFORM prior with PRIOR_3RD/4TH_PARAMETER (the ",
                    "support [", p3, ", ", p4, "]) must leave P1/P2 empty ",
                    "(Dynare uniform_specification() rejects it too): write ",
                    "`uniform_pdf, , , ", p3, ", ", p4, "`.")
        if (!(is.finite(p3) && is.finite(p4) && p3 < p4))
          abort_row(nm, "a UNIFORM prior's support [P3, P4] = [", p3, ", ",
                    p4, "] must be finite with P3 < P4.")
        lo <- p3; hi <- p4
        p1v[i] <- (lo + hi) / 2; p2v[i] <- (hi - lo) / sqrt(12)
      } else {
        if (!(is.finite(p1) && is.finite(p2)))
          abort_row(nm, "a UNIFORM prior needs a finite MEAN (P1) and ",
                    "STANDARD DEVIATION (P2) -- Dynare convention -- or its ",
                    "bounds in P3/P4 with P1/P2 empty: `uniform_pdf, , , lb, ",
                    "ub` (got P1 = ", p1, ", P2 = ", p2, ").")
        if (!(p2 > 0))
          abort_row(nm, "a UNIFORM prior's STANDARD DEVIATION P2 must be > 0 ",
                    "(got ", p2, "; Dynare convention: P1 = mean, P2 = sd). ",
                    "For bounds write `uniform_pdf, , , lb, ub`.")
        ab <- .uniform_ab(p1, p2)
        lo <- ab[1]; hi <- ab[2]
        p1v[i] <- p1; p2v[i] <- p2
        uni_meansd[i] <- TRUE
      }
      p3v[i] <- lo; p4v[i] <- hi
      nat_lo[i] <- lo; nat_hi[i] <- hi
    } else {
      p1v[i] <- p1; p2v[i] <- p2
      if (identical(dist_base, "beta")) {
        ## Generalised beta on [P3, P4] (Dynare set_prior.m: NaN -> 0 / 1).
        ab <- .prior_gen_ab(p3, p4)
        if (!all(is.finite(ab)) || ab[1] >= ab[2])
          abort_row(nm, "a BETA prior's PRIOR_3RD/4TH_PARAMETER must be a ",
                    "finite interval [p3, p4] with p3 < p4 (got [", ab[1],
                    ", ", ab[2], "]).")
        p3v[i] <- p3; p4v[i] <- p4
        nat_lo[i] <- ab[1]; nat_hi[i] <- ab[2]
      } else if (dist_base %in% c("gamma", "inv_gamma", "inv_gamma1",
                                  "inv_gamma2")) {
        ## P3 = shift (lower end of the support). Dynare never reads P4 for
        ## these shapes (prior_bounds.m sets ub = Inf), so neither does dynhr
        ## -- but not silently: a pre-0.9.4 dynhr spec meant it as an upper
        ## TRUNCATION bound, which now belongs in the long form's UB.
        if (!is.na(p4) && is.finite(p4))
          .dynhr_warn(
            "extract_prior_spec: prior '", nm, "' (", toupper(dist_base),
            "): PRIOR_4TH_PARAMETER = ", p4, " is ignored, as in Dynare (P3 ",
            "is a shift of the support, not a bound). For an upper bound use ",
            "the long form NAME, INITVAL, LB, UB, SHAPE, P1, P2.",
            once = TRUE, key = paste0("prior-p4-ignored:", nm, ":", p4),
            class = "dynhr_warning_prior_p4_ignored")
        if (!is.na(p3) && !is.finite(p3))
          abort_row(nm, "the shift PRIOR_3RD_PARAMETER must be finite (got ",
                    p3, ").")
        s <- .prior_gen_ab(p3, NA_real_)[1]
        p3v[i] <- p3
        nat_lo[i] <- s + .default_lower(dist_base); nat_hi[i] <- Inf
      } else if (identical(dist_base, "normal")) {
        ## P3/P4 are truncation bounds (prior_bounds.m case 3); the density
        ## is the untruncated normal kernel.
        nat_lo[i] <- if (is.na(p3)) -Inf else p3
        nat_hi[i] <- if (is.na(p4))  Inf else p4
      } else {
        ## Unsupported name: validate_prior_spec() below reports it.
        nat_lo[i] <- .default_lower(dist_base); nat_hi[i] <- .default_upper(dist_base)
      }
    }

    lb <- cell_num(i, col_lb); ub <- cell_num(i, col_ub)
    if (!is.na(lb)) lbv[i] <- lb
    if (!is.na(ub)) ubv[i] <- ub
    initv[i] <- cell_num(i, col_init)
  }

  ## estimated_params_bounds block: replaces LB/UB (Dynare semantics), the
  ## intersection with the prior support below still applies.
  if (!is.null(model$estimated_params_bounds)) {
    bnd <- model$estimated_params_bounds
    bnd_cnames   <- tolower(names(bnd))
    bnd_name_col <- if ("name" %in% bnd_cnames) which(bnd_cnames == "name")[1] else 1
    for (i in seq_len(nrow(bnd))) {
      bnm <- trimws(as.character(bnd[i, bnd_name_col]))
      idx <- match(bnm, key)
      if (!is.na(idx)) {
        bvals <- suppressWarnings(as.numeric(as.character(
          bnd[i, setdiff(seq_len(ncol(bnd)), bnd_name_col)])))
        bvals <- bvals[!is.na(bvals)]
        if (length(bvals) >= 1) lbv[idx] <- bvals[1]
        if (length(bvals) >= 2) ubv[idx] <- bvals[2]
      }
    }
  }

  ## corr LB/UB are clamped into [-1, 1] (Dynare set_prior.m).
  lbv[is_corr] <- pmin(pmax(lbv[is_corr], -1), 1)
  ubv[is_corr] <- pmax(pmin(ubv[is_corr],  1), -1)

  lower <- pmax(nat_lo, lbv)
  upper <- pmin(nat_hi, ubv)
  empty <- which(!(lower < upper))
  if (length(empty) > 0L) {
    i <- empty[1]
    abort_row(key[i], "the prior support [", nat_lo[i], ", ", nat_hi[i],
              "] and the estimated_params bounds [LB, UB] = [", lbv[i], ", ",
              ubv[i], "] do not overlap -- no admissible value.")
  }

  ## estimated_params_init overrides INITVAL (Dynare: do_parameter_initialization).
  epi <- model$estimated_params_init
  if (is.data.frame(epi) && nrow(epi) > 0L && all(c("name", "init") %in% names(epi))) {
    for (j in seq_len(nrow(epi))) {
      ty <- if ("type" %in% names(epi)) tolower(trimws(as.character(epi$type[j]))) else NA
      k  <- trimws(as.character(epi$name[j]))
      if (identical(ty, "corr") && "name2" %in% names(epi))
        k <- .corr_key(k, trimws(as.character(epi$name2[j])))
      else if (identical(ty, "skew"))
        k <- .skew_key(k)
      v  <- suppressWarnings(as.numeric(epi$init[j]))
      idx <- match(k, key)
      if (!is.na(idx) && length(v) == 1L && is.finite(v)) initv[idx] <- v
    }
  }

  ## Probable bounds-style uniform rows read under the Dynare mean/sd
  ## convention (never for the P3/P4 form).
  for (i in which(uni_meansd)) {
    cal <- initv[i]
    if (!is.finite(cal) && identical(row_kind[i], "parameter")) {
      pv <- model$param_values
      if (!is.null(pv) && key[i] %in% names(pv))
        cal <- suppressWarnings(as.numeric(pv[[key[i]]]))
    }
    .uniform_param_warn(key[i], row_kind[i], p1v[i], p2v[i],
                        nat_lo[i], nat_hi[i], cal)
  }

  priors <- data.frame(name = key, distribution = dist_v, p1 = p1v, p2 = p2v,
                       lower = lower, upper = upper, stringsAsFactors = FALSE)

  ## `stderr eps_X` -> `sig_X` rename, ONLY when the shocks block routes the
  ## shock's std through sig_X (its stderr_expr / variance_expr references
  ## it). The pre-0.9.4 rule renamed whenever a parameter sig_X existed; with
  ## `var eps_X; stderr 0.1;` the draw then went into an unused sig_X and the
  ## likelihood ignored it (it relied on the removed sig_<shock> variance
  ## heuristic). Without the rename the draw is injected under the shock name,
  ## which .get_shock_stderr() reads first (Priority 0).
  param_names <- names(model$param_values) %||% character(0)
  if (length(param_names) > 0L) {
    rename_map <- character(0)
    for (i in seq_len(nrow(priors))) {
      nm <- priors$name[i]
      if (!(nm %in% param_names) && grepl("^eps_", nm)) {
        sig_nm <- .shock_to_sig(nm)
        if (sig_nm %in% param_names && .shock_expr_references(model, nm, sig_nm)) {
          rename_map[nm] <- sig_nm
          priors$name[i] <- sig_nm
        }
      }
    }
    if (length(rename_map) > 0L && verbose) {
      .dynhr_cat(sprintf("  Renamed %d prior(s) from 'stderr eps_*' to model param 'sig_*':\n",
                  length(rename_map)))
      for (k in names(rename_map)) .dynhr_cat(sprintf("    %s -> %s\n", k, rename_map[[k]]))
    }
  }

  ## Dynare layout: p1/p2 ARE the mean/sd for every shape (uniform included).
  priors$mean <- priors$p1
  priors$std  <- priors$p2
  priors$p3   <- p3v
  priors$p4   <- p4v
  priors$init <- initv

  ## Fail loud at CONSTRUCTION on a malformed / typo'd spec (unknown distribution,
  ## duplicate names, degenerate support) rather than letting it surface as a
  ## silent flat prior or a NaN density deep inside an MCMC run.
  validate_prior_spec(priors, where = "extract_prior_spec")

  if (verbose) {
    .dynhr_cat(sprintf("  Extracted %d priors:\n", nrow(priors)))
    show <- c("name", "distribution", "p1", "p2", "lower", "upper")
    if (any(!is.na(priors$p3) | !is.na(priors$p4))) show <- c(show, "p3", "p4")
    if (any(is.finite(priors$init))) show <- c(show, "init")
    print(priors[, show], row.names = FALSE, right = FALSE)
    .dynhr_cat("\n")
  }
  priors
}

#' Warn when a UNIFORM row read as Dynare MEAN/SD looks like a bounds spec
#'
#' Called by extract_prior_spec() only for rows WITHOUT P3/P4. Triggers:
#' (a) the calibrated / INITVAL value \code{cal} lies in \code{[p1, p2]} but outside the
#' implied support \code{[lo, hi]}; (b) \code{p2 > p1 >= 0} (sd larger than the mean,
#' the shape of \code{uniform_pdf, 0, 1}); (c) the support crosses the hard
#' boundary of a \code{stderr} (< 0) or \code{corr} (outside [-1, 1]) row.
#' @noRd
.uniform_param_warn <- function(name, kind, p1, p2, lo, hi, cal = NA_real_) {
  why <- character(0)
  if (is.finite(cal) && p1 < p2 && cal >= p1 && cal <= p2 &&
      (cal < lo || cal > hi))
    why <- c(why, sprintf(paste0(
      "the calibrated/initial value %s lies in [P1, P2] = [%s, %s] but ",
      "OUTSIDE that support"), format(cal), format(p1), format(p2)))
  if (p2 > p1 && p1 >= 0)
    why <- c(why, sprintf(paste0(
      "P2 = %s exceeds P1 = %s (a standard deviation larger than the mean ",
      "is the signature of a bounds-style row)"), format(p2), format(p1)))
  if (identical(kind, "stderr") && lo < 0)
    why <- c(why, "the support extends below 0 for a standard deviation")
  if (identical(kind, "corr") && (lo < -1 || hi > 1))
    why <- c(why, "the support leaves [-1, 1] for a correlation")
  if (length(why) == 0L) return(invisible(FALSE))
  .dynhr_warn(
    "extract_prior_spec: UNIFORM prior '", name, "' was read with Dynare's ",
    "convention P1 = MEAN (", format(p1), "), P2 = STANDARD DEVIATION (",
    format(p2), "), i.e. the support [", format(lo), ", ", format(hi),
    "] = P1 -+ sqrt(3) P2; but ", paste(why, collapse = "; and "), ". ",
    "If [", format(p1), ", ", format(p2), "] were meant as the BOUNDS, put ",
    "them in P3/P4 with P1/P2 empty: `NAME, INIT, LB, UB, uniform_pdf, , , ",
    format(p1), ", ", format(p2), ";` (short form `NAME, uniform_pdf, , , ",
    format(p1), ", ", format(p2), ";`).",
    class = "dynhr_warning_uniform_parameterisation")
  invisible(TRUE)
}

#' Does the shocks block's stderr / variance expression for \code{shock} reference
#' the parameter \code{par}?
#'
#' Token match on the raw expression text (identifier boundaries on both
#' sides), no evaluation.
#' @noRd
.shock_expr_references <- function(model, shock, par) {
  sv <- model$shocks$variances
  if (!is.data.frame(sv) || !("name" %in% names(sv))) return(FALSE)
  idx <- which(sv$name == shock)
  if (length(idx) == 0L) return(FALSE)
  pat <- paste0("(^|[^A-Za-z0-9_.])", gsub(".", "\\.", par, fixed = TRUE),
                "($|[^A-Za-z0-9_.])")
  for (col in intersect(c("stderr_expr", "variance_expr"), names(sv))) {
    ex <- as.character(sv[[col]][idx])
    ex <- ex[!is.na(ex)]
    if (length(ex) && any(grepl(pat, ex))) return(TRUE)
  }
  FALSE
}

## Canonical set of prior distributions the evaluators understand. Single source
## of truth shared by log_prior() / log_prior_density() / the analytic gradient /
## the SMC prior sampler / the parameter-transform support resolver. `inv_gamma1`
## is an accepted alias of `inv_gamma`.
.dynhr_prior_dists <- function()
  c("beta", "gamma", "normal", "inv_gamma", "inv_gamma1", "inv_gamma2", "uniform")

## Normalise a distribution name to canonical form: lower-case, trimmed, with
## Dynare's `_pdf` suffix stripped ("beta_pdf" / "Beta" / "beta " all -> "beta").
## Shared by validate_prior_spec(), log_prior(), log_prior_density() and the
## analytic-gradient prior score so they accept the same forms.
## extract_prior_spec() emits already-canonical names, but specs hand-built and
## passed straight to the evaluators (tests, user code) may carry the raw
## Dynare name.
##
## HOT-PATH NOTE (2026-08-05, Windows benchmark diagnosis): this used to call
## trimws(), whose internal perl = TRUE sub() costs ~200us PER CALL on some
## Windows builds (vs single-digit us on macOS) -- and log_prior() called it
## once per parameter per posterior evaluation, dominating dynhr_benchmark()
## there. Two defenses now: (1) the fast path below returns already-canonical
## input untouched with no regex at all (the extract_prior_spec() case, i.e.
## every in-package caller); (2) the slow path avoids trimws in favor of a
## single non-perl gsub. Keep both properties if editing.
.normalize_dist <- function(d) {
  d <- as.character(d)
  if (all(d %in% .dynhr_prior_dists())) return(d)
  sub("_pdf$", "", tolower(gsub("^[ \t\r\n]+|[ \t\r\n]+$", "", d)))
}

#' Validate a prior_spec data.frame, failing loud on a malformed / typo'd spec
#'
#' Catches at CONSTRUCTION the silent-failure class the S7 scoping report flagged:
#' an unrecognised / misspelled distribution name (otherwise a silent improper
#' flat prior on the hot path), duplicate / empty names, missing columns, and
#' degenerate support / scale. Cheap -- called once per estimation setup, never in
#' the per-draw loop. Only flags values that can only be mistakes; specs that
#' legitimately defer a bound to a downstream default (NA p1/p2/lower/upper) are
#' left alone.
#'
#' @param prior_spec A prior_spec data.frame (from extract_prior_spec()).
#' @param where Character label for the error message's call site.
#' @return \code{prior_spec}, invisibly.
#' @noRd
validate_prior_spec <- function(prior_spec, where = "prior_spec") {
  if (!is.data.frame(prior_spec))
    stop(where, ": prior_spec must be a data.frame; got '",
         class(prior_spec)[1], "'.", call. = FALSE)
  req  <- c("name", "distribution", "p1", "p2", "lower", "upper")
  miss <- setdiff(req, names(prior_spec))
  if (length(miss) > 0L)
    stop(where, ": prior_spec is missing required column(s): ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  if (nrow(prior_spec) == 0L)
    stop(where, ": prior_spec has no rows (no estimated parameters).",
         call. = FALSE)

  nm <- as.character(prior_spec$name)
  if (anyNA(nm) || any(!nzchar(trimws(nm))))
    stop(where, ": every prior must have a non-empty name.", call. = FALSE)
  dup <- unique(nm[duplicated(nm)])
  if (length(dup) > 0L)
    stop(where, ": duplicate prior name(s): ", paste(dup, collapse = ", "),
         ".", call. = FALSE)

  ok  <- .dynhr_prior_dists()
  d   <- .normalize_dist(prior_spec$distribution)
  bad <- !(d %in% ok)
  if (any(bad))
    stop(where, ": unsupported prior distribution(s): ",
         paste0(nm[bad], " (\"", prior_spec$distribution[bad], "\")",
                collapse = "; "),
         ". Supported: ", paste(ok, collapse = ", "),
         ". A misspelled name would otherwise be given a silent flat (improper) ",
         "prior.", call. = FALSE)

  ## Per-row support / scale sanity. Each check fires only when the value is
  ## present (non-NA) and can only be a mistake -- so valid specs are untouched.
  p1 <- suppressWarnings(as.numeric(prior_spec$p1))
  p2 <- suppressWarnings(as.numeric(prior_spec$p2))
  lo <- suppressWarnings(as.numeric(prior_spec$lower))
  hi <- suppressWarnings(as.numeric(prior_spec$upper))
  g3 <- suppressWarnings(as.numeric(prior_spec$p3 %||% rep(NA_real_, nrow(prior_spec))))
  g4 <- suppressWarnings(as.numeric(prior_spec$p4 %||% rep(NA_real_, nrow(prior_spec))))
  for (i in seq_len(nrow(prior_spec))) {
    who <- nm[i]; di <- d[i]
    if (!is.na(lo[i]) && !is.na(hi[i]) && lo[i] > hi[i])
      stop(where, ": prior '", who, "' has lower (", lo[i], ") > upper (",
           hi[i], ").", call. = FALSE)
    if (di == "uniform") {
      ## Dynare convention: support [p3, p4] when both are given, else
      ## p1 -+ sqrt(3) p2 (p1 = mean, p2 = sd).
      if (!is.na(g3[i]) && !is.na(g4[i])) {
        if (!(g3[i] < g4[i]))
          stop(where, ": uniform prior '", who, "' needs support p3 < p4 ",
               "(got [", g3[i], ", ", g4[i], "]).", call. = FALSE)
        ## p1/p2 alongside the bounds must BE the mean/sd of [p3, p4] --
        ## catches a hand-built spec still carrying bounds in p1/p2.
        m <- (g3[i] + g4[i]) / 2; sdv <- (g4[i] - g3[i]) / sqrt(12)
        tol <- 1e-8 * max(1, abs(g3[i]), abs(g4[i]))
        if ((!is.na(p1[i]) && abs(p1[i] - m) > tol) ||
            (!is.na(p2[i]) && abs(p2[i] - sdv) > tol))
          stop(where, ": uniform prior '", who, "' on [p3, p4] = [", g3[i],
               ", ", g4[i], "] has p1/p2 = (", p1[i], ", ", p2[i], "), which ",
               "are not its mean/sd (", m, ", ", sdv, "). Dynare convention: ",
               "p1 = mean, p2 = sd -- leave them NA or set them to those ",
               "values.", call. = FALSE)
      } else if (!is.na(p2[i]) && p2[i] <= 0) {
        stop(where, ": uniform prior '", who, "' needs p2 (standard deviation, ",
             "Dynare convention) > 0 (got ", p2[i], "); for bounds use p3/p4.",
             call. = FALSE)
      }
      next
    }
    if (!is.na(p2[i]) && p2[i] <= 0)
      stop(where, ": prior '", who, "' (", di, ") needs p2 (sd/scale) > 0 (got ",
           p2[i], ").", call. = FALSE)
    ## Generalised beta on [p3, p4] / shift p3 (Dynare; NA = 0 / 1).
    ab <- .prior_gen_ab(g3[i], g4[i])
    if (di == "beta" && !is.na(p1[i]) && (p1[i] <= ab[1] || p1[i] >= ab[2]))
      stop(where, ": beta prior '", who, "' needs mean p1 in (", ab[1], ", ",
           ab[2], ") (got ", p1[i], ").", call. = FALSE)
    if (di %in% c("gamma", "inv_gamma", "inv_gamma1", "inv_gamma2") &&
        !is.na(p1[i]) && p1[i] <= ab[1])
      stop(where, ": ", di, " prior '", who, "' needs mean p1 > ",
           if (ab[1] == 0) "0" else paste0("the shift p3 = ", ab[1]),
           " (got ", p1[i], ").", call. = FALSE)
  }
  invisible(prior_spec)
}

.default_lower <- function(dist) {
  switch(dist,
    inv_gamma = 1e-8, inv_gamma2 = 1e-8,
    beta = 0, gamma = 0, normal = -Inf,
    0)
}

.default_upper <- function(dist) {
  switch(dist,
    beta = 1, normal = Inf, gamma = Inf,
    inv_gamma = Inf, inv_gamma2 = Inf,
    Inf)
}

.shock_to_sig <- function(shock_nm) {
  core <- sub("^eps_", "", shock_nm)
  core <- sub("^e_", "", core)
  core <- sub("_$", "", core)
  paste0("sig_", core)
}
