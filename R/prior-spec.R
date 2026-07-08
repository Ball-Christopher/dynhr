## R/prior-spec.R
## --------------------------------------------------------------------------
## Phase-2 split from estimate-monolith.R.
##
## extract_prior_spec(): build the prior_spec data.frame from a parsed model.
## Helper functions for shock-to-sigma renaming and distribution defaults.
## --------------------------------------------------------------------------

#' Extract prior specification from a parsed model
#'
#' Reads the `estimated_params` block from the parsed model and returns a
#' data.frame suitable for `make_log_posterior()`.
#'
#' Applies the eps_* -> sig_* rename for `stderr` priors (the model's shock
#' standard deviation parameters use the sig_ prefix, while the estimated_params
#' block uses the shock name eps_*).
#'
#' @param model   dynhr_mod (from parse_mod())
#' @param verbose Logical. Print extracted prior info to console (default TRUE).
#' @return data.frame with columns: name, distribution, p1, p2, lower, upper,
#'         mean, std
#' @noRd
extract_prior_spec <- function(model, verbose = TRUE) {
  ep <- model$estimated_params
  if (is.null(ep) || nrow(ep) == 0) {
    stop("No estimated_params found in model (add an estimated_params block ",
         "to the .mod file, or supply model$estimated_params directly)")
  }

  cnames <- tolower(names(ep))
  if ("name" %in% cnames && "prior" %in% cnames) {
    col_name <- which(cnames == "name")[1]
    col_dist <- which(cnames == "prior")[1]
    col_p1   <- which(cnames == "p1")[1]
    col_p2   <- which(cnames == "p2")[1]
    col_p3   <- which(cnames == "p3")[1]
    col_p4   <- which(cnames == "p4")[1]
  } else {
    col_name <- 1; col_dist <- 2; col_p1 <- 3; col_p2 <- 4
    col_p3 <- if (ncol(ep) >= 5) 5 else NA
    col_p4 <- if (ncol(ep) >= 6) 6 else NA
  }

  priors <- data.frame(
    name = character(), distribution = character(),
    p1 = numeric(), p2 = numeric(), lower = numeric(), upper = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(ep))) {
    nm        <- trimws(as.character(ep[i, col_name]))
    dist      <- tolower(trimws(as.character(ep[i, col_dist])))
    dist_base <- sub("_pdf$", "", dist)

    vals <- c(
      as.numeric(as.character(ep[i, col_p1])),
      as.numeric(as.character(ep[i, col_p2])),
      if (!is.na(col_p3)) as.numeric(as.character(ep[i, col_p3])) else NA,
      if (!is.na(col_p4)) as.numeric(as.character(ep[i, col_p4])) else NA
    )
    vals <- vals[!is.na(vals)]

    if (dist_base == "uniform") {
      lo <- if (length(vals) >= 1) vals[1] else 0
      hi <- if (length(vals) >= 2) vals[2] else 1
      priors[nrow(priors) + 1L, ] <- list(nm, dist_base, lo, hi, lo, hi)
    } else {
      p1 <- if (length(vals) >= 1) vals[1] else NA
      p2 <- if (length(vals) >= 2) vals[2] else NA
      lo <- if (length(vals) >= 3) vals[3] else .default_lower(dist_base)
      hi <- if (length(vals) >= 4) vals[4] else .default_upper(dist_base)
      priors[nrow(priors) + 1L, ] <- list(nm, dist_base, p1, p2, lo, hi)
    }
  }

  ## Merge bounds from estimated_params_bounds block.
  if (!is.null(model$estimated_params_bounds)) {
    bnd <- model$estimated_params_bounds
    bnd_cnames   <- tolower(names(bnd))
    bnd_name_col <- if ("name" %in% bnd_cnames) which(bnd_cnames == "name")[1] else 1
    for (i in seq_len(nrow(bnd))) {
      bnm <- trimws(as.character(bnd[i, bnd_name_col]))
      idx <- match(bnm, priors$name)
      if (!is.na(idx)) {
        bvals <- suppressWarnings(as.numeric(as.character(
          bnd[i, setdiff(seq_len(ncol(bnd)), bnd_name_col)])))
        bvals <- bvals[!is.na(bvals)]
        if (length(bvals) >= 1) priors$lower[idx] <- bvals[1]
        if (length(bvals) >= 2) priors$upper[idx] <- bvals[2]
      }
    }
  }

  ## Auto-rename: if 'eps_X' is not a model parameter but 'sig_X' is,
  ## swap the prior name to 'sig_X' so make_log_posterior() can look it up.
  param_names <- names(model$param_values) %||% character(0)
  if (length(param_names) > 0L) {
    rename_map <- character(0)
    for (i in seq_len(nrow(priors))) {
      nm <- priors$name[i]
      if (!(nm %in% param_names) && grepl("^eps_", nm)) {
        sig_nm <- .shock_to_sig(nm)
        if (sig_nm %in% param_names) {
          rename_map[nm] <- sig_nm
          priors$name[i] <- sig_nm
        }
      }
    }
    if (length(rename_map) > 0L && verbose) {
      cat(sprintf("  Renamed %d prior(s) from 'stderr eps_*' to model param 'sig_*':\n",
                  length(rename_map)))
      for (k in names(rename_map)) cat(sprintf("    %s -> %s\n", k, rename_map[[k]]))
    }
  }

  priors$mean <- ifelse(priors$distribution == "uniform",
                        (priors$p1 + priors$p2) / 2, priors$p1)
  priors$std  <- ifelse(priors$distribution == "uniform",
                        (priors$p2 - priors$p1) / sqrt(12), priors$p2)

  ## Fail loud at CONSTRUCTION on a malformed / typo'd spec (unknown distribution,
  ## duplicate names, degenerate support) rather than letting it surface as a
  ## silent flat prior or a NaN density deep inside an MCMC run.
  validate_prior_spec(priors, where = "extract_prior_spec")

  if (verbose) {
    cat(sprintf("  Extracted %d priors:\n", nrow(priors)))
    print(priors[, c("name", "distribution", "p1", "p2", "lower", "upper")],
          row.names = FALSE, right = FALSE)
    cat("\n")
  }
  priors
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
## analytic-gradient prior score so they accept the same forms. extract_prior_spec()
## already strips _pdf, but specs hand-built and passed straight to the evaluators
## (tests, user code) may carry the raw Dynare name.
.normalize_dist <- function(d) sub("_pdf$", "", tolower(trimws(as.character(d))))

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
#' @return `prior_spec`, invisibly.
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
  for (i in seq_len(nrow(prior_spec))) {
    who <- nm[i]; di <- d[i]
    if (!is.na(lo[i]) && !is.na(hi[i]) && lo[i] > hi[i])
      stop(where, ": prior '", who, "' has lower (", lo[i], ") > upper (",
           hi[i], ").", call. = FALSE)
    if (di == "uniform") {
      if (!is.na(p1[i]) && !is.na(p2[i]) && p1[i] >= p2[i])
        stop(where, ": uniform prior '", who, "' needs lower < upper (got ",
             p1[i], ", ", p2[i], ").", call. = FALSE)
      next
    }
    if (!is.na(p2[i]) && p2[i] <= 0)
      stop(where, ": prior '", who, "' (", di, ") needs p2 (sd/scale) > 0 (got ",
           p2[i], ").", call. = FALSE)
    if (di == "beta" && !is.na(p1[i]) && (p1[i] <= 0 || p1[i] >= 1))
      stop(where, ": beta prior '", who, "' needs mean p1 in (0, 1) (got ",
           p1[i], ").", call. = FALSE)
    if (di %in% c("gamma", "inv_gamma", "inv_gamma1", "inv_gamma2") &&
        !is.na(p1[i]) && p1[i] <= 0)
      stop(where, ": ", di, " prior '", who, "' needs mean p1 > 0 (got ",
           p1[i], ").", call. = FALSE)
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
