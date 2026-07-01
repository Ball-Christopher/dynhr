## R/model-registry.R
## --------------------------------------------------------------------------
## Model registry: parse @dynhr-model metadata from .mod files, validate
## against the parsed model structure, and generate the model table.
##
## parse_mod_metadata(path)          -- read @dynhr-model block from a .mod
## check_mod_metadata(model, meta)   -- compare metadata vs. parsed model
## model_registry(model_dir)         -- build data.frame of all registered models
## print_model_registry(model_dir)   -- pretty-print the registry table
## --------------------------------------------------------------------------

## ---- Metadata format -------------------------------------------------------
## Metadata lives inside a @dynhr-model ... @dynhr-model-end comment block:
##
##   // @dynhr-model
##   // name: Real Business Cycle (baseline)
##   // short: rbc
##   // tier: 1
##   // endo: 7
##   // exo: 1
##   // params: 8
##   // has_lead: true
##   // has_measurement: false
##   // has_estimation: false
##   // has_analytical_ss: false
##   // features: perturbation, kalman_filter, stoch_simul
##   // refs: ...
##   // notes: ...
##   // @dynhr-model-end
##
## Scalar boolean keys: has_lead, has_measurement, has_estimation, has_analytical_ss
## Scalar integer keys: tier, endo, exo, params
## Free text keys: name, short, features, refs, notes
## ---------------------------------------------------------------------------

#' Parse @dynhr-model metadata from a .mod file
#'
#' Returns NULL (with a message) if no metadata block is found.
#'
#' @param path  Path to a .mod file
#' @return Named list of metadata fields, or NULL
#' @export
parse_mod_metadata <- function(path) {
  if (!file.exists(path))
    stop(sprintf("File not found: %s", path))

  lines <- readLines(path, warn = FALSE)

  start_idx <- which(grepl("@dynhr-model\\b", lines, perl = TRUE) &
                     !grepl("@dynhr-model-end", lines))
  end_idx   <- which(grepl("@dynhr-model-end", lines))

  if (length(start_idx) == 0 || length(end_idx) == 0) {
    message(sprintf("  No @dynhr-model block in %s", basename(path)))
    return(NULL)
  }

  block <- lines[(start_idx[1] + 1):(end_idx[1] - 1)]

  meta <- list(path = path, file = basename(path))

  for (ln in block) {
    ln <- trimws(ln)
    ln <- sub("^//\\s*", "", ln)
    if (!grepl("^[a-z_]+\\s*:", ln)) next

    key <- trimws(sub(":.*", "", ln))
    val <- trimws(sub("^[^:]+:\\s*", "", ln))

    if (key %in% c("tier", "endo", "exo", "params")) {
      meta[[key]] <- as.integer(val)
    } else if (key %in% c("has_lead", "has_measurement", "has_estimation",
                           "has_analytical_ss")) {
      meta[[key]] <- tolower(val) %in% c("true", "yes", "1")
    } else {
      meta[[key]] <- val
    }
  }

  meta
}


#' Validate metadata against a parsed model object
#'
#' Compares the counts and flags in a @dynhr-model block against what
#' parse_mod() actually found.  Returns a data.frame of check results.
#'
#' @param model     dynhr_mod from parse_mod(), or the path to a .mod file
#'   (then parsed with parse_mod() and meta read from its @dynhr-model block)
#' @param meta      Metadata list from parse_mod_metadata(); ignored when
#'   model is a path
#' @param warn      If TRUE, emit warnings for mismatches (default TRUE)
#' @return Invisible data.frame: check, expected, found, ok
#' @export
check_mod_metadata <- function(model, meta, warn = TRUE) {
  if (is.character(model) && length(model) == 1L) {
    meta  <- parse_mod_metadata(model)
    model <- parse_mod(model, verbose = FALSE)
  }
  # When no @dynhr-model block is present, return an empty result rather
  # than erroring — the caller can inspect model$metadata for @dynhr: blocks.
  if (is.null(meta)) {
    return(invisible(data.frame(check = character(0), expected = character(0),
                                found = character(0), ok = logical(0),
                                stringsAsFactors = FALSE)))
  }
  checks <- list()

  .chk <- function(name, expected, found) {
    ok <- isTRUE(all.equal(expected, found))
    checks[[length(checks) + 1]] <<- data.frame(
      check    = name,
      expected = as.character(expected),
      found    = as.character(found),
      ok       = ok,
      stringsAsFactors = FALSE
    )
  }

  n_endo   <- length(model$var_names    %||% model$endogenous %||% character(0))
  n_exo    <- length(model$varexo_names %||% model$exogenous  %||% character(0))
  n_params <- length(model$param_values %||% numeric(0))

  if (!is.null(meta$endo))   .chk("endo",   meta$endo,   n_endo)
  if (!is.null(meta$exo))    .chk("exo",    meta$exo,    n_exo)
  if (!is.null(meta$params)) .chk("params", meta$params, n_params)

  eqs <- model$equations %||% model$model_equations
  if (is.null(eqs) || length(eqs) == 0) { has_lead_found <- NA } else {
    eq_strs <- vapply(eqs, function(e) paste(deparse(e), collapse = ""), character(1))
    has_lead_found <- any(grepl("(+1)", eq_strs, fixed = TRUE))
  }

  if (!is.null(meta$has_lead) && !is.na(has_lead_found))
    .chk("has_lead", meta$has_lead, has_lead_found)

  ep <- model$estimated_params
  has_est_found <- !is.null(ep) &&
                   ((is.data.frame(ep) && nrow(ep) > 0) ||
                    (is.list(ep) && !is.data.frame(ep) && length(ep) > 0))

  if (!is.null(meta$has_estimation))
    .chk("has_estimation", meta$has_estimation, has_est_found)

  result <- do.call(rbind, checks)

  if (warn && !is.null(result) && any(!result$ok)) {
    bad <- result[!result$ok, ]
    for (i in seq_len(nrow(bad))) {
      warning(sprintf(
        "Metadata mismatch in '%s': %s expected=%s found=%s",
        meta$file %||% "?", bad$check[i], bad$expected[i], bad$found[i]
      ))
    }
  }

  invisible(result)
}


#' Build a model registry data.frame from all .mod files in a directory
#'
#' Reads metadata from every .mod file that has a @dynhr-model block.
#'
#' @param model_dir  Directory to search (default: package extdata/models)
#' @param recursive  Search subdirectories (default FALSE)
#' @return data.frame with one row per model
#' @export
model_registry <- function(model_dir = NULL, recursive = FALSE) {
  if (is.null(model_dir))
    model_dir <- system.file("extdata", "models", package = "dynhr")

  mods <- list.files(model_dir, pattern = "\\.mod$",
                     full.names = TRUE, recursive = recursive)

  rows <- lapply(mods, function(p) {
    meta <- parse_mod_metadata(p)
    if (is.null(meta)) return(NULL)

    data.frame(
      tier             = meta$tier %||% NA_integer_,
      short            = meta$short %||% tools::file_path_sans_ext(basename(p)),
      name             = meta$name %||% NA_character_,
      endo             = meta$endo  %||% NA_integer_,
      exo              = meta$exo   %||% NA_integer_,
      params           = meta$params %||% NA_integer_,
      has_lead         = meta$has_lead         %||% NA,
      has_measurement  = meta$has_measurement  %||% NA,
      has_estimation   = meta$has_estimation   %||% NA,
      has_analytical_ss = meta$has_analytical_ss %||% NA,
      features         = meta$features %||% NA_character_,
      refs             = meta$refs  %||% NA_character_,
      notes            = meta$notes %||% NA_character_,
      file             = basename(p),
      stringsAsFactors = FALSE
    )
  })

  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0)
    return(data.frame())

  reg <- do.call(rbind, rows)
  reg[order(reg$tier, reg$short), ]
}


#' Print the model registry as a markdown table
#'
#' @param model_dir  Directory of .mod files (default: package extdata/models),
#'   or a registry data.frame from \code{model_registry()} (printed as-is,
#'   no re-scan)
#' @param output     "cat" (default) or "return" (character vector)
#' @export
print_model_registry <- function(model_dir = NULL, output = "cat") {
  reg <- if (is.data.frame(model_dir)) model_dir else model_registry(model_dir)
  if (nrow(reg) == 0) {
    cat("No models with @dynhr-model metadata found.\n")
    return(invisible(character(0)))
  }

  .yn <- function(x) ifelse(is.na(x), "?", ifelse(x, "yes", "no"))
  .n  <- function(x) ifelse(is.na(x), "?", as.character(x))

  hdr <- sprintf("| %-5s | %-20s | %5s | %4s | %6s | %5s | %5s | %5s |\n",
                 "Tier", "Short", "Endo", "Exo", "Params",
                 "Lead", "Meas", "Est")
  sep <- paste0("|", paste(rep("-", 7), collapse=""), "|",
                paste(rep("-", 22), collapse=""), "|",
                paste(rep("-", 7), collapse=""), "|",
                paste(rep("-", 6), collapse=""), "|",
                paste(rep("-", 8), collapse=""), "|",
                paste(rep("-", 7), collapse=""), "|",
                paste(rep("-", 7), collapse=""), "|",
                paste(rep("-", 7), collapse=""), "|\n")

  rows_txt <- vapply(seq_len(nrow(reg)), function(i) {
    r <- reg[i, ]
    sprintf("| %-5s | %-20s | %5s | %4s | %6s | %5s | %5s | %5s |\n",
            .n(r$tier), r$short, .n(r$endo), .n(r$exo), .n(r$params),
            .yn(r$has_lead), .yn(r$has_measurement), .yn(r$has_estimation))
  }, character(1))

  lines <- c(hdr, sep, rows_txt)

  if (output == "cat") {
    cat(lines)
    invisible(lines)
  } else {
    lines
  }
}
