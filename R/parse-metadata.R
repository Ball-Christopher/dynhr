## R/parse-metadata.R
## --------------------------------------------------------------------------
## @dynhr:* metadata block extractors: shock_categories, benchmarks,
## variable_labels, narratives. These are dynhr-specific extensions to the
## Dynare .mod format that drive the diagnostic suite.
##
## Phase-1 split from parser-monolith.R (no logic changes).
## --------------------------------------------------------------------------

#' @export
extract_mod_metadata <- function(mod_file_or_lines) {
  txt <- if (length(mod_file_or_lines) == 1 && file.exists(mod_file_or_lines))
    readLines(mod_file_or_lines) else mod_file_or_lines

  # -- Parse %(key='value', ...) annotations from var/varexo lines --
  pct_pattern <- "([a-zA-Z_][a-zA-Z0-9_]*)\\s+%\\((.+?)\\)"
  var_annotations <- list()

  for (line in txt) {
    ms <- gregexpr(pct_pattern, line, perl = TRUE)
    hits <- regmatches(line, ms)[[1]]
    for (hit in hits) {
      parts <- regmatches(hit, regexec(pct_pattern, hit, perl = TRUE))[[1]]
      varname <- parts[2]
      kvtext  <- parts[3]
      kv_re   <- "([a-zA-Z_][a-zA-Z0-9_]*)\\s*=\\s*'([^']*)'"
      kvs     <- regmatches(kvtext, gregexpr(kv_re, kvtext, perl = TRUE))[[1]]
      meta    <- list()
      for (kv in kvs) {
        kv_parts <- regmatches(kv, regexec(kv_re, kv, perl = TRUE))[[1]]
        meta[[kv_parts[2]]] <- kv_parts[3]
      }
      var_annotations[[varname]] <- meta
    }
  }

  # -- Parse // @dynhr:<block_name> ... // @dynhr:end blocks --
  dynhr_blocks <- list()
  in_block <- FALSE
  block_name <- NULL
  block_lines <- character()

  for (line in txt) {
    stripped <- trimws(line)
    tag_match <- regmatches(stripped, regexec("^//\\s*@dynhr:(\\w+)", stripped))[[1]]
    if (length(tag_match) == 2) {
      tag <- tag_match[2]
      if (tag == "end") {
        if (in_block) {
          dynhr_blocks[[block_name]] <- block_lines
          in_block <- FALSE
          block_lines <- character()
        }
      } else {
        in_block <- TRUE
        block_name <- tag
        block_lines <- character()
      }
    } else if (in_block && grepl("^//", stripped)) {
      content <- sub("^//\\s*", "", stripped)
      if (nzchar(content)) block_lines <- c(block_lines, content)
    }
  }

  # -- Parse shock_categories --
  shock_categories <- list()      # category -> c(shock_names)
  shock_to_category <- character() # shock_name -> category
  if ("shock_categories" %in% names(dynhr_blocks)) {
    for (line in dynhr_blocks$shock_categories) {
      if (grepl(":", line)) {
        parts    <- strsplit(line, ":\\s*", perl = TRUE)[[1]]
        category <- trimws(parts[1])
        shocks   <- trimws(strsplit(parts[2], ",\\s*")[[1]])
        shock_categories[[category]] <- shocks
        for (s in shocks) shock_to_category[s] <- category
      }
    }
  }

  # -- Parse benchmarks --
  benchmarks <- list()
  if ("benchmarks" %in% names(dynhr_blocks)) {
    for (line in .fold_continuation_lines(dynhr_blocks$benchmarks)) {
      colon_pos <- regexpr(":", line, fixed = TRUE)
      if (colon_pos < 1) next
      bname <- trimws(substr(line, 1, colon_pos - 1))
      rest  <- trimws(substr(line, colon_pos + 1, nchar(line)))
      # Parse key=value or key="value with spaces"
      kv_re <- '([a-zA-Z_]+)\\s*=\\s*"([^"]+)"|([a-zA-Z_]+)\\s*=\\s*([^,]+)'
      raw_kvs <- gregexpr(kv_re, rest, perl = TRUE)
      kv_hits <- regmatches(rest, raw_kvs)[[1]]
      bm <- list()
      for (kv in kv_hits) {
        # Try quoted form first
        m1 <- regmatches(kv, regexec('([a-zA-Z_]+)\\s*=\\s*"([^"]+)"', kv, perl = TRUE))[[1]]
        if (length(m1) == 3) { bm[[m1[2]]] <- m1[3]; next }
        m2 <- regmatches(kv, regexec('([a-zA-Z_]+)\\s*=\\s*([^,]+)', kv, perl = TRUE))[[1]]
        if (length(m2) == 3) bm[[m2[2]]] <- trimws(m2[3])
      }
      # L21: warn and skip if no recognised benchmark key=value pairs were parsed
      .bm_known_keys <- c("variable", "shock", "sign", "min", "max", "type",
                          "peak_horizon", "peak_magnitude", "description")
      if (length(bm) == 0 || !any(names(bm) %in% .bm_known_keys)) {
        warning(sprintf("@dynhr:benchmarks: could not parse '%s' as key=value fields; skipped.", bname),
                call. = FALSE)
        next
      }
      benchmarks[[bname]] <- bm
    }
  }

  # -- Parse variable_labels --
  var_labels <- list()
  if ("variable_labels" %in% names(dynhr_blocks)) {
    for (line in dynhr_blocks$variable_labels) {
      colon_pos <- regexpr(":", line, fixed = TRUE)
      if (colon_pos < 1) next
      vname <- trimws(substr(line, 1, colon_pos - 1))
      label <- trimws(substr(line, colon_pos + 1, nchar(line)))
      var_labels[[vname]] <- label
    }
  }
  # Merge in %(long_name=...) annotations (lower priority)
  for (vname in names(var_annotations)) {
    if (is.null(var_labels[[vname]]) &&
        "long_name" %in% names(var_annotations[[vname]])) {
      var_labels[[vname]] <- var_annotations[[vname]]$long_name
    }
  }

  # -- Parse expectations --
  expectations <- .parse_expectations_block(dynhr_blocks)

  # -- Parse deep-parameter taxonomy --
  deep <- .parse_deep_block(dynhr_blocks)

  structure(
    list(
      annotations       = var_annotations,
      shock_categories  = shock_categories,
      shock_to_category = shock_to_category,
      benchmarks        = benchmarks,
      var_labels        = var_labels,
      narratives        = .extract_narratives(mod_file_or_lines),
      expectations      = expectations,
      deep              = deep
    ),
    class = "dynhr_metadata"
  )
}

# Helper: get display label for a variable
var_label <- function(metadata, varname, with_code = TRUE) {
  lab <- metadata$var_labels[[varname]]
  if (is.null(lab)) return(varname)
  if (with_code) sprintf("%s (%s)", lab, varname) else lab
}


# ---------------------------------------------------------------------------
#' Extract @dynhr:narratives block from .mod file
#'
#' Parses narrative episode definitions. Each line has format:
#'   name: date_range=START:END, shock=SHOCK_OR_CATEGORY, variable=VAR,
#'         sign=SIGN, min_sd=N, description="TEXT"
#'
#' @param mod_file  Path to .mod file
#' @return list of narrative episode specs
#' @noRd
# ---------------------------------------------------------------------------
.extract_narratives <- function(mod_file) {
  lines <- readLines(mod_file, warn = FALSE)
  in_block <- FALSE
  raw <- character()

  for (line in lines) {
    stripped <- trimws(sub("^//\\s*", "", line))
    if (grepl("@dynhr:narratives", stripped, fixed = TRUE)) {
      in_block <- TRUE; next
    }
    if (in_block && grepl("@dynhr:end", stripped, fixed = TRUE)) {
      in_block <- FALSE; next
    }
    if (in_block && nchar(stripped) > 0 && !grepl("^Syntax:|^\\s*-\\s", stripped)) {
      raw <- c(raw, stripped)
    }
  }

  raw <- raw[grepl("date_range\\s*=", raw)]

  if (length(raw) == 0) return(list())

  lapply(raw, function(r) {
    # Split name from rest: "gfc_demand: date_range=..."
    colon_pos <- regexpr(":", r)
    if (colon_pos < 1) return(NULL)

    ep_name <- trimws(substr(r, 1, colon_pos - 1))
    rest    <- trimws(substr(r, colon_pos + 1, nchar(r)))

    # Parse key=value pairs (handle description="..." with quotes)
    # First extract description if present (it may contain commas)
    description <- ""
    desc_match <- regmatches(rest, regexpr('description\\s*=\\s*"[^"]*"', rest))
    if (length(desc_match) > 0) {
      description <- sub('.*description\\s*=\\s*"([^"]*)".*', "\\1", desc_match)
      rest <- sub(',?\\s*description\\s*=\\s*"[^"]*"', "", rest)
    }

    # Parse remaining key=value pairs
    pairs <- strsplit(rest, ",")[[1]]
    kv <- list()
    for (p in pairs) {
      p <- trimws(p)
      eq <- regexpr("=", p)
      if (eq > 0) {
        key <- trimws(substr(p, 1, eq - 1))
        val <- trimws(substr(p, eq + 1, nchar(p)))
        kv[[key]] <- val
      }
    }

    # Parse date_range: "1997-Q3:1998-Q2" -> c("1997-Q3", "1998-Q2")
    date_range <- NULL
    if (!is.null(kv$date_range)) {
      date_range <- trimws(strsplit(kv$date_range, ":")[[1]])
    }

    list(
      name        = ep_name,
      date_range  = date_range,
      shock       = trimws(kv$shock %||% ""),
      variable    = trimws(kv$variable %||% "y"),
      sign        = trimws(kv$sign %||% "negative"),
      min_sd      = as.numeric(kv$min_sd %||% "1.0"),
      description = description
    )
  }) |> Filter(Negate(is.null), x = _)
}


# ---------------------------------------------------------------------------
# Fold continuation lines in a @dynhr:* block.
#
# A line starts a new entry iff it has a "name:" header — i.e. a colon that
# appears before the first "=".  Any line that does NOT have such a header is
# a continuation of the previous entry and is appended (comma-joined) to it.
# If the very first line has no header it is left as-is (will be caught by the
# colon check in the caller and produce a warning).
# ---------------------------------------------------------------------------
.fold_continuation_lines <- function(lines) {
  if (length(lines) == 0) return(lines)
  is_new_entry <- function(line) {
    colon <- regexpr(":", line, fixed = TRUE)
    eq    <- regexpr("=", line, fixed = TRUE)
    # A new entry must have a colon AND that colon precedes the first "="
    colon > 0 && (eq < 1 || colon < eq)
  }
  folded <- character()
  for (line in lines) {
    if (is_new_entry(line) || length(folded) == 0) {
      folded <- c(folded, line)
    } else {
      folded[length(folded)] <- paste0(folded[length(folded)], ", ", trimws(line))
    }
  }
  folded
}

# ---------------------------------------------------------------------------
# @dynhr:expectations block parser
#
# Each line in the block has the form:
#   check_name: type="<type>", <key>=<value>, ..., description="<text>"
#
# Supported types:
#   data_ratio   -- ratio of two data-column means: numerator / denominator
#   data_mean    -- mean of a data column
#   data_sd      -- standard deviation of a data column
#   param_range  -- calibrated / estimated parameter falls within range
#   irf_sign     -- sign of IRF for a given variable/shock at a given horizon
#
# Every check must have min= and/or max= fields; description= is optional.
# ---------------------------------------------------------------------------
.parse_expectations_block <- function(dynhr_blocks) {
  if (!("expectations" %in% names(dynhr_blocks))) return(list())

  kv_re <- '([a-zA-Z_][a-zA-Z0-9_]*)\\s*=\\s*"([^"]+)"|([a-zA-Z_][a-zA-Z0-9_]*)\\s*=\\s*([^,\\s][^,]*)'

  lapply(.fold_continuation_lines(dynhr_blocks$expectations), function(line) {
    colon_pos <- regexpr(":", line, fixed = TRUE)
    if (colon_pos < 1) return(NULL)

    check_name <- trimws(substr(line, 1, colon_pos - 1))
    rest       <- trimws(substr(line, colon_pos + 1, nchar(line)))

    raw_kvs <- regmatches(rest, gregexpr(kv_re, rest, perl = TRUE))[[1]]
    kv <- list()
    for (hit in raw_kvs) {
      m1 <- regmatches(hit, regexec('([a-zA-Z_][a-zA-Z0-9_]*)\\s*=\\s*"([^"]+)"', hit, perl = TRUE))[[1]]
      if (length(m1) == 3) { kv[[m1[2]]] <- m1[3]; next }
      m2 <- regmatches(hit, regexec('([a-zA-Z_][a-zA-Z0-9_]*)\\s*=\\s*([^,\\s][^,]*)', hit, perl = TRUE))[[1]]
      if (length(m2) == 3) kv[[m2[2]]] <- trimws(m2[3])
    }

    # L21: warn and skip if type is missing or no key=value pairs were parsed
    if (length(kv) == 0 || is.null(kv$type)) {
      warning(sprintf("@dynhr:expectations: could not parse '%s' (missing type= or no key=value fields); skipped.",
                      check_name), call. = FALSE)
      return(NULL)
    }

    spec <- list(
      name        = check_name,
      type        = kv$type,
      description = kv$description %||% check_name
    )

    # Numeric fields
    if (!is.null(kv$min)) spec$min <- as.numeric(kv$min)
    if (!is.null(kv$max)) spec$max <- as.numeric(kv$max)

    # Variable / reference fields
    if (!is.null(kv$variable))    spec$variable    <- kv$variable
    if (!is.null(kv$numerator))   spec$numerator   <- kv$numerator
    if (!is.null(kv$denominator)) spec$denominator <- kv$denominator
    if (!is.null(kv$shock))       spec$shock       <- kv$shock
    if (!is.null(kv$horizon))     spec$horizon     <- as.integer(kv$horizon)
    if (!is.null(kv$sign))        spec$sign        <- kv$sign

    spec
  }) |> Filter(Negate(is.null), x = _)
}


# ---------------------------------------------------------------------------
# @dynhr:deep block parser  (deep-parameter taxonomy)
#
# Declares the structural class and (optionally) the reduced-form image of each
# model parameter, so the deep-parameter diagnostics (D33/D34/D35) and the
# Deep-Parameter Passport can (a) tell deep primitives from auxiliary shock
# parameters and (b) partition policy vs private parameters for the Lucas-
# critique invariance test.
#
# Each line has the form
#   <name>[ <name> ...]: class=<class>, role=<role>, reduced_form=<k>, from=<...>
#
#   class  : preference | technology | rigidity | policy | labour | open |
#            shock | unknown
#   role   : primitive | reduced_form | auxiliary
#   reduced_form : (primitives only) the data-facing coefficient this primitive
#                  maps to, e.g. theta_H -> kappa_H
#   from   : (reduced_form rows) the primitives it is built from, e.g.
#            from="theta_H,beta"
#
# Several parameters that share attributes may be listed on one line, e.g.
#   rho_a sig_a rho_pref sig_pref: class=shock, role=auxiliary
#
# Returns a list keyed by parameter name; each element is a spec list.
# ---------------------------------------------------------------------------
.parse_deep_block <- function(dynhr_blocks) {
  if (!("deep" %in% names(dynhr_blocks))) return(list())

  kv_re <- '([a-zA-Z_][a-zA-Z0-9_]*)\\s*=\\s*"([^"]+)"|([a-zA-Z_][a-zA-Z0-9_]*)\\s*=\\s*([^,\\s][^,]*)'

  specs <- lapply(dynhr_blocks$deep, function(line) {
    colon_pos <- regexpr(":", line, fixed = TRUE)
    if (colon_pos < 1) return(NULL)

    names_part <- trimws(substr(line, 1, colon_pos - 1))
    rest       <- trimws(substr(line, colon_pos + 1, nchar(line)))
    if (!nzchar(names_part)) return(NULL)

    raw_kvs <- regmatches(rest, gregexpr(kv_re, rest, perl = TRUE))[[1]]
    kv <- list()
    for (hit in raw_kvs) {
      m1 <- regmatches(hit, regexec('([a-zA-Z_][a-zA-Z0-9_]*)\\s*=\\s*"([^"]+)"', hit, perl = TRUE))[[1]]
      if (length(m1) == 3) { kv[[m1[2]]] <- m1[3]; next }
      m2 <- regmatches(hit, regexec('([a-zA-Z_][a-zA-Z0-9_]*)\\s*=\\s*([^,\\s][^,]*)', hit, perl = TRUE))[[1]]
      if (length(m2) == 3) kv[[m2[2]]] <- trimws(m2[3])
    }

    base <- list(
      class        = tolower(kv$class %||% "unknown"),
      role         = tolower(kv$role  %||% "primitive"),
      reduced_form = kv$reduced_form %||% NA_character_,
      from         = kv$from %||% NA_character_
    )

    pnames <- strsplit(names_part, "\\s+")[[1]]
    pnames <- pnames[nzchar(pnames)]
    lapply(pnames, function(pn) c(list(name = pn), base))
  })

  # Flatten (each line may yield several params) and key by name.
  flat <- unlist(Filter(Negate(is.null), specs), recursive = FALSE)
  out  <- list()
  for (sp in flat) if (!is.null(sp$name)) out[[sp$name]] <- sp
  out
}
