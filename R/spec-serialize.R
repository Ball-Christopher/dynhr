## R/spec-serialize.R
## ---------------------------------------------------------------------------
## write_spec() / read_spec(): estimation specs as YAML, JSON or RDS.
##
## RDS is canonical (the spec by value, model and data included). YAML and
## JSON go through a schema-driven DOCUMENT layer: the spec becomes a plain
## tree of scalars, sequences and mappings, and is rebuilt from it with every
## field's type taken from the schema (R/estimation-spec.R), not guessed from
## the text. Two encodings are used:
##
##   * a typed field whose value is a plain instance of its schema type is
##     written bare (`n_draws: 5000`, `freq_band: [0.0, 3.14...]`) and coerced
##     back to that type on read;
##   * any other value (free-form fields, NA, named vectors, matrices, S3
##     lists, empty vectors) is written TAGGED:
##       {.r_type: <typeof>, value: [...], attributes: {name: <value>, ...}}
##     which restores the exact type and every attribute. Inside a tagged
##     atomic value NA is null; NaN, Inf and -Inf are strings, and a
##     subnormal double m * 2^-1074 is the string "subnormal:<m>" (the YAML
##     reader cannot read subnormals, and R's decimal parser is not exact).
##     Free-form scalars and named lists are written bare when that is
##     unambiguous (a double always carries a decimal point or exponent, so it
##     reads back as a double).
##
## Doubles: JSON via jsonlite `digits = I(17)` (17 significant digits); YAML
## via .spec_fmt_double() (the shortest of 15-17 significant digits that the
## YAML reader maps back exactly -- yaml's own `precision = 17` is lossy below
## 1). Both round-trip every finite double exactly.
## Closures, environments and other language objects cannot be written as
## text: the writer lists their fields and points to RDS. The model and the
## data are referenced by path + hash (an RDS sidecar next to the spec file
## when they exist only in memory); read_spec() checks the hash.
## ---------------------------------------------------------------------------

.spec_rtag <- ".r_type"

.spec_plain_double <- function(x)
  is.finite(x) && (x == 0 || abs(x) >= .Machine$double.xmin)

.spec_good_names <- function(nm)
  !is.null(nm) && all(nzchar(nm)) && !anyNA(nm) && !anyDuplicated(nm) &&
    !(.spec_rtag %in% nm)

## One atomic vector's elements as a list of document scalars.
.spec_enc_atomic <- function(x) {
  tp <- typeof(x)
  attributes(x) <- NULL   # element access must not dispatch (Date, factor, ...)
  out <- vector("list", length(x))
  for (i in seq_along(x)) {
    v <- x[[i]]
    out[i] <- list(
      if (tp == "double") {
        if (is.nan(v)) "NaN"
        else if (is.na(v)) NULL
        else if (is.infinite(v)) if (v > 0) "Inf" else "-Inf"
        else if (.spec_plain_double(v)) v
        else sprintf("subnormal:%.0f", v / 2^-1074)   # exact: m * 2^-1074
      } else if (is.na(v)) NULL else v)
  }
  out
}

## Encode any R value. `bad` (an environment) collects the paths of values
## that have no text form.
.spec_enc_any <- function(x, path, bad) {
  if (is.null(x)) return(NULL)
  tp <- typeof(x)
  if (isS4(x) || !tp %in% c("logical", "integer", "double", "character", "list")) {
    bad$paths <- c(bad$paths, path)
    return(NULL)
  }
  at <- attributes(x)
  if (is.null(at) && length(x) == 1L && tp != "list") {
    if (!is.na(x) && (tp != "double" || .spec_plain_double(x))) return(x)
  }
  if (tp == "list" && length(x) && identical(names(at), "names") &&
      .spec_good_names(names(x))) {
    out <- lapply(names(x), function(nm)
      .spec_enc_any(x[[nm]], paste0(path, "$", nm), bad))
    names(out) <- names(x)
    return(out)
  }
  val <- if (tp == "list") {
    lapply(seq_along(x), function(i) {
      el <- x[[i]]
      if (is.null(el)) setNames(list("NULL"), .spec_rtag)
      else .spec_enc_any(el, sprintf("%s[[%d]]", path, i), bad)
    })
  } else .spec_enc_atomic(x)
  out <- setNames(list(tp, val), c(.spec_rtag, "value"))
  if (length(at)) {
    enc_at <- lapply(names(at), function(nm)
      .spec_enc_any(at[[nm]], paste0(path, "@", nm), bad))
    names(enc_at) <- names(at)
    out$attributes <- enc_at
  }
  out
}

.spec_is_tagged <- function(d)
  is.list(d) && !is.null(names(d)) && .spec_rtag %in% names(d)

.spec_dec_elem <- function(e, tp) {
  na <- is.null(e) || (is.atomic(e) && length(e) == 1L && is.na(e) &&
                         !(is.double(e) && is.nan(e)))
  switch(tp,
    double    = if (na) NA_real_
                else if (is.character(e) && startsWith(e, "subnormal:"))
                  as.numeric(substring(e, 11L)) * 2^-1074
                else as.numeric(e),
    integer   = if (na) NA_integer_ else as.integer(e),
    logical   = if (na) NA else as.logical(e),
    character = if (na) NA_character_ else as.character(e))
}

## Decode a value written by .spec_enc_any().
.spec_dec_any <- function(d) {
  if (is.null(d)) return(NULL)
  if (.spec_is_tagged(d)) {
    tp <- d[[.spec_rtag]]
    if (identical(tp, "NULL")) return(NULL)
    val <- d$value
    els <- if (is.null(val)) list() else if (is.list(val)) val else as.list(val)
    x <- switch(tp,
      list = lapply(els, .spec_dec_any),
      double = , integer = , logical = , character =
        vapply(els, .spec_dec_elem, switch(tp, double = 0, integer = 0L,
                                           logical = NA, character = ""),
               tp = tp, USE.NAMES = FALSE),
      .dynhr_abort("read_spec: unknown value type \"", format(tp), "\".",
                   class = "dynhr_error_spec_parse"))
    if (length(d$attributes)) attributes(x) <- lapply(d$attributes, .spec_dec_any)
    return(x)
  }
  if (is.list(d)) {
    if (is.null(names(d)))
      .dynhr_abort("read_spec: an untagged sequence has no R type; write ",
                   "vectors as {.r_type: ..., value: [...]}.",
                   class = "dynhr_error_spec_parse")
    return(lapply(d, .spec_dec_any))
  }
  if (length(d) != 1L)
    .dynhr_abort("read_spec: an untagged sequence has no R type; write ",
                 "vectors as {.r_type: ..., value: [...]}.",
                 class = "dynhr_error_spec_parse")
  d
}

## Typed field: bare when the value is a plain instance of its type.
.spec_enc_field <- function(v, def, path, bad) {
  if (is.null(v)) return(NULL)
  if (!isTRUE(def$text)) {
    bad$paths <- c(bad$paths, path)
    return(NULL)
  }
  plain <- is.null(attributes(v)) && length(v) >= 1L && !anyNA(v)
  if (plain && def$type %in% c("lgl1", "int1", "num1", "chr1", "num")) {
    if (is.double(v) && !all(vapply(v, .spec_plain_double, logical(1))))
      return(.spec_enc_any(v, path, bad))
    if (def$type == "num" && length(v) > 1L) return(as.list(v))
    return(v)
  }
  .spec_enc_any(v, path, bad)
}

.spec_dec_field <- function(d, def) {
  if (is.null(d)) return(NULL)
  if (.spec_is_tagged(d) || def$type %in% c("any", "list")) return(.spec_dec_any(d))
  u <- unlist(d, use.names = FALSE)
  switch(def$type,
    lgl1 = as.logical(u), int1 = as.integer(u), num1 = , num = as.numeric(u),
    chr1 = as.character(u), u)
}

.spec_enc_component <- function(x, component, bad) {
  sch <- .spec_component_schema(component, x$method)
  out <- lapply(names(sch), function(nm)
    .spec_enc_field(x[[nm]], sch[[nm]], paste0(component, "$", nm), bad))
  names(out) <- names(sch)
  out
}

.spec_dec_component <- function(d, component) {
  if (!is.list(d) || is.null(names(d)))
    .dynhr_abort("read_spec: `", component, "` must be a mapping of fields.",
                 class = "dynhr_error_spec_parse")
  method <- if (identical(component, "sampler")) tolower(as.character(d$method %||% "rwmh"))
  if (identical(component, "sampler") && !method %in% .spec_samplers)
    .dynhr_abort("read_spec: unknown sampler method \"", method, "\".",
                 class = "dynhr_error_spec_parse")
  sch <- .spec_component_schema(component, method)
  vals <- lapply(names(d), function(nm)
    if (nm %in% names(sch)) .spec_dec_field(d[[nm]], sch[[nm]]) else d[[nm]])
  names(vals) <- names(d)
  .spec_build(component, vals, method = method)
}

# ---------------------------------------------------------------------------
# Spec <-> document
# ---------------------------------------------------------------------------

.spec_sidecar <- function(path, what)
  paste0(sub("\\.[^./\\\\]*$", "", basename(path)), ".", what, ".rds")

## The document for `spec`, written to `path` (sidecars go next to it).
.spec_to_doc <- function(spec, path) {
  bad <- new.env(parent = emptyenv())
  bad$paths <- character(0)
  dir <- dirname(path)
  mp <- spec$model
  dropped <- c(if (!is.null(mp$compiled)) "model$compiled",
               if (!is.null(mp$solved)) "model$solved")

  ## ---- model: by path when the file still parses to this model ----------
  by_path <- !is.null(mp$path) && file.exists(mp$path) &&
    identical(parse_mod(mp$path, verbose = FALSE), mp$mod)
  model_doc <- if (by_path) {
    list(source = "file", path = mp$path, hash = .spec_file_hash(mp$path))
  } else {
    sc <- .spec_sidecar(path, "model")
    saveRDS(mp$mod, file.path(dir, sc))
    list(source = "sidecar", path = sc, hash = .spec_file_hash(file.path(dir, sc)))
  }
  model_doc$max_order  <- mp$max_order
  model_doc$prior_spec <- .spec_enc_any(mp$prior_spec, "model$prior_spec", bad)

  ## ---- data --------------------------------------------------------------
  data_doc <- if (!is.null(spec$data$path)) {
    list(source = "file", path = spec$data$path,
         hash = .spec_file_hash(spec$data$path))
  } else if (!is.null(spec$data$value)) {
    sc <- .spec_sidecar(path, "data")
    saveRDS(spec$data$value, file.path(dir, sc))
    list(source = "sidecar", path = sc, hash = .spec_file_hash(file.path(dir, sc)))
  }

  sl <- .spec_sampler_list(spec$sampler)
  samp_doc <- if (!length(sl)) NULL
    else if (length(sl) == 1L) .spec_enc_component(sl[[1L]], "sampler", bad)
    else lapply(sl, .spec_enc_component, component = "sampler", bad = bad)

  doc <- list(
    spec_version  = spec$spec_version,
    dynhr_version = unname(as.character(utils::packageVersion("dynhr"))),
    spec_hash     = spec$hashes$spec,
    model         = model_doc,
    data          = data_doc,
    obs_vars      = as.list(spec$obs_vars),
    likelihood    = .spec_enc_component(spec$likelihood, "likelihood", bad),
    mode          = .spec_enc_component(spec$mode, "mode", bad),
    sampler       = samp_doc,
    compute       = .spec_enc_component(spec$compute, "compute", bad),
    outputs       = .spec_enc_component(spec$outputs, "outputs", bad),
    options       = .spec_enc_any(spec$options, "options", bad))
  if (length(bad$paths))
    .dynhr_abort("write_spec: these fields hold values with no text form ",
                 "(functions, environments, a precomputed mode result, ...): ",
                 paste(bad$paths, collapse = ", "),
                 ". Write the spec with format = \"rds\" instead.",
                 class = "dynhr_error_spec_not_serialisable")
  if (length(dropped))
    .dynhr_inform("write_spec: ", paste(dropped, collapse = " and "),
                  " (cached, rebuilt from the model) not written to the text ",
                  "format.", class = "dynhr_message_spec_cache_dropped")
  doc
}

## Resolve a referenced file (absolute, or relative to the spec file) and
## check its hash.
.spec_ref_file <- function(ref, dir, what, force) {
  p <- as.character(ref$path %||% "")
  if (!nzchar(p))
    .dynhr_abort("read_spec: the ", what, " reference has no path.",
                 class = "dynhr_error_spec_parse")
  if (identical(ref$source, "sidecar") || !file.exists(p)) {
    alt <- file.path(dir, p)
    if (file.exists(alt)) p <- alt
  }
  if (!file.exists(p))
    .dynhr_abort("read_spec: the ", what, " file ", p, " does not exist.",
                 class = "dynhr_error_spec_missing_file")
  h <- .spec_file_hash(p)
  want <- as.character(ref$hash %||% NA_character_)
  same_algo <- identical(sub(":.*$", "", h), sub(":.*$", "", want))
  if (same_algo && !identical(h, want)) {
    msg <- paste0("the ", what, " file ", p, " has changed since the spec was ",
                  "written (hash ", want, ", now ", h, ")")
    if (!isTRUE(force))
      .dynhr_abort("read_spec: ", msg, ". Pass force = TRUE to read it anyway.",
                   class = "dynhr_error_spec_hash_mismatch")
    .dynhr_warn("read_spec: ", msg, "; read anyway (force = TRUE).",
                class = "dynhr_warning_spec_hash_mismatch")
  }
  p
}

.spec_from_doc <- function(doc, dir, force) {
  if (!is.list(doc) || is.null(names(doc)))
    .dynhr_abort("read_spec: not a spec document.", class = "dynhr_error_spec_parse")
  sv <- doc$spec_version
  if (!is.numeric(sv) || length(sv) != 1L || !identical(as.integer(sv), .spec_version))
    .dynhr_abort("read_spec: spec_version ", format(sv), " is not supported ",
                 "(this dynhr reads version ", .spec_version, ").",
                 class = "dynhr_error_spec_version")
  md <- doc$model
  mp_path <- .spec_ref_file(md, dir, "model", force)
  mod <- if (identical(md$source, "sidecar")) readRDS(mp_path)
         else parse_mod(mp_path, verbose = FALSE)
  model <- .spec_model_part(mod, max_order = as.integer(md$max_order %||% 1L),
                            prior_spec = .spec_dec_any(md$prior_spec))
  dd <- doc$data
  data <- if (is.null(dd)) list(value = NULL, path = NULL) else {
    dp <- .spec_ref_file(dd, dir, "data", force)
    if (identical(dd$source, "sidecar")) .spec_data_part(readRDS(dp))
    else .spec_data_part(dp)
  }
  sd <- doc$sampler
  sampler <- if (is.null(sd)) NULL
    else if (is.null(names(sd))) {
      sl <- lapply(sd, .spec_dec_component, component = "sampler")
      if (length(sl) == 1L) sl[[1L]] else structure(sl, class = "dynhr_sampler_sequence")
    } else .spec_dec_component(sd, "sampler")
  parts <- list(
    model      = model,
    data       = data,
    obs_vars   = .spec_resolve_obs_vars(as.character(unlist(doc$obs_vars)), model$mod),
    likelihood = .spec_dec_component(doc$likelihood %||% setNames(list(), character(0)), "likelihood"),
    mode       = .spec_dec_component(doc$mode %||% setNames(list(), character(0)), "mode"),
    sampler    = sampler,
    compute    = .spec_dec_component(doc$compute %||% setNames(list(), character(0)), "compute"),
    outputs    = .spec_dec_component(doc$outputs %||% setNames(list(), character(0)), "outputs"),
    options    = .spec_dec_any(doc$options) %||% list())
  spec <- .spec_assemble(parts)
  old <- as.character(doc$spec_hash %||% NA_character_)
  if (identical(sub(":.*$", "", old), sub(":.*$", "", spec$hashes$spec)) &&
      !identical(old, spec$hashes$spec))
    .dynhr_inform("read_spec: the spec differs from the one written (its hash ",
                  "changed): it was edited, or a default came from the ",
                  "current options.", class = "dynhr_message_spec_edited")
  spec
}

# ---------------------------------------------------------------------------
# YAML text with comments
# ---------------------------------------------------------------------------

## The shortest decimal form (15 to 17 significant digits) of a finite normal
## double that the YAML reader maps back to exactly that double, always with a
## decimal point so YAML resolves it as a float. Checked with the YAML reader
## itself: R's own as.numeric() is not correctly rounded (on arm64 it misreads
## most 17-digit strings with large exponents).
.spec_fmt_double <- function(v) {
  for (p in 15:17) {
    s <- sprintf("%.*g", p, v)
    if (!grepl("[.]", s)) s <- if (grepl("e", s)) sub("e", ".0e", s) else paste0(s, ".0")
    ## a shorter form can round past the largest double: the reader then
    ## warns and returns NA, which simply fails the check
    if (identical(suppressWarnings(yaml::yaml.load(s)), v)) return(s)
  }
  sprintf("%.16e", v)
}

## Doubles in a document written verbatim with .spec_fmt_double(): yaml's own
## `precision = 17` counts DECIMALS for |x| < 1 (0.00012864095435953938 is
## written 0.00012864095435954, which reads back as a different double).
.spec_yaml_prep <- function(d) {
  if (is.list(d)) {
    for (i in seq_along(d)) if (!is.null(d[[i]])) d[i] <- list(.spec_yaml_prep(d[[i]]))
    return(d)
  }
  if (is.double(d)) {
    out <- lapply(d, function(v) structure(.spec_fmt_double(v), class = "verbatim"))
    return(if (length(d) == 1L) out[[1L]] else out)
  }
  d
}

.spec_yaml_lines <- function(x) {
  txt <- yaml::as.yaml(.spec_yaml_prep(x), precision = 17L)
  strsplit(sub("\n$", "", txt), "\n", fixed = TRUE)[[1L]]
}

.spec_yaml_component <- function(d, component, method = NULL, indent = "  ") {
  sch <- .spec_component_schema(component, method)
  lines <- character(0)
  for (nm in names(d)) {
    fl <- .spec_yaml_lines(setNames(list(d[[nm]]), nm))
    lines <- c(lines, paste0(indent, "# ", sch[[nm]]$doc), paste0(indent, fl))
  }
  lines
}

.spec_yaml_text <- function(doc) {
  top <- function(nm) .spec_yaml_lines(doc[nm])
  out <- c(
    sprintf("# dynhr estimation spec, spec_version %d, written by dynhr %s",
            doc$spec_version, doc$dynhr_version),
    paste0("# spec hash: ", doc$spec_hash),
    "# Read it back with read_spec(). Numbers carry 17 significant digits so",
    "# they read back exactly. A value written as {.r_type: ..., value: ...}",
    "# is an R value stored with its exact type and attributes.",
    top("spec_version"), top("dynhr_version"), top("spec_hash"),
    "# Model: a .mod file (or an RDS sidecar) checked against its hash",
    top("model"),
    "# Data: a CSV file (or an RDS sidecar) checked against its hash",
    if (is.null(doc$data)) "data: ~" else top("data"),
    "# Observable names, in data-column order",
    top("obs_vars"))
  for (cp in .spec_components) {
    out <- c(out, paste0("# ", .spec_component_doc[[cp]]))
    d <- doc[[cp]]
    if (identical(cp, "sampler")) {
      if (is.null(d)) {
        out <- c(out, "sampler: ~")
      } else if (!is.null(names(d))) {
        out <- c(out, "sampler:", .spec_yaml_component(d, "sampler", d$method))
      } else {
        out <- c(out, "sampler:")
        for (item in d) {
          sch <- .spec_component_schema("sampler", item$method)
          first <- TRUE
          for (nm in names(item)) {
            fl <- .spec_yaml_lines(setNames(list(item[[nm]]), nm))
            cm <- paste0("# ", sch[[nm]]$doc)
            if (first) {
              out <- c(out, paste0("  ", cm), paste0("  - ", fl[1L]),
                       if (length(fl) > 1L) paste0("    ", fl[-1L]))
              first <- FALSE
            } else {
              out <- c(out, paste0("    ", cm), paste0("    ", fl))
            }
          }
        }
      }
    } else {
      out <- c(out, paste0(cp, ":"), .spec_yaml_component(d, cp))
    }
  }
  c(out, "# Package options that do not change results (snapshot)",
    if (is.null(doc$options)) "options: ~" else top("options"))
}

.spec_need <- function(pkg, fmt) {
  if (!requireNamespace(pkg, quietly = TRUE))
    .dynhr_abort("write_spec / read_spec: the ", fmt, " format needs the ",
                 "suggested package '", pkg, "'; install it or use ",
                 "format = \"rds\".", class = "dynhr_error_missing_package")
}

.spec_format_of <- function(path) {
  ext <- tolower(sub("^.*\\.", "", basename(path)))
  if (ext %in% c("yaml", "yml")) "yaml" else if (ext == "json") "json"
  else if (ext == "rds") "rds" else "yaml"
}

#' Write and read estimation specs (YAML, JSON, RDS)
#'
#' \code{write_spec()} saves a \code{\link{dynhr_estimation_spec}};
#' \code{read_spec()} reads it back. RDS is the canonical format: the spec by
#' value, model and data included. YAML (via the suggested package
#' \pkg{yaml}) and JSON (via \pkg{jsonlite}) are human-readable and editable:
#' YAML carries a header (spec version, dynhr version, spec hash) and a comment
#' line for every field.
#'
#' Reading a text spec rebuilds every field with the type the spec schema
#' gives it and validates the result, so a spec written and read back is
#' \code{identical()} to the original: doubles are written with 17
#' significant digits, and values that plain YAML/JSON cannot represent
#' exactly (integer versus double, \code{NA}, \code{NaN}, \code{Inf}, named
#' vectors, matrices, empty vectors, S3 lists such as data.frames) are written
#' in a tagged form \code{{.r_type: ..., value: [...], attributes: {...}}}.
#' A field left out of an edited file takes its default (options-sourced
#' defaults from the current options).
#'
#' In the text formats the model and the data are referenced by path and
#' content hash; a model or data set that exists only in memory (or a model
#' that no longer matches its file) is saved as an RDS sidecar next to the
#' spec file (\code{<name>.model.rds}, \code{<name>.data.rds}).
#' \code{read_spec()} refuses a referenced file whose hash changed (error class
#' \code{dynhr_error_spec_hash_mismatch}) unless \code{force = TRUE}. Values
#' with no text form -- functions (e.g. system-prior densities), environments,
#' a precomputed mode result -- make \code{write_spec()} fail with class
#' \code{dynhr_error_spec_not_serialisable}, listing the fields; use RDS for
#' such specs. Cached compiled models are not written (they are rebuilt).
#'
#' @param spec A \code{dynhr_estimation_spec}.
#' @param path File path.
#' @param format \code{"yaml"}, \code{"json"} or \code{"rds"}; by default from
#'   the file extension (\code{.yaml}/\code{.yml}, \code{.json}, \code{.rds};
#'   anything else is YAML).
#' @param force \code{read_spec()}: read even when a referenced model or data
#'   file changed since the spec was written (with a warning).
#' @return \code{write_spec()}: \code{path}, invisibly. \code{read_spec()}: the
#'   \code{dynhr_estimation_spec}.
#' @seealso \code{\link{dynhr_estimation_spec}}
#' @examples
#' \donttest{
#' if (requireNamespace("yaml", quietly = TRUE)) {
#'   mod <- system.file("extdata/models/nk_demo.mod", package = "dynhr")
#'   dat <- system.file("extdata/models/nk_demo_data.csv", package = "dynhr")
#'   spec <- dynhr_estimation_spec(mod, data = dat,
#'                                 obs_vars = c("ygr", "infl", "intr"))
#'   f <- tempfile(fileext = ".yaml")
#'   write_spec(spec, f)
#'   identical(read_spec(f), spec)
#' }
#' }
#' @export
write_spec <- function(spec, path, format = c("yaml", "json", "rds")) {
  if (!inherits(spec, "dynhr_estimation_spec"))
    .dynhr_abort("write_spec: `spec` must be a dynhr_estimation_spec.",
                 class = "dynhr_error_bad_argument")
  format <- if (missing(format)) .spec_format_of(path) else match.arg(format)
  if (identical(format, "rds")) {
    saveRDS(spec, path)
    return(invisible(path))
  }
  .spec_need(if (identical(format, "yaml")) "yaml" else "jsonlite", format)
  doc <- .spec_to_doc(spec, path)
  txt <- if (identical(format, "yaml")) {
    .spec_yaml_text(doc)
  } else {
    as.character(jsonlite::toJSON(doc, auto_unbox = TRUE, digits = I(17),
                                  always_decimal = TRUE, null = "null",
                                  pretty = TRUE))
  }
  writeLines(txt, path, useBytes = TRUE)
  invisible(path)
}

#' @rdname write_spec
#' @export
read_spec <- function(path, format = NULL, force = FALSE) {
  if (!is.character(path) || length(path) != 1L || !file.exists(path))
    .dynhr_abort("read_spec: file not found: ", format(path),
                 class = "dynhr_error_spec_missing_file")
  format <- format %||% .spec_format_of(path)
  format <- match.arg(format, c("yaml", "json", "rds"))
  if (identical(format, "rds")) {
    spec <- readRDS(path)
    if (!inherits(spec, "dynhr_estimation_spec"))
      .dynhr_abort("read_spec: ", path, " does not hold a dynhr_estimation_spec.",
                   class = "dynhr_error_spec_parse")
    if (!identical(spec$spec_version, .spec_version))
      .dynhr_abort("read_spec: spec_version ", format(spec$spec_version),
                   " is not supported (this dynhr reads version ",
                   .spec_version, ").", class = "dynhr_error_spec_version")
    return(spec)
  }
  .spec_need(if (identical(format, "yaml")) "yaml" else "jsonlite", format)
  doc <- if (identical(format, "yaml")) {
    yaml::read_yaml(path, eval.expr = FALSE)
  } else {
    jsonlite::fromJSON(path, simplifyVector = FALSE)
  }
  .spec_from_doc(doc, dirname(path), force)
}
