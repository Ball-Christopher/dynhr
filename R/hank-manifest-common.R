## R/hank-manifest-common.R
## --------------------------------------------------------------------------
## Family-agnostic machinery shared by the one-, two- and three-asset run
## manifests and fingerprints.
##
## The three-asset provenance tooling (R/hank-manifest3.R,
## hank_het3_fingerprint()) came first and its helpers were NAMED with a "3"
## even though most of them have nothing three-asset about them: resolving the
## git commit, reading peak RSS, hashing a fixed field list, flattening a
## manifest to one row, printing it. Those live here now, once, and the three
## per-family wrappers supply only what genuinely differs:
##
##   * the accepted block / bare-solve classes,
##   * the asset-axis names and their grid sizes,
##   * the grid-mass axis set (a | b,a | d,f,a),
##   * the aggregate field names (A,C | B,A,C,CHI | D_agg,F_agg,A_agg,C,CHI,PHI),
##   * the literal, order-fixed fingerprint field list.
##
## CONTRACT: the three-asset outputs are on an external paper's reproducibility
## path. Every field, the as.data.frame() column set AND order, the printed
## text and the fingerprint hash string must remain bit-identical to the
## pre-refactor implementation. Anything added here must be additive.
## --------------------------------------------------------------------------


#' Resolve the git commit a manifest should record, best effort
#'
#' Never \code{stop()}s -- a manifest call must not be able to fail a long
#' production run, so every step is wrapped and a failure at any stage just
#' falls through to the next source.
#'
#' Resolution order, first hit wins:
#' \enumerate{
#'   \item a LIVE \code{git rev-parse HEAD}, but only when a \code{.git} entry
#'     is findable from the package root -- a directory in a normal checkout,
#'     or a FILE in a git worktree (this repo's own dev setup), so testing for
#'     a directory specifically would be wrong here. This is the
#'     \code{devtools::load_all()} case and must never be stale relative to
#'     HEAD.
#'   \item \code{system.file("GIT_COMMIT", package = "dynhr")}, written at
#'     release time (never by this function): the commit the package was BUILT
#'     from, not necessarily HEAD of whatever tree it now sits in.
#'   \item \code{packageDescription()$RemoteSha} / \code{$GitCommit}, if a
#'     remote-install workflow (e.g. \code{remotes::install_github()})
#'     recorded one.
#'   \item \code{NA_character_}.
#' }
#'
#' @return A length-1 character: a 40-hex-char SHA, a remote-recorded string,
#'   or \code{NA_character_}.
#' @keywords internal
.hank_manifest_git_commit <- function() {
  is_sha1 <- function(z) grepl("^[0-9a-f]{40}$", z)

  root <- tryCatch(system.file(package = "dynhr"), error = function(e) "")
  ## pkgload::load_all() (what devtools::load_all() uses) shims system.file()
  ## to resolve into <source>/inst, not the source root itself -- a real
  ## R CMD INSTALL tree has no separate "inst" (its contents are merged up
  ## a level at install time), so this correction is safe in both cases: it
  ## only fires when the path really does end in "inst".
  if (is.character(root) && length(root) == 1L && nzchar(root) &&
      basename(root) == "inst" && dir.exists(dirname(root)))
    root <- dirname(root)
  if (is.character(root) && length(root) == 1L && nzchar(root) &&
      file.exists(file.path(root, ".git"))) {
    commit <- tryCatch({
      ## system2() shells out (system() under the hood), so an unquoted path
      ## containing spaces -- this very repo's own path does -- silently
      ## truncates at the first space. shQuote() is required, not optional.
      out <- suppressWarnings(system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"),
                                       stdout = TRUE, stderr = FALSE))
      st <- attr(out, "status")
      if (!is.null(st) && !identical(st, 0L)) character(0) else out
    }, error = function(e) character(0))
    if (length(commit) == 1L && is_sha1(commit)) return(commit)
  }

  gc_path <- tryCatch(system.file("GIT_COMMIT", package = "dynhr"),
                       error = function(e) "")
  if (is.character(gc_path) && length(gc_path) == 1L && nzchar(gc_path) &&
      file.exists(gc_path)) {
    txt <- tryCatch(suppressWarnings(readLines(gc_path, n = 1L, warn = FALSE)),
                     error = function(e) character(0))
    txt <- trimws(txt)
    if (length(txt) == 1L && is_sha1(txt)) return(txt)
  }

  desc <- tryCatch(utils::packageDescription("dynhr"), error = function(e) NULL)
  if (is.list(desc)) {
    for (fld in c("RemoteSha", "GitCommit")) {
      v <- desc[[fld]]
      if (is.character(v) && length(v) == 1L && !is.na(v) && nzchar(v))
        return(v)
    }
  }

  NA_character_
}


#' Boundary and interior masses of a stationary distribution, any asset count
#'
#' The family-agnostic generalisation of the three-asset
#' \code{.hank3_grid_mass}: for each ASSET axis it reports the mass sitting on
#' that axis's lower and upper gridpoint, plus the mass interior to every
#' asset axis at once. The income axis is never a boundary axis -- it is an
#' exogenous Markov state with no "constraint" to sit at -- so it is summed
#' over throughout.
#'
#' ONE implementation, so the cell-order convention is written down (and can
#' go stale) in exactly one place.
#'
#' @param D Numeric \code{prod(dims)} stationary distribution in package cell
#'   order (income slowest, the last asset axis fastest).
#' @param dims Integer vector of array extents in FASTEST-FIRST order, i.e.
#'   the shape that \code{array(D, dims)} reproduces the vector with:
#'   \code{c(n_a, n_e)} one-asset, \code{c(n_a, n_b, n_e)} two-asset,
#'   \code{c(n_a, n_f, n_d, n_e)} three-asset.
#' @param axis_names Character, same length as \code{dims}: the reported name
#'   of each asset axis, and \code{NA_character_} for the income axis (which
#'   contributes no boundary masses). Output names run in REVERSE \code{dims}
#'   order -- slowest asset axis first -- so the three-asset result reads
#'   \code{d, f, a}, the two-asset one \code{b, a}.
#'
#' @return A named numeric vector: \code{<axis>_lo}, \code{<axis>_hi} for each
#'   asset axis (slowest first), then \code{interior}. \code{interior} is
#'   \code{NA_real_} when any asset axis has fewer than 3 points (no distinct
#'   interior to sum).
#' @keywords internal
.hank_grid_mass_axes <- function(D, dims, axis_names) {
  dims <- as.integer(dims)
  Darr <- array(D, dims)
  k <- length(dims)
  asset <- which(!is.na(axis_names))

  ## Slowest-varying asset axis first, so the reported order is the economic
  ## one (d, f, a) rather than the storage one (a, f, d).
  ord <- rev(asset)

  full <- lapply(dims, seq_len)
  interior_mass <- if (all(dims[asset] >= 3L)) {
    idx <- full
    for (j in asset) idx[[j]] <- 2:(dims[j] - 1L)
    sum(do.call(`[`, c(list(Darr), idx, list(drop = FALSE))))
  } else NA_real_

  out <- numeric(0)
  for (j in ord) {
    for (side in c("lo", "hi")) {
      idx <- full
      idx[[j]] <- if (side == "lo") 1L else dims[j]
      out[paste0(axis_names[j], "_", side)] <-
        sum(do.call(`[`, c(list(Darr), idx, list(drop = FALSE))))
    }
  }
  c(out, interior = interior_mass)
}


#' One atomic, length-1 label for a block's \code{Pi_inputs}
#'
#' Keeps the manifest a single data-frame row. Values are formatted to full
#' double precision -- this is provenance, and a rounded transition rate would
#' make two genuinely different households look identical.
#'
#' @param Pi_inputs The block's \code{Pi_inputs} list, or \code{NULL}.
#' @return A length-1 character, or \code{NA_character_} for a price-only
#'   block.
#' @keywords internal
.hank_transition_inputs_label <- function(Pi_inputs) {
  if (is.null(Pi_inputs) || length(Pi_inputs) == 0L) return(NA_character_)
  nm <- names(Pi_inputs)
  paste(paste0(nm, "=", vapply(Pi_inputs, format, character(1), digits = 17)),
        collapse = ", ")
}


#' Assemble a run manifest from an already-classified household object
#'
#' The shared body of \code{\link{hank_het_manifest}},
#' \code{\link{hank_het2_manifest}} and \code{\link{hank_het3_manifest}}: it
#' reads the fields every family carries under the same names (backend,
#' threads, iterations, converged, gaps, timings, \code{Pi_inputs}) and splices
#' in the per-family pieces its caller has already computed.
#'
#' Field ORDER is fixed here and is part of the contract -- it is the column
#' order of \code{as.data.frame()} on the result.
#'
#' @param x The household object (a block, or a bare solve), already
#'   validated and unwrapped by the caller.
#' @param is_block Logical: whether \code{x} carries a stationary distribution
#'   and block-level aggregates. Bare solves report the block-only fields as
#'   \code{NA} rather than erroring.
#' @param peak_rss Logical: when \code{FALSE}, skip the
#'   \code{\link{.hank_peak_rss}} call and report \code{NA_real_}.
#' @param sizes Named list of the family's grid sizes, in reporting order
#'   (e.g. \code{list(n_e =, n_d =, n_f =, n_a =)}).
#' @param state_count Integer total cell count.
#' @param grid_mass The family's grid-mass vector, or \code{NULL} when no
#'   distribution is available.
#' @param agg_fields Character vector of the family's aggregate field names,
#'   read off \code{x} when \code{is_block} and \code{NA_real_} otherwise.
#' @param class Character: the S3 class to stamp on the result.
#'
#' @return A named list of class \code{class}.
#' @keywords internal
.hank_manifest_build <- function(x, is_block, peak_rss, sizes, state_count,
                                 grid_mass, agg_fields, class) {
  elapsed_solve <- if (!is.null(x$elapsed_solve)) x$elapsed_solve else
    if (!is.null(x$elapsed)) x$elapsed else NA_real_
  elapsed_dist <- if (!is.null(x$elapsed_dist)) x$elapsed_dist else NA_real_

  ## Exact indexing (not `$`): a block that predates one of these fields must
  ## report NA, not partial-match a longer sibling name.
  aggs <- lapply(agg_fields, function(f) {
    v <- x[[f, exact = TRUE]]
    if (is_block && !is.null(v)) v else NA_real_
  })
  names(aggs) <- agg_fields

  structure(
    c(list(package_version = as.character(utils::packageVersion("dynhr")),
           git_commit = .hank_manifest_git_commit(),
           backend = if (!is.null(x$backend)) x$backend else NA_character_,
           threads = if (!is.null(x$threads)) x$threads else NA_integer_),
      sizes,
      list(state_count = state_count,
           iterations = if (!is.null(x$iterations)) x$iterations else NA_integer_,
           converged = if (!is.null(x$converged)) x$converged else NA,
           last_value_gap = if (!is.null(x$last_value_gap)) x$last_value_gap else NA_real_,
           last_policy_gap = if (!is.null(x$last_policy_gap)) x$last_policy_gap else NA_real_,
           dist_converged = if (is_block && !is.null(x$dist_converged)) x$dist_converged else NA,
           elapsed_solve = elapsed_solve,
           elapsed_dist = elapsed_dist,
           peak_rss_bytes = if (isTRUE(peak_rss)) .hank_peak_rss() else NA_real_,
           grid_mass = grid_mass,
           ## Which transition-probability inputs this block carries, and at
           ## what steady-state values. Provenance, not diagnostics: two runs
           ## whose manifests agree on everything else but differ here were
           ## solved against DIFFERENT households, and a Jacobian cache keyed
           ## on the rest of this row would silently conflate them. NA when
           ## the block is price-only.
           transition_inputs = .hank_transition_inputs_label(x$Pi_inputs)),
      aggs),
    class = class)
}


#' Validate the \code{peak_rss} flag of a manifest call
#'
#' @param peak_rss The value supplied by the user.
#' @param caller Name of the calling manifest function, for the message.
#' @return \code{invisible(NULL)}; called for its error.
#' @keywords internal
.hank_manifest_check_peak <- function(peak_rss, caller) {
  if (!(is.logical(peak_rss) && length(peak_rss) == 1L && !is.na(peak_rss)))
    stop(caller, ": 'peak_rss' must be a single TRUE/FALSE.")
  invisible(NULL)
}


#' Flatten any HANK run manifest to exactly one data-frame row
#'
#' Shared body of the three \code{as.data.frame} methods. The only non-scalar
#' field, \code{grid_mass}, is unnested into \code{grid_mass_*} columns
#' (filled with \code{NA_real_} when the source manifest carries no
#' distribution).
#'
#' @param x A HANK run manifest.
#' @param gm_names Character vector of the family's grid-mass axis names, in
#'   the order the columns must appear.
#' @return A one-row \code{data.frame}.
#' @keywords internal
.hank_manifest_as_df <- function(x, gm_names) {
  gm <- x$grid_mass
  gm_vals <- if (is.null(gm)) stats::setNames(rep(NA_real_, length(gm_names)), gm_names)
             else gm[gm_names]
  names(gm_vals) <- paste0("grid_mass_", gm_names)

  scalar <- unclass(x)[setdiff(names(x), "grid_mass")]
  out <- c(scalar, as.list(gm_vals))
  as.data.frame(out, stringsAsFactors = FALSE)
}


#' Print any HANK run manifest
#'
#' Shared body of the three \code{print} methods.
#'
#' @param x A HANK run manifest.
#' @param cls Character: the class name to show in the header.
#' @param size_detail Character: the parenthesised per-axis size list for the
#'   backend/threads/states line (e.g. \code{"e=2 d=3 f=3 a=3"}).
#' @param agg_fields Character vector of the family's aggregate field names.
#' @return \code{x}, invisibly.
#' @keywords internal
.hank_manifest_print <- function(x, cls, size_detail, agg_fields) {
  cat(sprintf("<%s: dynhr %s @ %s>\n", cls, x$package_version,
              if (is.na(x$git_commit)) "<unknown commit>" else
                substr(x$git_commit, 1L, 10L)))
  cat(sprintf("  backend = %s, threads = %s, states = %d (%s)\n",
              x$backend, x$threads, x$state_count, size_detail))
  cat(sprintf("  iterations = %s, converged = %s, dist_converged = %s\n",
              x$iterations, x$converged, x$dist_converged))
  cat(sprintf("  transition_inputs = %s\n",
              if (is.na(x$transition_inputs)) "<none: price-only block>"
              else x$transition_inputs))
  cat(sprintf("  last_value_gap = %s, last_policy_gap = %s\n",
              format(x$last_value_gap, digits = 4),
              format(x$last_policy_gap, digits = 4)))
  cat(sprintf("  elapsed_solve = %ss, elapsed_dist = %ss, peak_rss = %s bytes\n",
              format(x$elapsed_solve, digits = 4),
              format(x$elapsed_dist, digits = 4),
              format(x$peak_rss_bytes, digits = 6)))
  if (!is.null(x$grid_mass)) {
    cat("  grid_mass:\n")
    print(x$grid_mass)
  } else {
    cat("  grid_mass: <no distribution available>\n")
  }
  cat(sprintf("  %s\n",
              paste(vapply(agg_fields,
                           function(f) sprintf("%s = %s", f,
                                               format(x[[f]], digits = 4)),
                           character(1)),
                    collapse = ", ")))
  invisible(x)
}


#' Hash an explicitly-ordered field list into a content fingerprint
#'
#' Shared tail of \code{\link{hank_het_fingerprint}},
#' \code{\link{hank_het2_fingerprint}} and
#' \code{\link{hank_het3_fingerprint}}. Takes the field list ALREADY built by
#' the caller in a literal, fixed order -- never \code{names(block)}, so a
#' field appended to a block in future cannot silently change every existing
#' fingerprint.
#'
#' @param parts A list of the fields to hash, in the family's fixed order.
#' @param n Number of hex characters to return (1-16).
#' @param caller Name of the calling fingerprint function, for the error
#'   message.
#' @return A length-1 character string of \code{n} lowercase hex digits.
#' @keywords internal
.hank_fingerprint_hash <- function(parts, n, caller) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 || n > 16)
    stop(caller, ": 'n' must be a single integer in 1:16.")
  raw <- serialize(parts, connection = NULL, xdr = TRUE, version = 3L)
  substr(hank_fnv1a64_cpp(raw), 1L, as.integer(n))
}
