## R/hank-manifest.R
## --------------------------------------------------------------------------
## ONE-asset provenance tooling: the run manifest and the content fingerprint,
## the direct counterparts of hank_het3_manifest() / hank_het3_fingerprint().
##
## Everything family-agnostic lives in R/hank-manifest-common.R; this file
## carries only what is genuinely one-asset: which objects are accepted, the
## single asset axis, the (A, C) aggregates, and the literal fingerprint field
## list.
## --------------------------------------------------------------------------


#' Is this a bare one-asset household solve?
#'
#' \code{\link{hank_egm_solve}} returns a PLAIN LIST (unlike
#' \code{\link{hank_egm3_solve}}, which carries an S3 class), so the manifest
#' has to duck-type it. The signature is the one-asset policy set on a single
#' asset grid, with no liquid/illiquid or three-asset grid in sight.
#'
#' @param x The object to classify.
#' @return \code{TRUE} for a bare \code{\link{hank_egm_solve}} return value,
#'   \code{FALSE} otherwise.
#' @keywords internal
.hank_is_egm_solve <- function(x) {
  is.list(x) &&
    !is.null(x[["a_grid", exact = TRUE]]) &&
    !is.null(x[["a", exact = TRUE]]) &&
    !is.null(x[["c", exact = TRUE]]) &&
    !is.null(x[["Va", exact = TRUE]]) &&
    is.null(x[["b_grid", exact = TRUE]]) &&
    is.null(x[["d_grid", exact = TRUE]])
}


#' Run-manifest metadata for a one-asset solve or block
#'
#' The one-asset counterpart of \code{\link{hank_het3_manifest}}: collects, as
#' DATA rather than printed text, everything a run manifest needs to record
#' about one household solve -- backend, worker count, state count,
#' iterations, wall time, peak memory, convergence flags, boundary/interior
#' grid masses, package version and git commit. It is a thin READER over what
#' \code{\link{hank_het_block}} and \code{\link{hank_egm_solve}} already carry,
#' plus the handful neither carries on its own.
#'
#' \code{x} may be either a \code{\link{hank_het_block}} (a household solve
#' plus its stationary distribution) or a bare \code{\link{hank_egm_solve}}
#' return value (no distribution). Fields that only make sense with a
#' distribution (\code{dist_converged}, the aggregates, and \code{grid_mass}
#' when \code{x} carries no \code{D}) are \code{NA} (or \code{NULL} for
#' \code{grid_mass}) on a bare solve, never an error.
#'
#' \strong{Fields the one-asset family does not carry.} Both
#' \code{\link{hank_het_block}} and \code{\link{hank_egm_solve}} now record
#' \code{backend}, \code{threads}, \code{iterations}, \code{converged} and
#' their timings, so those are REAL values here, not \code{NA}. Two
#' exceptions, both deliberate:
#' \itemize{
#'   \item \code{threads} is always \code{1} -- the one-asset kernel is
#'     serial by measurement, not by omission (see
#'     \code{\link{hank_egm_solve}}), so \code{1} is the true worker count
#'     rather than a placeholder.
#'   \item \code{last_value_gap} / \code{last_policy_gap} are \code{NA_real_}
#'     whenever the solve ran the compiled backend or the borrowing-wedge
#'     fixed point, neither of which returns them; only the R backend of
#'     \code{\link{hank_egm_solve}} measures both.
#' }
#' \code{elapsed_dist} is \code{NA_real_} on a bare solve, which computes no
#' distribution.
#'
#' \strong{git_commit} resolution never fails the call, and
#' \strong{peak_rss_bytes} is a peak-since-process-start figure: see
#' \code{\link{hank_het3_manifest}} for both.
#'
#' @param x A \code{\link{hank_het_block}}, a bare \code{\link{hank_egm_solve}}
#'   return value, or a \code{hank_het_block_failed} wrapper (unwrapped to the
#'   household it carries, so the manifest reports the iterations and gaps that
#'   explain the failure while leaving the block-only aggregate fields
#'   \code{NA}).
#' @param peak_rss Logical (default \code{TRUE}). When \code{FALSE}, skip the
#'   \code{\link{.hank_peak_rss}} call and report \code{NA_real_} for
#'   \code{peak_rss_bytes} -- for callers who solve many small blocks in a loop.
#'
#' @return An object of class \code{hank_het_manifest}: a named list with
#'   \code{package_version}, \code{git_commit}, \code{backend},
#'   \code{threads}, \code{n_e}, \code{n_a}, \code{state_count},
#'   \code{iterations}, \code{converged}, \code{last_value_gap},
#'   \code{last_policy_gap}, \code{dist_converged}, \code{elapsed_solve},
#'   \code{elapsed_dist}, \code{peak_rss_bytes}, \code{grid_mass} (the named
#'   \code{a_lo}, \code{a_hi}, \code{interior} vector, or \code{NULL} with no
#'   distribution), \code{transition_inputs}, and the aggregates \code{A} and
#'   \code{C}.
#' @seealso \code{\link{hank_het_fingerprint}} (identity, not provenance),
#'   \code{\link{hank_het2_manifest}}, \code{\link{hank_het3_manifest}},
#'   \code{\link{as.data.frame.hank_het_manifest}},
#'   \code{\link{print.hank_het_manifest}}
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' ag  <- hank_asset_grid(20, 10, 0)
#' blk <- hank_het_block(ag, inc$Pi, inc$e, beta = 0.96, eis = 1,
#'                       r = 0.01, w = 1)
#' mf <- hank_het_manifest(blk)
#' mf$state_count
#' as.data.frame(mf)
#' @export
hank_het_manifest <- function(x, peak_rss = TRUE) {
  if (inherits(x, "hank_het_block_failed")) x <- x$hh
  if (inherits(x, "hank_het2_block") || inherits(x, "hank_het2d_block") ||
      inherits(x, "hank_het3_block") ||
      !(inherits(x, "hank_het_block") || .hank_is_egm_solve(x)))
    stop("hank_het_manifest: 'x' must be a hank_het_block, a bare ",
         "hank_egm_solve return value, or a hank_het_block_failed object.")
  .hank_manifest_check_peak(peak_rss, "hank_het_manifest")

  is_block <- inherits(x, "hank_het_block")

  n_e <- if (!is.null(x$n_e)) x$n_e else length(x$y)
  n_a <- if (!is.null(x$n_a)) x$n_a else length(x$a_grid)
  state_count <- n_e * n_a

  D <- x[["D", exact = TRUE]]
  have_dist <- !is.null(D) && is.numeric(D) && length(D) == n_e * n_a
  ## Cell order (hank_forward_operator / .hank_mat_to_vec): e slowest, a
  ## fastest -- array(D, c(n_a, n_e)) reproduces the vector directly.
  grid_mass <- if (have_dist)
    .hank_grid_mass_axes(D, c(n_a, n_e), c("a", NA_character_)) else NULL

  .hank_manifest_build(
    x, is_block = is_block, peak_rss = peak_rss,
    sizes = list(n_e = n_e, n_a = n_a),
    state_count = state_count, grid_mass = grid_mass,
    agg_fields = c("A", "C"),
    class = "hank_het_manifest")
}


#' Coerce a one-asset run manifest to one data-frame row
#'
#' Flattens a \code{\link{hank_het_manifest}} into exactly one row of atomic,
#' length-1 columns -- the shape a run-manifest table needs. The only
#' non-scalar field, \code{grid_mass}, is unnested into three
#' \code{grid_mass_*} columns (filled with \code{NA_real_} when the source
#' manifest carries no distribution).
#'
#' @param x A \code{\link{hank_het_manifest}}.
#' @param row.names,optional,... Ignored; present for S3 signature
#'   compatibility with \code{\link[base]{as.data.frame}}.
#' @return A one-row \code{data.frame}.
#' @seealso \code{\link{hank_het_manifest}}
#' @export
as.data.frame.hank_het_manifest <- function(x, row.names = NULL,
                                            optional = FALSE, ...) {
  .hank_manifest_as_df(x, c("a_lo", "a_hi", "interior"))
}


#' Print a one-asset run manifest
#'
#' @param x A \code{\link{hank_het_manifest}}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{hank_het_manifest}}
#' @export
print.hank_het_manifest <- function(x, ...) {
  .hank_manifest_print(x, "hank_het_manifest",
                       sprintf("e=%d a=%d", x$n_e, x$n_a),
                       c("A", "C"))
}


#' Content fingerprint of a one-asset household block
#'
#' A short, stable hash of everything that makes a \code{\link{hank_het_block}}
#' the block it is: the asset grid and borrowing constraint, the income
#' process (\code{Pi}, \code{e}), preferences, the prices (including the
#' transfer and its incidence, and any borrowing-rate wedge), the transition
#' inputs, the converged policies and the stationary distribution. The
#' one-asset counterpart of \code{\link{hank_het3_fingerprint}}.
#'
#' \strong{What this is for.} Keying a cache of expensive derived objects -- a
#' sequence-space Jacobian, a GE solve -- on the block that produced them.
#' Hashing the serialised block itself is the obvious alternative and is a
#' trap: it folds in fields with no bearing on the mathematics, so two
#' identical households hash differently and the cache never hits.
#'
#' \strong{What it deliberately excludes}: the sparse forward operator
#' \code{Lambda} (a deterministic function of the policy and \code{Pi} that is
#' already hashed), the derived aggregates \code{A} and \code{C}, and
#' \code{dist_converged}. If you need to distinguish HOW a block was produced,
#' that is what \code{\link{hank_het_manifest}} is for -- fingerprint answers
#' "is this the same household?", the manifest answers "how did this run go?".
#'
#' The field list is fixed and enumerated literally, never
#' \code{names(block)}, so a field appended to the block in future cannot
#' silently change every existing fingerprint.
#'
#' @param block A \code{\link{hank_het_block}}.
#' @param n Number of hex characters to return (1-16, default 16).
#' @return A length-1 character string of \code{n} lowercase hex digits.
#' @seealso \code{\link{hank_het_manifest}} (run provenance, not identity),
#'   \code{\link{hank_het2_fingerprint}}, \code{\link{hank_het3_fingerprint}}
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' ag  <- hank_asset_grid(20, 10, 0)
#' blk <- hank_het_block(ag, inc$Pi, inc$e, beta = 0.96, eis = 1,
#'                       r = 0.01, w = 1)
#' hank_het_fingerprint(blk)
#' @export
hank_het_fingerprint <- function(block, n = 16L) {
  if (!inherits(block, "hank_het_block"))
    stop("hank_het_fingerprint: 'block' must be a hank_het_block.")
  ## Order is fixed and explicit, NOT names(block).
  parts <- list(
    block$a_grid, block$Pi, block$e,
    block$beta, block$eis,
    block$r, block$w,
    .hank_block_tr(block), .hank_block_omega(block),
    block$r_minus, block$amin,
    block$a, block$c, block$D,
    ## Transition inputs by NAME as well as value: two blocks whose Pi_fn
    ## arguments differ only in name are different wirings of the DAG.
    names(block$Pi_inputs), unlist(block$Pi_inputs, use.names = FALSE)
  )
  .hank_fingerprint_hash(parts, n, "hank_het_fingerprint")
}
