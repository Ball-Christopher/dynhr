## R/hank-manifest2.R
## --------------------------------------------------------------------------
## TWO-asset (liquid/illiquid, convex adjustment cost) provenance tooling: the
## run manifest and the content fingerprint, the direct counterparts of
## hank_het3_manifest() / hank_het3_fingerprint().
##
## Everything family-agnostic lives in R/hank-manifest-common.R; this file
## carries only what is genuinely two-asset: which objects are accepted, the
## (b, a) asset axes, the (B, A, C, CHI) aggregates, and the literal
## fingerprint field list.
##
## The DISCRETE-adjustment block (hank_het2d_block) is a different household
## with different fields and is deliberately REJECTED here rather than
## silently mis-read -- same policy as .hank_reject_het2().
## --------------------------------------------------------------------------


#' Is this a bare two-asset household solve?
#'
#' \code{\link{hank_egm2_solve}} returns a PLAIN LIST (unlike
#' \code{\link{hank_egm3_solve}}, which carries an S3 class), so the manifest
#' has to duck-type it: two asset grids, two marginal values, and no
#' three-asset grid.
#'
#' @param x The object to classify.
#' @return \code{TRUE} for a bare \code{\link{hank_egm2_solve}} return value,
#'   \code{FALSE} otherwise.
#' @keywords internal
.hank_is_egm2_solve <- function(x) {
  is.list(x) &&
    !is.null(x[["b_grid", exact = TRUE]]) &&
    !is.null(x[["a_grid", exact = TRUE]]) &&
    !is.null(x[["Vb", exact = TRUE]]) &&
    !is.null(x[["Va", exact = TRUE]]) &&
    !is.null(x[["b", exact = TRUE]]) &&
    !is.null(x[["a", exact = TRUE]]) &&
    !is.null(x[["c", exact = TRUE]]) &&
    is.null(x[["d_grid", exact = TRUE]])
}


#' Run-manifest metadata for a two-asset solve or block
#'
#' The two-asset counterpart of \code{\link{hank_het3_manifest}}: collects, as
#' DATA rather than printed text, everything a run manifest needs to record
#' about one household solve -- backend, worker count, state count,
#' iterations, wall time, peak memory, convergence flags, boundary/interior
#' grid masses on BOTH asset axes, package version and git commit. It is a
#' thin READER over what \code{\link{hank_het2_block}} and
#' \code{\link{hank_egm2_solve}} already carry, plus the handful neither
#' carries on its own.
#'
#' \code{x} may be either a \code{\link{hank_het2_block}} (a household solve
#' plus its stationary joint distribution) or a bare
#' \code{\link{hank_egm2_solve}} return value (no distribution). Fields that
#' only make sense with a distribution (\code{dist_converged}, the aggregates,
#' and \code{grid_mass} when \code{x} carries no \code{D}) are \code{NA} (or
#' \code{NULL} for \code{grid_mass}) on a bare solve, never an error.
#'
#' A \code{\link{hank_het2d_block}} (the DISCRETE-adjustment two-asset
#' household) is rejected: it is a different problem with different fields,
#' and reading it through this collector would report an aggregate set that
#' does not describe it.
#'
#' \strong{Fields the two-asset family does not carry.}
#' \code{\link{hank_het2_block}} does not record \code{backend},
#' \code{threads}, \code{iterations}, \code{converged}, the per-iteration gaps
#' or any timing, so those are \code{NA} for a block (a bare
#' \code{\link{hank_egm2_solve}} does carry \code{iterations} and
#' \code{converged}, and they are reported).
#'
#' \strong{git_commit} resolution never fails the call, and
#' \strong{peak_rss_bytes} is a peak-since-process-start figure: see
#' \code{\link{hank_het3_manifest}} for both.
#'
#' @param x A \code{\link{hank_het2_block}}, a bare
#'   \code{\link{hank_egm2_solve}} return value, or a
#'   \code{hank_het2_block_failed} wrapper (unwrapped to the household it
#'   carries, so the manifest reports the iterations and gaps that explain the
#'   failure while leaving the block-only aggregate fields \code{NA}).
#' @param peak_rss Logical (default \code{TRUE}). When \code{FALSE}, skip the
#'   \code{\link{.hank_peak_rss}} call and report \code{NA_real_} for
#'   \code{peak_rss_bytes} -- for callers who solve many small blocks in a loop.
#'
#' @return An object of class \code{hank_het2_manifest}: a named list with
#'   \code{package_version}, \code{git_commit}, \code{backend},
#'   \code{threads}, \code{n_e}, \code{n_b}, \code{n_a}, \code{state_count},
#'   \code{iterations}, \code{converged}, \code{last_value_gap},
#'   \code{last_policy_gap}, \code{dist_converged}, \code{elapsed_solve},
#'   \code{elapsed_dist}, \code{peak_rss_bytes}, \code{grid_mass} (the named
#'   \code{b_lo}, \code{b_hi}, \code{a_lo}, \code{a_hi}, \code{interior}
#'   vector, or \code{NULL} with no distribution), \code{transition_inputs},
#'   and the aggregates \code{B}, \code{A}, \code{C}, \code{CHI}.
#' @seealso \code{\link{hank_het2_fingerprint}} (identity, not provenance),
#'   \code{\link{hank_het_manifest}}, \code{\link{hank_het3_manifest}},
#'   \code{\link{as.data.frame.hank_het2_manifest}},
#'   \code{\link{print.hank_het2_manifest}}
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' blk <- hank_het2_block(hank_asset_grid(20, 8L, 0),
#'                        hank_asset_grid(30, 8L, 0), inc$Pi, inc$e,
#'                        beta = 0.95, eis = 0.5, rb = 0.005, ra = 0.02,
#'                        w = 1, n_k = 8L)
#' mf <- hank_het2_manifest(blk)
#' mf$state_count
#' as.data.frame(mf)
#' @export
hank_het2_manifest <- function(x, peak_rss = TRUE) {
  if (inherits(x, "hank_het2_block_failed")) x <- x$hh
  if (inherits(x, "hank_het2d_block"))
    stop("hank_het2_manifest: this is a DISCRETE-ADJUSTMENT two-asset block ",
         "(hank_het2d_block), whose fields and aggregates differ; there is ",
         "no het2d manifest.")
  if (inherits(x, "hank_het_block") || inherits(x, "hank_het3_block") ||
      !(inherits(x, "hank_het2_block") || .hank_is_egm2_solve(x)))
    stop("hank_het2_manifest: 'x' must be a hank_het2_block, a bare ",
         "hank_egm2_solve return value, or a hank_het2_block_failed object.")
  .hank_manifest_check_peak(peak_rss, "hank_het2_manifest")

  is_block <- inherits(x, "hank_het2_block")

  n_e <- if (!is.null(x$n_e)) x$n_e else length(x$y)
  n_b <- if (!is.null(x$n_b)) x$n_b else length(x$b_grid)
  n_a <- if (!is.null(x$n_a)) x$n_a else length(x$a_grid)
  state_count <- n_e * n_b * n_a

  D <- x[["D", exact = TRUE]]
  have_dist <- !is.null(D) && is.numeric(D) && length(D) == n_e * n_b * n_a
  ## Cell order (hank_forward_operator2 / .hank2_arr_to_vec): e slowest, then
  ## b, with a fastest -- array(D, c(n_a, n_b, n_e)) reproduces the vector
  ## directly.
  grid_mass <- if (have_dist)
    .hank_grid_mass_axes(D, c(n_a, n_b, n_e), c("a", "b", NA_character_))
    else NULL

  .hank_manifest_build(
    x, is_block = is_block, peak_rss = peak_rss,
    sizes = list(n_e = n_e, n_b = n_b, n_a = n_a),
    state_count = state_count, grid_mass = grid_mass,
    agg_fields = c("B", "A", "C", "CHI"),
    class = "hank_het2_manifest")
}


#' Coerce a two-asset run manifest to one data-frame row
#'
#' Flattens a \code{\link{hank_het2_manifest}} into exactly one row of atomic,
#' length-1 columns. The only non-scalar field, \code{grid_mass}, is unnested
#' into five \code{grid_mass_*} columns (filled with \code{NA_real_} when the
#' source manifest carries no distribution).
#'
#' @param x A \code{\link{hank_het2_manifest}}.
#' @param row.names,optional,... Ignored; present for S3 signature
#'   compatibility with \code{\link[base]{as.data.frame}}.
#' @return A one-row \code{data.frame}.
#' @seealso \code{\link{hank_het2_manifest}}
#' @export
as.data.frame.hank_het2_manifest <- function(x, row.names = NULL,
                                             optional = FALSE, ...) {
  .hank_manifest_as_df(x, c("b_lo", "b_hi", "a_lo", "a_hi", "interior"))
}


#' Print a two-asset run manifest
#'
#' @param x A \code{\link{hank_het2_manifest}}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{hank_het2_manifest}}
#' @export
print.hank_het2_manifest <- function(x, ...) {
  .hank_manifest_print(x, "hank_het2_manifest",
                       sprintf("e=%d b=%d a=%d", x$n_e, x$n_b, x$n_a),
                       c("B", "A", "C", "CHI"))
}


#' Content fingerprint of a two-asset household block
#'
#' A short, stable hash of everything that makes a
#' \code{\link{hank_het2_block}} the block it is: both asset grids and the
#' multiplier grid, the income process (\code{Pi}, \code{e}), preferences, the
#' adjustment-cost parameters, the prices (including the transfer and the LTV
#' \code{theta_coll}), the transition inputs, the converged policies and the
#' stationary joint distribution. The two-asset counterpart of
#' \code{\link{hank_het3_fingerprint}}.
#'
#' \strong{What this is for.} Keying a cache of expensive derived objects -- a
#' long-horizon Jacobian, a GE solve -- on the block that produced them, at a
#' scale where hashing the serialised block itself is a trap (it folds in
#' fields with no bearing on the mathematics, so two identical households
#' solved on different days hash differently and the cache never hits).
#'
#' \strong{What it deliberately excludes}: the sparse forward operator
#' \code{Lambda} (a deterministic function of the policies and \code{Pi} that
#' are already hashed), the derived aggregates, and \code{dist_converged}. Use
#' \code{\link{hank_het2_manifest}} to distinguish HOW a block was produced --
#' fingerprint answers "is this the same household?", the manifest answers
#' "how did this run go?".
#'
#' Both the gap policy \code{b} and the TRUE liquid position \code{b_liq} are
#' hashed: under collateral (\code{theta_coll > 0}) they differ, and a cache
#' keyed on only one of them would conflate two genuinely different
#' households.
#'
#' The field list is fixed and enumerated literally, never
#' \code{names(block)}, so a field appended to the block in future cannot
#' silently change every existing fingerprint.
#'
#' @param block A \code{\link{hank_het2_block}}.
#' @param n Number of hex characters to return (1-16, default 16).
#' @return A length-1 character string of \code{n} lowercase hex digits.
#' @seealso \code{\link{hank_het2_manifest}} (run provenance, not identity),
#'   \code{\link{hank_het_fingerprint}}, \code{\link{hank_het3_fingerprint}}
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 3)
#' blk <- hank_het2_block(hank_asset_grid(20, 8L, 0),
#'                        hank_asset_grid(30, 8L, 0), inc$Pi, inc$e,
#'                        beta = 0.95, eis = 0.5, rb = 0.005, ra = 0.02,
#'                        w = 1, n_k = 8L)
#' hank_het2_fingerprint(blk)
#' @export
hank_het2_fingerprint <- function(block, n = 16L) {
  if (!inherits(block, "hank_het2_block"))
    stop("hank_het2_fingerprint: 'block' must be a hank_het2_block.")
  ## Order is fixed and explicit, NOT names(block).
  parts <- list(
    block$b_grid, block$a_grid, block$k_grid, block$Pi, block$e,
    block$beta, block$eis,
    block$rb, block$ra, block$w,
    .hank_block_tr(block),
    if (is.null(block$theta_coll)) 0 else block$theta_coll,
    block$chi0, block$chi1, block$chi2,
    block$b, block$a, block$c, block$chi, block$b_liq, block$D,
    ## Transition inputs by NAME as well as value: two blocks whose Pi_fn
    ## arguments differ only in name are different wirings of the DAG.
    names(block$Pi_inputs), unlist(block$Pi_inputs, use.names = FALSE)
  )
  .hank_fingerprint_hash(parts, n, "hank_het2_fingerprint")
}
