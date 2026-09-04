## R/hank-manifest3.R
## --------------------------------------------------------------------------
## A8: three-asset run-manifest metadata. The paper's release-request
## contract needs enough EXPOSED DATA (not printed text) per solve/block to
## record, in an external run manifest: backend, worker count, state count,
## iterations, wall time, peak RSS, convergence gaps, grid masses, package
## version and git commit. hank_egm3_solve()/hank_het3_block() already carry
## most of this; hank_het3_manifest() below is the single collector that
## reads it off either object and adds what neither carries on its own
## (package_version, git_commit, peak_rss_bytes).
## --------------------------------------------------------------------------


## Best-effort git commit resolution. Historic name, kept because it is
## referenced by the three-asset tests and docs; the implementation is the
## family-agnostic `.hank_manifest_git_commit()` in R/hank-manifest-common.R
## (nothing about resolving a commit is three-asset).
.hank3_resolve_git_commit <- function() .hank_manifest_git_commit()


#' Run-manifest metadata for a three-asset solve or block
#'
#' Collects, as DATA rather than printed text, everything the three-asset
#' paper's run manifest needs to record about one solve: backend, worker
#' count, state counts, iterations, wall time, peak memory, convergence
#' gaps, boundary/interior grid masses, package version and git commit. This
#' is a thin READER over fields \code{\link{hank_egm3_solve}} and
#' \code{\link{hank_het3_block}} already carry (iterations, converged,
#' gaps, state counts, backend, threads, elapsed timings), plus the handful
#' neither carries on its own.
#'
#' \code{x} may be either a \code{\link{hank_het3_block}} (a household
#' solve plus its stationary distribution) or a bare
#' \code{\link{hank_egm3_solve}} (no distribution). Fields that only make
#' sense with a distribution (\code{dist_converged}, the aggregates, and
#' \code{grid_mass} when \code{x} carries no \code{D}) are \code{NA} (or
#' \code{NULL} for \code{grid_mass}) on a bare solve, never an error.
#'
#' \strong{git_commit} resolution never fails the call: see
#' \code{.hank3_resolve_git_commit} for the fallback chain (live
#' \code{git rev-parse HEAD} in a source tree, then a release-time
#' \code{GIT_COMMIT} file, then \code{packageDescription()} remote fields,
#' then \code{NA_character_}). A missing or unreachable \code{git} binary is
#' an \code{NA}, not an error.
#'
#' \strong{peak_rss_bytes} is a PEAK-SINCE-PROCESS-START figure (see
#' \code{.hank_peak_rss}): monotone non-decreasing for the life of the R
#' session, not the memory this call or this solve alone used.
#'
#' @param x A \code{\link{hank_het3_block}}, a bare
#'   \code{\link{hank_egm3_solve}}, or a \code{hank_het3_block_failed} from
#'   \code{hank_het3_block(strict = FALSE)}. A failed object is unwrapped to
#'   the household it carries, so the manifest reports
#'   \code{converged = FALSE} together with the iterations and gaps that
#'   explain the failure, and leaves the block-only aggregate fields
#'   \code{NA} -- a calibration driver hitting the convergence frontier wants
#'   exactly that row.
#' @param peak_rss Logical (default \code{TRUE}). When \code{FALSE}, skip
#'   the \code{\link{.hank_peak_rss}} call and report \code{NA_real_} for
#'   \code{peak_rss_bytes} -- for callers who solve many small blocks in a
#'   loop and want to time only the OS call they actually need.
#'
#' @return An object of class \code{hank_het3_manifest}: a named list with
#'   \describe{
#'     \item{\code{package_version}}{Character:
#'       \code{utils::packageVersion("dynhr")}.}
#'     \item{\code{git_commit}}{Character, 40-hex-char SHA or
#'       \code{NA_character_}; see Details.}
#'     \item{\code{backend}}{Character \code{x$backend} (\code{"cpp"} or
#'       \code{"R"}), or \code{NA_character_} if \code{x} predates A8.}
#'     \item{\code{threads}}{Integer \code{x$threads}, or \code{NA_integer_}.}
#'     \item{\code{n_e}, \code{n_d}, \code{n_f}, \code{n_a}}{Grid sizes, read
#'       from \code{x} directly (blocks) or derived from \code{x$y} /
#'       \code{x$d_grid} / \code{x$f_grid} / \code{x$a_grid} (bare solves).}
#'     \item{\code{state_count}}{\code{x$state_count}.}
#'     \item{\code{iterations}, \code{converged}}{\code{x$iterations},
#'       \code{x$converged}.}
#'     \item{\code{last_value_gap}, \code{last_policy_gap}}{As carried by
#'       \code{x} (\code{NA} on the singleton-\code{f_grid} reduction path,
#'       same as \code{x} itself).}
#'     \item{\code{dist_converged}}{\code{x$dist_converged} for a block;
#'       \code{NA} for a bare solve.}
#'     \item{\code{elapsed_solve}}{\code{x$elapsed_solve} (blocks) or
#'       \code{x$elapsed} (bare solves); \code{NA_real_} if neither is
#'       present.}
#'     \item{\code{elapsed_dist}}{\code{x$elapsed_dist} (blocks only);
#'       \code{NA_real_} otherwise.}
#'     \item{\code{peak_rss_bytes}}{\code{\link{.hank_peak_rss}}\code{()}
#'       when \code{peak_rss = TRUE} and the platform can report it;
#'       \code{NA_real_} otherwise.}
#'     \item{\code{grid_mass}}{The named length-7 vector from
#'       \code{.hank3_grid_mass} (\code{d_lo}, \code{d_hi}, \code{f_lo},
#'       \code{f_hi}, \code{a_lo}, \code{a_hi}, \code{interior}) when a
#'       stationary distribution is available on \code{x}; \code{NULL}
#'       otherwise. \code{as.data.frame()} flattens this into seven
#'       \code{grid_mass_*} columns, filled with \code{NA_real_} when this
#'       field is \code{NULL}.}
#'     \item{\code{transition_inputs}}{Character: the block's
#'       \code{Pi_inputs} names and steady-state values (e.g.
#'       \code{"f=0.6, s=0.1"}), or \code{NA_character_} for a price-only
#'       block. Provenance rather than diagnostics -- two runs whose manifests
#'       agree on every other field but differ here were solved against
#'       DIFFERENT households, so a Jacobian cache keyed on the rest of the
#'       row would silently conflate them. Values are formatted at full
#'       double precision for the same reason.}
#'     \item{\code{D_agg}, \code{F_agg}, \code{A_agg}, \code{C}, \code{CHI},
#'       \code{PHI}}{Aggregate holdings/consumption/adjustment-cost
#'       resources, for a block; \code{NA_real_} for a bare solve.}
#'   }
#' @seealso \code{\link{hank_egm3_solve}}, \code{\link{hank_het3_block}},
#'   \code{\link{hank_euler3_residual}} (the other caller of
#'   \code{.hank3_grid_mass}), \code{\link{as.data.frame.hank_het3_manifest}},
#'   \code{\link{print.hank_het3_manifest}}
#' @examples
#' Pi <- matrix(c(.9, .1, .1, .9), 2)
#' dg <- seq(-.1, .5, length.out = 3); fg <- seq(0, .5, length.out = 3)
#' ag <- seq(0, .6, length.out = 3);   e  <- c(.8, 1.1)
#' blk <- hank_het3_block(dg, fg, ag, Pi, e, beta = .97, eis = .5,
#'                        rd = .01, rf = .015, ra = .02, w = 1,
#'                        chi1 = .2, phi1 = .1, maxit = 250)
#' mf <- hank_het3_manifest(blk)
#' mf$state_count
#' as.data.frame(mf)
#' @export
hank_het3_manifest <- function(x, peak_rss = TRUE) {
  ## A calibration driver that hits the convergence frontier is exactly when a
  ## manifest row is most useful: iterations against the cap, both gaps and the
  ## wall time ARE the evidence for what to change next. So accept the failed
  ## object rather than making the caller know to reach inside it, and unwrap to
  ## the household it carries -- whose converged = FALSE and gaps then flow
  ## through the ordinary path unchanged. Block-only fields (aggregates,
  ## dist_converged) stay NA, which is honest: the block never completed.
  if (inherits(x, "hank_het3_block_failed")) x <- x$hh
  if (!(inherits(x, "hank_het3_block") || inherits(x, "hank_egm3_solve")))
    stop("hank_het3_manifest: 'x' must be a hank_het3_block, a ",
         "hank_egm3_solve, or a hank_het3_block_failed object.")
  .hank_manifest_check_peak(peak_rss, "hank_het3_manifest")

  is_block <- inherits(x, "hank_het3_block")

  n_e <- if (!is.null(x$n_e)) x$n_e else length(x$y)
  n_d <- if (!is.null(x$n_d)) x$n_d else length(x$d_grid)
  n_f <- if (!is.null(x$n_f)) x$n_f else length(x$f_grid)
  n_a <- if (!is.null(x$n_a)) x$n_a else length(x$a_grid)
  state_count <- if (!is.null(x$state_count)) x$state_count else n_e * n_d * n_f * n_a

  D <- x$D
  have_dist <- !is.null(D) && is.numeric(D) &&
    length(D) == n_e * n_d * n_f * n_a
  grid_mass <- if (have_dist) .hank3_grid_mass(D, n_e, n_d, n_f, n_a) else NULL

  .hank_manifest_build(
    x, is_block = is_block, peak_rss = peak_rss,
    sizes = list(n_e = n_e, n_d = n_d, n_f = n_f, n_a = n_a),
    state_count = state_count, grid_mass = grid_mass,
    agg_fields = c("D_agg", "F_agg", "A_agg", "C", "CHI", "PHI"),
    class = "hank_het3_manifest")
}


#' Coerce a three-asset run manifest to one data-frame row
#'
#' Flattens a \code{\link{hank_het3_manifest}} into exactly one row of
#' atomic, length-1 columns -- the shape a run-manifest table (the paper
#' \code{rbind}s or \code{dplyr::bind_rows}s one of these per solve) needs.
#' The only non-scalar field, \code{grid_mass}, is unnested into seven
#' \code{grid_mass_*} columns (filled with \code{NA_real_} when the source
#' manifest carries no distribution).
#'
#' @param x A \code{\link{hank_het3_manifest}}.
#' @param row.names,optional,... Ignored; present for S3 signature
#'   compatibility with \code{\link[base]{as.data.frame}}.
#' @return A one-row \code{data.frame}.
#' @seealso \code{\link{hank_het3_manifest}}
#' @export
as.data.frame.hank_het3_manifest <- function(x, row.names = NULL,
                                              optional = FALSE, ...) {
  .hank_manifest_as_df(
    x, c("d_lo", "d_hi", "f_lo", "f_hi", "a_lo", "a_hi", "interior"))
}


#' Print a three-asset run manifest
#'
#' @param x A \code{\link{hank_het3_manifest}}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{hank_het3_manifest}}
#' @export
print.hank_het3_manifest <- function(x, ...) {
  .hank_manifest_print(
    x, "hank_het3_manifest",
    sprintf("e=%d d=%d f=%d a=%d", x$n_e, x$n_d, x$n_f, x$n_a),
    c("D_agg", "F_agg", "A_agg", "C", "CHI", "PHI"))
}
