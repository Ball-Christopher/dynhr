## R/hank-print.R
## --------------------------------------------------------------------------
## ONE compact print method for the HANK object family.
##
## Every HANK constructor stamps `hank_block` as the LAST element of its class
## vector (`class = c("<specific>", "hank_block")`), so this single method is
## the default printer for all of them. Before it existed the ~20 method-less
## HANK classes fell through to print.default(), which dumps every field --
## including the stationary distribution (an n_e x n_a, or n_e x n_d x n_f x
## n_a, array) -- across thousands of terminal lines.
##
## The method reads only fields that are PRESENT (via [[..., exact = TRUE]]),
## never invents one, and never prints an array. Constructors were not
## changed apart from their class vectors.
## --------------------------------------------------------------------------


## Grid-size fields, in the order a reader expects them (income states first,
## then the asset dimensions from the fastest-moving one outward).
.hank_print_grid_fields <- c("n_e", "n_z", "n_a", "n_b", "n_d", "n_f", "n_k")

## Aggregate / price fields worth a header line. Only length-1 finite numerics
## are shown; anything else is skipped silently.
.hank_print_agg_fields <- c("A", "B", "C", "D_agg", "F_agg", "A_agg",
                            "CHI", "PHI", "K", "Y", "N", "L",
                            "r", "rd", "rf", "ra", "rb", "w", "Tr")

## Run-metadata fields.
.hank_print_meta_fields <- c("backend", "threads", "iterations", "converged",
                             "dist_converged")

.hank_print_get <- function(x, nm) {
  if (!is.list(x)) return(NULL)
  x[[nm, exact = TRUE]]
}

## A scalar rendered for one header line, or NULL if it is not a scalar.
.hank_print_scalar <- function(v) {
  if (is.null(v) || length(v) != 1L) return(NULL)
  if (is.logical(v)) return(if (is.na(v)) "NA" else if (v) "TRUE" else "FALSE")
  if (is.character(v)) return(v)
  if (is.numeric(v)) {
    if (!is.finite(v)) return("NA")
    if (is.integer(v) || v == round(v)) return(format(v, scientific = FALSE))
    return(trimws(formatC(v, digits = 4, format = "g")))
  }
  NULL
}

## "k = v" pairs for the fields of `nms` that are present and scalar.
.hank_print_pairs <- function(x, nms) {
  out <- character(0)
  for (nm in nms) {
    s <- .hank_print_scalar(.hank_print_get(x, nm))
    if (!is.null(s)) out <- c(out, paste0(nm, " = ", s))
  }
  out
}

## Print `pairs` under `label`, wrapping at `width` and capping at `max_lines`
## so no single section can run away.
.hank_print_section <- function(label, pairs, width = 72L, max_lines = 4L) {
  if (!length(pairs)) return(invisible(NULL))
  ind <- paste0("  ", formatC(label, width = -10L), ": ")
  cont <- strrep(" ", nchar(ind))
  line <- ""
  lines <- character(0)
  for (p in pairs) {
    cand <- if (nzchar(line)) paste(line, p, sep = "  ") else p
    if (nchar(cand) > width && nzchar(line)) {
      lines <- c(lines, line); line <- p
    } else {
      line <- cand
    }
  }
  if (nzchar(line)) lines <- c(lines, line)
  if (length(lines) > max_lines)
    lines <- c(lines[seq_len(max_lines)], "...")
  cat(ind, lines[1], "\n", sep = "")
  for (l in lines[-1]) cat(cont, l, "\n", sep = "")
  invisible(NULL)
}


#' Compact print method for HANK objects
#'
#' The single default printer for the HANK object family: heterogeneous-agent
#' blocks (one-, two- and three-asset), block specs, models, steady states,
#' linearisations and emulators. Every HANK constructor carries
#' \code{"hank_block"} as the last element of its class vector, so this method
#' catches all of them.
#'
#' It prints a fixed-size header -- class, grid sizes, aggregates and prices,
#' backend / convergence metadata -- and \strong{never} prints the stationary
#' distribution, the policy arrays or the transition matrix, however large
#' they are. Use \code{names(x)} to list the fields and \code{x$D},
#' \code{x$Lambda}, ... to reach them.
#'
#' @param x   A HANK object (any class inheriting from \code{"hank_block"}).
#' @param ... Ignored, for S3 consistency.
#' @return \code{x}, invisibly.
#' @seealso \code{\link{hank_het_block}}, \code{\link{hank_het2_block}},
#'   \code{\link{hank_het3_block}}, \code{\link{hank_model}}
#' @examples
#' inc <- hank_income_rouwenhorst(0.9, 0.7, 5)
#' ag  <- hank_asset_grid(50, 100, 0)
#' blk <- hank_het_block(ag, inc$Pi, inc$e, beta = 0.96, eis = 1,
#'                       r = 0.01, w = 1.0)
#' print(blk)   # header only -- the distribution is NOT dumped
#' @export
print.hank_block <- function(x, ...) {
  cls <- class(x)[1]
  cat(sprintf("<%s>\n", cls))

  if (!is.list(x)) {
    cat("  (not a list; ", length(x), " element(s))\n", sep = "")
    return(invisible(x))
  }

  ## --- block SPECs (name / kind / inputs / outputs / block) -----------------
  kind <- .hank_print_get(x, "kind")
  inner <- .hank_print_get(x, "block")
  if (!is.null(kind) && is.character(kind)) {
    nm <- .hank_print_scalar(.hank_print_get(x, "name"))
    cat(sprintf("  %-10s: %s%s\n", "spec", kind,
                if (is.null(nm)) "" else paste0("  (name = ", nm, ")")))
    inp <- .hank_print_get(x, "inputs")
    outp <- .hank_print_get(x, "outputs")
    if (is.character(inp))
      .hank_print_section("inputs", inp)
    if (is.character(outp))
      .hank_print_section("outputs", outp)
  }

  ## Grid sizes / aggregates are read from the inner block when this is a spec.
  src <- if (!is.null(kind) && is.list(inner)) inner else x

  .hank_print_section("grid", .hank_print_pairs(src, .hank_print_grid_fields))
  .hank_print_section("aggregates", .hank_print_pairs(src, .hank_print_agg_fields))
  .hank_print_section("run", .hank_print_pairs(src, .hank_print_meta_fields))

  el <- c(.hank_print_get(src, "elapsed_solve"), .hank_print_get(src, "elapsed_dist"),
          .hank_print_get(src, "elapsed"))
  el <- suppressWarnings(as.numeric(el))
  el <- el[is.finite(el)]
  if (length(el))
    cat(sprintf("  %-10s: %.2f sec\n", "elapsed", sum(el)))

  ## --- model containers: list the blocks, not their contents ---------------
  blks <- .hank_print_get(x, "blocks")
  if (is.list(blks) && length(blks)) {
    bn <- vapply(blks, function(b) {
      n <- if (is.list(b)) b[["name", exact = TRUE]] else NULL
      k <- if (is.list(b)) b[["kind", exact = TRUE]] else NULL
      if (is.character(n) && length(n) == 1L)
        paste0(n, if (is.character(k) && length(k) == 1L) paste0(" [", k, "]") else "")
      else class(b)[1]
    }, character(1))
    .hank_print_section("blocks", bn)
  }

  nms <- names(x)
  cat(sprintf("  %-10s: %d (names(x) to list; nothing is printed here)\n",
              "fields", length(nms)))
  invisible(x)
}
