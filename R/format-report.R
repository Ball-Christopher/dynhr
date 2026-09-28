## R/format-report.R
## Shared line-builders for dynhr's textual reports.
##
## The house pattern, established by R/diag-summary.R and adopted here:
##
##   format_<x>()  -- pure; returns a character vector of LINES, no I/O
##   print_<x>()   -- ONE emit call, returns the lines invisibly
##   write_<x>()   -- the same lines to a file
##
## The reason this matters beyond tidiness: a report that emits once works on
## ANY stream. The hand-rolled alternative -- `cat(sprintf(...))` per cell
## inside a nested loop -- is bound to `cat()` forever, because every call is a
## separate write and a message stream turns each one into its own line. That
## is exactly how print_moments() rendered a correlation table as one line per
## cell when it was migrated.
##
## Build lines, emit once, and the choice of stream stops being a property of
## the formatting code.

#' Format a numeric matrix as fixed-width report lines
#'
#' @param m Numeric matrix. Column labels come from `colnames(m)` when
#'   `col_labels` is not given.
#' @param row_labels Row labels; defaults to `rownames(m)`.
#' @param col_labels Column labels; defaults to `colnames(m)`.
#' @param row_width Width of the leading label column.
#' @param col_width Width of each value column.
#' @param digits Decimal places for the values.
#' @param max_label Truncate labels to this many characters (`NULL` to keep
#'   them whole).
#' @param corner Text for the header's top-left cell.
#' @return Character vector: one header line followed by one line per row.
#' @noRd
.fmt_matrix_lines <- function(m, row_labels = rownames(m),
                              col_labels = colnames(m),
                              row_width = 15L, col_width = 10L,
                              digits = 4L, max_label = NULL, corner = "") {
  m <- as.matrix(m)
  if (is.null(row_labels)) row_labels <- as.character(seq_len(nrow(m)))
  if (is.null(col_labels)) col_labels <- as.character(seq_len(ncol(m)))
  trunc1 <- function(x) if (is.null(max_label)) x else substr(x, 1L, max_label)
  hdr <- paste0(formatC(corner, width = -row_width),
                paste(formatC(trunc1(col_labels), width = col_width),
                      collapse = ""))
  body <- vapply(seq_len(nrow(m)), function(i) {
    paste0(formatC(trunc1(row_labels[i]), width = -row_width),
           paste(formatC(m[i, ], format = "f", digits = digits,
                         width = col_width),
                 collapse = ""))
  }, character(1))
  c(hdr, body)
}

#' A section rule of `n` dashes
#' @noRd
.fmt_rule <- function(n = 70L) strrep("-", n)

#' Format a two-column "label: value" table as report lines
#'
#' @param labels Character vector of row labels.
#' @param ... Named numeric/character vectors, one per value column; the names
#'   become the column headers.
#' @param label_width,col_width Field widths.
#' @param digits Decimal places for numeric columns.
#' @return Character vector: header, then one line per label.
#' @noRd
.fmt_cols_lines <- function(labels, ..., label_width = 20L, col_width = 12L,
                            digits = 6L) {
  cols <- list(...)
  hdr <- paste0(formatC("Variable", width = -label_width), " ",
                paste(formatC(names(cols), width = col_width), collapse = " "))
  body <- vapply(seq_along(labels), function(i) {
    vals <- vapply(cols, function(cl) {
      v <- cl[[i]]
      if (is.numeric(v)) formatC(v, format = "f", digits = digits,
                                 width = col_width)
      else formatC(as.character(v), width = col_width)
    }, character(1))
    paste0(formatC(labels[i], width = -label_width), " ",
           paste(vals, collapse = " "))
  }, character(1))
  c(hdr, body)
}


# ---------------------------------------------------------------------------
# The rendering classification
# ---------------------------------------------------------------------------
# dynhr writes to the console for two different reasons, and they belong on
# different streams:
#
#   RENDERING -- an object or report showing itself because the user asked it
#     to. Its output IS the return value in visible form, so it goes to stdout
#     via cat()/print(), where capture.output() and sink() can take it.
#
#   MESSAGING -- progress, diagnostics and warnings emitted while doing
#     something else. These go through .dynhr_cat()/.dynhr_inform() so
#     dynhr_set_verbosity() can silence them.
#
# S3 methods (print./summary./format./show.) are rendering by definition, and
# so is anything they delegate to -- the rule is a CALL-GRAPH CLOSURE, because
# print.hank_het3_manifest() delegates its whole body to .hank_manifest_print().
#
# This vector names the rendering entry points that are NOT S3 methods and so
# cannot be found by that pattern. `print_moments` is the cautionary one: its
# name uses an UNDERSCORE, so `^print[.]` misses it entirely.
#
# Adding a name here is a deliberate act meaning "this function's printed
# output is its product". If instead the output is progress inside a longer
# operation, do NOT add it -- use .dynhr_cat() and let the verbosity level
# govern it. (stoch_simul() does exactly that with format_moments().)
#
# tests/testthat/test-conditions.R reads this vector and asserts that no other
# function in R/ still calls cat().
.dynhr_report_printers <- c(
  "verify_steady_state",    # dynhr:::verify_steady_state() -- a report
  "hmc_summary",            # dynhr:::hmc_summary()         -- a summary table
  "smc_summary",            # dynhr:::smc_summary()         -- a summary table
  "print_moments",          # the print half of format_moments()
  "print_model_registry"    # exported registry listing
)


# ---------------------------------------------------------------------------
# Escaping at the sink boundary
# ---------------------------------------------------------------------------
# Diagnostic `$summary` / `$model_name` text is user- and model-derived, and it
# is interpolated into three different sinks: raw HTML, Pandoc markdown, and a
# markdown pipe-table cell (which Typst and HTML both go through). Before 0.9.4
# each sink had its own partial escaper in the Quarto templates, and each one
# missed a different character: the HTML one escaped only `[` and `]` (so `<`
# opened a tag and ate the rest of the cell), the PDF one escaped `@ < >` but
# not the pipe (so a multi-line summary was silently truncated at the first
# line). `.escape_for()` is the single boundary function; the templates call
# nothing else.
#
# Mechanism: HTML gets entity escaping; every markdown-flavoured target gets
# Pandoc's `escaped_char` (a backslash before an ASCII punctuation character
# makes it literal in EVERY writer -- html, typst and latex alike).

#' Escape text for one of the report sinks
#'
#' @param x      Character vector (coerced). `NA` is preserved as `NA`.
#' @param target One of `"html"` (entity escaping for raw HTML),
#'   `"markdown"` (Pandoc markdown body text), `"md_table"` (a markdown
#'   pipe-table cell -- newlines become `"; "`), `"typst"` (markdown headed for
#'   the Typst writer; same rules as `"markdown"`), or `"title"` (a document
#'   title -- markdown escaping with whitespace collapsed).
#' @return Character vector the same length as `x`.
#' @noRd
.escape_for <- function(x, target = c("html", "markdown", "md_table",
                                      "typst", "title")) {
  target <- match.arg(target)
  if (length(x) == 0L) return(character(0))
  x  <- as.character(x)
  na <- is.na(x)
  x[na] <- ""

  out <- if (identical(target, "html")) {
    ## `&` MUST go first, or the ampersands introduced below are re-escaped.
    y <- gsub("&", "&amp;",   x, fixed = TRUE)
    y <- gsub("<", "&lt;",    y, fixed = TRUE)
    y <- gsub(">", "&gt;",    y, fixed = TRUE)
    y <- gsub("\"", "&quot;", y, fixed = TRUE)
    y <- gsub("'", "&#39;", y, fixed = TRUE)
    ## The raw-HTML sinks in the templates (dashboard cells, <li> summary
    ## lines, <p> action/explanation text) are INLINE HTML inside a Pandoc
    ## document, so Pandoc still parses markdown between the tags: `*star*`
    ## became <em>, `[PASS]` an implicit link, `$x$` a MathJax span (review
    ## A5/A6). Backslash-escape the markdown-active set as well; Pandoc
    ## consumes the backslashes, so the rendered text is the original.
    y <- gsub("\\", "\\\\", y, fixed = TRUE)
    gsub("([][*`$|])", "\\\\\\1", y)
  } else {
    ## Backslash first, for the same reason.
    y <- gsub("\\", "\\\\", x, fixed = TRUE)
    gsub("([][*_`$<>@#~^|{}])", "\\\\\\1", y)
  }

  if (identical(target, "md_table"))
    out <- gsub("[\r\n]+", "; ", out)
  if (identical(target, "title"))
    out <- trimws(gsub("[[:space:]]+", " ", out))

  out[na] <- NA_character_
  out
}


#' Format numbers to a fixed number of significant digits, without signed zero
#'
#' `sprintf("%.6f", x)` gives a state root of 0.85 six decimal places of
#' spurious precision and turns a numerically-zero root into the string
#' `"-0.000000"`. This formats to `digits` significant figures and snaps any
#' value with `abs(x) <= zap` to a plain `"0"`.
#'
#' @param x      Numeric vector.
#' @param digits Significant digits (default 4).
#' @param zap    Absolute magnitude at or below which `x` prints as `"0"`.
#' @return Character vector.
#' @noRd
.fmt_sig <- function(x, digits = 4L, zap = 0) {
  v <- as.numeric(x)
  vapply(v, function(z) {
    if (is.na(z))       return("NA")
    if (!is.finite(z))  return(if (z > 0) "Inf" else "-Inf")
    if (abs(z) <= zap)  return("0")
    s <- trimws(formatC(z, format = "g", digits = digits))
    ## formatC() can still emit "-0" for a value that rounds to zero.
    if (grepl("^-0(\\.0*)?$", s)) s <- sub("^-", "", s)
    s
  }, character(1))
}


## ---------------------------------------------------------------------------
## Report figure sizing and pagination
##
## Both report templates (inst/templates/report.qmd, report-pdf.qmd) used to
## carry their own copy of the height rule; the copies had already drifted in
## how they clamped the result.  The rule lives HERE now and the templates
## call the spooled copy (see .report_plotutil_source()), so there is one.
##
## The second half of the story is that a height rule alone cannot save a
## 68-parameter facet plot: the templates cap the figure at 40 in (HTML) /
## 9 in (PDF), so 12-18 facet rows arrive as unreadable 0.5 in strips.
## .paginate_plot() splits such a plot into several plots of readable height
## instead of squashing one.
## ---------------------------------------------------------------------------

#' Figure height for a plot that carries no \code{dynhr_fig_height}
#'
#' Grows with the number of facet ROWS and with the number of discrete y
#' levels (a 20-observable ridge plot or a 20 x 20 correlation heatmap at the
#' flat default was unreadable at NZSIM scale), clamped to
#' \code{[default, maximum]}.
#'
#' @param p       A ggplot.
#' @param default Height used when nothing suggests a taller figure.
#' @param maximum Hard ceiling (the template's page / ggsave limit).
#' @return A single numeric height in inches.
#' @noRd
.fig_height_for <- function(p, default = 6.5, maximum = 40) {
  if (!inherits(p, "gg") || !requireNamespace("ggplot2", quietly = TRUE))
    return(default)
  b   <- ggplot2::ggplot_build(p)
  lay <- b$layout$layout
  n_panel <- if (is.data.frame(lay) && nrow(lay) > 0L) nrow(lay) else 1L
  n_col   <- if (is.data.frame(lay) && !is.null(lay$COL) && length(lay$COL) > 0L)
    max(lay$COL) else 1L
  n_row   <- ceiling(max(n_panel, 1L) / max(n_col, 1L))
  pp  <- b$layout$panel_params
  ybr <- if (length(pp) >= 1L && !is.null(pp[[1L]]$y)) pp[[1L]]$y$breaks else NULL
  n_y <- if (is.null(ybr) || !is.atomic(ybr)) 0L else sum(!is.na(ybr))
  h <- default
  if (n_row > 1L) h <- max(h, 1.6 + 1.25 * n_row)
  if (n_y > 8L)   h <- max(h, 1.6 + 0.28 * n_y)
  min(max(h, default), maximum)
}


#' Facet variable names of a ggplot
#'
#' Handles \code{facet_wrap} (\code{params$facets}) and \code{facet_grid}
#' (\code{params$rows} / \code{params$cols}).
#' @noRd
.facet_var_names <- function(p) {
  f <- p$facet
  if (is.null(f) || inherits(f, "FacetNull")) return(character(0))
  pr <- f$params
  if (is.null(pr)) return(character(0))
  nms <- c(names(pr$facets), names(pr$rows), names(pr$cols))
  nms <- nms[!is.na(nms) & nzchar(nms)]
  unique(nms)
}


#' Name of the variable mapped to one aesthetic, or NA
#'
#' @param q A quosure / formula / symbol taken from a ggplot mapping.
#' @noRd
.aes_var_name <- function(q) {
  if (is.null(q)) return(NA_character_)
  e <- q
  if (inherits(e, "quosure") || inherits(e, "formula")) {
    ## A quosure is a one-sided formula `~expr` (plus an env attribute); its
    ## expression is the last element, as for a formula. unclass() first so
    ## rlang's deprecated `[[.quosure` method is not dispatched -- this keeps
    ## the helper free of an rlang dependency (it also runs, as source, in
    ## the Quarto subprocess; see .report_plotutil_source()).
    e <- unclass(e)
    e <- e[[length(e)]]
  }
  if (is.name(e)) as.character(e) else NA_character_
}


#' Row key for a set of facet / level columns
#' @noRd
.paginate_key <- function(df) {
  if (!is.data.frame(df) || ncol(df) == 0L) return(character(0))
  do.call(paste, c(lapply(df, as.character), list(sep = "\r")))
}


#' Copy one ggplot layer, giving the copy its own data
#'
#' A layer is a ggproto object, i.e. an environment: assigning to
#' \code{layer$data} would mutate the layer shared by every page AND by the
#' caller's original plot.  \code{ggproto(NULL, layer, data = )} is not an
#' option either (it sends ggplot2 4.x's S7 layer machinery into infinite
#' recursion), so the environment is copied by hand.
#' @noRd
.paginate_clone_layer <- function(l, dat) {
  e <- new.env(parent = parent.env(l))
  for (nm in ls(l, all.names = TRUE)) assign(nm, get(nm, envir = l), envir = e)
  class(e) <- class(l)
  e$data <- dat
  e
}


#' One page: the plot with its data -- and any layer data carrying the same
#' columns -- restricted to \code{keys}
#' @noRd
.paginate_page <- function(p, sub, vars, keys) {
  pg   <- p + sub
  lyrs <- pg$layers
  if (length(lyrs) > 0L) {
    for (i in seq_along(lyrs)) {
      ld <- lyrs[[i]]$data
      if (!is.data.frame(ld) || nrow(ld) == 0L) next
      if (!all(vars %in% names(ld))) next
      keep <- .paginate_key(ld[, vars, drop = FALSE]) %in% keys
      lyrs[[i]] <- .paginate_clone_layer(lyrs[[i]], ld[keep, , drop = FALSE])
    }
    pg$layers <- lyrs
  }
  pg
}


#' Title the pages and give each its own height attribute
#' @noRd
.paginate_finish <- function(p, pages, default, maximum) {
  n   <- length(pages)
  ttl <- p$labels$title
  ttl <- if (is.null(ttl) || !is.character(ttl) || !nzchar(ttl[[1L]])) ""
         else ttl[[1L]]
  lapply(seq_len(n), function(k) {
    tag <- sprintf("(page %d of %d)", k, n)
    pk  <- pages[[k]] +
      ggplot2::labs(title = if (nzchar(ttl)) paste(ttl, tag) else tag)
    attr(pk, "dynhr_fig_height") <- .fig_height_for(pk, default, maximum)
    pk
  })
}


#' Split an over-tall report figure into several readable pages
#'
#' Two cases are split; everything else is returned unchanged (the helper
#' never errors, so a template may call it on every plot):
#'
#' \itemize{
#'   \item a faceted plot with more than \code{max_rows} panel rows -- the
#'     facet levels are chunked, in their existing order, into pages of
#'     \code{max_rows} rows each;
#'   \item a single-panel plot whose discrete axis (y, or x when the bars are
#'     flipped) carries more than \code{max_y} levels -- the levels are
#'     chunked the same way and the other axis is left whole, so a heatmap
#'     keeps all its columns on every page.
#' }
#'
#' @param p        A ggplot.
#' @param max_rows Facet rows per page.
#' @param max_y    Discrete-axis levels per page.
#' @param max_pages Most pages a figure may be split into. A figure that would
#'   need more (the NZSIM 1554-row D1 Jacobian heatmap asked for 60 PDF pages)
#'   is returned unchanged: it is a lookup table rather than a picture, and
#'   sixty pages of it help nobody.
#' @param default,maximum Passed to \code{.fig_height_for()} for the per-page
#'   \code{dynhr_fig_height} attribute.
#' @return A list of ggplots; length 1 (holding \code{p} itself) when the plot
#'   was not split.
#' @noRd
.paginate_plot <- function(p, max_rows = 12L, max_y = 60L,
                           default = 6.5, maximum = 40, max_pages = 8L) {
  one <- list(p)
  if (!inherits(p, "gg") || !requireNamespace("ggplot2", quietly = TRUE))
    return(one)
  dat <- p$data
  if (!is.data.frame(dat) || nrow(dat) == 0L) return(one)
  b   <- ggplot2::ggplot_build(p)
  lay <- b$layout$layout
  if (!is.data.frame(lay) || nrow(lay) == 0L) return(one)
  if (length(max_pages) != 1L || !is.finite(max_pages) || max_pages < 1L)
    return(one)
  if (nrow(lay) > 1L)
    return(.paginate_facets(p, dat, lay, max_rows, default, maximum, max_pages))
  .paginate_discrete(p, dat, b, max_y, default, maximum, max_pages)
}


#' Facet case of .paginate_plot()
#' @noRd
.paginate_facets <- function(p, dat, lay, max_rows, default, maximum,
                             max_pages = 8L) {
  one <- list(p)
  if (length(max_rows) != 1L || !is.finite(max_rows) || max_rows < 1L) return(one)
  fv <- .facet_var_names(p)
  if (length(fv) == 0L) return(one)
  if (!all(fv %in% names(dat)) || !all(fv %in% names(lay))) return(one)
  drop <- p$facet$params$drop
  if (!is.null(drop) && !isTRUE(drop)) return(one)
  if (is.null(lay$ROW) || length(lay$ROW) == 0L) return(one)
  rows <- sort(unique(lay$ROW))
  if (length(rows) <= max_rows) return(one)
  if (length(rows) > max_rows * max_pages) return(one)
  grp     <- split(rows, ceiling(seq_along(rows) / max_rows))
  key_dat <- .paginate_key(dat[, fv, drop = FALSE])
  pages   <- lapply(grp, function(rg) {
    keys <- unique(.paginate_key(lay[lay$ROW %in% rg, fv, drop = FALSE]))
    .paginate_page(p, dat[key_dat %in% keys, , drop = FALSE], fv, keys)
  })
  .paginate_finish(p, pages, default, maximum)
}


#' Discrete-axis case of .paginate_plot()
#' @noRd
.paginate_discrete <- function(p, dat, b, max_y, default, maximum,
                               max_pages = 8L) {
  one <- list(p)
  if (length(max_y) != 1L || !is.finite(max_y) || max_y < 1L) return(one)
  ax   <- "y"
  sc   <- b$layout$panel_scales_y
  disc <- length(sc) >= 1L && isTRUE(sc[[1L]]$is_discrete())
  if (!disc) {
    ax   <- "x"
    sc   <- b$layout$panel_scales_x
    disc <- length(sc) >= 1L && isTRUE(sc[[1L]]$is_discrete())
  }
  if (!disc) return(one)
  v <- .aes_var_name(p$mapping[[ax]])
  if (is.na(v) && length(p$layers) > 0L) {
    for (l in p$layers) {
      cand <- .aes_var_name(l$mapping[[ax]])
      if (!is.na(cand)) { v <- cand; break }
    }
  }
  if (is.na(v) || !(v %in% names(dat))) return(one)
  col <- dat[[v]]
  lv  <- if (is.factor(col)) levels(droplevels(col))
         else sort(unique(as.character(col)))
  if (length(lv) <= max_y) return(one)
  if (length(lv) > max_y * max_pages) return(one)
  grp   <- split(lv, ceiling(seq_along(lv) / max_y))
  chr   <- as.character(col)
  pages <- lapply(grp, function(keys)
    .paginate_page(p, dat[chr %in% keys, , drop = FALSE], v, keys))
  .paginate_finish(p, pages, default, maximum)
}
