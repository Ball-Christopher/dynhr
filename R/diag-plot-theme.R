## R/diag-plot-theme.R
## --------------------------------------------------------------------------
## dynhr plot theme + Paul Tol colour palettes.
##
## Primary colour brand: Paul Tol "vibrant" palette in the order
##   blue -> teal -> magenta -> orange -> cyan -> red, with grey reserved for
##   missing / invalid data.  Area fills use Paul Tol "light" in the same
##   conceptual order.  Diverging plots use Paul Tol "sunset"; sequential use
##   cividis (viridis option "E").
##
## The visual style follows the RBNZ v3 ggplot tools (clean minimal panels,
## bottom legend, horizontal x-axis line, no y axis line, dotted horizontal
## gridlines) but with Source Sans 3 / Source Serif Pro fonts and the Tol
## palettes above.
##
## Backwards compatibility: the old names dynhr_palette,
## theme_dynhr_diagnostic(), scale_colour_dynhr(), scale_fill_dynhr() remain
## exported and resolve to the new vibrant/light scales and theme.
## --------------------------------------------------------------------------


# ============================================================================
# Plot-provenance infrastructure
# ============================================================================

#' Build a provenance descriptor for diagnostic plots
#'
#' Pass the returned object as \code{meta} to any diagnostic function that
#' supports it.  When \code{show_caption = TRUE} (the default), every plot
#' produced by the diagnostic will carry a footer caption showing the model
#' name, an 8-character data fingerprint, and the run date.
#'
#' @param model_name  Short human-readable model identifier (e.g. "03j_v4").
#' @param data        Observation matrix used in estimation.  Its fingerprint
#'   is embedded in the caption so plots remain traceable after the data file
#'   is updated.
#' @param show_caption Logical -- append the provenance caption to every plot
#'   (default \code{TRUE}).
#' @return A \code{dynhr_diag_meta} list.
#' @export
diag_meta <- function(model_name = NULL, data = NULL, show_caption = TRUE) {
  structure(
    list(
      model_name   = model_name,
      data_hash    = if (!is.null(data)) .short_hash(data) else NULL,
      show_caption = isTRUE(show_caption),
      timestamp    = format(Sys.Date(), "%Y-%m-%d")
    ),
    class = "dynhr_diag_meta"
  )
}


#' 8-character deterministic fingerprint (no external dependencies)
#'
#' @param x Any R object.
#' @return 8-character lowercase hexadecimal string.
#' @noRd
.short_hash <- function(x) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(substr(digest::digest(x, algo = "md5"), 1L, 8L))
  }
    vals <- as.numeric(x)
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) return("00000000")
    vals <- head(vals, 200L)
    idx  <- seq_along(vals) - 1L
    h    <- sum(vals * (31^(idx %% 15L))) * 1e4
    sprintf("%08x", bitwAnd(as.integer(abs(h)) %% (2L^30L - 1L), 0x3FFFFFFF))
}


#' Attach a provenance caption to a ggplot
#'
#' @param plot  A ggplot object.
#' @param meta  A \code{dynhr_diag_meta} from \code{diag_meta()}, or NULL.
#' @return The plot with \code{labs(caption = ...)} appended.
#' @noRd
.apply_meta <- function(plot, meta) {
  if (is.null(meta) || !isTRUE(meta$show_caption)) return(plot)
  if (!requireNamespace("ggplot2", quietly = TRUE))  return(plot)
  parts <- character(0)
  if (!is.null(meta$model_name)) parts <- c(parts, meta$model_name)
  if (!is.null(meta$data_hash))  parts <- c(parts, sprintf("data[%s]", meta$data_hash))
  if (!is.null(meta$timestamp))  parts <- c(parts, meta$timestamp)
  if (length(parts) == 0L) return(plot)
  caption <- paste(parts, collapse = " * ")
  plot + ggplot2::labs(caption = caption)
}

## Lazy ggplot2 namespace accessor â€” use inside functions, not at source time.
.ensure_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for diagnostic plots. ",
         "Install it with install.packages('ggplot2').")
  }
  getNamespace("ggplot2")
}


# ============================================================================
# Paul Tol colour palettes
# ============================================================================
# Reference: Paul Tol, "Colour Schemes" technical note v3.2 (SRON).
# Vibrant + light are colour-blind-friendly qualitative palettes; sunset is the
# diverging scheme; cividis is the perceptually-uniform sequential scheme used
# via viridis::scale_*_viridis_*(option = "E").

#' Paul Tol vibrant palette, named (qualitative, high-contrast)
#' @export
tol_vibrant <- c(
  blue    = "#0077BB",
  cyan    = "#33BBEE",
  teal    = "#009988",
  orange  = "#EE7733",
  red     = "#CC3311",
  magenta = "#EE3377",
  grey    = "#BBBBBB"
)

#' Paul Tol light palette, named (qualitative, suitable for area fills)
#' @export
tol_light <- c(
  light_blue   = "#77AADD",
  orange       = "#EE8866",
  light_yellow = "#EEDD88",
  pink         = "#FFAABB",
  light_cyan   = "#99DDFF",
  mint         = "#44BB99",
  pear         = "#BBCC33",
  olive        = "#AAAA00",
  pale_grey    = "#DDDDDD"
)

#' Paul Tol sunset diverging palette (11 stops, low -> high)
#' @export
tol_sunset <- c(
  "#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1", "#C2E4EF",
  "#EAECCC",
  "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"
)

# ----------------------------------------------------------------------------
# dynhr line/fill palettes in the user-specified order
#   line (vibrant): blue, teal, magenta, orange, cyan, red, (grey = NA)
#   fill (light) :  light_blue, mint, pink, orange, light_cyan, light_yellow,
#                   (pale_grey = NA)
# "Light fills in roughly the vibrant order, where matches exist": for red
# there is no light-red in Tol's light palette so we substitute light_yellow
# as a warm contrasting fill.
# ----------------------------------------------------------------------------

#' Default dynhr line colour palette (Tol vibrant in dynhr order)
#' @export
dynhr_palette_vibrant <- unname(tol_vibrant[c("blue", "teal", "magenta",
                                              "orange", "cyan", "red")])

#' Default dynhr fill palette (Tol light, aligned with vibrant order)
#' @export
dynhr_palette_light <- unname(c(
  tol_light["light_blue"],
  tol_light["mint"],
  tol_light["pink"],
  tol_light["orange"],
  tol_light["light_cyan"],
  tol_light["light_yellow"]
))

#' Colour used for missing / invalid / out-of-sample data
#' @export
dynhr_na_colour <- unname(tol_vibrant["grey"])

#' Colour used for missing / invalid fill regions
#' @export
dynhr_na_fill <- unname(tol_light["pale_grey"])

#' Primary brand colour for headings and accents (vibrant blue)
#' @export
dynhr_primary_colour <- unname(tol_vibrant["blue"])


# ============================================================================
# Backwards-compatible palette names used throughout R/diag-*.R
# ============================================================================

#' Legacy colour list (kept for callers; values now Paul Tol-flavoured)
#' @export
dynhr_colours <- list(
  dark_blue   = "#003B73",
  mid_blue    = dynhr_primary_colour,            # vibrant blue
  light_blue  = unname(tol_light["light_blue"]),
  teal        = unname(tol_vibrant["teal"]),
  orange      = unname(tol_vibrant["orange"]),
  red         = unname(tol_vibrant["red"]),
  green       = "#117733",                       # Tol muted green
  grey        = dynhr_na_colour,
  light_grey  = dynhr_na_fill,
  white       = "#FFFFFF"
)

#' Legacy discrete palette (Tol vibrant + filler), length >= 9 so existing
#' diagnostics that index into it cannot trip an out-of-bounds error
#' @export
dynhr_palette <- c(
  dynhr_palette_vibrant,                   # 6 vibrant
  dynhr_na_colour,                         # 7 grey
  dynhr_palette_light[c(2, 3)]             # 8-9 mint, pink for variety
)


# ============================================================================
# Theme
# ============================================================================

#' dynhr ggplot2 theme (clean minimal panels, Tol-coloured, Source-Sans fonts)
#'
#' Clean minimal panel with horizontal y gridlines, x axis line and ticks
#' only, bottom legend, and large readable Source Sans 3 typography.
#' If \code{base_family} is not installed on the system, falls back to
#' "sans" silently.
#'
#' @param base_size   Numeric base font size (default 14).
#' @param base_family Base sans-serif family (default "Source Sans 3").
#' @export
theme_dynhr <- function(base_size = 14,
                        base_family = "Source Sans 3") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for theme_dynhr(). ",
         "Install it with install.packages('ggplot2').")
  }

  fam <- .dynhr_resolve_family(base_family)

  text_colour    <- "#444444"
  heading_colour <- "#1A1A1A"
  grid_colour    <- "#DDDDDD"
  accent_colour  <- dynhr_primary_colour   # vibrant blue

  `%+replace%` <- ggplot2::`%+replace%`
  ggplot2::theme_minimal(base_size = base_size, base_family = fam) %+replace%
    ggplot2::theme(
      text             = ggplot2::element_text(family = fam, colour = text_colour),

      plot.title       = ggplot2::element_text(size = ggplot2::rel(1.25), face = "bold",
                                               colour = heading_colour, hjust = 0,
                                               margin = ggplot2::margin(b = 6)),
      plot.title.position = "plot",
      plot.subtitle    = ggplot2::element_text(size = ggplot2::rel(1.00), hjust = 0,
                                               colour = text_colour,
                                               margin = ggplot2::margin(b = 10)),
      plot.caption     = ggplot2::element_text(size = ggplot2::rel(0.80), hjust = 0,
                                               colour = "#777777",
                                               margin = ggplot2::margin(t = 8)),
      plot.caption.position = "plot",
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin      = ggplot2::margin(10, 16, 10, 10),

      panel.background   = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border       = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = grid_colour, linewidth = 0.30),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.spacing.x    = ggplot2::unit(2.0, "lines"),
      panel.spacing.y    = ggplot2::unit(1.8, "lines"),

      axis.line.x  = ggplot2::element_line(colour = "#222222", linewidth = 0.5),
      axis.line.y  = ggplot2::element_blank(),
      axis.title   = ggplot2::element_text(size = ggplot2::rel(0.95), face = "bold",
                                           colour = heading_colour),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8), hjust = 1),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8), angle = 90),
      axis.text    = ggplot2::element_text(size = ggplot2::rel(0.90), colour = text_colour),
      axis.text.x  = ggplot2::element_text(margin = ggplot2::margin(t = 4)),
      axis.text.y  = ggplot2::element_text(margin = ggplot2::margin(r = 4)),
      axis.ticks.x = ggplot2::element_line(colour = "#222222", linewidth = 0.3),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.length.x = ggplot2::unit(3, "pt"),

      legend.position      = "bottom",
      legend.justification = "center",
      legend.direction     = "horizontal",
      legend.title         = ggplot2::element_text(size = ggplot2::rel(0.90), face = "bold",
                                                   colour = heading_colour),
      legend.text          = ggplot2::element_text(size = ggplot2::rel(0.90), colour = text_colour),
      legend.key           = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.background    = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.margin        = ggplot2::margin(t = 10),

      strip.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      # family = "sans": the display font's bold "fi"/"ig" ligatures render as
      # "fl"/"lg" in strips ("phi" -> "phl", "sig" -> "slg")
      strip.text       = ggplot2::element_text(size = ggplot2::rel(0.95), face = "bold",
                                               family = "sans",
                                               colour = heading_colour,
                                               margin = ggplot2::margin(t = 6, b = 6))
    )
}

#' Legacy alias retained so existing diagnostics keep rendering
#' @param base_size Numeric base font size (default 14). Forwarded to
#'   \code{\link{theme_dynhr}}.
#' @export
theme_dynhr_diagnostic <- function(base_size = 14) theme_dynhr(base_size = base_size)

#' Compact dynhr theme for many-facet plots
#'
#' Builds on \code{theme_dynhr(base_size)} but maximises panel area: tighter
#' panel spacing, no y-axis title/text/ticks/line, and a smaller strip label.
#' The x-axis line and ticks are retained so the time axis remains readable.
#' Use together with \code{geom_dynhr_zero()} to add a visible zero baseline
#' when the y-axis labels are suppressed.
#'
#' @param base_size Numeric base font size (default 12).
#' @return A ggplot2 theme object.
#' @export
theme_dynhr_compact <- function(base_size = 12) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for theme_dynhr_compact(). ",
         "Install it with install.packages('ggplot2').")
  }
  theme_dynhr(base_size = base_size) +
    ggplot2::theme(
      panel.spacing.x  = ggplot2::unit(0.3, "lines"),
      panel.spacing.y  = ggplot2::unit(0.3, "lines"),
      axis.title.y     = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank(),
      axis.line.y      = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(
        size   = ggplot2::rel(0.8),
        face   = "bold",
        family = "sans",
        colour = "#1A1A1A",
        margin = ggplot2::margin(t = 3, b = 3)
      ),
      plot.margin      = ggplot2::margin(6, 10, 6, 6)
    )
}

#' Zero-baseline geom for compact facet plots
#'
#' Returns a \code{geom_hline} at \code{yintercept = 0} styled with the dynhr
#' grey colour and a thin linewidth.  Drop this into any plot that suppresses
#' y-axis labels (e.g. when using \code{theme_dynhr_compact()}) so viewers can
#' still locate zero without a printed tick.
#'
#' @return A ggplot2 layer.
#' @export
geom_dynhr_zero <- function() {
  ggplot2::geom_hline(yintercept = 0,
                      colour     = dynhr_colours$grey,
                      linewidth  = 0.3)
}

# Resolve a font family to "sans" silently if not installed.
.dynhr_resolve_family <- function(family) {
  if (!is.character(family) || length(family) != 1L || !nzchar(family)) return("sans")
  ok <- FALSE
  if (requireNamespace("systemfonts", quietly = TRUE)) {
    sf <- systemfonts::system_fonts()
    if (!is.null(sf) && any(grepl(family, sf$family, ignore.case = TRUE, fixed = FALSE))) {
      ok <- TRUE
    }
  } else {
    # Best-effort: trust caller; ggplot will fall back at render time.
    ok <- TRUE
  }
  if (ok) family else "sans"
}


# ============================================================================
# Colour and fill scales
# ============================================================================

## Internal: build a palette function that returns N colours, interpolating
## via grDevices::colorRampPalette when N exceeds the base palette length.
.dynhr_palette_fun <- function(base) {
  force(base)
  function(n) {
    if (n <= length(base)) base[seq_len(n)]
    else grDevices::colorRampPalette(base, space = "Lab")(n)
  }
}

#' Discrete line/colour scale using Paul Tol vibrant in dynhr order
#'
#' If there are more levels than the 6 vibrant colours, the palette is
#' interpolated in Lab space.
#' @param ... Passed to \code{ggplot2::discrete_scale}.
#' @param na.value Colour for NA values (default: Tol grey).
#' @export
scale_colour_dynhr_vibrant <- function(..., na.value = dynhr_na_colour) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  ggplot2::discrete_scale(
    aesthetics = "colour",
    palette    = .dynhr_palette_fun(dynhr_palette_vibrant),
    na.value   = na.value,
    ...
  )
}

#' @rdname scale_colour_dynhr_vibrant
#' @export
scale_color_dynhr_vibrant <- scale_colour_dynhr_vibrant

#' Discrete fill scale using Paul Tol light in dynhr order
#'
#' If there are more levels than the 6 light colours, the palette is
#' interpolated in Lab space.
#' @param ... Passed to \code{ggplot2::discrete_scale}.
#' @param na.value Fill for NA values (default: Tol pale grey).
#' @export
scale_fill_dynhr_light <- function(..., na.value = dynhr_na_fill) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  ggplot2::discrete_scale(
    aesthetics = "fill",
    palette    = .dynhr_palette_fun(dynhr_palette_light),
    na.value   = na.value,
    ...
  )
}

#' Diverging colour scale using Paul Tol sunset (continuous)
#' @param midpoint Numeric midpoint (default 0).
#' @param ... Passed to \code{ggplot2::scale_colour_gradientn}.
#' @param na.value Colour for NA values.
#' @export
scale_colour_dynhr_sunset <- function(..., midpoint = 0, na.value = dynhr_na_colour) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  ggplot2::scale_colour_gradientn(colours = tol_sunset, na.value = na.value, ...)
}

#' @rdname scale_colour_dynhr_sunset
#' @export
scale_color_dynhr_sunset <- scale_colour_dynhr_sunset

#' Diverging fill scale using Paul Tol sunset (continuous)
#' @param midpoint Numeric midpoint of the diverging scale (default 0).
#' @param ... Additional arguments passed to
#'   \code{ggplot2::scale_fill_gradientn}.
#' @param na.value Fill colour for missing values (default: dynhr NA fill).
#' @export
scale_fill_dynhr_sunset <- function(..., midpoint = 0, na.value = dynhr_na_fill) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  ggplot2::scale_fill_gradientn(colours = tol_sunset, na.value = na.value, ...)
}

#' Sequential cividis colour scale (viridis option "E")
#' @param ... Additional arguments passed to
#'   \code{ggplot2::scale_colour_viridis_c}.
#' @param na.value Colour for NA values (default: dynhr NA colour).
#' @export
scale_colour_dynhr_cividis <- function(..., na.value = dynhr_na_colour) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  ggplot2::scale_colour_viridis_c(option = "E", na.value = na.value, ...)
}

#' @rdname scale_colour_dynhr_cividis
#' @export
scale_color_dynhr_cividis <- scale_colour_dynhr_cividis

#' Sequential cividis fill scale (viridis option "E")
#' @param ... Additional arguments passed to
#'   \code{ggplot2::scale_fill_viridis_c}.
#' @param na.value Fill colour for NA values (default: dynhr NA fill).
#' @export
scale_fill_dynhr_cividis <- function(..., na.value = dynhr_na_fill) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  ggplot2::scale_fill_viridis_c(option = "E", na.value = na.value, ...)
}


# ----------------------------------------------------------------------------
# Legacy aliases (existing diagnostics call these names directly)
# ----------------------------------------------------------------------------

#' Legacy: colour scale -- now Tol vibrant in dynhr order
#' @param ... Arguments passed to \code{\link{scale_colour_dynhr_vibrant}}.
#' @export
scale_colour_dynhr <- function(...) scale_colour_dynhr_vibrant(...)

#' Legacy: fill scale -- now Tol light in dynhr order
#' @param ... Arguments passed to \code{\link{scale_fill_dynhr_light}}.
#' @export
scale_fill_dynhr <- function(...) scale_fill_dynhr_light(...)
