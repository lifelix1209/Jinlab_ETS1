# ============================================================
# Shared figure font setup (Arial via showtext)
#
# cairo_pdf rounds glyph advances to whole pixels, which gives uneven
# letter spacing at small sizes. showtext draws glyphs with exact
# metrics from the font file, so spacing is consistent. Text is
# embedded as glyph outlines (not editable as text).
# ============================================================

FONT_FAMILY <- "ArialFig"

register_arial <- function(family = FONT_FAMILY, dpi = 300) {
  font_dirs <- c("/System/Library/Fonts/Supplemental", "/Library/Fonts",
                 "C:/Windows/Fonts", "/usr/share/fonts/truetype/msttcorefonts")
  files <- c(regular = "Arial.ttf", bold = "Arial Bold.ttf",
             italic = "Arial Italic.ttf", bolditalic = "Arial Bold Italic.ttf")
  find_font <- function(f) {
    hits <- file.path(font_dirs, f)
    hits[file.exists(hits)][1]
  }
  paths <- vapply(files, find_font, character(1))
  if (anyNA(paths)) {
    stop("Arial font files not found: ", paste(files[is.na(paths)], collapse = ", "))
  }
  sysfonts::font_add(family, regular = paths[["regular"]], bold = paths[["bold"]],
                     italic = paths[["italic"]], bolditalic = paths[["bolditalic"]])
  showtext::showtext_opts(dpi = dpi)
  showtext::showtext_auto(enable = TRUE)
  invisible(family)
}

save_figure_pdf <- function(plot, path, width, height) {
  ggplot2::ggsave(path, plot, width = width, height = height, device = grDevices::pdf)
}
