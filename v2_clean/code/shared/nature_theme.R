# ============================================================
# Nature-style theme, palette and export helpers for Figure 3 / S3.
#
# Final sizes follow Nature: 183 mm double column, text 5-7 pt, panel
# labels 8 pt bold lowercase, Arial, editable text (cairo_pdf, no
# outlined glyphs), axis lines and ticks, no gridlines.
#
# The multi-panel alignment gate (nature-figure skill) runs when the
# environment variable NATURE_FIGURE_DIR points at the skill folder and
# NATURE_QA_PYTHON at a Python interpreter; otherwise it is skipped with
# a message.
# ============================================================

FONT     <- "Arial"
MM       <- 1 / 25.4
PT_TEXT  <- 6     # ordinary text
PT_SMALL <- 5.5   # ticks / dense labels (>= 5 pt floor)
PT_LABEL <- 8     # panel letters

# The alignment gate measures layouts on grDevices::pdf(), which only knows
# PostScript fonts; Helvetica metrics match Arial. Exports use cairo_pdf.
if (!FONT %in% names(grDevices::pdfFonts())) {
  grDevices::pdfFonts(Arial = grDevices::pdfFonts()$Helvetica)
}

# Convert a font size in pt to ggplot2 text size (mm).
pt <- function(x) x / ggplot2::.pt

# Rich-text superscript held at the 5 pt glyph floor (ggtext would
# otherwise shrink <sup> to 80% of the surrounding size).
sup <- function(x, size = 5.2) sprintf("<sup><span style='font-size:%.1fpt'>%s</span></sup>", size, x)

GROUP_COLS <- c(Control = "#3B75AF", "β-cat cKO" = "#E8892B")
NEUTRAL    <- "#BDBDBD"

# Figure annotation = v1 labels with two marker-based corrections
# (00_preprocessing/03_reassess_unknown_immune.R): cluster 7 "Unknown" -> Immune cell, cluster 14
# "Immune cell" -> Cholangiocyte. The two pericentral states carry the
# genotype colours; other populations use muted, distinguishable tones.
CLUSTER_RELABEL <- c("7" = "Immune cell", "14" = "Cholangiocyte")
CELLTYPE_V1_ORDER <- c("Hepatocyte", "Glul+ Hep", "eGFP+ Hep", "PeriportalHep", "Cholangiocyte",
                       "Endothelial", "Kupffer cell", "Stellate cell", "Immune cell", "Blood cell")
CELLTYPE_V1_COLS <- setNames(
  c("#AEBCCD", "#3B75AF", "#E8892B", "#7FA7C9", "#4DA3A8",
    "#5E9E7E", "#A68A6D", "#C2B280", "#8E7CB0", "#C25B56"),
  CELLTYPE_V1_ORDER)

theme_nature <- function(base = PT_TEXT) {
  ggplot2::theme_classic(base_size = base, base_family = FONT) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black"),
      axis.text = ggplot2::element_text(colour = "black", size = PT_SMALL),
      axis.title = ggplot2::element_text(size = base),
      axis.line = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = grid::unit(0.8, "mm"),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(size = base, face = "italic"),
      legend.title = ggplot2::element_text(size = PT_SMALL),
      legend.text = ggplot2::element_text(size = PT_SMALL),
      legend.key.size = grid::unit(2.5, "mm"),
      legend.background = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(2, 2, 2, 2, "mm"),
      plot.tag = ggplot2::element_text(size = PT_LABEL, face = "bold", family = FONT),
      plot.tag.position = c(0, 1)
    )
}

theme_umap <- function() {
  theme_nature() +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.line = ggplot2::element_line(linewidth = 0.25,
                                                     arrow = grid::arrow(length = grid::unit(1.2, "mm"))))
}

# Resolve gtable track sizes to pt: absolute units are converted on a
# device of the final size; "null" units share the remaining space.
resolve_tracks_pt <- function(sizes, total_pt) {
  is_null <- vapply(seq_along(sizes), function(i) grid::unitType(sizes[i]) == "null", logical(1))
  abs_pt <- ifelse(is_null, 0, grid::convertWidth(sizes, "pt", valueOnly = TRUE))
  null_val <- ifelse(is_null, as.numeric(sizes), 0)
  free <- max(total_pt - sum(abs_pt), 0)
  abs_pt + if (sum(null_val) > 0) free * null_val / sum(null_val) else 0
}

# Measure the top-level panel cells of a flat patchwork at final size and
# write the manifest read by the skill's audit_panel_alignment.py.
write_alignment_manifest <- function(plot, path, w_mm, h_mm, exemptions = list()) {
  w_pt <- w_mm * MM * 72
  h_pt <- h_mm * MM * 72
  grDevices::pdf(NULL, width = w_mm * MM, height = h_mm * MM)
  on.exit(grDevices::dev.off(), add = TRUE)
  gt <- patchwork::patchworkGrob(plot)
  rows <- gt$layout[grepl("^panel", gt$layout$name) & gt$layout$name != "panel-area", , drop = FALSE]
  rows <- rows[order(rows$t, rows$l), , drop = FALSE]
  x <- c(0, cumsum(resolve_tracks_pt(gt$widths, w_pt)))
  y <- c(0, cumsum(resolve_tracks_pt(gt$heights, h_pt)))
  panels <- lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    list(id = letters[i], bbox_pt = c(x[r$l], h_pt - y[r$b + 1], x[r$r + 1], h_pt - y[r$t]),
         grid_id = "patchwork-grid-1", row_start = r$t - 1L, row_stop = r$b,
         col_start = r$l - 1L, col_stop = r$r)
  })
  manifest <- list(schema_version = 1L, backend = "r-patchwork",
                   figure = list(width_pt = w_pt, height_pt = h_pt), panels = panels, exemptions = exemptions)
  jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  invisible(path)
}

# Optional blocking alignment gate using the nature-figure skill auditor.
run_alignment_gate <- function(plot, stem, w_mm, h_mm, exemptions = list()) {
  skill <- Sys.getenv("NATURE_FIGURE_DIR")
  py    <- Sys.getenv("NATURE_QA_PYTHON")
  if (!nzchar(skill) || !nzchar(py) || !inherits(plot, "patchwork")) {
    message("Alignment gate skipped for ", basename(stem))
    return(invisible(NULL))
  }
  manifest <- write_alignment_manifest(plot, paste0(stem, ".alignment-layout.json"), w_mm, h_mm, exemptions)
  args <- c(file.path(skill, "scripts", "audit_panel_alignment.py"), manifest,
            "--json-out", paste0(stem, ".alignment.json"), "--tolerance-pt", "1.5",
            "--gutter-tolerance-pt", "1.5", "--strict")
  out <- suppressWarnings(system2(py, shQuote(args), stdout = TRUE, stderr = TRUE))
  message(paste(out, collapse = "\n"))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) stop("Alignment gate failed (exit ", status, ") for ", basename(stem))
  invisible(manifest)
}

# Vector PDF device. cairo_pdf (via pango) rounds glyph advances to whole
# 1/72-inch units, giving uneven letter spacing at 5-7 pt; the macOS
# Quartz device lays text out with CoreText (exact spacing) and embeds
# Arial as editable TrueType. cairo_pdf is only a fallback off macOS.
open_pdf_device <- function(path, w, h) {
  if (isTRUE(capabilities("aqua"))) {
    grDevices::quartz(type = "pdf", file = path, width = w, height = h, family = FONT)
  } else {
    warning("Quartz unavailable; cairo_pdf letter spacing may be uneven.")
    grDevices::cairo_pdf(path, width = w, height = h, family = FONT)
  }
}

# Export editable PDF + SVG, and a raster file (PNG preview for main
# figures, 300-dpi TIFF for supplementary figures).
save_nature <- function(plot, stem, w_mm, h_mm, raster = c("png", "tiff"), dpi = 600, ...) {
  raster <- match.arg(raster)
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  run_alignment_gate(plot, stem, w_mm, h_mm, ...)
  w <- w_mm * MM
  h <- h_mm * MM
  tryCatch({
    open_pdf_device(paste0(stem, ".pdf"), w, h)
    print(plot)
    grDevices::dev.off()
    svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
    print(plot)
    grDevices::dev.off()
    if (raster == "png") {
      ragg::agg_png(paste0(stem, ".png"), width = w, height = h, units = "in", res = dpi)
    } else {
      ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in", res = 300,
                     compression = "lzw")
    }
    print(plot)
    grDevices::dev.off()
  }, error = function(e) {
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
    stop("Export failed for ", stem, ": ", conditionMessage(e))
  })
  invisible(stem)
}
