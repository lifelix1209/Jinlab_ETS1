# ============================================================
# Figure 3a-c panel builders (Nature style).
# Requires v2_clean/code/shared/nature_theme.R and fig3_analysis.R.
# ============================================================

fmt_n <- function(x) format(x, big.mark = ",", trim = TRUE)

# ---------------- a: workflow schematic ----------------
schematic_boxes <- function(n_nuclei) {
  data.frame(
    id = c("nuc", "mult", "rna", "atac", "int"),
    xmin = c(60, 86, 120, 120, 146), xmax = c(79, 112, 138, 138, 178),
    ymin = c(12, 11, 20, 4, 9), ymax = c(22, 23, 28, 12, 25),
    label = c("Nuclei<br>isolation", "10x Chromium<br>Multiome<br>(ATAC + GEX)",
              "snRNA-seq", "snATAC-seq",
              paste0("QC, doublet removal<br>WNN integration<br><b>", fmt_n(n_nuclei), " nuclei</b>"))
  )
}

schematic_arrows <- function() {
  data.frame(x = c(55, 79, 112, 112, 138, 138), xend = c(59.5, 85.5, 119.5, 119.5, 145.5, 145.5),
             y = c(17, 17, 17, 17, 24, 8), yend = c(17, 17, 24, 8, 19, 15))
}

genotype_rows <- function(n_mice) {
  data.frame(
    y = c(23, 11), genotype = factor(names(GROUP_COLS), levels = names(GROUP_COLS)),
    label = c(
      paste0("<b>Control</b> (n = ", n_mice[["Control"]], " mice)<br>",
             "<i>Axin2-rtTA</i>; <i>mTmG</i>; <i>β-cat</i>", sup("f/f")),
      paste0("<b>β-cat cKO</b> (n = ", n_mice[["β-cat cKO"]], " mice)<br>",
             "<i>Axin2-rtTA</i>; <i>TetO-Cre</i>; <i>mTmG</i>; <i>β-cat</i>", sup("f/del"))))
}

plot_fig3a <- function(n_nuclei, n_mice) {
  boxes <- schematic_boxes(n_nuclei)
  geno <- genotype_rows(n_mice)
  ggplot2::ggplot() +
    ggplot2::geom_point(data = geno, ggplot2::aes(2, y, colour = genotype), size = 2.2) +
    ggtext::geom_richtext(data = geno, ggplot2::aes(5, y, label = label), hjust = 0, size = pt(PT_TEXT),
                          family = FONT, fill = NA, label.colour = NA, lineheight = 1.15) +
    ggplot2::geom_segment(data = geno, ggplot2::aes(x = 50, xend = 55, y = y, yend = 17), linewidth = 0.3) +
    ggplot2::geom_rect(data = boxes, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = "#F2F4F7", colour = "black", linewidth = 0.25) +
    ggtext::geom_richtext(data = boxes, ggplot2::aes((xmin + xmax) / 2, (ymin + ymax) / 2, label = label),
                          size = pt(PT_TEXT), family = FONT, fill = NA, label.colour = NA, lineheight = 1.1) +
    ggplot2::geom_segment(data = schematic_arrows(), ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                          linewidth = 0.3, arrow = grid::arrow(length = grid::unit(1, "mm"), type = "closed")) +
    ggplot2::scale_colour_manual(values = GROUP_COLS, guide = "none") +
    ggplot2::coord_cartesian(xlim = c(0, 184), ylim = c(1, 31), expand = FALSE, clip = "off") +
    ggplot2::theme_void(base_family = FONT) +
    ggplot2::theme(plot.margin = ggplot2::margin(2, 1, 1, 1, "mm"),
                   plot.tag = ggplot2::element_text(size = PT_LABEL, face = "bold", family = FONT),
                   plot.tag.position = c(0, 1))
}

# ---------------- b: UMAP ----------------
# Curated direct-label positions (UMAP units), placed in empty space next
# to each population of the fixed UMAP; checked by the collision audit.
UMAP_LABELS <- data.frame(
  label = c("Hepatocyte", "PeriportalHep", "Cholangiocyte", "Blood cell", "eGFP+ Hep", "Glul+ Hep",
            "Endothelial", "Stellate cell", "Kupffer cell", "Immune cell"),
  UMAP1 = c(-2.3, -2.5, 1.9, 3.1, 3.0, -0.3, 6.6, 8.3, 9.5, 6.1),
  UMAP2 = c(0.6, 5.65, 5.65, 0.9, -3.0, -5.3, 1.9, 3.4, -3.0, -1.5),
  hjust = c(0.5, 0.5, 0.5, 0, 0, 0, 0.5, 1, 0.5, 1)
)

plot_fig3b_umap <- function(cells) {
  set.seed(42)
  d <- cells[sample(nrow(cells)), ]
  ggplot2::ggplot(d, ggplot2::aes(UMAP1, UMAP2)) +
    ggrastr::rasterise(ggplot2::geom_point(ggplot2::aes(colour = celltype), size = 0.35, stroke = 0),
                       dpi = 600, dev = "ragg") +
    ggplot2::geom_text(data = UMAP_LABELS, ggplot2::aes(label = label, hjust = hjust),
                       size = pt(PT_TEXT), family = FONT) +
    ggplot2::scale_colour_manual(values = CELLTYPE_V1_COLS, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(add = c(0.4, 1.4))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = c(0.4, 0.9))) +
    ggplot2::labs(x = "UMAP 1", y = "UMAP 2") +
    theme_umap()
}

plot_fig3b_split <- function(cells) {
  focus <- c("Glul+ Hep", "eGFP+ Hep")
  d <- transform(cells, focus = ifelse(as.character(celltype) %in% focus, as.character(celltype), "Other"))
  d <- d[order(d$focus != "Other"), ]
  n <- dplyr::count(cells, genotype)
  n$label <- paste0("n = ", fmt_n(n$n))
  ggplot2::ggplot(d, ggplot2::aes(UMAP1, UMAP2)) +
    ggrastr::rasterise(ggplot2::geom_point(ggplot2::aes(colour = focus), size = 0.25, stroke = 0),
                       dpi = 600, dev = "ragg") +
    ggplot2::geom_text(data = n, ggplot2::aes(x = Inf, y = -Inf, label = label), hjust = 1.05, vjust = -0.5,
                       size = pt(PT_SMALL), family = FONT) +
    ggplot2::facet_wrap(~genotype, ncol = 1) +
    ggplot2::scale_colour_manual(values = c(CELLTYPE_V1_COLS[focus], Other = "#D6D6D6"),
                                 breaks = focus, name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 1.6))) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_umap() +
    ggplot2::theme(axis.line = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "plain"),
                   legend.position = "bottom", legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.box.spacing = grid::unit(0.5, "mm"))
}

# ---------------- c: candidate bubble plot ----------------
# Same encoding as the original Fig. 3C: genes on x (sorted by %), mean UMI per nucleus
# on y (zeros included), bubble area = % nuclei with >= 1 UMI, fill = mean
# UMI (sequential ramp of the cKO colour), % printed above each bubble
# because area alone is hard to read.
BUBBLE_MAX_MM <- 6
BUBBLE_PCT_LIMIT <- 60

bubble_label_y <- function(d, mm_per_unit) {
  radius_mm <- 0.5 * BUBBLE_MAX_MM * sqrt(pmin(d$pct, BUBBLE_PCT_LIMIT) / BUBBLE_PCT_LIMIT)
  d$mean_umi + (radius_mm + 0.8) / mm_per_unit
}

FIG3C_STYLE <- list(
  "β-cat cKO" = list(low = "#FCE8D5", high = "#B4561C", title = paste0("β-cat cKO eGFP", sup("+"), " Hep")),
  "Control"   = list(low = "#DCE8F3", high = "#1F4E79", title = paste0("Control Glul", sup("+"), " Hep"))
)

# One genotype group; genes ordered left to right by % nuclei with >= 1 UMI
# (descending). The y axis is shared across groups so the panels compare.
plot_fig3c <- function(group, which = "β-cat cKO", panel_height_mm = 42) {
  st <- FIG3C_STYLE[[which]]
  if (is.null(st)) stop("Unknown group for Fig. 3c: ", which)
  d <- group[group$group == which, ]
  d$gene <- factor(d$gene, levels = d$gene[order(-d$pct, -d$mean_umi)])
  y_top <- max(group$mean_umi) * 1.25
  d$label_y <- bubble_label_y(d, panel_height_mm / y_top)
  d$label <- sprintf("%.1f%%", d$pct)
  ggplot2::ggplot(d, ggplot2::aes(gene, mean_umi)) +
    ggplot2::geom_point(ggplot2::aes(size = pct, fill = mean_umi), shape = 21, colour = "black", stroke = 0.25) +
    ggplot2::geom_text(ggplot2::aes(y = label_y, label = label), vjust = 0, size = pt(PT_SMALL), family = FONT) +
    ggplot2::scale_size_area(max_size = BUBBLE_MAX_MM, limits = c(0, BUBBLE_PCT_LIMIT), breaks = c(10, 30, 50),
                             name = "Nuclei with\n\u2265 1 UMI (%)") +
    ggplot2::scale_fill_gradient(low = st$low, high = st$high, limits = c(0, max(group$mean_umi)),
                                 breaks = c(0, 0.5, 1), name = "Mean UMI\nper nucleus") +
    ggplot2::scale_y_continuous(limits = c(0, y_top), expand = ggplot2::expansion(mult = c(0.03, 0))) +
    ggplot2::guides(size = ggplot2::guide_legend(order = 1, override.aes = list(fill = GROUP_COLS[[which]])),
                    fill = ggplot2::guide_colourbar(order = 2, barwidth = grid::unit(2, "mm"),
                                                    barheight = grid::unit(10, "mm"))) +
    ggplot2::labs(x = NULL, y = "Mean UMI per nucleus", subtitle = st$title) +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic", angle = 45, hjust = 1),
                   legend.position = "right", legend.box.spacing = grid::unit(1, "mm"),
                   legend.spacing.y = grid::unit(2, "mm"),
                   plot.subtitle = ggtext::element_markdown(size = PT_TEXT))
}

# Both groups side by side per gene, genes ordered by cKO % (descending).
# Fill = genotype; bubble area = % nuclei with >= 1 UMI; y = mean UMI.
plot_fig3c_combined <- function(group, panel_height_mm = 42, dodge = 0.8) {
  d <- group[group$group %in% names(GROUP_COLS), ]
  cko <- d[d$group == "β-cat cKO", ]
  d$gene <- factor(d$gene, levels = cko$gene[order(-cko$pct, -cko$mean_umi)])
  d$group <- factor(d$group, levels = names(GROUP_COLS))
  y_top <- max(d$mean_umi) * 1.25
  d$label_y <- bubble_label_y(d, panel_height_mm / y_top)
  d$label <- sprintf("%.1f%%", d$pct)
  pos <- ggplot2::position_dodge(width = dodge)
  ggplot2::ggplot(d, ggplot2::aes(gene, mean_umi, group = group)) +
    ggplot2::geom_point(ggplot2::aes(size = pct, fill = group), shape = 21, colour = "black", stroke = 0.25,
                        position = pos) +
    ggplot2::geom_text(ggplot2::aes(y = label_y, label = label), vjust = 0, size = pt(PT_SMALL - 0.5),
                       family = FONT, position = pos) +
    ggplot2::scale_size_area(max_size = BUBBLE_MAX_MM, limits = c(0, BUBBLE_PCT_LIMIT), breaks = c(10, 30, 50),
                             name = "Nuclei with\n\u2265 1 UMI (%)") +
    ggplot2::scale_fill_manual(values = GROUP_COLS, labels = c(Control = paste0("Control Glul", sup("+"), " Hep"),
                                                               "β-cat cKO" = paste0("β-cat cKO eGFP", sup("+"), " Hep")),
                               name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, y_top), expand = ggplot2::expansion(mult = c(0.03, 0))) +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1, override.aes = list(size = 2.5)),
                    size = ggplot2::guide_legend(order = 2, override.aes = list(fill = "grey70"))) +
    ggplot2::labs(x = NULL, y = "Mean UMI per nucleus") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic", angle = 45, hjust = 1),
                   legend.position = "right", legend.box.spacing = grid::unit(1, "mm"),
                   legend.spacing.y = grid::unit(2, "mm"), legend.text = ggtext::element_markdown(size = PT_SMALL))
}

# ---------------- composite a-c ----------------
compose_fig3 <- function(res, n_mice) {
  a <- plot_fig3a(nrow(res$cells), n_mice)
  b1 <- plot_fig3b_umap(res$cells)
  b2 <- plot_fig3b_split(res$cells)
  c <- plot_fig3c(res$group)
  # Flat design (no nesting) so the alignment gate can measure every panel.
  design <- c(patchwork::area(1, 1, 1, 3), patchwork::area(2, 1), patchwork::area(2, 2), patchwork::area(2, 3))
  patchwork::wrap_plots(a, b1, b2, c) +
    patchwork::plot_layout(design = design, widths = c(62, 28, 70), heights = c(32, 60)) +
    patchwork::plot_annotation(tag_levels = list(c("a", "b", "", "c")))
}
