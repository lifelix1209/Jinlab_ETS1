# ============================================================
# Panel builders for the RNA-ATAC figure (rna_atac_run.R).
# ============================================================

BROAD_ORDER <- c("Endothelial", "Hepatocyte", "Kupffer cell", "Stellate cell", "Immune cell", "Blood cell",
                 "Cholangiocyte")
BROAD_COLS <- c(Endothelial = "#2E7D5B", Hepatocyte = "#AEBCCD", "Kupffer cell" = "#C4A98C",
                "Stellate cell" = "#D8CB9A", "Immune cell" = "#B3A6CB", "Blood cell" = "#D9A29F",
                Cholangiocyte = "#93C4C7")
track_cols <- function(ligand) {
  stats::setNames(c("#1B5E20", "#81B29A", "#7F8FA6", "#A68A6D", "#B8A15E"),
                  c(paste0(ligand, "+ EC"), paste0(ligand, "− EC"), TRACK_NON_EC))
}
LINK_COL <- "#6A3D9A"
ETS_COL  <- "#C0392B"

# ---------------- a: ATAC UMAP ----------------
plot_atac_umap <- function(um) {
  set.seed(SEED)
  um <- um[sample(nrow(um)), ]
  um <- um[order(um$broad == "Endothelial"), ]
  um$broad <- factor(um$broad, levels = BROAD_ORDER)
  lab <- stats::aggregate(cbind(UMAP1, UMAP2) ~ broad, um[um$broad %in% BROAD_ORDER[1:4], ], stats::median)
  n_ec <- sum(um$broad == "Endothelial")
  lab$label <- ifelse(lab$broad == "Endothelial", paste0("Endothelial\n(n = ", fmt_n(n_ec), ")"),
                      as.character(lab$broad))
  ggplot2::ggplot(um, ggplot2::aes(UMAP1, UMAP2)) +
    ggrastr::rasterise(ggplot2::geom_point(ggplot2::aes(colour = broad), size = 0.2, stroke = 0), dpi = 600,
                       dev = "ragg") +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = label), size = pt(PT_SMALL), family = FONT,
                             fontface = ifelse(lab$broad == "Endothelial", "bold", "plain"), seed = SEED,
                             box.padding = 0.2, min.segment.length = 0, segment.size = 0.2) +
    ggplot2::scale_colour_manual(values = BROAD_COLS, guide = "none") +
    ggplot2::labs(x = paste0("ATAC UMAP 1 (", fmt_n(nrow(um)), " nuclei)"), y = "ATAC UMAP 2") +
    theme_umap() + ggplot2::theme(aspect.ratio = 1)
}

# ---------------- b: RNA expression vs gene activity ----------------
plot_rna_vs_ga <- function(dot) {
  d <- dot[dot$celltype %in% BROAD_ORDER & dot$gene %in% DOT_GENES, ]
  d$celltype <- factor(d$celltype, levels = rev(BROAD_ORDER))
  d$gene <- factor(d$gene, levels = DOT_GENES)
  d$modality <- factor(d$modality, levels = c("RNA expression", "Gene activity"))
  ggplot2::ggplot(d, ggplot2::aes(gene, celltype)) +
    ggplot2::geom_point(ggplot2::aes(size = pct, fill = scaled), shape = 21, stroke = 0.2) +
    ggplot2::facet_wrap(~modality, nrow = 1) +
    ggplot2::scale_size_area(max_size = 3, limits = c(0, 100), breaks = c(10, 50, 90), name = "% nuclei") +
    ggplot2::scale_fill_gradient(low = "#F2F2F2", high = "#1B5E20", limits = c(0, 1), breaks = c(0, 1),
                                 name = "Scaled mean") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = grid::unit(1.6, "mm"), barheight = grid::unit(8, "mm"))) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic", angle = 90, hjust = 1, vjust = 0.5),
                   strip.text = ggplot2::element_text(face = "plain"), legend.box.spacing = grid::unit(1, "mm"),
                   panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.25),
                   axis.line = ggplot2::element_blank())
}

# ---------------- c: locus peak specificity heatmap ----------------
peak_heat_long <- function(x) {
  m <- x$cpm
  z <- t(scale(t(m)))
  z[is.na(z)] <- 0
  nets <- table(factor(x$ets$peak, levels = rownames(m)))
  d <- data.frame(locus = x$ligand, peak = rep(rownames(m), times = ncol(m)),
                  group = rep(sub("^Wnt2|^Wnt9b", "Ligand", colnames(m)), each = nrow(m)),
                  cpm = as.vector(m), z = as.vector(z), row.names = NULL)
  d$ets <- as.integer(nets[d$peak])
  d$linked <- d$peak %in% x$links$peak
  d$pos <- as.numeric(sub("^chr[^-]+-(\\d+)-.*$", "\\1", d$peak))
  d
}

plot_peak_heatmap <- function(loci) {
  d <- do.call(rbind, lapply(loci, peak_heat_long))
  d <- d[stats::ave(d$cpm, d$peak, FUN = max) >= 5, ]
  d$group <- factor(d$group, levels = c("Ligand+ EC", "Ligand− EC", TRACK_NON_EC),
                    labels = c("Ligand+ EC", "Ligand− EC", "Hepatocyte", "Kupffer", "Stellate"))
  d$label <- paste0(sub("-.*$", "", d$peak), ":", format(round(d$pos / 1e3), big.mark = ","), " kb",
                    ifelse(d$linked, " †", ""), ifelse(d$ets > 0, paste0(" (", d$ets, " ETS)"), ""))
  d$label <- factor(d$label, levels = unique(d$label[order(d$locus, -d$pos)]))
  d$locus <- factor(d$locus, levels = names(LOCI))
  ggplot2::ggplot(d, ggplot2::aes(group, label, fill = pmax(pmin(z, 2), -2))) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::facet_grid(locus ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"), limits = c(-2, 2),
                                  name = "Accessibility\n(row z-score)") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = grid::unit(1.6, "mm"), barheight = grid::unit(8, "mm"))) +
    ggplot2::scale_x_discrete(expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature() +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
                   axis.text.y = ggplot2::element_text(size = PT_SMALL - 0.5),
                   strip.text.y = ggplot2::element_text(face = "italic", angle = 0),
                   panel.spacing = grid::unit(1, "mm"), legend.box.spacing = grid::unit(1, "mm"))
}

# ---------------- d, e: locus tracks ----------------
track_xscale <- function(region) {
  s <- as.numeric(strsplit(region, "-")[[1]][2:3])
  ggplot2::scale_x_continuous(limits = s, expand = c(0, 0), labels = function(x) paste0(round(x / 1e6, 3)))
}

blank_x <- function() ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                                      axis.line.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank())

coverage_panel <- function(x) {
  cv <- x$coverage
  ymax <- max(cv$signal)
  ncount <- stats::setNames(x$n$Freq, x$n$group)
  strip_lab <- stats::setNames(paste0(levels(cv$group), " (", fmt_n(ncount[levels(cv$group)]), ")"), levels(cv$group))
  ggplot2::ggplot(cv, ggplot2::aes(pos, signal, fill = group)) +
    ggplot2::geom_rect(data = x$peaks[x$peaks$linked, ], ggplot2::aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
                       fill = "#EFE7F6", inherit.aes = FALSE) +
    ggplot2::geom_area(linewidth = 0, outline.type = "full") +
    ggplot2::facet_wrap(~group, ncol = 1, strip.position = "left", labeller = ggplot2::as_labeller(strip_lab)) +
    ggplot2::scale_fill_manual(values = track_cols(x$ligand), guide = "none") +
    ggplot2::scale_y_continuous(breaks = NULL, limits = c(0, ymax), expand = c(0, 0)) +
    track_xscale(x$region) +
    ggplot2::labs(y = NULL, subtitle = paste0("Normalized Tn5 insertions (range 0\u2013", signif(ymax, 2), ")")) +
    theme_nature() + blank_x() +
    ggplot2::theme(strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, face = "plain", size = PT_SMALL),
                   strip.placement = "outside", panel.spacing = grid::unit(0.4, "mm"),
                   plot.margin = ggplot2::margin(4, 2, 1, 2, "mm"),
                   axis.line.y = ggplot2::element_blank(), plot.subtitle = ggplot2::element_text(size = PT_SMALL, hjust = 1),
                   panel.border = ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill = NA, colour = NA))
}

peaks_panel <- function(x) {
  pk <- x$peaks
  pk$type <- ifelse(pk$linked, "Linked to gene", "Peak")
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = pk, ggplot2::aes(xmin = start, xmax = end, ymin = 0.55, ymax = 1, fill = type)) +
    ggplot2::geom_segment(data = x$ets, ggplot2::aes(x = (start + end) / 2, xend = (start + end) / 2, y = 0, yend = 0.42),
                          colour = ETS_COL, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c("Linked to gene" = LINK_COL, Peak = "#8C8C8C"), name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0.2, 0.78), labels = c("ETS motif", "Peaks")) +
    track_xscale(x$region) + ggplot2::labs(y = NULL) +
    theme_nature() + blank_x() +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "none")
}

links_panel <- function(x, tss) {
  lk <- x$links
  if (nrow(lk) == 0) return(ggplot2::ggplot() + ggplot2::theme_void())
  pk_mid <- (as.numeric(sub("^[^-]+-(\\d+)-.*", "\\1", lk$peak)) + as.numeric(sub(".*-(\\d+)$", "\\1", lk$peak))) / 2
  d <- data.frame(x = pk_mid, xend = tss, score = lk$score)
  ggplot2::ggplot(d) +
    ggplot2::geom_curve(ggplot2::aes(x = x, xend = xend, y = 0, yend = 0, colour = score), curvature = 0.35,
                        linewidth = 0.45, ncp = 10) +
    ggplot2::scale_colour_gradient(low = "#CDB8E0", high = LINK_COL, name = "Link score",
                                   breaks = function(l) signif(l, 2)) +
    ggplot2::scale_y_continuous(limits = c(-1, 0.05), breaks = -0.5, labels = "Links") +
    ggplot2::guides(colour = ggplot2::guide_colourbar(barwidth = grid::unit(1.4, "mm"), barheight = grid::unit(6, "mm"))) +
    track_xscale(x$region) + ggplot2::labs(y = NULL) +
    theme_nature() + blank_x() +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "none")
}

gene_panel <- function(x) {
  g <- x$genes[x$genes$type == "exon", ]
  ex <- g
  body <- stats::aggregate(cbind(start, end) ~ gene_name + strand, g, function(v) v[1])
  body$start <- tapply(g$start, g$gene_name, min)[body$gene_name]
  body$end <- tapply(g$end, g$gene_name, max)[body$gene_name]
  s <- as.numeric(strsplit(x$region, "-")[[1]][2:3])
  body$lab_x <- pmin(pmax((pmax(body$start, s[1]) + pmin(body$end, s[2])) / 2, s[1]), s[2])
  body$label <- paste0(body$gene_name, ifelse(body$strand == "-", " ←", " →"))
  ggplot2::ggplot() +
    ggplot2::geom_segment(data = body, ggplot2::aes(x = start, xend = end, y = 0, yend = 0), linewidth = 0.3) +
    ggplot2::geom_rect(data = ex, ggplot2::aes(xmin = start, xmax = end, ymin = -0.3, ymax = 0.3), fill = "#1F3A5F") +
    ggplot2::geom_text(data = body, ggplot2::aes(x = lab_x, y = 0.95, label = label), size = pt(PT_SMALL),
                       family = FONT, fontface = "italic") +
    ggplot2::scale_y_continuous(limits = c(-0.5, 1.35), breaks = NULL) +
    ggplot2::coord_cartesian(clip = "on") +
    track_xscale(x$region) +
    ggplot2::labs(x = paste0(sub("-.*$", "", x$region), " position (Mb)"), y = NULL) +
    theme_nature() + ggplot2::theme(axis.line.y = ggplot2::element_blank())
}

ligand_tss <- function(x) {
  g <- x$genes[x$genes$gene_name == x$ligand, ]
  if (unique(g$strand)[1] == "-") max(g$end) else min(g$start)
}

plot_locus <- function(x) {
  p <- patchwork::wrap_plots(coverage_panel(x), peaks_panel(x), links_panel(x, ligand_tss(x)), gene_panel(x),
                             ncol = 1, heights = c(5, 0.7, 0.8, 0.8))
  patchwork::wrap_elements(full = p)
}
