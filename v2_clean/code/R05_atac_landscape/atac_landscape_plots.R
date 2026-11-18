# ============================================================
# Panel builders for the chromatin-landscape figure (atac_landscape_run.R).
# Uses BROAD_ORDER / BROAD_COLS from R04 rna_atac_plots.R.
# ============================================================

CAT_COLS <- c(Distal = "#5DB761", Exonic = "#86C5E8", Intronic = "#5B2C8F", Promoter = "#F2C14E")
LANDMARKS_OLD <- c("Sds", "Cyp2f2", "Hal", "Pck1", "Cps1", "Gls2", "Hamp", "Igfbp2", "Aldh1b1", "Mfsd2a",
               "Glul", "Cyp2e1", "Axin2", "Lgr5", "Oat", "Cyp1a2", "Lect2", "Tbx3", "Rhbg", "Cyp27a1")
LAYER_GAP <- 0.45

# ---------------- a: transcriptome <-> epigenome linked UMAPs ----------------
plot_linked_umaps <- function(u) {
  set.seed(SEED)
  u <- u[u$celltype %in% BROAD_ORDER, ]
  u$celltype <- factor(u$celltype, levels = BROAD_ORDER)
  top <- 1 + LAYER_GAP
  pts <- rbind(data.frame(x = u$rna_x, y = u$rna_y + top, celltype = u$celltype),
               data.frame(x = u$atac_x, y = u$atac_y, celltype = u$celltype))
  pts <- pts[sample(nrow(pts)), ]
  lk <- u[sample(nrow(u), min(N_LINKS, nrow(u))), ]
  lab <- rbind(
    transform(stats::aggregate(cbind(x = rna_x, y = rna_y) ~ celltype, u, stats::median), y = y + top),
    stats::aggregate(cbind(x = atac_x, y = atac_y) ~ celltype, u, stats::median))
  lab <- lab[lab$celltype %in% BROAD_ORDER[1:4], ]
  ggplot2::ggplot() +
    ggrastr::rasterise(ggplot2::geom_segment(data = lk, ggplot2::aes(x = atac_x, y = atac_y, xend = rna_x,
                                                                      yend = rna_y + top, colour = celltype),
                                             linewidth = 0.08, alpha = 0.08), dpi = 600, dev = "ragg") +
    ggrastr::rasterise(ggplot2::geom_point(data = pts, ggplot2::aes(x, y, colour = celltype), size = 0.15,
                                           stroke = 0), dpi = 600, dev = "ragg") +
    ggrepel::geom_text_repel(data = lab, ggplot2::aes(x, y, label = celltype), size = pt(PT_SMALL), family = FONT,
                             seed = SEED, box.padding = 0.15, min.segment.length = 0, segment.size = 0.2) +
    ggplot2::geom_hline(yintercept = 1 + LAYER_GAP / 2, linewidth = 0.25, linetype = "22", colour = "grey50") +
    ggplot2::annotate("text", x = -0.06, y = c(top + 0.5, 0.5), label = c("Transcriptome", "Epigenome"),
                      angle = 90, size = pt(PT_TEXT), family = FONT) +
    ggplot2::scale_colour_manual(values = BROAD_COLS, name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 1.6, alpha = 1, linewidth = 0), nrow = 3)) +
    ggplot2::coord_cartesian(xlim = c(-0.1, 1.02), clip = "off") +
    ggplot2::labs(subtitle = paste0(fmt_n(nrow(u)), " paired nuclei; ", fmt_n(nrow(lk)), " links shown")) +
    ggplot2::theme_void(base_family = FONT) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = PT_SMALL),
                   legend.key.size = grid::unit(2.2, "mm"), legend.spacing = grid::unit(0, "mm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0),
                   plot.subtitle = ggplot2::element_text(size = PT_TEXT, hjust = 0.5),
                   plot.margin = ggplot2::margin(2, 2, 2, 2, "mm"),
                   plot.tag = ggplot2::element_text(size = PT_LABEL, face = "bold", family = FONT))
}

composition_df <- function(comp) {
  lv <- c(intersect(BROAD_ORDER, unique(comp$celltype)), "Union")
  comp$celltype <- factor(comp$celltype, levels = lv)
  comp$category <- factor(comp$category, levels = c("Distal", "Exonic", "Intronic", "Promoter"))
  comp
}

# ---------------- b: peaks per cell type (vertical, Fig. ref 4B) ----------------
plot_peaks_vertical <- function(comp) {
  d <- composition_df(comp)
  ggplot2::ggplot(d, ggplot2::aes(celltype, n / 1e3, fill = category)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_fill_manual(values = CAT_COLS, name = NULL) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = NULL, y = "Peaks (thousands)") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
                   legend.position = c(0.02, 1), legend.justification = c(0, 1),
                   legend.key.size = grid::unit(2, "mm"))
}

# ---------------- c: peaks per cell type (horizontal with N, Fig. ref 5C) ----------------
plot_peaks_horizontal <- function(comp) {
  d <- composition_df(comp)
  d$celltype <- factor(d$celltype, levels = rev(levels(d$celltype)))
  tot <- d[!duplicated(d$celltype), c("celltype", "total")]
  tot$label <- paste0("N = ", sprintf("%.1fk", tot$total / 1e3))
  ggplot2::ggplot(d, ggplot2::aes(n / 1e3, celltype, fill = category)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(data = tot, ggplot2::aes(x = total / 1e3, y = celltype, label = label), hjust = -0.1,
                       size = pt(PT_SMALL), family = FONT, inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = CAT_COLS, name = NULL, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.25)), position = "top") +
    ggplot2::labs(x = "Number of peaks (thousands)", y = NULL) +
    theme_nature()
}

# ---------------- d: marker-locus coverage matrix (Fig. ref 4E) ----------------
plot_marker_tracks <- function(tr) {
  cov <- tr$coverage[tr$coverage$gene %in% TRACK_GENES, ]
  tr$models <- tr$models[tr$models$gene %in% TRACK_GENES, ]
  types <- intersect(BROAD_ORDER, unique(cov$group))
  genes <- unique(cov$gene)
  cov$group <- factor(cov$group, levels = types)
  cov$gene <- factor(cov$gene, levels = genes)
  cov$scaled <- stats::ave(cov$signal, cov$gene, FUN = function(v) v / max(v))
  rng <- stats::aggregate(signal ~ gene, cov, max)
  ncell <- stats::setNames(tr$n$Freq, tr$n$celltype)
  strip_y <- stats::setNames(paste0(types, " (", fmt_n(ncell[types]), ")"), types)
  tracks <- ggplot2::ggplot(cov, ggplot2::aes(pos, scaled, fill = group)) +
    ggplot2::geom_area(outline.type = "full", linewidth = 0) +
    ggplot2::facet_grid(group ~ gene, scales = "free_x", switch = "y",
                        labeller = ggplot2::labeller(group = strip_y)) +
    ggplot2::scale_fill_manual(values = BROAD_COLS, guide = "none") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = NULL, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   axis.line.x = ggplot2::element_line(linewidth = 0.2, colour = "grey60"),
                   axis.line.y = ggplot2::element_blank(),
                   strip.text.x = ggplot2::element_text(face = "italic"),
                   strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, face = "plain", size = PT_SMALL),
                   strip.placement = "outside", panel.spacing.y = grid::unit(0.3, "mm"),
                   panel.spacing.x = grid::unit(1.5, "mm"),
                   panel.border = ggplot2::element_blank())
  mod <- tr$models
  mod$gene <- factor(mod$gene, levels = genes)
  span <- stats::aggregate(pos ~ gene, cov, range)
  body <- data.frame(gene = factor(unique(mod$gene), levels = genes),
                     start = tapply(mod$start, mod$gene, min)[unique(as.character(mod$gene))],
                     end = tapply(mod$end, mod$gene, max)[unique(as.character(mod$gene))],
                     strand = tapply(mod$strand, mod$gene, `[`, 1)[unique(as.character(mod$gene))])
  body$arrow <- ifelse(body$strand == "-", "←", "→")
  lim <- data.frame(gene = factor(span$gene, levels = genes), x = c(span$pos[, 1], span$pos[, 2]))
  models <- ggplot2::ggplot() +
    ggplot2::geom_blank(data = lim, ggplot2::aes(x = x, y = 0)) +
    ggplot2::geom_segment(data = body, ggplot2::aes(x = start, xend = end, y = 0, yend = 0), linewidth = 0.3,
                          colour = "#1F3A5F") +
    ggplot2::geom_rect(data = mod, ggplot2::aes(xmin = start, xmax = end, ymin = -0.35, ymax = 0.35), fill = "#1F3A5F") +
    ggplot2::geom_text(data = body, ggplot2::aes(x = (start + end) / 2, y = 0.9, label = arrow), size = pt(PT_SMALL),
                       family = FONT) +
    ggplot2::facet_grid(. ~ gene, scales = "free_x") +
    ggplot2::scale_y_continuous(limits = c(-0.6, 1.3), breaks = NULL) +
    ggplot2::scale_x_continuous(breaks = NULL) +
    ggplot2::labs(x = NULL, y = "Gene") +
    theme_nature() +
    ggplot2::theme(strip.text = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
                   panel.spacing.x = grid::unit(1.5, "mm"), axis.title.y = ggplot2::element_text(size = PT_SMALL))
  patchwork::wrap_elements(full = patchwork::wrap_plots(tracks, models, ncol = 1, heights = c(10, 1)))
}

# ---------------- e, f: zonation heatmaps (Fig. ref 6e) ----------------
# Rows: zonated features ordered PV-high -> CV-high by the bin of their
# smoothed maximum. Values: loess-smoothed (span ZONE_SPAN) log2 CPM along
# the 40 PV->CV bins, row z-scored; reference pseudobulks (control Glul+
# Hep, cKO eGFP+ Hep) are scaled with the same row mean / s.d.
ZONE_SPAN <- 0.35
ZONE_CLIP <- 1.5
REF_W <- 3  # display width of each reference column, in bins
HEAT_COLS <- c("#3B4CC0", "#7B9FF9", "#C0D4F5", "#F2F2F2", "#F5C4AC", "#EE8468", "#B40426")
PP_COL <- "#C9D65B"
PC_COL <- "#2A8C82"
LABEL_GENES <- c("Sds", "Cyp2f2", "Hal", "Pck1", "Cps1", "Gls2", "Glul", "Cyp2e1", "Axin2", "Lgr5", "Oat", "Tbx3")

smooth_rows <- function(m) {
  idx <- seq_len(ncol(m))
  t(apply(m, 1, function(v) stats::predict(stats::loess(v ~ idx, span = ZONE_SPAN, degree = 2))))
}

zone_long <- function(z) {
  feats <- z$features$feature
  m <- log2(z$cpm[feats, , drop = FALSE] + 1)
  sm <- smooth_rows(m)
  mu <- rowMeans(sm)
  s <- pmax(apply(sm, 1, stats::sd), 1e-9)
  zm <- (sm - mu) / s
  zr <- (log2(z$ref[feats, , drop = FALSE] + 1) - mu) / s
  nb <- ncol(m)
  nf <- length(feats)
  heat <- rbind(
    data.frame(block = "zone", x = rep(seq_len(nb), each = nf), row = rep(rev(seq_len(nf)), times = nb),
               z = as.vector(zm)),
    data.frame(block = "ref", x = rep(seq_len(2 * REF_W), each = nf), row = rep(rev(seq_len(nf)), times = 2 * REF_W),
               z = as.vector(zr[, rep(1:2, each = REF_W)])))
  heat$block <- factor(heat$block, levels = c("zone", "ref"))
  heat$zc <- pmax(pmin(heat$z, ZONE_CLIP), -ZONE_CLIP)
  list(heat = heat, nb = nb, nf = nf, feats = feats, rho = z$features$rho)
}

zone_heat <- function(z, what, fill_lab, show_genes = FALSE) {
  L <- zone_long(z)
  nf <- L$nf
  nb <- L$nb
  bar_lo <- nf + nf * 0.025
  bar_hi <- nf + nf * 0.075
  bar <- data.frame(block = factor("zone", levels = c("zone", "ref")), xmin = seq_len(nb) - 0.5,
                    xmax = seq_len(nb) + 0.5, ymin = bar_lo, ymax = bar_hi)
  bar_fill <- grDevices::colorRampPalette(c(PP_COL, PC_COL))(nb)
  bar_txt <- data.frame(block = factor("zone", levels = c("zone", "ref")), x = c(1.2, nb - 0.2),
                        y = (bar_lo + bar_hi) / 2, label = c("PV", "CV"), hjust = c(0, 1),
                        col = c("black", "white"))
  ref_txt <- data.frame(block = factor("ref", levels = c("zone", "ref")), x = c(REF_W, 2 * REF_W) - (REF_W - 1) / 2, y = bar_lo + nf * 0.005,
                        label = c("Ctrl Glul+ Hep", "cKO eGFP+ Hep"))
  n_pv <- sum(L$rho < 0)
  n_cv <- sum(L$rho > 0)
  ybreaks <- c(nf - n_pv / 2 + 0.5, n_cv / 2 + 0.5)
  ylabs <- c(paste0("PV-high\n(n = ", fmt_n(n_pv), ")"), paste0("CV-high\n(n = ", fmt_n(n_cv), ")"))
  g <- ggplot2::ggplot(L$heat, ggplot2::aes(x, row, fill = zc)) +
    ggplot2::geom_raster(interpolate = FALSE) +
    ggplot2::geom_hline(data = data.frame(block = factor(c("zone", "ref"), levels = c("zone", "ref"))),
                        ggplot2::aes(yintercept = n_cv + 0.5), linewidth = 0.25, colour = "white") +
    ggplot2::geom_rect(data = bar, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = bar_fill, colour = NA, inherit.aes = FALSE) +
    ggplot2::geom_text(data = bar_txt, ggplot2::aes(x, y, label = label, hjust = hjust), colour = bar_txt$col,
                       size = pt(PT_SMALL), family = FONT, fontface = "bold", inherit.aes = FALSE) +
    ggplot2::geom_text(data = ref_txt, ggplot2::aes(x, y, label = label), angle = 90, hjust = 0,
                       size = pt(PT_SMALL), family = FONT, inherit.aes = FALSE) +
    ggplot2::facet_grid(. ~ block, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradientn(colours = HEAT_COLS, limits = c(-ZONE_CLIP, ZONE_CLIP),
                                  breaks = c(-ZONE_CLIP, 0, ZONE_CLIP), name = fill_lab) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = ybreaks, labels = ylabs) +
    ggplot2::coord_cartesian(ylim = c(0.5, bar_hi), clip = "off") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = grid::unit(14, "mm"), barheight = grid::unit(1.6, "mm"),
                                                    title.position = "top", title.hjust = 0.5)) +
    ggplot2::labs(x = paste0(fmt_n(nf), " ", what), y = NULL) +
    theme_nature() +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = PT_SMALL, lineheight = 0.9),
                   axis.title.x = ggplot2::element_text(size = PT_SMALL),
                   strip.text = ggplot2::element_blank(), panel.spacing.x = grid::unit(1.2, "mm"),
                   legend.position = "bottom", legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.box.spacing = grid::unit(1, "mm"),
                   plot.margin = ggplot2::margin(13, 2, 2, 2, "mm"))
  if (show_genes) {
    lab <- data.frame(block = factor("ref", levels = c("zone", "ref")), feature = intersect(L$feats, LABEL_GENES))
    lab$row <- nf + 1 - match(lab$feature, L$feats)
    g <- g + ggrepel::geom_text_repel(data = lab, ggplot2::aes(x = 2 * REF_W + 0.5, y = row, label = feature), inherit.aes = FALSE,
                                      direction = "y", hjust = 0, nudge_x = 2, xlim = c(2 * REF_W + 2, Inf),
                                      size = pt(PT_SMALL), family = FONT, fontface = "italic", force = 1.5,
                                      box.padding = 0.1, segment.size = 0.15, segment.colour = "grey45",
                                      min.segment.length = 0, seed = SEED) +
      ggplot2::theme(plot.margin = ggplot2::margin(13, 11, 2, 2, "mm"))
  }
  g
}

plot_zonation_heatmaps <- function(zh) {
  (zone_heat(zh$peaks, "zonated regions", "Normalized region accessibility") |
     zone_heat(zh$RNA, "zonated genes", "Normalized gene expression", show_genes = TRUE)) +
    patchwork::plot_layout(widths = c(1, 1))
}
