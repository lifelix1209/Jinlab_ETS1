# ============================================================
# Panel builders for the EC response figure (ec_response_run.R).
# ============================================================

AMB_COLS <- c("Ambient-dominated" = "#C9A27E", "EC-intrinsic" = "#2F6F73")
NPC_TYPES <- c("Hepatocyte", "Endothelial", "Kupffer cell", "Stellate cell", "Immune cell")
GSEA_ORDER <- c("Naive pseudobulk", "Corrected pseudobulk", "Corrected nucleus-level")

short_term <- function(x, width = 44) {
  x <- sub("^(HALLMARK|REACTOME|GOBP)_", "", x)
  x <- gsub("_", " ", tolower(x))
  x <- gsub("\\bmapk\\b", "MAPK", x)
  x <- sub("^(.)", "\\U\\1", x, perl = TRUE)
  x <- gsub("\\berk1\\b", "ERK1", x)
  x <- gsub("\\berk2\\b", "ERK2", x)
  x <- gsub("\\bvegf\\b", "VEGF", x)
  ifelse(nchar(x) > width, paste0(substr(x, 1, width - 1), "…"), x)
}

# ---------------- a: ambient contamination per nucleus ----------------
plot_contamination <- function(meta) {
  d <- meta[meta$celltype %in% NPC_TYPES, ]
  d$celltype <- factor(d$celltype, levels = NPC_TYPES, labels = c("Hep", "EC", "Kupffer", "Stellate", "Immune"))
  pos <- ggplot2::position_dodge(width = 0.75)
  ggplot2::ggplot(d, ggplot2::aes(celltype, 100 * contamination, fill = genotype)) +
    ggplot2::geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.25, position = pos) +
    ggplot2::scale_fill_manual(values = GROUP_COLS, name = NULL) +
    ggplot2::scale_y_continuous(limits = c(0, 100), expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = NULL, y = "Ambient RNA per nucleus\n(DecontX, %)") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top",
                   legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.spacing = grid::unit(0.5, "mm"))
}

volcano_base <- function(d, x, p, colour_by, cols, labels = NULL, ylab = "−log10 P") {
  d$y <- -log10(d[[p]])
  d$x <- d[[x]]
  d$col <- d[[colour_by]]
  d <- d[order(d$col != names(cols)[1]), ]
  g <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggrastr::rasterise(ggplot2::geom_point(ggplot2::aes(colour = col), size = 0.45, stroke = 0, alpha = 0.8),
                       dpi = 600, dev = "ragg") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 1.6, alpha = 1))) +
    ggplot2::labs(x = "log2 fold change (cKO / control)", y = ylab) +
    theme_nature() +
    ggplot2::theme(legend.position = "top", legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.box.spacing = grid::unit(0.5, "mm"))
  if (!is.null(labels)) {
    lab <- d[d$gene %in% labels, ]
    g <- g + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = gene), size = pt(PT_SMALL), family = FONT,
                                      fontface = "italic", segment.size = 0.2, min.segment.length = 0,
                                      max.overlaps = Inf, box.padding = 0.4, force = 3, seed = SEED)
  }
  g
}

# ---------------- b: naive volcano, coloured by ambient status ----------------
plot_naive_volcano <- function(naive) {
  naive$status <- factor(ifelse(naive$intrinsic, "EC-intrinsic", "Ambient-dominated"), levels = names(AMB_COLS))
  dn <- naive[!naive$intrinsic & naive$logFC < 0, ]
  lab <- utils::head(dn$gene[order(dn$PValue)], 6)
  volcano_base(naive, "logFC", "PValue", "status", AMB_COLS, labels = lab) +
    ggplot2::labs(subtitle = "All genes") +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = PT_TEXT))
}

# ---------------- d: corrected volcano ----------------
plot_corrected_volcano <- function(corr) {
  corr$status <- factor(ifelse(corr$PValue < 0.05, ifelse(corr$logFC > 0, "Up (P < 0.05)", "Down (P < 0.05)"), "NS"),
                        levels = c("Up (P < 0.05)", "Down (P < 0.05)", "NS"))
  cols <- c("Up (P < 0.05)" = GROUP_COLS[["β-cat cKO"]], "Down (P < 0.05)" = GROUP_COLS[["Control"]], NS = "#C8C8C8")
  top <- utils::head(grep("^Gm\\d|Rik$", corr$gene[order(corr$PValue)], value = TRUE, invert = TRUE), 5)
  volcano_base(corr, "logFC", "PValue", "status", cols, labels = unique(c(top, setdiff(SHOW_GENES, c("Wnt2", "Wnt9b"))))) +
    ggplot2::labs(subtitle = paste0("EC-intrinsic genes; min FDR = ", signif(min(corr$FDR), 2))) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = PT_TEXT))
}

# ---------------- c: GSEA, naive vs corrected ----------------
select_terms <- function(gsea, n_naive = 5, n_corr = 4) {
  nv <- gsea[gsea$analysis == "Naive pseudobulk", ]
  cn <- gsea[gsea$analysis == "Corrected nucleus-level", ]
  pick <- c(utils::head(nv$pathway[nv$NES < 0], n_naive), utils::head(nv$pathway[nv$NES > 0], n_naive),
            utils::head(cn$pathway[cn$NES > 0], n_corr), utils::head(cn$pathway[cn$NES < 0], n_corr),
            "GOBP_MAPK_CASCADE", "GOBP_ERK1_AND_ERK2_CASCADE", "REACTOME_SIGNALING_BY_VEGF",
            "REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES")
  unique(pick)
}

plot_gsea_compare <- function(gsea) {
  terms <- select_terms(gsea)
  d <- gsea[gsea$pathway %in% terms, ]
  d$analysis <- factor(d$analysis, levels = GSEA_ORDER,
                       labels = c("Naive\npseudobulk", "Corrected\npseudobulk", "Corrected\nnucleus-level"))
  ord <- gsea[gsea$analysis == "Naive pseudobulk" & gsea$pathway %in% terms, ]
  d$term <- factor(short_term(d$pathway), levels = rev(unique(short_term(terms))))
  d$sig <- -log10(pmax(d$p_BH, 1e-10))
  ggplot2::ggplot(d, ggplot2::aes(analysis, term)) +
    ggplot2::geom_point(ggplot2::aes(fill = pmax(pmin(NES, 2.5), -2.5), size = sig), shape = 21, stroke = 0.2) +
    ggplot2::geom_point(data = d[d$p_BH < 0.05, ], ggplot2::aes(size = sig), shape = 21, fill = NA, stroke = 0.5) +
    ggplot2::scale_fill_gradientn(colours = c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"),
                                  limits = c(-2.5, 2.5), name = "NES") +
    ggplot2::scale_size_area(max_size = 3.2, limits = c(0, 10), breaks = c(1.3, 5, 10), labels = c("1.3", "5", "≥10"),
                             name = "−log10 FDR") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = grid::unit(1.8, "mm"), barheight = grid::unit(10, "mm"),
                                                    order = 1)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = PT_SMALL), legend.box.spacing = grid::unit(1, "mm"),
                   panel.grid.major = ggplot2::element_line(linewidth = 0.15, colour = "grey90"))
}

# ---------------- e: per-library expression ----------------
plot_library_genes <- function(lib) {
  lib$gene <- factor(lib$gene, levels = SHOW_GENES)
  m <- stats::aggregate(cpm ~ gene + genotype, lib, mean)
  ggplot2::ggplot(lib, ggplot2::aes(genotype, cpm, colour = genotype)) +
    ggplot2::geom_col(data = m, ggplot2::aes(fill = genotype), colour = NA, width = 0.65, alpha = 0.35) +
    ggplot2::geom_point(size = 0.9) +
    ggplot2::facet_wrap(~gene, nrow = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(values = GROUP_COLS, name = NULL) +
    ggplot2::scale_fill_manual(values = GROUP_COLS, guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(x = NULL, y = "Pseudobulk CPM\n(EC-intrinsic genes)") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
                   legend.position = "bottom", legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.box.spacing = grid::unit(0.5, "mm"), panel.spacing = grid::unit(1.5, "mm"))
}
