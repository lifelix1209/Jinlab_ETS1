# ============================================================
# Panel builders for the NicheNet figure (nichenet_run.R).
# Rows (ligands) share one order across panels: cKO-response activity.
# ============================================================

ligand_levels <- function(res) rev(res$shown)

tile_theme <- function() {
  theme_nature() +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(face = "italic"),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic"),
                   legend.position = "bottom", legend.margin = ggplot2::margin(0, 0, 0, 0),
                   legend.box.spacing = grid::unit(0.5, "mm"),
                   panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.25),
                   plot.subtitle = ggplot2::element_text(size = PT_TEXT, hjust = 0.5))
}

bar_guide <- function() ggplot2::guide_colourbar(barwidth = grid::unit(11, "mm"), barheight = grid::unit(1.6, "mm"),
                                                  title.position = "top")
three_breaks <- function(lim) signif(seq(lim[1], lim[2], length.out = 3), 2)

# ---------------- a: ligand activity ----------------
plot_ligand_activity <- function(res) {
  d <- res$activity[res$activity$ligand %in% res$shown, ]
  d$ligand <- factor(d$ligand, levels = ligand_levels(res))
  n_set <- stats::setNames(res$activity$n_genes[!duplicated(res$activity$gene_set)],
                           res$activity$gene_set[!duplicated(res$activity$gene_set)])
  d$gene_set <- factor(d$gene_set, levels = names(n_set), labels = paste0(c("Response", "Wnt-EC"), "\n(n = ", n_set, ")"))
  ggplot2::ggplot(d, ggplot2::aes(gene_set, ligand, fill = aupr_corrected)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = rank), size = pt(PT_SMALL - 0.5), family = FONT) +
    ggplot2::scale_fill_gradient2(low = "#4393C3", mid = "#F7F7F7", high = "#E8892B", midpoint = 0,
                                  name = "AUPR\n(corrected)", breaks = function(l) signif(l * 0.9, 1)) +
    ggplot2::guides(fill = bar_guide()) +
    ggplot2::labs(x = NULL, y = "Hepatocyte ligands", subtitle = "Activity (rank)") +
    tile_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = "plain", size = PT_SMALL))
}

# ---------------- b: ligand -> target regulatory potential ----------------
plot_ligand_targets <- function(res) {
  d <- res$targets
  if (is.null(d) || nrow(d) == 0) stop("No ligand-target links among shown ligands")
  d$ligand <- factor(d$ligand, levels = ligand_levels(res))
  ord <- stats::aggregate(weight ~ target, d, sum)
  d$target <- factor(d$target, levels = ord$target[order(-ord$weight)])
  ggplot2::ggplot(d, ggplot2::aes(target, ligand, fill = weight)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_fill_gradient(low = "#EFE7F6", high = "#6A3D9A", name = "Regulatory\npotential", breaks = three_breaks) +
    ggplot2::guides(fill = bar_guide()) +
    ggplot2::labs(x = NULL, y = NULL, subtitle = "Targets") +
    tile_theme() + ggplot2::theme(axis.text.y = ggplot2::element_blank())
}

# ---------------- c: ligand -> receptor prior ----------------
plot_ligand_receptors <- function(res) {
  d <- res$receptors
  d$ligand <- factor(d$ligand, levels = ligand_levels(res))
  ggplot2::ggplot(d, ggplot2::aes(receptor, ligand, fill = weight)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::scale_fill_gradient(low = "#FCE4EF", high = "#C2185B", name = "Prior interaction\npotential", breaks = three_breaks) +
    ggplot2::guides(fill = bar_guide()) +
    ggplot2::labs(x = NULL, y = NULL, subtitle = "EC receptors") +
    tile_theme() + ggplot2::theme(axis.text.y = ggplot2::element_blank())
}

# ---------------- d: sender expression ----------------
plot_sender <- function(res) {
  d <- res$sender
  d$ligand <- factor(d$ligand, levels = ligand_levels(res))
  d$group <- factor(d$group, levels = HEP_GROUPS, labels = c("Ctrl Glul+", "cKO eGFP+", "cKO other"))
  d$scaled <- stats::ave(d$mean, d$ligand, FUN = function(v) (v - mean(v)) / max(stats::sd(v), 1e-9))
  ggplot2::ggplot(d, ggplot2::aes(group, ligand)) +
    ggplot2::geom_point(ggplot2::aes(size = 100 * pct, fill = scaled), shape = 21, stroke = 0.2) +
    ggplot2::scale_size_area(max_size = 3, limits = c(0, 100), breaks = c(25, 75), name = "% nuclei") +
    ggplot2::scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"), limits = c(-1.2, 1.2),
                                  oob = scales::squish, name = "Scaled\nmean", breaks = c(-1, 0, 1)) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ggplot2::guides(fill = bar_guide(), size = ggplot2::guide_legend(title.position = "top", ncol = 1)) +
    ggplot2::labs(x = NULL, y = NULL, subtitle = "Sender") +
    tile_theme() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), legend.box = "vertical",
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, face = "plain"))
}

# ---------------- e: prior regulatory potential toward niche genes ----------------
plot_niche_prior <- function(res) {
  d <- res$prior[res$prior$ligand %in% res$shown, ]
  d$ligand <- factor(d$ligand, levels = ligand_levels(res))
  d$target <- factor(d$target, levels = NICHE_GENES)
  n_lig <- res$n$n_potential_ligands
  ggplot2::ggplot(d, ggplot2::aes(target, ligand, fill = percentile)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = round(percentile)), size = pt(PT_SMALL - 0.5), family = FONT,
                       colour = ifelse(d$percentile > 75, "white", "black")) +
    ggplot2::scale_fill_gradient(low = "#F7F7F7", high = "#1B5E20", limits = c(0, 100),
                                 name = paste0("Percentile among\n", n_lig, " ligands"), breaks = c(0, 50, 100)) +
    ggplot2::guides(fill = bar_guide()) +
    ggplot2::labs(x = NULL, y = NULL, subtitle = "Niche-gene prior") +
    tile_theme() + ggplot2::theme(axis.text.y = ggplot2::element_blank())
}
