# ============================================================
# Plot helpers for v2_clean/code/Figure4/fig4_run.R (theme from shared/nature_theme.R).
# Layout follows manuscript Fig. 4B-C; every panel is saved on its own.
#   Fig4B_TF_venn          TF overlap (Wnt2 vs Wnt9b TF sets), Venn diagram
#   Fig4B_TF_filter_matrix motif / up-in-ligand+ criteria for candidate TFs
#   Fig4C_Ets1_expression  Ets1, ligand- vs ligand+ EC (Wilcoxon P)
#   FigS4A_TF_venn         motif TFs vs up genes (Wnt2, Wnt9b) and TF overlap
#   FigS4B_Tead1_expression        as Fig4C for Tead1
#   FigS4_Nr2f2_expression         Nr2f2 (reference, not in the figure)
#   FigS4_other_TF_expression      Klf2 and Prdm1 (reference)
#   FigS4_volcano_<ligand> genome-wide ligand+ vs ligand- EC
# All P values are nominal two-sided Wilcoxon P (FindMarkers); TF screen uses P < P_CUT.
# ============================================================

LIGAND_COLS <- c("#BDBDBD", "#3B75AF")
UP_COL      <- "#C25B56"
VENN_COLS   <- c(Wnt2 = "#3B75AF", Wnt9b = "#E8892B")

format_sci <- function(p) {
  e <- floor(log10(p))
  sprintf("'%.1f' %%*%% 10^%d", p / 10^e, e)
}

fmt_p <- function(p) {
  ifelse(is.na(p), "'n.t.'",
         ifelse(p < 1e-3, sprintf("italic(P) == %s", format_sci(p)), sprintf("italic(P) == '%.2g'", p)))
}

# ---------------- volcano ----------------
plot_volcano <- function(de, ligand, label_genes) {
  d <- transform(de[de$gene != ligand, ], y = -log10(pmax(p_val, 1e-300)),
                 cls = ifelse(gene %in% label_genes, "focal", ifelse(up, "up", "ns")))
  d <- d[order(match(d$cls, c("ns", "up", "focal"))), ]
  lab <- d[d$cls == "focal" | (d$up & rank(-d$y) <= 6), ]
  ggplot(d, aes(avg_log2FC, y)) +
    geom_hline(yintercept = -log10(P_CUT), linewidth = 0.25, linetype = "dashed", colour = "grey40") +
    geom_point(aes(colour = cls, size = cls == "focal"), stroke = 0) +
    scale_size_manual(values = c(`FALSE` = 0.5, `TRUE` = 1.1), guide = "none") +
    ggrepel::geom_label_repel(data = lab, aes(label = gene, colour = cls), size = pt(PT_SMALL), family = FONT,
                              fontface = "italic", min.segment.length = 0, segment.size = 0.2,
                              max.overlaps = Inf, seed = SEED, show.legend = FALSE,
                              fill = scales::alpha("white", 0.85), label.size = 0, label.padding = 0.08,
                              box.padding = 0.4) +
    scale_colour_manual(values = c(ns = "grey75", up = UP_COL, focal = "black"), guide = "none") +
    labs(x = bquote(log[2]~"fold change ("*italic(.(ligand))^"+"~"vs"~italic(.(ligand))^"–"~"EC)"),
         y = bquote(-log[10]~italic(P))) +
    theme_nature()
}

# ---------------- Venn ----------------
tf_sets <- function(res) lapply(res$tfs, function(t) t$gene[t$up])

# Two-set Venn (equal circles). `titles` are plotmath strings; `note` sits in the overlap.
plot_venn2 <- function(a, b, titles, cols, bottom = NULL, note = NULL) {
  shared <- intersect(a, b)
  n <- c(length(setdiff(a, shared)), length(shared), length(setdiff(b, shared)))
  circles <- data.frame(set = c("a", "b"), x0 = c(-0.55, 0.55), y0 = 0, r = 1)
  p <- ggplot() +
    ggforce::geom_circle(data = circles, aes(x0 = x0, y0 = y0, r = r, fill = set, colour = set),
                         alpha = 0.3, linewidth = 0.4) +
    annotate("text", x = c(-0.95, 0, 0.95), y = 0.08, label = n, size = pt(PT_LABEL), family = FONT) +
    annotate("text", x = c(-0.85, 0.85), y = 1.12, label = titles, parse = TRUE, size = pt(PT_TEXT),
             family = FONT, colour = unname(cols), vjust = 0, lineheight = 0.9) +
    scale_fill_manual(values = c(a = unname(cols[1]), b = unname(cols[2])), guide = "none") +
    scale_colour_manual(values = c(a = unname(cols[1]), b = unname(cols[2])), guide = "none") +
    coord_fixed(xlim = c(-1.75, 1.75), ylim = c(-1.35, 1.75), clip = "off") +
    theme_void(base_family = FONT)
  if (!is.null(note)) {
    p <- p + annotate("text", x = 0, y = -0.22, label = note, size = pt(PT_SMALL), family = FONT,
                      fontface = "bold", lineheight = 0.9)
  }
  if (!is.null(bottom)) {
    p <- p + annotate("text", x = 0, y = -1.12, label = bottom, parse = TRUE, size = pt(PT_TEXT),
                      family = FONT, vjust = 1)
  }
  p
}

plot_tf_venn <- function(res) {
  s <- tf_sets(res)
  plot_venn2(s$Wnt2, s$Wnt9b, VENN_COLS,
             titles = c(sprintf("atop('TFs in', italic(Wnt2)^'+'~EC~'(%d)')", length(s$Wnt2)),
                        sprintf("atop('TFs in', italic(Wnt9b)^'+'~EC~'(%d)')", length(s$Wnt9b))),
             bottom = "'TFs in'~italic(Wnt)^'+'~EC",
             note = paste(toupper(intersect(s$Wnt2, s$Wnt9b)), collapse = "\n"))
}

# Motif-mapped TF genes at the ligand locus vs genes up in ligand+ EC; overlap = TFs in ligand+ EC.
plot_motif_de_venn <- function(res, ligand, cols) {
  plot_venn2(res$tfs[[ligand]]$gene, res$de[[ligand]]$gene[res$de[[ligand]]$up], cols,
             titles = c(sprintf("atop('Motifs in accessible', 'chromatin at'~italic(%s))", ligand),
                        sprintf("atop('Up in', italic(%s)^'+'~EC)", ligand)),
             bottom = sprintf("'TFs in'~italic(%s)^'+'~EC", ligand))
}

plot_s4a <- function(res) {
  (plot_motif_de_venn(res, "Wnt2", c("#6FA36B", "#8FB4D9")) |
     plot_motif_de_venn(res, "Wnt9b", c("#E3A35F", "#8CC97A")) |
     plot_tf_venn(res))
}

# ---------------- filter matrix ----------------
filter_matrix <- function(res) {
  s <- tf_sets(res)
  genes <- unique(c(FOCAL_TFS, intersect(s$Wnt2, s$Wnt9b)))
  crit <- c("Wnt2 motif", "Wnt9b motif", "Up in Wnt2+ EC", "Up in Wnt9b+ EC")
  hit <- function(g) c(g %in% res$tfs$Wnt2$gene, g %in% res$tfs$Wnt9b$gene,
                       isTRUE(res$de$Wnt2$up[match(g, res$de$Wnt2$gene)]),
                       isTRUE(res$de$Wnt9b$up[match(g, res$de$Wnt9b$gene)]))
  d <- do.call(rbind, lapply(genes, function(g) data.frame(gene = g, criterion = crit, met = hit(g))))
  d$all <- stats::ave(d$met, d$gene, FUN = all)
  d$gene <- factor(toupper(d$gene), levels = rev(toupper(genes)))
  d$criterion <- factor(d$criterion, levels = crit)
  d
}

plot_tf_filter <- function(res) {
  d <- filter_matrix(res)
  sel <- match(unique(as.character(d$gene[d$all == 1])), levels(d$gene))
  hl <- data.frame(ymin = sel - 0.4, ymax = sel + 0.4)
  ggplot(d, aes(criterion, gene)) +
    geom_rect(data = hl, aes(xmin = 0.5, xmax = 4.5, ymin = ymin, ymax = ymax), inherit.aes = FALSE,
              fill = "#F2D4D4") +
    geom_point(aes(shape = met), size = 1.6, stroke = 0.35, fill = "black") +
    scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), breaks = c(TRUE, FALSE),
                       labels = c("Met", "Not met"), name = NULL) +
    scale_x_discrete(labels = c("Wnt2 motif", "Wnt9b motif", expression("Up in"~italic(Wnt2)^"+"~EC),
                                expression("Up in"~italic(Wnt9b)^"+"~EC))) +
    labs(x = NULL, y = NULL) +
    theme_nature() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.line = element_blank(),
          axis.ticks = element_blank(), legend.position = "bottom")
}

# ---------------- expression boxplots ----------------
expr_panel <- function(res, gene, ligand) {
  d <- res$long[res$long$gene == gene & res$long$ligand == ligand, ]
  d$grp <- factor(d$group, levels = paste0(ligand, c("-", "+")))
  p <- res$focal$p_val[res$focal$gene == gene & res$focal$ligand == ligand]
  top <- max(d$expr, 0.5) * 1.08
  ggplot(d, aes(grp, expr)) +
    geom_jitter(width = 0.18, height = 0, size = 0.25, stroke = 0, colour = "grey55") +
    geom_boxplot(aes(fill = grp), width = 0.55, outlier.shape = NA, linewidth = 0.25, alpha = 0.85) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 1, fill = "white", stroke = 0.3) +
    annotate("segment", x = 1, xend = 2, y = top, yend = top, linewidth = 0.25) +
    annotate("text", x = 1.5, y = top * 1.03, label = fmt_p(p), parse = TRUE, vjust = 0,
             size = pt(PT_SMALL), family = FONT) +
    scale_fill_manual(values = stats::setNames(LIGAND_COLS, levels(d$grp)), guide = "none") +
    scale_x_discrete(labels = c(bquote(italic(.(ligand))^"–"~EC), bquote(italic(.(ligand))^"+"~EC))) +
    scale_y_continuous(limits = c(0, top * 1.15), expand = expansion(mult = c(0.02, 0))) +
    labs(x = NULL, y = bquote(italic(.(gene))~"expression")) +
    theme_nature()
}

plot_gene_expression <- function(res, gene) {
  expr_panel(res, gene, "Wnt9b") | expr_panel(res, gene, "Wnt2")
}

plot_tf_boxes <- function(res, genes) {
  panels <- unlist(lapply(genes, function(g) list(expr_panel(res, g, "Wnt9b"), expr_panel(res, g, "Wnt2"))),
                   recursive = FALSE)
  wrap_plots(panels, ncol = 4)
}

# ---------------- export ----------------
export_figures <- function(res, dir) {
  focal <- intersect(FOCAL_TFS, unique(c(res$de$Wnt2$gene, res$de$Wnt9b$gene)))
  save_nature(plot_tf_venn(res), file.path(dir, "Fig4B_TF_venn"), 70, 48)
  save_nature(plot_tf_filter(res), file.path(dir, "Fig4B_TF_filter_matrix"), 55, 62)
  save_nature(plot_gene_expression(res, "Ets1"), file.path(dir, "Fig4C_Ets1_expression"), 70, 50)
  save_nature(plot_s4a(res), file.path(dir, "FigS4A_TF_venn"), 183, 62)
  save_nature(plot_gene_expression(res, "Tead1"), file.path(dir, "FigS4B_Tead1_expression"), 70, 50)
  save_nature(plot_tf_boxes(res, "Nr2f2"), file.path(dir, "FigS4_Nr2f2_expression"), 70, 50)
  save_nature(plot_tf_boxes(res, c("Klf2", "Prdm1")), file.path(dir, "FigS4_other_TF_expression"), 160, 50)
  for (l in LIGANDS) {
    save_nature(plot_volcano(res$de[[l]], l, focal), file.path(dir, paste0("FigS4_volcano_", l)), 75, 65)
  }
}
