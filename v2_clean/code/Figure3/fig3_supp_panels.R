# ============================================================
# Supplementary panels for Figure 3 (Nature style).
# Requires v2_clean/code/shared/nature_theme.R and fig3_analysis.R.
#   S3a  candidate expression, control Glul+ vs cKO eGFP+ hepatocytes
#   S3-markers  marker dot plot supporting the v1 annotation
# ============================================================

# Exact P as rich text, e.g. "P = 3.2 x 10^-5".
format_p_rich <- function(p) {
  vapply(p, function(v) {
    if (v >= 0.001) return(sprintf("*P* = %.3f", v))
    e <- floor(log10(v))
    sprintf("*P* = %.1f × 10%s", v / 10^e, sup(paste0("−", -e)))
  }, character(1))
}

p_annotations <- function(long, stats) {
  top <- long |>
    dplyr::group_by(gene) |>
    dplyr::summarise(ymax = max(expr), .groups = "drop")
  d <- merge(top, stats[, c("gene", "p")])
  transform(d, gene = factor(gene, levels = CANDIDATES), y = ymax * 1.08, ylab = ymax * 1.2,
            label = format_p_rich(p))
}

group_axis_labels <- function(stats) {
  c(Control = paste0("Control<br>Glul", sup("+"), " Hep<br>(n = ", stats$n_control[1], ")"),
    "β-cat cKO" = paste0("β-cat cKO<br>eGFP", sup("+"), " Hep<br>(n = ", stats$n_cko[1], ")"))
}

plot_s3a <- function(long, stats, libs) {
  ann <- p_annotations(long, stats)
  set.seed(42)
  ggplot2::ggplot(long, ggplot2::aes(group, expr)) +
    ggrastr::rasterise(
      ggplot2::geom_jitter(ggplot2::aes(colour = group), width = 0.18, height = 0, size = 0.55,
                           stroke = 0, alpha = 0.6), dpi = 600, dev = "ragg") +
    ggplot2::geom_boxplot(ggplot2::aes(fill = group), width = 0.45, outlier.shape = NA, linewidth = 0.25,
                          alpha = 0.55, colour = "black") +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 23, size = 1.3, fill = "white",
                          colour = "black", stroke = 0.3) +
    ggplot2::geom_point(data = libs, ggplot2::aes(y = mean_expr), shape = 21, fill = "black", colour = "white",
                        stroke = 0.2, size = 1.1, position = ggplot2::position_nudge(x = 0.33)) +
    ggplot2::geom_segment(data = ann, ggplot2::aes(x = 1, xend = 2, y = y, yend = y), linewidth = 0.25) +
    ggtext::geom_richtext(data = ann, ggplot2::aes(x = 1.5, y = ylab, label = label), size = pt(PT_SMALL),
                          family = FONT, fill = NA, label.colour = NA, label.padding = grid::unit(0, "pt")) +
    ggplot2::facet_wrap(~gene, nrow = 2, scales = "free_y") +
    ggplot2::scale_colour_manual(values = GROUP_COLS, guide = "none") +
    ggplot2::scale_fill_manual(values = GROUP_COLS, guide = "none") +
    ggplot2::scale_x_discrete(labels = group_axis_labels(stats)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.1))) +
    ggplot2::labs(x = NULL, y = "ln(UMI per 10,000 + 1)") +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggtext::element_markdown(size = PT_SMALL, lineheight = 1.35),
                   panel.spacing.x = grid::unit(3, "mm"), panel.spacing.y = grid::unit(2, "mm"))
}

plot_markers_v1 <- function(markers) {
  n <- unique(markers[, c("celltype", "n")])
  lab <- setNames(paste0(n$celltype, " (", fmt_n(n$n), ")"), n$celltype)
  ggplot2::ggplot(markers, ggplot2::aes(gene, celltype)) +
    ggplot2::geom_point(ggplot2::aes(size = pct, colour = scaled)) +
    ggplot2::scale_size_area(max_size = 3.2, breaks = c(25, 50, 75), name = "Nuclei\ndetecting (%)") +
    ggplot2::scale_colour_gradient(low = "#E6EBF2", high = "#1F3B63", breaks = c(0, 1),
                                   labels = c("Min", "Max"), name = "Scaled\nmean") +
    ggplot2::scale_y_discrete(labels = lab) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(face = "italic", angle = 90, hjust = 1, vjust = 0.5),
                   axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA, linewidth = 0.25),
                   legend.key.height = grid::unit(2.5, "mm"))
}
