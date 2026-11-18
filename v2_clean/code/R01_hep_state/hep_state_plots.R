# ============================================================
# Panel builders for the hepatocyte state map (hep_state_run.R).
# Requires v2_clean/code/shared/nature_theme.R and hep_state_analysis.R.
# ============================================================

STATE_COLS <- c("Control Glul+ Hep" = "#3B75AF", "cKO eGFP+ Hep" = "#E8892B",
                "cKO Glul+ Hep" = "#7B5EA7", "Other" = "#D9D9D9")
STATE_LABELS <- c("Control Glul+ Hep" = paste0("Control Glul", sup("+"), " Hep"),
                  "cKO eGFP+ Hep" = paste0("cKO eGFP", sup("+"), " Hep"),
                  "cKO Glul+ Hep" = paste0("cKO Glul", sup("+"), " Hep"),
                  "Other" = "Other hepatocytes")
RESID_COLS <- c("cKO eGFP+ Hep" = "#E8892B", "cKO Glul+ Hep" = "#7B5EA7", "cKO other Hep" = "#8C8C8C")
SEQ_COLS <- c("#FFFFD9", "#C7E9B4", "#41B6C4", "#225EA8", "#081D58")
DIV_COLS <- c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B")

fmt_p <- function(p) {
  ifelse(p < 1e-3, sprintf("P = %.0e", p), sprintf("P = %.2g", p)) |>
    sub(pattern = "e-0?(\\d+)$", replacement = paste0(" × 10", sup("−\\1")))
}

genotype_strip <- function() {
  ggplot2::theme(strip.text = ggplot2::element_text(face = "plain", size = PT_TEXT),
                 panel.spacing = grid::unit(2, "mm"))
}

# ---------------- a: identity-loss state map ----------------
plot_state_map <- function(cells) {
  set.seed(SEED)
  bg <- cells[cells$state == "Other", ]
  bg <- bg[sample(nrow(bg)), ]
  fg <- cells[cells$state != "Other", ]
  cen <- fg |>
    dplyr::group_by(genotype, state) |>
    dplyr::summarise(x = stats::median(pp_score), y = stats::median(pc_score), .groups = "drop")
  ctrl <- cen[cen$state == "Control Glul+ Hep", ]
  ko <- cen[cen$state == "cKO eGFP+ Hep", ]
  arrow_df <- data.frame(genotype = factor("β-cat cKO", levels = levels(cells$genotype)),
                         x = ctrl$x, y = ctrl$y, xend = ko$x, yend = ko$y)
  ghost <- transform(ctrl, genotype = factor("β-cat cKO", levels = levels(cells$genotype)))
  ggplot2::ggplot(mapping = ggplot2::aes(pp_score, pc_score)) +
    ggrastr::rasterise(ggplot2::geom_point(data = bg, colour = STATE_COLS[["Other"]], size = 0.15, stroke = 0),
                       dpi = 600, dev = "ragg") +
    ggrastr::rasterise(ggplot2::geom_point(data = fg, ggplot2::aes(colour = state), size = 0.35, stroke = 0,
                                           alpha = 0.8), dpi = 600, dev = "ragg") +
    ggplot2::geom_point(data = ghost, ggplot2::aes(x, y), shape = 21, fill = "white",
                        colour = STATE_COLS[["Control Glul+ Hep"]], size = 1.6, stroke = 0.5) +
    ggplot2::geom_segment(data = arrow_df, ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                          linewidth = 0.35, arrow = grid::arrow(length = grid::unit(1.2, "mm"), type = "closed")) +
    ggplot2::facet_wrap(~genotype, nrow = 1) +
    ggplot2::scale_colour_manual(values = STATE_COLS, labels = STATE_LABELS[1:3], breaks = STATE_LEVELS[1:3],
                                 name = NULL) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 1.6, alpha = 1), nrow = 2)) +
    ggplot2::labs(x = "Periportal signature score", y = "Pericentral WNT-responsive score") +
    theme_nature() + genotype_strip() +
    ggplot2::theme(legend.position = "bottom", legend.text = ggtext::element_markdown(size = PT_SMALL),
                   legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.spacing = grid::unit(0.5, "mm"))
}

# ---------------- b: candidate-factor score on the same map ----------------
mean_if_dense <- function(v) if (length(v) >= 8) mean(v) else NA_real_

plot_factor_map <- function(cells, bins = 18) {
  lim <- stats::quantile(cells$factor_score, c(0.02, 0.98))
  contour_d <- cells[cells$state %in% STATE_LEVELS[1:2], ]
  ggplot2::ggplot(cells, ggplot2::aes(pp_score, pc_score)) +
    ggplot2::stat_summary_2d(ggplot2::aes(z = factor_score), fun = mean_if_dense, bins = bins, na.rm = TRUE) +
    ggplot2::stat_density_2d(data = contour_d, ggplot2::aes(colour = state), breaks = 0.4, linewidth = 0.35,
                             contour_var = "ndensity") +
    ggplot2::facet_wrap(~genotype, nrow = 1) +
    ggplot2::scale_fill_gradientn(colours = SEQ_COLS, limits = lim, oob = scales::squish, na.value = NA,
                                  name = "Igf1/Nampt/Vegfa\nmodule score") +
    ggplot2::scale_colour_manual(values = STATE_COLS, guide = "none") +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = grid::unit(1.8, "mm"),
                                                    barheight = grid::unit(10, "mm"))) +
    ggplot2::labs(x = "Periportal signature score", y = "Pericentral WNT-responsive score") +
    theme_nature() + genotype_strip() +
    ggplot2::theme(legend.position = "right", legend.box.spacing = grid::unit(1, "mm"))
}

# ---------------- c: occupancy of control-defined zonation deciles ----------------
plot_zone_occupancy <- function(occ) {
  m <- occ |>
    dplyr::group_by(genotype, zone_bin) |>
    dplyr::summarise(pct = mean(pct), .groups = "drop")
  pos <- ggplot2::position_dodge(width = 0.8)
  ggplot2::ggplot(m, ggplot2::aes(zone_bin, pct, fill = genotype)) +
    ggplot2::geom_col(width = 0.75, position = pos) +
    ggplot2::geom_point(data = occ, ggplot2::aes(group = genotype), position = pos, size = 0.5,
                        colour = "black", stroke = 0) +
    ggplot2::geom_hline(yintercept = 100 / N_BINS, linewidth = 0.25, linetype = "22") +
    ggplot2::scale_fill_manual(values = GROUP_COLS, name = NULL) +
    ggplot2::scale_x_discrete(labels = function(x) ifelse(x == "1", "CV", ifelse(x == as.character(N_BINS), "PV", ""))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = "Zonation decile (control-defined)", y = "Hepatocyte nuclei (%)") +
    theme_nature() +
    ggplot2::theme(legend.position = c(0.02, 1), legend.justification = c(0, 1),
                   axis.ticks.x = ggplot2::element_blank())
}

# ---------------- d: heatmap along the zonation coordinate ----------------
plot_zonation_heatmap <- function(hm) {
  ref_lab <- c("Control Glul+ Hep" = "C", "cKO eGFP+ Hep" = "K")
  hm$block <- factor(hm$block, levels = c("Control", "β-cat cKO", "Pericentral groups"),
                     labels = c("Control", "β-cat cKO", "Ref."))
  hm$xlab <- ifelse(is.na(hm$bin), unname(ref_lab[hm$column]),
                    ifelse(hm$bin == "1", "CV", ifelse(hm$bin == as.character(N_BINS), "PV", "")))
  hm$xkey <- factor(hm$column, levels = unique(hm$column[order(hm$block, suppressWarnings(as.integer(hm$bin)))]))
  hm$gene <- factor(hm$gene, levels = rev(unlist(HEATMAP_GENES, use.names = FALSE)))
  labs <- stats::setNames(hm$xlab, hm$xkey)
  ggplot2::ggplot(hm, ggplot2::aes(xkey, gene, fill = pmax(pmin(z, 2), -2))) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.15) +
    ggplot2::facet_grid(module ~ block, scales = "free", space = "free") +
    ggplot2::scale_fill_gradientn(colours = DIV_COLS, limits = c(-2, 2), name = "Row z-score") +
    ggplot2::scale_x_discrete(labels = labs, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(barwidth = grid::unit(1.8, "mm"),
                                                    barheight = grid::unit(10, "mm"))) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_nature() +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(face = "italic", size = 5),
                   axis.text.x = ggtext::element_markdown(size = PT_SMALL, angle = 0),
                   strip.text.x = ggplot2::element_text(face = "plain", size = PT_TEXT),
                   strip.text.y = ggplot2::element_blank(), panel.spacing = grid::unit(0.8, "mm"),
                   legend.box.spacing = grid::unit(1, "mm"))
}

# ---------------- e: candidate factors along the zonation coordinate ----------------
plot_factor_profiles <- function(prof, ref) {
  ref$colour <- ref$state
  ggplot2::ggplot(prof, ggplot2::aes(zonation, mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean - sem, ymax = mean + sem, fill = genotype), alpha = 0.25,
                         colour = NA) +
    ggplot2::geom_line(ggplot2::aes(colour = genotype), linewidth = 0.4) +
    ggplot2::geom_errorbar(data = ref, ggplot2::aes(y = mean, xmin = z_lo, xmax = z_hi, colour = colour),
                           width = 0, linewidth = 0.3, orientation = "y", inherit.aes = FALSE) +
    ggplot2::geom_point(data = ref, ggplot2::aes(zonation, mean, fill = colour), shape = 23, size = 1.4,
                        stroke = 0.25, colour = "black", inherit.aes = FALSE) +
    ggplot2::facet_wrap(~gene, nrow = 1, scales = "free_y") +
    ggplot2::scale_colour_manual(values = c(GROUP_COLS, STATE_COLS), breaks = names(GROUP_COLS), name = NULL) +
    ggplot2::scale_fill_manual(values = c(GROUP_COLS, STATE_COLS), breaks = STATE_LEVELS[1:2],
                               labels = STATE_LABELS[1:2], name = NULL) +
    ggplot2::scale_x_continuous(breaks = c(0, 1), labels = c("CV", "PV"), expand = ggplot2::expansion(add = 0.03)) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1, override.aes = list(linewidth = 0.6)),
                    fill = ggplot2::guide_legend(order = 2, override.aes = list(size = 1.6))) +
    ggplot2::labs(x = "Inferred zonation coordinate", y = "Mean expression (ln)") +
    theme_nature() +
    ggplot2::theme(legend.position = "bottom", legend.text = ggtext::element_markdown(size = PT_SMALL),
                   legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.spacing = grid::unit(0.5, "mm"),
                   panel.spacing = grid::unit(2.5, "mm"))
}

# ---------------- f: zone-matched residuals ----------------
plot_zone_residuals <- function(resid) {
  pos <- ggplot2::position_dodge(width = 0.65)
  lab <- resid[resid$group == "cKO eGFP+ Hep", ]
  lab$mean_residual <- lab$beta
  lab$label <- fmt_p(lab$p)
  resid <- dplyr::rename(resid, mean_residual = beta)
  xr <- max(abs(c(resid$ci_lo, resid$ci_hi)), na.rm = TRUE)
  ggplot2::ggplot(resid, ggplot2::aes(mean_residual, gene, colour = group)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, linetype = "22") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_lo, xmax = ci_hi), width = 0, linewidth = 0.35, position = pos,
                           orientation = "y") +
    ggplot2::geom_point(size = 1.1, position = pos) +
    ggtext::geom_richtext(data = lab, ggplot2::aes(x = xr * 1.12, y = gene, label = label), hjust = 0,
                          size = pt(PT_SMALL), family = FONT, fill = NA, label.colour = NA, colour = "black",
                          label.padding = grid::unit(0, "pt"), inherit.aes = FALSE) +
    ggplot2::scale_colour_manual(values = RESID_COLS, name = NULL,
                                 labels = c("cKO eGFP+ Hep" = paste0("cKO eGFP", sup("+"), " Hep"),
                                            "cKO Glul+ Hep" = paste0("cKO Glul", sup("+"), " Hep"),
                                            "cKO other Hep" = "cKO other Hep")) +
    ggplot2::scale_x_continuous(limits = c(-xr * 1.05, xr * 1.95), breaks = scales::breaks_pretty(4)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(x = "Zone-adjusted difference vs control\n(ln expression, β ± 95% CI)", y = NULL) +
    theme_nature() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(face = "italic"), legend.position = "bottom",
                   legend.direction = "vertical", legend.text = ggtext::element_markdown(size = PT_SMALL),
                   legend.margin = ggplot2::margin(0, 0, 0, 0), legend.box.spacing = grid::unit(0.5, "mm"))
}
