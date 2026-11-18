# ============================================================
# R04 — RNA-ATAC integration: ATAC-defined populations, Wnt2/Wnt9b gene
# activity, and cell-type-specific accessibility at the Wnt2 / Wnt9b loci
# with ETS-class motif sites and peak-to-gene links (pre-Fig. 4A panel).
#
# Input : data/processed/seu_bko_lib1245_singlets.rds (+ data/raw/atac_fragments.tsv.gz)
# Output: v2_clean/results/R04_rna_atac/  composite, single panels, source_data/
# Usage : Rscript v2_clean/code/R04_rna_atac/rna_atac_run.R (FORCE_ANALYSIS=1 to recompute)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

source("v2_clean/code/shared/nature_theme.R")
source("v2_clean/code/Figure3/fig3_analysis.R")
source("v2_clean/code/Figure3/fig3_main_panels.R")
source("v2_clean/code/R04_rna_atac/rna_atac_analysis.R")
source("v2_clean/code/R04_rna_atac/rna_atac_plots.R")

INPUT_RDS <- "data/processed/seu_bko_lib1245_singlets.rds"
OUT_DIR   <- "v2_clean/results/R04_rna_atac"
SRC_DIR   <- file.path(OUT_DIR, "source_data")
CACHE_RDS <- file.path(SRC_DIR, "analysis_results.rds")

load_results <- function() {
  if (file.exists(CACHE_RDS) && !identical(Sys.getenv("FORCE_ANALYSIS"), "1")) return(readRDS(CACHE_RDS))
  obj <- suppressPackageStartupMessages({ library(Signac); load_multiome(INPUT_RDS) })
  res <- run_rna_atac(obj)
  write_atac_source_data(res, SRC_DIR)
  saveRDS(res, CACHE_RDS)
  res
}

compose_r04 <- function(res) {
  top <- (plot_atac_umap(res$umap) | plot_rna_vs_ga(res$dot) | plot_peak_heatmap(res$loci)) +
    plot_layout(widths = c(52, 62, 30))
  bottom <- (plot_locus(res$loci$Wnt2) | plot_locus(res$loci$Wnt9b))
  (top / bottom) + plot_layout(heights = c(58, 70)) + plot_annotation(tag_levels = "a") & ggplot2::theme(plot.tag = ggplot2::element_text(size = PT_LABEL, face = "bold", family = FONT))
}

main <- function() {
  res <- load_results()
  save_nature(compose_r04(res), file.path(OUT_DIR, "R04_rna_atac_composite"), 183, 140)
  save_nature(plot_atac_umap(res$umap), file.path(OUT_DIR, "R04a_atac_umap"), 60, 58)
  save_nature(plot_rna_vs_ga(res$dot), file.path(OUT_DIR, "R04b_rna_vs_gene_activity"), 75, 50)
  save_nature(plot_peak_heatmap(res$loci), file.path(OUT_DIR, "R04c_peak_specificity"), 70, 70)
  save_nature(plot_locus(res$loci$Wnt2), file.path(OUT_DIR, "R04d_Wnt2_tracks"), 90, 75)
  save_nature(plot_locus(res$loci$Wnt9b), file.path(OUT_DIR, "R04e_Wnt9b_tracks"), 90, 75)
  print(res$links, digits = 3)
  message("Done. Figures in ", OUT_DIR)
}

tryCatch(main(), error = function(e) {
  message("R04 failed: ", conditionMessage(e))
  quit(status = 1)
})
