# ============================================================
# R05 — Chromatin landscape panels modelled on the reference figures:
#   linked transcriptome/epigenome UMAPs, peaks per cell type by genomic
#   category, marker-locus coverage matrix (incl. Wnt2/Wnt9b), and
#   hepatocyte zonation heatmaps of accessible regions and genes.
#
# Input : data/processed/seu_bko_lib1245_singlets.rds, data/raw/atac_fragments.tsv.gz
#         v2_clean/results/R01_hep_state/source_data/analysis_results.rds (zonation)
#         v2_clean/results/R04_rna_atac/source_data/analysis_results.rds (ATAC UMAP)
#         MACS3 executable (env MACS3_PATH)
# Output: v2_clean/results/R05_atac_landscape/  composite, single panels, source_data/
# Usage : MACS3_PATH=/path/to/macs3 Rscript v2_clean/code/R05_atac_landscape/atac_landscape_run.R
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
source("v2_clean/code/R05_atac_landscape/atac_landscape_analysis.R")
source("v2_clean/code/R05_atac_landscape/atac_landscape_plots.R")

INPUT_RDS <- "data/processed/seu_bko_lib1245_singlets.rds"
R01_RDS   <- "v2_clean/results/R01_hep_state/source_data/analysis_results.rds"
R04_RDS   <- "v2_clean/results/R04_rna_atac/source_data/analysis_results.rds"
OUT_DIR   <- "v2_clean/results/R05_atac_landscape"
SRC_DIR   <- file.path(OUT_DIR, "source_data")
CACHE_RDS <- file.path(SRC_DIR, "analysis_results.rds")

run_analysis <- function() {
  for (p in c(INPUT_RDS, R01_RDS, R04_RDS)) if (!file.exists(p)) stop("Missing ", p)
  macs <- Sys.getenv("MACS3_PATH")
  if (!nzchar(macs)) stop("Set MACS3_PATH to the macs3 executable")
  obj <- suppressPackageStartupMessages({ library(Signac); load_multiome(INPUT_RDS) })
  pk_rds <- file.path(SRC_DIR, "celltype_peaks.rds")
  peak_sets <- if (file.exists(pk_rds)) readRDS(pk_rds) else call_celltype_peaks(obj, macs, file.path(tempdir(), "macs3"))
  saveRDS(peak_sets, pk_rds)
  ann <- Signac::Annotation(obj[["peaks"]])
  res <- list(peak_sets = peak_sets, composition = peak_composition(peak_sets, ann),
              umaps = linked_umaps(obj, readRDS(R04_RDS)$umap), tracks = marker_tracks(obj),
              zonation = zonation_heatmaps(obj, readRDS(R01_RDS)$cells))
  write_landscape_source_data(res, SRC_DIR)
  saveRDS(res, CACHE_RDS)
  res
}

load_results <- function() {
  if (file.exists(CACHE_RDS) && !identical(Sys.getenv("FORCE_ANALYSIS"), "1")) return(readRDS(CACHE_RDS))
  dir.create(SRC_DIR, recursive = TRUE, showWarnings = FALSE)
  run_analysis()
}

compose_r05 <- function(res) {
  a <- plot_linked_umaps(res$umaps)
  b <- plot_peaks_vertical(res$composition)
  c <- plot_peaks_horizontal(res$composition)
  d <- plot_marker_tracks(res$tracks)
  e <- plot_zonation_heatmaps(res$zonation)
  top <- (a | (b / c)) + plot_layout(widths = c(80, 103))
  (top / d / e) + plot_layout(heights = c(88, 72, 80)) + plot_annotation(tag_levels = "a") & ggplot2::theme(plot.tag = ggplot2::element_text(size = PT_LABEL, face = "bold", family = FONT))
}

main <- function() {
  res <- load_results()
  save_nature(compose_r05(res), file.path(OUT_DIR, "R05_atac_landscape_composite"), 183, 247)
  save_nature(plot_linked_umaps(res$umaps), file.path(OUT_DIR, "R05a_linked_umaps"), 80, 95)
  save_nature(plot_peaks_vertical(res$composition), file.path(OUT_DIR, "R05b_peaks_per_celltype"), 80, 50)
  save_nature(plot_peaks_horizontal(res$composition), file.path(OUT_DIR, "R05c_peaks_per_celltype_h"), 90, 45)
  save_nature(plot_marker_tracks(res$tracks), file.path(OUT_DIR, "R05d_marker_tracks"), 183, 72)
  save_nature(plot_zonation_heatmaps(res$zonation), file.path(OUT_DIR, "R05e_zonation_heatmaps"), 120, 85)
  print(stats::aggregate(n ~ celltype, res$composition, sum))
  print(vapply(c("peaks", "RNA"), function(m) nrow(res$zonation[[m]]$features), 1L))
  message("Done. Figures in ", OUT_DIR)
}

if (identical(Sys.getenv("R05_ANALYSIS_ONLY"), "1")) {
  tryCatch(invisible(load_results()), error = function(e) { message("R05 analysis failed: ", conditionMessage(e)); quit(status = 1) })
} else {
  tryCatch(main(), error = function(e) {
    message("R05 failed: ", conditionMessage(e))
    quit(status = 1)
  })
}
