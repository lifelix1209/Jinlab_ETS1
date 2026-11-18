# ============================================================
# Figure 3a-c and related supplementary panels, re-analysed on the
# cleaned multiome object and drawn in Nature style.
#
# Input : data/processed/seu_bko_lib1245_singlets.rds  (v2_clean/code/00_preprocessing/02_subset_clean_object.R;
#         libraries 1, 2 control and 4, 5 cKO; scDblFinder singlets;
#         original UMAP; v1 annotation in `celltype`)
# Output: v2_clean/results/Figure3/
#           Fig3_abc.{pdf,svg,png}       composite a-c, 183 mm wide
#           Fig3a/b/c.{pdf,svg,png}      single panels for assembly
#           Fig3c_control_Glul.{pdf,svg,png}  Fig. 3c layout for control Glul+ Hep
#           Fig3c_combined.{pdf,svg,png}      control and cKO side by side per gene
#           FigS3a_candidates.{pdf,svg,tiff}
#           FigS3_markers_v1.{pdf,svg,tiff}
#           source_data/*.csv            values behind every panel
#
# Analyses are cached in source_data/analysis_results.rds; set
# FORCE_ANALYSIS=1 to recompute from the Seurat object.
#
# Usage (from project root):
#   Rscript v2_clean/code/Figure3/fig3_run.R
#   NATURE_FIGURE_DIR=<skill dir> NATURE_QA_PYTHON=<python> Rscript ...
#     (also runs the blocking panel-alignment gate)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source("v2_clean/code/shared/nature_theme.R")
source("v2_clean/code/Figure3/fig3_analysis.R")
source("v2_clean/code/Figure3/fig3_main_panels.R")
source("v2_clean/code/Figure3/fig3_supp_panels.R")

INPUT_RDS <- "data/processed/seu_bko_lib1245_singlets.rds"
OUT_DIR   <- "v2_clean/results/Figure3"
SRC_DIR   <- file.path(OUT_DIR, "source_data")
CACHE_RDS <- file.path(SRC_DIR, "analysis_results.rds")

load_results <- function() {
  if (file.exists(CACHE_RDS) && !identical(Sys.getenv("FORCE_ANALYSIS"), "1")) {
    return(readRDS(CACHE_RDS))
  }
  if (!file.exists(INPUT_RDS)) stop("Missing ", INPUT_RDS, ". Run v2_clean/code/00_preprocessing/02_subset_clean_object.R first.")
  obj <- tryCatch(suppressPackageStartupMessages({ library(Signac); readRDS(INPUT_RDS) }),
                  error = function(e) stop("Failed to read ", INPUT_RDS, ": ", conditionMessage(e)))
  res <- run_fig3_analysis(obj)
  write_source_data(res, SRC_DIR)
  saveRDS(res, CACHE_RDS)
  res
}

mice_per_group <- function(cells) {
  n <- tapply(cells$Library, cells$genotype, function(v) length(unique(v)))
  setNames(as.integer(n), names(n))
}

export_main <- function(res) {
  n_mice <- mice_per_group(res$cells)
  save_nature(compose_fig3(res, n_mice), file.path(OUT_DIR, "Fig3_abc"), 183, 92,
              exemptions = list(list(panels = c("b", "c", "d"), checks = c("panel-width", "horizontal-gutter"),
                                     reason = paste("Row 2 holds a hero UMAP, a narrow genotype split and a",
                                                    "scatter with its own y axis; widths differ by design"))))
  save_nature(plot_fig3a(nrow(res$cells), n_mice), file.path(OUT_DIR, "Fig3a"), 183, 32)
  b <- (plot_fig3b_umap(res$cells) + plot_fig3b_split(res$cells)) +
    patchwork::plot_layout(widths = c(62, 28))
  save_nature(b, file.path(OUT_DIR, "Fig3b"), 96, 60,
              exemptions = list(list(panels = c("a", "b"), checks = "panel-width",
                                     reason = "Main UMAP and genotype split differ in width by design")))
  save_nature(plot_fig3c(res$group), file.path(OUT_DIR, "Fig3c"), 72, 60)
  save_nature(plot_fig3c(res$group, "Control"), file.path(OUT_DIR, "Fig3c_control_Glul"), 72, 60)
  save_nature(plot_fig3c_combined(res$group), file.path(OUT_DIR, "Fig3c_combined"), 125, 60)
}

export_supp <- function(res) {
  save_nature(plot_s3a(res$long, res$stats, res$library), file.path(OUT_DIR, "FigS3a_candidates"),
              183, 82, raster = "tiff")
  save_nature(plot_markers_v1(res$markers), file.path(OUT_DIR, "FigS3_markers_v1"), 150, 62, raster = "tiff")
}

main <- function() {
  res <- load_results()
  export_main(res)
  export_supp(res)
  print(res$stats, digits = 3)
  print(as.data.frame(res$group), digits = 3)
  message("Done. Figures in ", OUT_DIR)
}

tryCatch(main(), error = function(e) {
  message("Figure 3 failed: ", conditionMessage(e))
  quit(status = 1)
})
