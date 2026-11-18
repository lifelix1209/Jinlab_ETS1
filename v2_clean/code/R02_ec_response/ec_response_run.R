# ============================================================
# R02 — Endothelial transcriptional response to pericentral beta-cat
# deletion, with hepatocyte ambient-RNA correction.
#
# Input : data/processed/seu_bko_lib1245_singlets.rds
#         data/processed/rna_decontx_counts.rds, rna_decontx_meta.csv
#           (v2_clean/code/00_preprocessing/04_decontx_ambient.R)
# Output: v2_clean/results/R02_ec_response/
#   R02_ec_response_composite.{pdf,svg,png}, single panels, source_data/
# Usage : Rscript v2_clean/code/R02_ec_response/ec_response_run.R
#         (FORCE_ANALYSIS=1 to recompute the cached analysis)
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source("v2_clean/code/shared/nature_theme.R")
source("v2_clean/code/Figure3/fig3_analysis.R")
source("v2_clean/code/R02_ec_response/ec_response_analysis.R")
source("v2_clean/code/R02_ec_response/ec_response_plots.R")

INPUT_RDS <- "data/processed/seu_bko_lib1245_singlets.rds"
DX_COUNTS <- "data/processed/rna_decontx_counts.rds"
DX_META   <- "data/processed/rna_decontx_meta.csv"
OUT_DIR   <- "v2_clean/results/R02_ec_response"
SRC_DIR   <- file.path(OUT_DIR, "source_data")
CACHE_RDS <- file.path(SRC_DIR, "analysis_results.rds")

load_results <- function() {
  if (file.exists(CACHE_RDS) && !identical(Sys.getenv("FORCE_ANALYSIS"), "1")) return(readRDS(CACHE_RDS))
  ec <- suppressPackageStartupMessages({ library(Signac); load_ec_counts(INPUT_RDS, DX_COUNTS) })
  res <- run_ec_response(ec)
  write_ec_source_data(res, SRC_DIR)
  saveRDS(res, CACHE_RDS)
  res
}

load_dx_meta <- function() {
  if (!file.exists(DX_META)) stop("Missing ", DX_META)
  m <- utils::read.csv(DX_META)
  m$genotype <- factor(m$genotype, levels = names(GROUP_COLS))
  m
}

compose_r02 <- function(res, dx_meta) {
  a <- plot_contamination(dx_meta)
  b <- plot_naive_volcano(res$naive)
  c <- plot_gsea_compare(res$gsea)
  d <- plot_corrected_volcano(res$corrected)
  e <- plot_library_genes(res$library)
  left <- (a / b / d) + plot_layout(heights = c(40, 48, 48))
  top <- (left | c) + plot_layout(widths = c(72, 111))
  (top / e) + plot_layout(heights = c(136, 34)) + plot_annotation(tag_levels = "a")
}

export_all <- function(res, dx_meta) {
  save_nature(compose_r02(res, dx_meta), file.path(OUT_DIR, "R02_ec_response_composite"), 183, 175)
  save_nature(plot_contamination(dx_meta), file.path(OUT_DIR, "R02a_ambient_contamination"), 62, 45)
  save_nature(plot_naive_volcano(res$naive), file.path(OUT_DIR, "R02b_naive_volcano"), 62, 55)
  save_nature(plot_gsea_compare(res$gsea), file.path(OUT_DIR, "R02c_gsea_compare"), 121, 120)
  save_nature(plot_corrected_volcano(res$corrected), file.path(OUT_DIR, "R02d_corrected_volcano"), 70, 60)
  save_nature(plot_library_genes(res$library), file.path(OUT_DIR, "R02e_library_genes"), 183, 38)
}

main <- function() {
  res <- load_results()
  export_all(res, load_dx_meta())
  print(res$n)
  print(utils::head(res$corrected, 15), digits = 3)
  g <- res$gsea
  print(utils::head(g[g$analysis == "Corrected nucleus-level", c("pathway", "NES", "p_BH")], 15), digits = 3)
  message("Done. Figures in ", OUT_DIR)
}

tryCatch(main(), error = function(e) {
  message("R02 failed: ", conditionMessage(e))
  quit(status = 1)
})
