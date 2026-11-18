# ============================================================
# R01 — Hepatocyte identity loss and signalling-state redistribution map.
#
# Input : data/processed/seu_bko_lib1245_singlets.rds
#         v2_clean/results/Figure3/source_data/analysis_results.rds (bubble plot)
# Output: v2_clean/results/R01_hep_state/
#   R01_hep_state_composite.{pdf,svg,png}   183 mm composite (a-g)
#   R01a_state_map ... R01f_zone_residuals  single panels
#   source_data/*.csv                       values behind every panel
#
# Analyses are cached in source_data/analysis_results.rds; set
# FORCE_ANALYSIS=1 to recompute.
# Usage (from project root):
#   Rscript v2_clean/code/R01_hep_state/hep_state_run.R
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

source("v2_clean/code/shared/nature_theme.R")
source("v2_clean/code/Figure3/fig3_analysis.R")
source("v2_clean/code/Figure3/fig3_main_panels.R")
BUBBLE_MAX_MM <- 4.6  # slightly smaller bubbles than Fig. 3c so labels clear neighbours in the composite
source("v2_clean/code/R01_hep_state/hep_state_analysis.R")
source("v2_clean/code/R01_hep_state/hep_state_plots.R")

INPUT_RDS <- "data/processed/seu_bko_lib1245_singlets.rds"
FIG3_RDS  <- "v2_clean/results/Figure3/source_data/analysis_results.rds"
OUT_DIR   <- "v2_clean/results/R01_hep_state"
SRC_DIR   <- file.path(OUT_DIR, "source_data")
CACHE_RDS <- file.path(SRC_DIR, "analysis_results.rds")

load_results <- function() {
  if (file.exists(CACHE_RDS) && !identical(Sys.getenv("FORCE_ANALYSIS"), "1")) return(readRDS(CACHE_RDS))
  h <- suppressPackageStartupMessages({ library(Signac); load_hepatocytes(INPUT_RDS) })
  res <- run_hep_state_analysis(h)
  write_hep_source_data(res, SRC_DIR)
  saveRDS(res, CACHE_RDS)
  res
}

load_bubble_data <- function() {
  if (!file.exists(FIG3_RDS)) stop("Missing ", FIG3_RDS, ". Run v2_clean/code/Figure3/fig3_run.R first.")
  readRDS(FIG3_RDS)$group
}

compose_r01 <- function(res, bubble) {
  a <- plot_state_map(res$cells)
  b <- plot_factor_map(res$cells)
  c <- plot_zone_occupancy(res$occupancy)
  d <- plot_zonation_heatmap(res$heatmap)
  e <- plot_factor_profiles(res$profiles, res$reference)
  f <- plot_fig3c_combined(bubble)
  g <- plot_zone_residuals(res$residuals)
  row1 <- (a | b | c) + plot_layout(widths = c(64, 64, 40))
  row2 <- (d | e) + plot_layout(widths = c(80, 100))
  row3 <- (f | g) + plot_layout(widths = c(128, 55))
  (row1 / row2 / row3) +
    plot_layout(heights = c(52, 62, 48)) +
    plot_annotation(tag_levels = "a")
}

export_all <- function(res, bubble) {
  save_nature(compose_r01(res, bubble), file.path(OUT_DIR, "R01_hep_state_composite"), 183, 182)
  save_nature(plot_state_map(res$cells), file.path(OUT_DIR, "R01a_state_map"), 80, 52)
  save_nature(plot_factor_map(res$cells), file.path(OUT_DIR, "R01b_factor_map"), 80, 48)
  save_nature(plot_zone_occupancy(res$occupancy), file.path(OUT_DIR, "R01c_zone_occupancy"), 50, 45)
  save_nature(plot_zonation_heatmap(res$heatmap), file.path(OUT_DIR, "R01d_zonation_heatmap"), 90, 55)
  save_nature(plot_factor_profiles(res$profiles, res$reference), file.path(OUT_DIR, "R01e_factor_profiles"), 100, 50)
  save_nature(plot_zone_residuals(res$residuals), file.path(OUT_DIR, "R01f_zone_residuals"), 65, 55)
}

main <- function() {
  res <- load_results()
  bubble <- load_bubble_data()
  export_all(res, bubble)
  print(res$tests, digits = 3)
  print(as.data.frame(res$residuals), digits = 3)
  message("Done. Figures in ", OUT_DIR)
}

tryCatch(main(), error = function(e) {
  message("R01 failed: ", conditionMessage(e))
  quit(status = 1)
})
