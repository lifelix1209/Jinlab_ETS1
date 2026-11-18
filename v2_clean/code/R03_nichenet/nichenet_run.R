# ============================================================
# R03 — Hepatocyte -> EC communication (NicheNet v2, mouse priors).
#
# Input : data/processed/seu_bko_lib1245_singlets.rds
#         v2_clean/results/R02_ec_response/source_data/analysis_results.rds
#         data/external/nichenet/*.rds  (Zenodo 7074291)
# Output: v2_clean/results/R03_nichenet/  composite, single panels, source_data/
# Usage : Rscript v2_clean/code/R03_nichenet/nichenet_run.R
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

source("v2_clean/code/shared/nature_theme.R")
source("v2_clean/code/Figure3/fig3_analysis.R")
source("v2_clean/code/R03_nichenet/nichenet_analysis.R")
source("v2_clean/code/R03_nichenet/nichenet_plots.R")

INPUT_RDS <- "data/processed/seu_bko_lib1245_singlets.rds"
EC_RDS    <- "v2_clean/results/R02_ec_response/source_data/analysis_results.rds"
OUT_DIR   <- "v2_clean/results/R03_nichenet"
SRC_DIR   <- file.path(OUT_DIR, "source_data")
CACHE_RDS <- file.path(SRC_DIR, "analysis_results.rds")

load_results <- function() {
  if (file.exists(CACHE_RDS) && !identical(Sys.getenv("FORCE_ANALYSIS"), "1")) return(readRDS(CACHE_RDS))
  if (!file.exists(EC_RDS)) stop("Missing ", EC_RDS, ". Run R02 first.")
  ec <- readRDS(EC_RDS)
  pri <- load_priors()
  nc <- suppressPackageStartupMessages({ library(Signac); load_nn_counts(INPUT_RDS) })
  res <- run_nichenet(nc, ec, pri)
  write_nn_source_data(res, SRC_DIR)
  saveRDS(res, CACHE_RDS)
  res
}

compose_r03 <- function(res) {
  n_t <- length(unique(res$targets$target))
  n_r <- length(unique(res$receptors$receptor))
  (plot_ligand_activity(res) | plot_ligand_targets(res) | plot_ligand_receptors(res) |
     plot_sender(res) | plot_niche_prior(res)) +
    plot_layout(widths = c(3.2, max(n_t, 3), max(n_r, 3), 3.2, length(NICHE_GENES))) +
    plot_annotation(tag_levels = "a")
}

main <- function() {
  res <- load_results()
  save_nature(compose_r03(res), file.path(OUT_DIR, "R03_nichenet_composite"), 183, 92,
              exemptions = list(list(panels = letters[1:5], checks = "panel-width",
                                     reason = "Heatmap widths follow the number of columns")))
  save_nature(plot_niche_prior(res), file.path(OUT_DIR, "R03e_niche_prior"), 50, 70)
  print(res$n)
  a <- res$activity
  print(utils::head(a[a$gene_set == "cKO response", ], 12), digits = 3)
  print(a[a$ligand %in% CANDIDATES_NN, ], digits = 3)
  message("Done. Figures in ", OUT_DIR)
}

tryCatch(main(), error = function(e) {
  message("R03 failed: ", conditionMessage(e))
  quit(status = 1)
})
