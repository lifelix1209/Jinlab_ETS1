# ============================================================
# Wnt2+/- and Wnt9b+/- endothelial comparisons and TF intersection
# (Figure 4B-C / Supplementary Figure S4A-B), re-run on the cleaned
# multiome object with a standard Seurat workflow.
#
# Input : data/processed/seu_bko_lib1245_singlets.rds
#         (v2_clean/code/00_preprocessing/02_subset_clean_object.R)
# Method: nominal Wilcoxon P < 0.1 (no multiple-testing correction for
#         selection); see header of v2_clean/code/Figure4/fig4_analysis.R. The same
#         thresholds are applied to Wnt2 and Wnt9b.
# Output: v2_clean/results/Figure4/
#           Fig4B_*, Fig4C_*, FigS4A_*, FigS4B_*, FigS4_*  one file per panel
#             (.pdf editable, .svg, .png preview); see fig4_plots.R header
#           tables/*.csv, Wnt_EC_TF_clean.xlsx, sessionInfo.txt
#
# Analyses are cached in tables/analysis_results.rds; set FORCE_ANALYSIS=1
# to recompute from the Seurat object.
# Usage (from project root): Rscript v2_clean/code/Figure4/fig4_run.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
  library(patchwork)
})

for (f in c("v2_clean/code/shared/nature_theme.R", "v2_clean/code/Figure3/fig3_analysis.R",
            "v2_clean/code/Figure4/fig4_analysis.R", "v2_clean/code/Figure4/fig4_plots.R")) {
  if (!file.exists(f)) stop("Helper not found: ", f, ". Run from the project root.")
  source(f)
}

INPUT_RDS <- "data/processed/seu_bko_lib1245_singlets.rds"
OUT_DIR   <- "v2_clean/results/Figure4"
TAB_DIR   <- file.path(OUT_DIR, "tables")
CACHE_RDS <- file.path(TAB_DIR, "analysis_results.rds")

load_results <- function() {
  if (file.exists(CACHE_RDS) && !identical(Sys.getenv("FORCE_ANALYSIS"), "1")) {
    res <- readRDS(CACHE_RDS)
    res$rules <- tf_rule_sensitivity(res$tfs)
    return(res)
  }
  message("Loading endothelial nuclei from ", INPUT_RDS)
  res <- run_wnt_tf_analysis(load_ec(INPUT_RDS))
  dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
  saveRDS(res, CACHE_RDS)
  res
}

write_tables <- function(res) {
  dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
  sheets <- list(group_sizes = res$groups,
                 DE_Wnt2pos_vs_neg = res$de$Wnt2, DE_Wnt9bpos_vs_neg = res$de$Wnt9b,
                 locus_peaks = do.call(rbind, res$peaks),
                 TF_Wnt2 = res$tfs$Wnt2, TF_Wnt9b = res$tfs$Wnt9b, TF_shared = res$shared, TF_by_P_rule = res$rules,
                 focal_TFs = res$focal,
                 motif_enrich_Wnt2 = res$enrichment$Wnt2, motif_enrich_Wnt9b = res$enrichment$Wnt9b)
  for (nm in names(sheets)) utils::write.csv(sheets[[nm]], file.path(TAB_DIR, paste0(nm, ".csv")), row.names = FALSE)
  wb <- openxlsx::createWorkbook()
  for (nm in names(sheets)) {
    openxlsx::addWorksheet(wb, nm)
    openxlsx::writeData(wb, nm, sheets[[nm]])
  }
  openxlsx::saveWorkbook(wb, file.path(OUT_DIR, "Wnt_EC_TF_clean.xlsx"), overwrite = TRUE)
  writeLines(utils::capture.output(utils::sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
}

print_summary <- function(res) {
  print(res$groups)
  for (l in LIGANDS) {
    t <- res$tfs[[l]]
    message(sprintf("%s: %d up genes (P < %.2f, log2FC > 0); %d locus peaks; %d motif-present TFs; %d TFs up",
                    l, sum(res$de[[l]]$up), P_CUT, nrow(res$peaks[[l]]), nrow(t), sum(t$up)))
    message("  TFs: ", paste(t$gene[t$up], collapse = ", "))
  }
  print(res$rules[, c("rule", "n_Wnt2", "n_Wnt9b", "n_shared", "shared_TFs")])
  message("Shared: ", if (nrow(res$shared)) paste(res$shared$gene, collapse = ", ") else "none")
  print(res$focal[, c("ligand", "gene", "avg_log2FC", "pct_pos", "pct_neg", "p_val", "p_val_adj_bonferroni",
                      "p_val_adj_BH")], digits = 3)
}

main <- function() {
  res <- load_results()
  write_tables(res)
  print_summary(res)
  if (!identical(Sys.getenv("SKIP_FIGURES"), "1")) export_figures(res, OUT_DIR)
  message("Done. Output in ", OUT_DIR)
}

tryCatch(main(), error = function(e) {
  message("ERROR in 21-wnt_ec_tf_clean.R: ", conditionMessage(e))
  quit(status = 1)
})
