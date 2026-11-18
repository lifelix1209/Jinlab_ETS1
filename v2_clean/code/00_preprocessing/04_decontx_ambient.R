# ============================================================
# Ambient-RNA correction of the snRNA counts with DecontX (celda).
#
# Why: endothelial (and other non-parenchymal) nuclei carry hepatocyte
# ambient RNA (Mup*, Saa*, Orm*, Apo*). Hepatocyte expression itself changes
# with genotype, so the ambient load differs between control and cKO and
# shows up as spurious "metabolic" differences in non-hepatocyte contrasts.
#
# Method: celda::decontX on the raw RNA counts of the clean object,
#   z     = figure cell-type labels (v1 + cluster 7/14 corrections),
#   batch = library (contamination profile estimated per library).
# Output (data/processed/):
#   rna_decontx_counts.rds  dgCMatrix, decontaminated counts (rounded)
#   rna_decontx_meta.csv    per-nucleus contamination estimate
# Usage (from project root):
#   Rscript v2_clean/code/00_preprocessing/04_decontx_ambient.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

source("v2_clean/code/shared/nature_theme.R")
source("v2_clean/code/Figure3/fig3_analysis.R")

INPUT_RDS  <- "data/processed/seu_bko_lib1245_singlets.rds"
OUT_COUNTS <- "data/processed/rna_decontx_counts.rds"
OUT_META   <- "data/processed/rna_decontx_meta.csv"
QC_DIR     <- "v2_clean/results/00_preprocessing/decontx"
SEED       <- 42

run_decontx <- function(counts, z, batch) {
  set.seed(SEED)
  res <- celda::decontX(x = counts, z = z, batch = batch, seed = SEED, verbose = TRUE)
  list(counts = round(res$decontXcounts), contamination = res$contamination)
}

summarise_contamination <- function(meta) {
  meta |>
    dplyr::group_by(celltype, genotype) |>
    dplyr::summarise(n = dplyr::n(), median_contamination = stats::median(contamination),
                     q25 = stats::quantile(contamination, 0.25), q75 = stats::quantile(contamination, 0.75),
                     .groups = "drop")
}

main <- function() {
  if (!file.exists(INPUT_RDS)) stop("Missing ", INPUT_RDS)
  obj <- readRDS(INPUT_RDS)
  counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
  ct <- figure_celltype(obj)
  if (anyNA(ct)) stop("Unlabelled nuclei in figure_celltype()")
  res <- run_decontx(counts, z = as.character(ct), batch = paste0("L", obj$Library))
  rownames(res$counts) <- rownames(counts)
  colnames(res$counts) <- colnames(counts)
  meta <- data.frame(barcode = colnames(obj), Library = paste0("L", obj$Library),
                     genotype = genotype_label(obj$Group), celltype = ct,
                     nCount_RNA = obj$nCount_RNA, contamination = res$contamination,
                     nCount_decontX = Matrix::colSums(res$counts), row.names = NULL)
  saveRDS(methods::as(res$counts, "CsparseMatrix"), OUT_COUNTS)
  utils::write.csv(meta, OUT_META, row.names = FALSE)
  dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
  summ <- summarise_contamination(meta)
  utils::write.csv(summ, file.path(QC_DIR, "contamination_by_celltype.csv"), row.names = FALSE)
  print(as.data.frame(summ), digits = 3)
  message("DecontX done: ", OUT_COUNTS)
}

tryCatch(main(), error = function(e) {
  message("DecontX failed: ", conditionMessage(e))
  quit(status = 1)
})
