# ============================================================
# Re-assess the v1 labels "Unknown" (cluster 7) and "Immune cell"
# (cluster 14) on the cleaned object (libraries 1, 2, 4, 5; singlets).
#
# Evidence collected
#   1. Top markers of each cluster vs all other nuclei (Wilcoxon,
#      log-normalized RNA, pct.1 >= 0.2, log2FC >= 1).
#   2. Canonical marker panel (immune lineages, cholangiocyte,
#      mesothelial, parenchymal reference) per cluster.
#   3. RNA/ATAC QC per cluster (low-quality / ambient check).
#   4. Lineage composition inside cluster 7 (B, T/NK, myeloid,
#      plasma) from marker co-detection.
#
# Output: v2_clean/results/00_preprocessing/annotation_reassessment/*.csv (+ printed summary)
# Usage (from project root): Rscript v2_clean/code/00_preprocessing/03_reassess_unknown_immune.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
})

INPUT_RDS <- "data/processed/seu_bko_lib1245_singlets.rds"
OUT_DIR   <- "v2_clean/results/00_preprocessing/annotation_reassessment"
TARGETS   <- c("7", "14")
REFERENCE <- c("0", "3", "5", "6", "12")

PANEL <- list(
  Pan_immune   = c("Ptprc", "Coro1a", "Cd52"),
  B_cell       = c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Ebf1", "Pax5"),
  Plasma       = c("Jchain", "Mzb1", "Xbp1", "Sdc1"),
  T_cell       = c("Cd3e", "Cd3g", "Cd2", "Il7r", "Lck"),
  NK_ILC       = c("Nkg7", "Gzma", "Klrb1c", "Ncr1"),
  Myeloid      = c("Lyz2", "Cd74", "H2-Ab1", "Itgax", "Ccr2", "Csf1r"),
  Neutrophil   = c("S100a8", "S100a9", "Retnlg"),
  pDC          = c("Siglech", "Bst2"),
  Kupffer      = c("Clec4f", "Vsig4", "Marco"),
  Cholangiocyte = c("Krt19", "Krt7", "Epcam", "Sox9", "Spp1", "Hnf1b", "Cftr"),
  Mesothelial  = c("Msln", "Upk3b", "Wt1"),
  Hepatocyte   = c("Alb", "Apoa1", "Serpina1a"),
  Ribosomal    = c("Rps27", "Rpl41", "Rps29")
)

LINEAGE_RULES <- list(
  B_cell  = c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Ebf1", "Pax5"),
  T_NK    = c("Cd3e", "Cd3g", "Cd2", "Nkg7", "Gzma", "Lck"),
  Myeloid = c("Lyz2", "Itgax", "Ccr2", "Csf1r", "S100a8", "S100a9"),
  Plasma  = c("Jchain", "Mzb1")
)

load_object <- function(path) {
  if (!file.exists(path)) stop("Missing ", path, ". Run v2_clean/code/00_preprocessing/02_subset_clean_object.R first.")
  obj <- tryCatch(readRDS(path), error = function(e) stop("Failed to read ", path, ": ", conditionMessage(e)))
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, scale.factor = 1e4, verbose = FALSE)
  Idents(obj) <- "seurat_clusters"
  obj
}

top_markers <- function(obj, cluster, n = 30) {
  FindMarkers(obj, ident.1 = cluster, min.pct = 0.2, logfc.threshold = 1, only.pos = TRUE,
              verbose = FALSE) |>
    tibble::rownames_to_column("gene") |>
    filter(p_val_adj < 0.05) |>
    arrange(desc(avg_log2FC)) |>
    head(n) |>
    mutate(cluster = cluster, .before = 1)
}

panel_table <- function(obj, clusters) {
  genes <- intersect(unlist(PANEL, use.names = FALSE), rownames(obj))
  cnt <- LayerData(obj, layer = "counts")[genes, ]
  cl <- as.character(obj$seurat_clusters)
  lineage <- rep(names(PANEL), lengths(PANEL))[match(genes, unlist(PANEL, use.names = FALSE))]
  do.call(rbind, lapply(clusters, function(k) {
    idx <- which(cl == k)
    data.frame(cluster = k, label_v1 = as.character(obj$celltype[idx[1]]), n = length(idx),
               lineage = lineage, gene = genes, pct = round(100 * Matrix::rowMeans(cnt[, idx] > 0), 1))
  }))
}

qc_table <- function(obj, annot) {
  dbl <- setNames(annot$scDblFinder_score, annot$barcode)
  obj@meta.data |>
    mutate(cluster = as.character(seurat_clusters), dbl_score = dbl[rownames(obj@meta.data)]) |>
    group_by(cluster, label_v1 = as.character(celltype)) |>
    summarise(n = n(), pct_cKO = round(100 * mean(Group == "BCatKO"), 1),
              med_UMI = median(nCount_RNA), med_genes = median(nFeature_RNA),
              med_ATAC = median(nCount_ATAC), med_TSS = round(median(TSS.enrichment), 2),
              med_dbl_score = round(median(dbl_score), 3), .groups = "drop") |>
    arrange(as.integer(cluster))
}

lineage_composition <- function(obj, cluster) {
  cells <- WhichCells(obj, idents = cluster)
  cnt <- LayerData(obj, layer = "counts")[, cells]
  hit <- sapply(LINEAGE_RULES, function(g) Matrix::colSums(cnt[intersect(g, rownames(cnt)), , drop = FALSE] > 0) >= 2)
  n_hit <- rowSums(hit)
  call <- ifelse(n_hit == 0, "None (<2 markers)", ifelse(n_hit > 1, "Mixed", colnames(hit)[max.col(hit)]))
  as.data.frame(table(call = call)) |> mutate(pct = round(100 * Freq / sum(Freq), 1), cluster = cluster)
}

main <- function() {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  obj <- load_object(INPUT_RDS)
  annot <- read.csv("data/processed/cell_annotation_v2.csv")
  res <- list(
    markers = bind_rows(lapply(TARGETS, top_markers, obj = obj)),
    panel = panel_table(obj, c(TARGETS, REFERENCE)),
    qc = qc_table(obj, annot),
    composition = bind_rows(lapply(TARGETS, lineage_composition, obj = obj))
  )
  for (k in names(res)) write.csv(res[[k]], file.path(OUT_DIR, paste0(k, ".csv")), row.names = FALSE)
  options(width = 200)
  print(tibble::as_tibble(res$markers[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2")]), n = 60)
  print(tidyr::pivot_wider(res$panel[, c("cluster", "lineage", "gene", "pct")], names_from = cluster,
                           values_from = pct), n = 80)
  print(res$qc, n = 20)
  print(res$composition)
}

tryCatch(main(), error = function(e) {
  message("Re-assessment failed: ", conditionMessage(e))
  quit(status = 1)
})
