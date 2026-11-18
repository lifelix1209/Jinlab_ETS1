# ============================================================
# Corrected cell-type annotation (v2), doublet calling, Figure 3B
#
# Annotation corrections (marker evidence, % of nuclei detecting):
#   cluster 7  "Unknown"        -> Immune cells (Ptprc 38%, Cd79a 14%, Cd3e 6%)
#   cluster 9  "Blood cell"     -> Erythroid (Hba-a1 90%, Hbb-bs 93%)
#   cluster 10 "eGFP+ Hep"      -> Cycling hepatocytes (Top2a 47%, Mki67 20%)
#   cluster 12 "Stellate"       -> Mesothelial cells (Msln 58%, Upk3b 44%)
#   cluster 13 "eGFP+ Hep"      -> IFN-responsive hepatocytes (Rsad2 62%, Isg15 51%)
#   cluster 14 "Immune cell"    -> Cholangiocytes (Spp1 90%, Epcam 22%, Ptprc 2%)
# Clusters are unchanged (WNN, resolution 0.2); only labels change.
#
# Doublets: scDblFinder on RNA counts, per library (samples = Library),
# computed once and cached in ANNOT_CSV (set FORCE_DOUBLETS=1 to recompute),
# cluster-aware (clusters = seurat_clusters). Doublets are flagged here and
# excluded in downstream comparisons; clustering/UMAP are not recomputed.
#
# Outputs
#   data/processed/cell_annotation_v2.csv   (one row per nucleus)
#   v2_clean/results/00_preprocessing/annotation_doublets/Fig3B_UMAP_celltype_v2.pdf
#   v2_clean/results/00_preprocessing/annotation_doublets/FigS_annotation_markers_v2.pdf
#   v2_clean/results/00_preprocessing/annotation_doublets/FigS_doublets_by_cluster.pdf
#
# Usage (from project root):
#   Rscript v2_clean/code/00_preprocessing/01_annotation_doublets.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(dplyr)
  library(ggplot2)
})

source("R/utils_fonts.R")
register_arial()

# ---------------- CONFIG ----------------
INPUT_RDS  <- "data/processed/BCatKO_multiome_clustered.rds"
ANNOT_CSV  <- "data/processed/cell_annotation_v2.csv"
FIG_DIR    <- "v2_clean/results/00_preprocessing/annotation_doublets"
SEED       <- 42

CLUSTER_TO_CELLTYPE_V2 <- c(
  "0"  = "Hepatocytes",
  "1"  = "Hepatocytes",
  "2"  = "Hepatocytes",
  "11" = "Hepatocytes",
  "4"  = "Pericentral Glul+ hepatocytes",
  "8"  = "Periportal hepatocytes",
  "10" = "Cycling hepatocytes",
  "13" = "IFN-responsive hepatocytes",
  "3"  = "Endothelial cells",
  "5"  = "Kupffer cells",
  "6"  = "Stellate cells",
  "12" = "Mesothelial cells",
  "7"  = "Immune cells",
  "9"  = "Erythroid cells",
  "14" = "Cholangiocytes"
)
# Labels used by scripts 06/07 (v1), kept for traceability.
CLUSTER_TO_CELLTYPE_V1 <- c(
  "0" = "Hepatocyte", "1" = "Hepatocyte", "2" = "Hepatocyte", "11" = "Hepatocyte",
  "3" = "Endothelial", "4" = "PericentralHep_GLUL+", "5" = "Kupffer cell",
  "6" = "Stellate cell", "12" = "Stellate cell", "7" = "Unknown", "8" = "PeriportalHep",
  "9" = "Blood cell", "10" = "PericentralHep_eGFP+", "13" = "PericentralHep_eGFP+",
  "14" = "Immune cell"
)
CELLTYPE_ORDER <- unique(CLUSTER_TO_CELLTYPE_V2[c("0", "4", "8", "10", "13", "3", "5", "6",
                                                  "12", "7", "9", "14")])
CELLTYPE_COLS <- setNames(
  c("#AFC3DC", "#2F4B7C", "#8FB3D9", "#E1604D", "#F2A541",
    "#00A087", "#7E6148", "#B09C85", "#91D1C2", "#8491B4", "#DC0000", "#4DBBD5"),
  CELLTYPE_ORDER)

MARKERS <- c("Alb", "Glul", "Cyp2e1", "Cyp2f2", "Sds", "Top2a", "Mki67", "Rsad2", "Isg15",
             "Stab2", "Kdr", "Clec4f", "Vsig4", "Reln", "Lrat", "Msln", "Upk3b",
             "Ptprc", "Cd79a", "Hba-a1", "Hbb-bs", "Spp1", "Epcam")

# ---------------- helpers ----------------
load_object <- function(path) {
  if (!file.exists(path)) stop("Input not found: ", path, ". Run from project root.")
  tryCatch(readRDS(path), error = function(e) stop("Failed to read ", path, ": ", conditionMessage(e)))
}

theme_fig <- function(base = 7) {
  theme_classic(base_size = base, base_family = FONT_FAMILY) +
    theme(text = element_text(color = "black"), axis.text = element_text(color = "black"),
          plot.title = element_text(face = "bold", size = base + 1.5),
          plot.subtitle = element_text(size = base - 0.5, color = "grey35"),
          legend.key.size = unit(3, "mm"))
}

save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = grDevices::pdf)
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

call_doublets <- function(obj) {
  sce <- SingleCellExperiment(assays = list(counts = LayerData(obj, assay = "RNA", layer = "counts")))
  sce$Library <- as.character(obj$Library)
  sce$cluster <- as.character(obj$seurat_clusters)
  set.seed(SEED)
  sce <- tryCatch(
    scDblFinder(sce, samples = "Library", clusters = "cluster", BPPARAM = BiocParallel::SerialParam(RNGseed = SEED)),
    error = function(e) stop("scDblFinder failed: ", conditionMessage(e)))
  data.frame(barcode = colnames(sce), scDblFinder_class = sce$scDblFinder.class,
             scDblFinder_score = sce$scDblFinder.score)
}

build_annotation <- function(obj, dbl) {
  m  <- obj@meta.data
  cl <- as.character(m$seurat_clusters)
  mt <- as.numeric(unlist(PercentageFeatureSet(obj, pattern = "^mt-", assay = "RNA")))
  data.frame(barcode = rownames(m), Library = as.character(m$Library), Group = as.character(m$Group),
             seurat_clusters = cl,
             celltype_v1 = unname(CLUSTER_TO_CELLTYPE_V1[cl]),
             celltype_v2 = unname(CLUSTER_TO_CELLTYPE_V2[cl]),
             nCount_RNA = m$nCount_RNA, nFeature_RNA = m$nFeature_RNA, percent_mt = mt) |>
    left_join(dbl, by = "barcode")
}

# ---------------- plots ----------------
plot_umap <- function(obj, annot) {
  emb <- Embeddings(obj, "umap")
  d <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2],
                  celltype = factor(annot$celltype_v2, levels = CELLTYPE_ORDER))
  n <- table(d$celltype)
  lab <- setNames(paste0(names(n), " (", trimws(format(as.integer(n), big.mark = ",")), ")"), names(n))
  set.seed(SEED)
  d <- d[sample(nrow(d)), ]
  ggplot(d, aes(UMAP1, UMAP2, color = celltype)) +
    geom_point(size = 0.35, stroke = 0, alpha = 0.9) +
    scale_color_manual(values = CELLTYPE_COLS, labels = lab, name = NULL) +
    guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
    labs(title = "Liver cell populations (snRNA + snATAC, WNN)",
         subtitle = paste0(format(nrow(d), big.mark = ","), " nuclei; 3 control and 3 β-cat cKO livers")) +
    coord_fixed() +
    theme_fig() +
    theme(axis.text = element_blank(), axis.ticks = element_blank())
}

plot_markers <- function(obj, annot) {
  cnt <- LayerData(obj, assay = "RNA", layer = "counts")[MARKERS, ]
  ln <- log1p(t(t(as.matrix(cnt)) / obj$nCount_RNA * 1e4))
  ct <- factor(annot$celltype_v2, levels = CELLTYPE_ORDER)
  d <- do.call(rbind, lapply(levels(ct), function(k) {
    idx <- which(ct == k)
    data.frame(celltype = k, gene = MARKERS, pct = 100 * rowMeans(cnt[, idx] > 0),
               mean = rowMeans(ln[, idx]))
  })) |>
    group_by(gene) |> mutate(scaled = (mean - min(mean)) / max(max(mean) - min(mean), 1e-9)) |> ungroup() |>
    mutate(gene = factor(gene, levels = MARKERS), celltype = factor(celltype, levels = rev(CELLTYPE_ORDER)))
  ggplot(d, aes(gene, celltype)) +
    geom_point(aes(size = pct, color = scaled)) +
    scale_size_area(max_size = 4, name = "% nuclei") +
    scale_color_gradient(low = "#E8ECF2", high = "#2F4B7C", name = "Scaled\nmean") +
    labs(x = NULL, y = NULL, title = "Marker genes supporting the v2 annotation") +
    theme_fig() +
    theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1))
}

plot_doublets <- function(annot) {
  d <- annot |>
    mutate(celltype = factor(celltype_v2, levels = CELLTYPE_ORDER)) |>
    group_by(celltype) |>
    summarise(n = n(), pct_doublet = 100 * mean(scDblFinder_class == "doublet"), .groups = "drop")
  ggplot(d, aes(pct_doublet, celltype)) +
    geom_col(fill = "#8C8C8C", width = 0.7) +
    geom_text(aes(label = paste0(round(pct_doublet, 1), "% (n=", n, ")")), hjust = -0.1, size = 2.2,
              family = FONT_FAMILY) +
    scale_y_discrete(limits = rev(CELLTYPE_ORDER)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(x = "Predicted doublets (%)", y = NULL, title = "scDblFinder doublet calls by cell type",
         subtitle = "Per-library, cluster-aware; doublets excluded from downstream comparisons") +
    theme_fig()
}

# ---------------- run ----------------
main <- function() {
  dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
  obj <- load_object(INPUT_RDS)
  unmapped <- setdiff(unique(as.character(obj$seurat_clusters)), names(CLUSTER_TO_CELLTYPE_V2))
  if (length(unmapped) > 0) stop("Clusters without v2 label: ", paste(unmapped, collapse = ", "))

  # Reuse existing doublet calls unless FORCE_DOUBLETS=1, so downstream
  # scripts always see the same calls.
  force <- identical(Sys.getenv("FORCE_DOUBLETS"), "1")
  dbl <- if (file.exists(ANNOT_CSV) && !force) {
    read.csv(ANNOT_CSV)[, c("barcode", "scDblFinder_class", "scDblFinder_score")]
  } else {
    call_doublets(obj)
  }
  annot <- build_annotation(obj, dbl)
  if (anyNA(annot$scDblFinder_class)) stop("Doublet calls missing for some nuclei.")
  write.csv(annot, ANNOT_CSV, row.names = FALSE)

  save_fig(plot_umap(obj, annot), "Fig3B_UMAP_celltype_v2", 5.6, 3.8)
  save_fig(plot_markers(obj, annot), "FigS_annotation_markers_v2", 6.4, 3.4)
  save_fig(plot_doublets(annot), "FigS_doublets_by_cluster", 4.6, 3.2)

  print(table(annot$celltype_v2, annot$scDblFinder_class))
  message("Done. Annotation: ", ANNOT_CSV, "; figures in ", FIG_DIR)
}

tryCatch(main(), error = function(e) {
  message("Annotation v2 failed: ", conditionMessage(e))
  quit(status = 1)
})
