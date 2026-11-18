# ============================================================
# Clean object: best two libraries per genotype, singlets only.
#
# Library choice (per-library medians, all nuclei):
#   Lib  Group   nuclei  UMI   genes  %mt>5  TSS   nucl.signal
#   1    wt      4620    1788  1176   0.2    3.89  0.88   keep
#   2    wt      3929    1578  1039   0.1    3.74  0.81   keep
#   3    wt       321     677   508   5.9    3.35  1.21   drop (PCA outlier)
#   4    BCatKO  3832    1617  1094   0.1    3.36  1.13   keep
#   5    BCatKO  4459    1339   947   0.5    3.68  0.99   keep
#   6    BCatKO  1276     784   586   0.7    3.34  1.28   drop
#
# Doublets: scDblFinder calls from data/processed/cell_annotation_v2.csv
# (00_preprocessing/01_annotation_doublets.R), so the calls match every downstream script.
#
# Kept unchanged: all assays (RNA, ATAC, peaks, SCT), existing annotations
# (celltype, detailed_celltype, seurat_clusters) plus celltype_v2, and all
# reductions incl. UMAP. No re-normalization, re-clustering or RunUMAP.
# The ATAC/peaks fragment path is repointed to data/raw/atac_fragments.tsv.gz.
#
# Output: data/processed/seu_bko_lib1245_singlets.rds
# Usage (from project root): Rscript v2_clean/code/00_preprocessing/02_subset_clean_object.R
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(ggplot2)
})

source("R/utils_annotation.R")

# ---------------- CONFIG ----------------
INPUT_RDS      <- "data/processed/seu_bko_annotated_with_ATAC_fragments.rds"
OUTPUT_RDS     <- "data/processed/seu_bko_lib1245_singlets.rds"
FRAGMENT_PATH  <- "data/raw/atac_fragments.tsv.gz"
KEEP_LIBRARIES <- c("1", "2", "4", "5")
QC_DIR         <- "v2_clean/results/00_preprocessing/clean_object_qc"
ANNOT_COLS     <- c("celltype_v2", "scDblFinder_class", "scDblFinder_score")

# ---------------- steps ----------------
load_object <- function(path) {
  if (!file.exists(path)) stop("Input not found: ", path)
  tryCatch(readRDS(path), error = function(e) stop("Failed to read ", path, ": ", conditionMessage(e)))
}

add_annotation <- function(obj, annot) {
  rownames(annot) <- annot$barcode
  missing <- setdiff(colnames(obj), annot$barcode)
  if (length(missing) > 0) stop(length(missing), " nuclei lack v2 annotation.")
  AddMetaData(obj, annot[colnames(obj), ANNOT_COLS])
}

select_cells <- function(obj) {
  m <- obj@meta.data
  keep <- as.character(m$Library) %in% KEEP_LIBRARIES & m$scDblFinder_class == "singlet"
  rownames(m)[keep]
}

repoint_fragments <- function(obj, path) {
  path <- normalizePath(path, mustWork = TRUE)
  if (!file.exists(paste0(path, ".tbi"))) stop("Missing tabix index: ", path, ".tbi")
  chrom_assays <- Filter(function(a) inherits(obj[[a]], "ChromatinAssay"), Assays(obj))
  Reduce(function(o, a) {
    frags <- Fragments(o[[a]])
    if (length(frags) == 0) return(o)
    # UpdatePath() refuses because the .tbi was rebuilt (MD5 differs), so
    # rebuild the fragment object keeping the original cell-name mapping.
    new_frags <- lapply(frags, function(f) {
      CreateFragmentObject(path = path, cells = GetFragmentData(f, "cells"), validate.fragments = FALSE)
    })
    Fragments(o[[a]]) <- NULL
    Fragments(o[[a]]) <- new_frags
    o
  }, chrom_assays, obj)
}

# Read fragments in a gene-dense region and check barcodes match the object.
check_fragment_barcodes <- function(obj, path, region = "chr7:30000000-31000000") {
  hits <- Rsamtools::scanTabix(path, param = GenomicRanges::GRanges(region))[[1]]
  if (length(hits) == 0) stop("No fragments read from ", region)
  bc <- unique(vapply(strsplit(hits, "\t", fixed = TRUE), `[`, character(1), 4))
  covered <- mean(colnames(obj) %in% bc)
  message(sprintf("Fragment check %s: %d fragments, %.1f%% of kept nuclei seen", region,
                  length(hits), 100 * covered))
  if (covered < 0.5) stop("Fragment barcodes do not match object cells (", round(100 * covered, 1), "%).")
  invisible(covered)
}

verify_unchanged <- function(orig, clean) {
  cells <- colnames(clean)
  same_umap <- identical(Embeddings(orig, "umap")[cells, ], Embeddings(clean, "umap"))
  cols <- c("celltype", "detailed_celltype", "seurat_clusters", "celltype_v2")
  same_annot <- all(vapply(cols, function(k) identical(orig@meta.data[cells, k], clean@meta.data[cells, k]),
                           logical(1)))
  same_atac <- identical(dim(clean[["ATAC"]]), c(nrow(orig[["ATAC"]]), length(cells)))
  if (!same_umap) stop("UMAP coordinates changed after subsetting.")
  if (!same_annot) stop("Annotations changed after subsetting.")
  if (!same_atac) stop("ATAC assay dimensions unexpected after subsetting.")
  invisible(TRUE)
}

summarise_cells <- function(orig, clean) {
  m <- orig@meta.data
  data.frame(
    Library = sort(unique(as.character(m$Library))),
    stringsAsFactors = FALSE
  ) |>
    transform(
      Group    = vapply(Library, function(l) as.character(m$Group[m$Library == l][1]), character(1)),
      n_all    = vapply(Library, function(l) sum(m$Library == l), integer(1)),
      n_dbl    = vapply(Library, function(l) sum(m$Library == l & m$scDblFinder_class == "doublet"), integer(1)),
      n_kept   = vapply(Library, function(l) sum(clean$Library == l), integer(1))
    )
}

plot_umap_check <- function(clean) {
  emb <- Embeddings(clean, "umap")
  d <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2], celltype = clean$celltype_v2,
                  Library = paste0("L", clean$Library))
  ggplot(d, aes(UMAP1, UMAP2, color = celltype)) +
    geom_point(size = 0.2, stroke = 0) +
    facet_wrap(~Library, nrow = 1) +
    guides(color = guide_legend(override.aes = list(size = 2))) +
    coord_fixed() +
    theme_classic(base_size = 8) +
    labs(title = "Libraries 1, 2 (WT) and 4, 5 (KO), singlets; original UMAP coordinates")
}

# ---------------- run ----------------
main <- function() {
  dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)
  orig <- add_annotation(load_object(INPUT_RDS), read_annotation_v2())
  clean <- subset(orig, cells = select_cells(orig)) |> repoint_fragments(FRAGMENT_PATH)
  verify_unchanged(orig, clean)
  check_fragment_barcodes(clean, FRAGMENT_PATH)

  summary_tbl <- summarise_cells(orig, clean)
  write.csv(summary_tbl, file.path(QC_DIR, "cells_per_library.csv"), row.names = FALSE)
  ggsave(file.path(QC_DIR, "UMAP_by_library_check.png"), plot_umap_check(clean),
         width = 10, height = 3.2, dpi = 200)

  tryCatch(saveRDS(clean, OUTPUT_RDS),
           error = function(e) stop("Failed to save ", OUTPUT_RDS, ": ", conditionMessage(e)))
  print(summary_tbl)
  print(table(clean$celltype_v2, clean$Group))
  message("Fragments: ", Fragments(clean[["ATAC"]])[[1]]@path)
  message("Saved ", ncol(clean), " nuclei to ", OUTPUT_RDS)
}

tryCatch(main(), error = function(e) {
  message("Cleaning failed: ", conditionMessage(e))
  quit(status = 1)
})
