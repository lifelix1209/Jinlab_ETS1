# ============================================================
# Figure 3 / S3 analyses on the cleaned multiome object
# (libraries 1, 2, 4, 5; scDblFinder singlets; v1 annotation).
#
# Same analysis as Methods (manuscript V13):
#   Groups: control nuclei in the Glul+ Hep cluster vs beta-cat cKO nuclei
#           in the eGFP+ Hep clusters (genotype AND v1 cluster label).
#   Fig. 3C: in the cKO eGFP+ group, mean UMI count per nucleus (zeros
#           included) and % nuclei with >= 1 UMI for each candidate.
#   Fig. S3A: LogNormalize (scale 1e4, ln(x + 1)); two-sided Wilcoxon
#           rank-sum test, normal approximation, per gene (uncorrected P,
#           BH-adjusted P also reported in the source data).
#   Markers: % nuclei detecting and mean ln-normalized expression per
#           v1 cell type.
# ============================================================

CANDIDATES <- c("Igf1", "Nrg4", "Vegfa", "Nampt", "Fgl1", "Cxcl1", "Cxcl10", "Angpt2")
VALIDATED  <- c("Igf1", "Nrg4", "Nampt", "Vegfa")
MARKERS_V1 <- c("Alb", "Cyp2e1", "Glul", "Oat", "eGFP", "Top2a", "Mki67", "Cyp2f2", "Sds",
                "Spp1", "Epcam", "Sox9", "Stab2", "Kdr", "Clec4f", "Vsig4", "Reln", "Lrat", "Msln",
                "Ptprc", "Cd52", "Ebf1", "Cd79a", "Hba-a1", "Hbb-bs")
GROUP_LEVELS <- c("Control", "β-cat cKO")
SCALE_FACTOR <- 1e4

# v1 labels with the marker-based cluster corrections (CLUSTER_RELABEL).
figure_celltype <- function(obj) {
  lab <- as.character(obj$celltype)
  cl <- as.character(obj$seurat_clusters)
  hit <- cl %in% names(CLUSTER_RELABEL)
  lab[hit] <- unname(CLUSTER_RELABEL[cl[hit]])
  factor(lab, levels = CELLTYPE_V1_ORDER)
}

genotype_label <- function(group) {
  out <- c(wt = "Control", BCatKO = "β-cat cKO")[as.character(group)]
  if (anyNA(out)) stop("Unexpected Group values: ", paste(unique(group[is.na(out)]), collapse = ", "))
  factor(unname(out), levels = GROUP_LEVELS)
}

# Per-nucleus table for UMAP panels.
cell_table <- function(obj) {
  emb <- SeuratObject::Embeddings(obj, "umap")
  egfp <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")["eGFP", colnames(obj)]
  data.frame(barcode = colnames(obj), UMAP1 = emb[, 1], UMAP2 = emb[, 2],
             celltype = figure_celltype(obj), seurat_clusters = as.character(obj$seurat_clusters),
             genotype = genotype_label(obj$Group), Library = paste0("L", obj$Library),
             eGFP_umi = as.numeric(egfp), row.names = NULL)
}

# Log-normalized candidate expression in the two pericentral groups.
pericentral_long <- function(obj) {
  geno <- genotype_label(obj$Group)
  ct <- as.character(obj$celltype)
  grp <- ifelse(geno == "Control" & ct == "Glul+ Hep", "Control",
                ifelse(geno == "β-cat cKO" & ct == "eGFP+ Hep", "β-cat cKO", NA))
  keep <- colnames(obj)[!is.na(grp)]
  counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")[, keep]
  lib_size <- Matrix::colSums(counts)
  raw <- as.matrix(counts[CANDIDATES, , drop = FALSE])
  ln <- log1p(t(t(raw) / lib_size * SCALE_FACTOR))
  data.frame(barcode = rep(keep, each = length(CANDIDATES)),
             gene = factor(rep(CANDIDATES, times = length(keep)), levels = CANDIDATES),
             group = factor(rep(grp[!is.na(grp)], each = length(CANDIDATES)), levels = GROUP_LEVELS),
             Library = rep(paste0("L", obj$Library[keep]), each = length(CANDIDATES)),
             umi = as.vector(raw), expr = as.vector(ln))
}

wilcoxon_table <- function(long) {
  res <- lapply(CANDIDATES, function(g) {
    d <- long[long$gene == g, ]
    x <- d$expr[d$group == "β-cat cKO"]
    y <- d$expr[d$group == "Control"]
    t <- stats::wilcox.test(x, y, alternative = "two.sided", exact = FALSE)
    data.frame(gene = g, n_control = length(y), n_cko = length(x),
               mean_control = mean(y), mean_cko = mean(x), W = unname(t$statistic), p = t$p.value)
  })
  out <- do.call(rbind, res)
  out$p_BH <- stats::p.adjust(out$p, "BH")
  out
}

group_summary <- function(long) {
  long |>
    dplyr::group_by(gene, group) |>
    dplyr::summarise(n = dplyr::n(), mean_umi = mean(umi), pct = 100 * mean(umi > 0),
                     mean_expr = mean(expr), .groups = "drop") |>
    dplyr::mutate(validated = gene %in% VALIDATED)
}

library_summary <- function(long) {
  long |>
    dplyr::group_by(gene, group, Library) |>
    dplyr::summarise(n = dplyr::n(), mean_expr = mean(expr), .groups = "drop")
}

marker_summary <- function(obj) {
  missing <- setdiff(MARKERS_V1, rownames(obj[["RNA"]]))
  if (length(missing) > 0) stop("Markers absent from RNA assay: ", paste(missing, collapse = ", "))
  counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
  ln <- log1p(t(t(as.matrix(counts[MARKERS_V1, ])) / Matrix::colSums(counts) * SCALE_FACTOR))
  ct <- figure_celltype(obj)
  res <- lapply(levels(ct), function(k) {
    idx <- which(ct == k)
    data.frame(celltype = k, gene = MARKERS_V1, n = length(idx),
               pct = 100 * rowMeans(ln[, idx, drop = FALSE] > 0), mean = rowMeans(ln[, idx, drop = FALSE]))
  })
  out <- do.call(rbind, res)
  rng <- tapply(out$mean, out$gene, function(v) c(min(v), max(v)))
  lo <- vapply(out$gene, function(g) rng[[g]][1], numeric(1))
  hi <- vapply(out$gene, function(g) rng[[g]][2], numeric(1))
  transform(out, scaled = (mean - lo) / pmax(hi - lo, 1e-9),
            gene = factor(gene, levels = MARKERS_V1),
            celltype = factor(celltype, levels = rev(CELLTYPE_V1_ORDER)))
}

# Run every analysis once; returns a named list of data frames.
run_fig3_analysis <- function(obj) {
  if (anyNA(figure_celltype(obj))) {
    stop("Labels missing from CELLTYPE_V1_ORDER: ",
         paste(setdiff(unique(as.character(obj$celltype)), CELLTYPE_V1_ORDER), collapse = ", "))
  }
  long <- pericentral_long(obj)
  list(cells = cell_table(obj), long = long, stats = wilcoxon_table(long),
       group = group_summary(long), library = library_summary(long), markers = marker_summary(obj))
}

write_source_data <- function(res, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  files <- c(cells = "Fig3b_umap_cells.csv", group = "Fig3c_S3a_group_summary.csv",
             stats = "FigS3a_wilcoxon.csv", library = "FigS3a_library_means.csv",
             long = "FigS3a_per_nucleus_expression.csv", markers = "FigS3_markers_v1.csv")
  for (k in names(files)) {
    utils::write.csv(res[[k]], file.path(dir, files[[k]]), row.names = FALSE)
  }
  invisible(file.path(dir, files))
}
