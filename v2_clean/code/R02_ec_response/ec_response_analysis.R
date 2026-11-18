# ============================================================
# Endothelial transcriptional response, control vs beta-cat cKO
# (analysis helpers for ec_response_run.R).
#
# Nuclei : v1 "Endothelial" (cluster 3); libraries L1, L2 control, L4, L5 cKO.
# Ambient: EC nuclei carry hepatocyte ambient RNA (DecontX median ~0.7-0.8,
#          lower in cKO). Per gene, ambient fraction in EC =
#          1 - sum(DecontX counts) / sum(raw counts) over EC nuclei
#          (00_preprocessing/04_decontx_ambient.R). EC-intrinsic genes =
#          ambient fraction < AMBIENT_MAX and >= MIN_COUNTS raw counts.
#          Ambient genes are removed BEFORE normalisation so that the
#          genotype difference in ambient load does not rescale EC genes.
# DE     : (1) naive pseudobulk, all genes (reference, shows the artefact);
#          (2) primary: pseudobulk per library on EC-intrinsic genes, edgeR
#              quasi-likelihood (robust), ~ genotype, n = 2 vs 2;
#          (3) exploratory: nucleus-level two-sided Wilcoxon on EC-intrinsic,
#              log-normalised data (nuclei as replicates; anti-conservative).
# GSEA   : fgsea (multilevel), ranking = sign(logFC) x -log10(P), MSigDB
#          mouse-native Hallmark (MH), Reactome (M2 CP:REACTOME), GO:BP (M5),
#          set size 15-500, BH across all sets.
# Program: Wnt-producing EC program = Wnt2+ or Wnt9b+ EC vs other EC
#          (nucleus-level Wilcoxon, EC-intrinsic data); used by R03.
# ============================================================

EC_LABEL     <- "Endothelial"
AMBIENT_MAX  <- 0.5
MIN_COUNTS   <- 10
MIN_PCT_NUC  <- 0.05
GSEA_MIN     <- 15
GSEA_MAX     <- 500
SEED         <- 42
CTRL_LIBS    <- c("L1", "L2")
SHOW_GENES   <- c("Hgf", "Rspo3", "Gab2", "Nos3", "Ets1", "Erg", "Fli1", "Wnt2", "Wnt9b", "Nrp1")

load_ec_counts <- function(obj_path, dx_counts_path) {
  for (p in c(obj_path, dx_counts_path)) if (!file.exists(p)) stop("Missing ", p)
  obj <- readRDS(obj_path)
  keep <- colnames(obj)[figure_celltype(obj) == EC_LABEL]
  raw <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")[, keep]
  md <- obj@meta.data[keep, c("Library", "Group", "nCount_RNA")]
  rm(obj); invisible(gc(verbose = FALSE))
  dx <- readRDS(dx_counts_path)
  list(raw = raw, dx = dx[rownames(raw), keep], meta = data.frame(barcode = keep, Library = paste0("L", md$Library),
       genotype = genotype_label(md$Group), nCount_RNA = md$nCount_RNA, row.names = NULL))
}

ambient_table <- function(raw, dx) {
  tot <- Matrix::rowSums(raw)
  frac <- 1 - Matrix::rowSums(dx) / pmax(tot, 1)
  data.frame(gene = rownames(raw), raw_counts = tot, ambient_fraction = frac,
             intrinsic = frac < AMBIENT_MAX & tot >= MIN_COUNTS, row.names = NULL)
}

pseudobulk <- function(raw, meta, genes) {
  sapply(split(meta$barcode, meta$Library), function(b) Matrix::rowSums(raw[genes, b, drop = FALSE]))
}

edger_pseudobulk <- function(pb) {
  grp <- factor(ifelse(colnames(pb) %in% CTRL_LIBS, "Control", "cKO"), levels = c("Control", "cKO"))
  y <- edgeR::DGEList(pb, group = grp)
  keep <- edgeR::filterByExpr(y, group = grp)
  y <- edgeR::normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  des <- stats::model.matrix(~grp)
  y <- edgeR::estimateDisp(y, des)
  fit <- edgeR::glmQLFit(y, des, robust = TRUE)
  tt <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = 2), n = Inf)$table
  cpm <- edgeR::cpm(y, log = FALSE)
  list(table = data.frame(gene = rownames(tt), tt, row.names = NULL), cpm = cpm)
}

lognorm <- function(m) {
  m <- methods::as(m, "CsparseMatrix")
  m@x <- log1p(m@x / rep.int(Matrix::colSums(m), diff(m@p)) * 1e4)
  m
}

wilcox_genes <- function(x, in_group, min_pct = MIN_PCT_NUC) {
  genes <- rownames(x)[Matrix::rowMeans(x > 0) >= min_pct]
  pr <- presto::wilcoxauc(as.matrix(x[genes, ]), ifelse(in_group, "a", "b"))
  pr <- pr[pr$group == "a", ]
  m_a <- Matrix::rowMeans(expm1(x[genes, in_group, drop = FALSE]))
  m_b <- Matrix::rowMeans(expm1(x[genes, !in_group, drop = FALSE]))
  out <- data.frame(gene = pr$feature, log2FC = log2((m_a[pr$feature] + 1e-3) / (m_b[pr$feature] + 1e-3)),
                    pct_in = pr$pct_in / 100, pct_out = pr$pct_out / 100, auc = pr$auc, p = pr$pval,
                    row.names = NULL)
  out$p_BH <- stats::p.adjust(out$p, "BH")
  out[order(out$p), ]
}

msig_sets <- function() {
  s <- rbind(msigdbr::msigdbr(species = "Mus musculus", db_species = "MM", collection = "MH"),
             msigdbr::msigdbr(species = "Mus musculus", db_species = "MM", collection = "M2",
                              subcollection = "CP:REACTOME"),
             msigdbr::msigdbr(species = "Mus musculus", db_species = "MM", collection = "M5",
                              subcollection = "GO:BP"))
  split(s$gene_symbol, s$gs_name)
}

run_gsea <- function(gene, logfc, p, sets, tag) {
  ok <- !is.na(p) & !is.na(logfc)
  rk <- sign(logfc[ok]) * -log10(pmax(p[ok], 1e-300))
  names(rk) <- gene[ok]
  set.seed(SEED)
  rk <- sort(rk + stats::rnorm(length(rk), 0, 1e-6))
  set.seed(SEED)
  g <- suppressWarnings(fgsea::fgsea(sets, rk, minSize = GSEA_MIN, maxSize = GSEA_MAX, nPermSimple = 10000))
  g <- g[order(g$pval), ]
  data.frame(analysis = tag, pathway = g$pathway, NES = g$NES, p = g$pval, p_BH = g$padj, size = g$size,
             leading_edge = vapply(g$leadingEdge, function(v) paste(utils::head(v, 15), collapse = ";"), ""),
             row.names = NULL)
}

library_expression <- function(cpm, genes) {
  genes <- intersect(genes, rownames(cpm))
  data.frame(gene = rep(genes, times = ncol(cpm)), Library = rep(colnames(cpm), each = length(genes)),
             genotype = factor(ifelse(rep(colnames(cpm), each = length(genes)) %in% CTRL_LIBS, "Control", "β-cat cKO"),
                               levels = names(GROUP_COLS)),
             cpm = as.vector(cpm[genes, , drop = FALSE]), row.names = NULL)
}

wnt_program <- function(raw_intr, x_intr) {
  pos <- raw_intr["Wnt2", ] > 0 | raw_intr["Wnt9b", ] > 0
  out <- wilcox_genes(x_intr, pos)
  out[!out$gene %in% c("Wnt2", "Wnt9b"), ]
}

run_ec_response <- function(ec) {
  amb <- ambient_table(ec$raw, ec$dx)
  intr <- amb$gene[amb$intrinsic]
  naive <- edger_pseudobulk(pseudobulk(ec$raw, ec$meta, rownames(ec$raw)))
  corr <- edger_pseudobulk(pseudobulk(ec$raw, ec$meta, intr))
  x_intr <- lognorm(ec$raw[intr, ])
  nuc <- wilcox_genes(x_intr, ec$meta$genotype == "β-cat cKO")
  naive$table$ambient_fraction <- amb$ambient_fraction[match(naive$table$gene, amb$gene)]
  naive$table$intrinsic <- naive$table$gene %in% intr
  sets <- msig_sets()
  gsea <- rbind(run_gsea(naive$table$gene, naive$table$logFC, naive$table$PValue, sets, "Naive pseudobulk"),
                run_gsea(corr$table$gene, corr$table$logFC, corr$table$PValue, sets, "Corrected pseudobulk"),
                run_gsea(nuc$gene, nuc$log2FC, nuc$p, sets, "Corrected nucleus-level"))
  list(ambient = amb, naive = naive$table, corrected = corr$table, nucleus = nuc, gsea = gsea,
       library = library_expression(corr$cpm, SHOW_GENES), program = wnt_program(ec$raw[intr, ], x_intr),
       intrinsic_genes = intr,
       n = data.frame(n_ec = ncol(ec$raw), n_intrinsic = length(intr),
                      n_control = sum(ec$meta$genotype == "Control"), n_cko = sum(ec$meta$genotype == "β-cat cKO")))
}

write_ec_source_data <- function(res, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (k in c("ambient", "naive", "corrected", "nucleus", "gsea", "library", "program", "n")) {
    utils::write.csv(res[[k]], file.path(dir, paste0("ec_", k, ".csv")), row.names = FALSE)
  }
  invisible(dir)
}
