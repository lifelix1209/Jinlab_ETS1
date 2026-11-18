# ============================================================
# Analysis helpers for v2_clean/code/Figure4/fig4_run.R
#
# Wnt2+/- and Wnt9b+/- endothelial comparisons and TF intersection,
# re-run on the cleaned object (libraries 1, 2, 4, 5; singlets) with a
# standard Seurat workflow and ONE rule applied to both ligands:
#
#   EC        : v1 "Endothelial" nuclei (cluster 3), both genotypes pooled.
#   Groups    : ligand+ = >= 1 UMI of Wnt2 (or Wnt9b); ligand- = 0 UMI.
#   RNA       : NormalizeData (LogNormalize, scale 1e4).
#   DE        : FindMarkers, two-sided Wilcoxon rank-sum (Seurat default),
#               min.pct 0.01, logfc.threshold 0. Up = nominal (unadjusted)
#               P < 0.1 and avg_log2FC > 0 (candidate screen, same for both ligands). Bonferroni and BH adjusted P are
#               reported in the tables only, not used for selection.
#   Locus     : MACS2 'peaks' assay peaks inside the manuscript locus spans.
#   Motifs    : JASPAR2020 CORE vertebrates, Signac AddMotifs (motifmatchr,
#               p < 5e-5, mm10). Enrichment: FindMotifs hypergeometric vs
#               40,000 GC-matched peaks, BH (reported, not used to filter).
#   TF sets   : single-factor motif present in >= 1 locus peak AND the
#               mapped mouse gene is Up. Shared = intersection of both sets.
# ============================================================

LIGANDS      <- c("Wnt2", "Wnt9b")
LOCUS_SPANS  <- c(Wnt2 = "chr6-17925014-18032048", Wnt9b = "chr11-103728955-103759011")
FOCAL_TFS    <- c("Ets1", "Klf2", "Nr2f2", "Prdm1", "Tead1")
EC_LABEL     <- "Endothelial"
P_CUT        <- 0.1
MIN_PCT      <- 0.01
N_BACKGROUND <- 40000
SEED         <- 42

# ---------------- cells ----------------
load_ec <- function(path) {
  if (!file.exists(path)) stop("Missing ", path, ". Run v2_clean/code/00_preprocessing/02_subset_clean_object.R first.")
  obj <- tryCatch(readRDS(path), error = function(e) stop("Failed to read ", path, ": ", conditionMessage(e)))
  keep <- colnames(obj)[figure_celltype(obj) == EC_LABEL]
  if (length(keep) < 100) stop("Too few endothelial nuclei: ", length(keep))
  SeuratObject::DefaultAssay(obj) <- "RNA"
  ec <- subset(obj, cells = keep)
  rm(obj); invisible(gc(verbose = FALSE))
  for (a in intersect(c("SCT", "ATAC"), SeuratObject::Assays(ec))) ec[[a]] <- NULL
  Seurat::NormalizeData(ec, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
}

add_ligand_groups <- function(ec) {
  cnt <- SeuratObject::LayerData(ec, assay = "RNA", layer = "counts")
  Reduce(function(o, l) {
    umi <- as.numeric(cnt[l, colnames(o)])
    o[[paste0(l, "_umi")]] <- umi
    o[[paste0(l, "_group")]] <- factor(ifelse(umi >= 1, "pos", "neg"), levels = c("neg", "pos"))
    o
  }, LIGANDS, ec)
}

group_summary <- function(ec) {
  md <- ec@meta.data
  do.call(rbind, lapply(LIGANDS, function(l) {
    g <- md[[paste0(l, "_group")]]
    do.call(rbind, lapply(levels(g), function(lv) {
      s <- md[g == lv, ]
      libs <- table(factor(paste0("n_L", s$Library), levels = paste0("n_L", c(1, 2, 4, 5))))
      cbind(data.frame(ligand = l, group = paste0(l, ifelse(lv == "pos", "+", "-")), n = nrow(s),
                       n_control = sum(s$Group == "wt"), n_cKO = sum(s$Group == "BCatKO"),
                       median_nCount_RNA = stats::median(s$nCount_RNA),
                       median_nFeature_RNA = stats::median(s$nFeature_RNA)),
            as.data.frame.list(as.list(libs)))
    }))
  }))
}

# ---------------- differential expression ----------------
run_de <- function(ec, ligand) {
  col <- paste0(ligand, "_group")
  res <- Seurat::FindMarkers(ec, group.by = col, ident.1 = "pos", ident.2 = "neg", test.use = "wilcox",
                             logfc.threshold = 0, min.pct = MIN_PCT, verbose = FALSE)
  g <- ec[[col, drop = TRUE]]
  out <- data.frame(ligand = ligand, gene = rownames(res), n_pos = sum(g == "pos"), n_neg = sum(g == "neg"),
                    avg_log2FC = res$avg_log2FC, pct_pos = res$pct.1, pct_neg = res$pct.2,
                    p_val = res$p_val, p_val_adj_bonferroni = res$p_val_adj,
                    p_val_adj_BH = stats::p.adjust(res$p_val, "BH"), row.names = NULL)
  out$up <- out$p_val < P_CUT & out$avg_log2FC > 0
  out[order(out$p_val, -out$avg_log2FC), ]
}

# ---------------- locus peaks & motifs ----------------
locus_peaks <- function(ec, ligand) {
  gr <- Signac::granges(ec[["peaks"]])
  span <- Signac::StringToGRanges(LOCUS_SPANS[[ligand]], sep = c("-", "-"))
  pk <- rownames(ec[["peaks"]])[IRanges::overlapsAny(gr, span)]
  if (length(pk) == 0) stop("No peaks in ", ligand, " locus span")
  cnt <- SeuratObject::LayerData(ec, assay = "peaks", layer = "counts")[pk, , drop = FALSE]
  g <- ec[[paste0(ligand, "_group"), drop = TRUE]]
  data.frame(ligand = ligand, peak = pk,
             pct_ligand_pos = Matrix::rowMeans(cnt[, g == "pos", drop = FALSE] > 0),
             pct_ligand_neg = Matrix::rowMeans(cnt[, g == "neg", drop = FALSE] > 0), row.names = NULL)
}

jaspar_pfm <- function() {
  TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020,
                          opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))
}

add_motifs <- function(ec, pfm) {
  genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  out <- Signac::AddMotifs(ec, genome = genome, pfm = pfm, assay = "peaks", verbose = FALSE)
  if (!"GC.percent" %in% colnames(out[["peaks"]][[]])) {
    out <- Signac::RegionStats(out, genome = genome, assay = "peaks", verbose = FALSE)
  }
  out
}

# Hypergeometric enrichment vs GC-matched background (BH in column p.adjust).
motif_enrichment <- function(ec, peaks, ligand) {
  meta <- ec[["peaks"]][[]]
  set.seed(SEED)
  bg <- Signac::MatchRegionStats(meta.feature = meta[setdiff(rownames(meta), peaks), ],
                                 query.feature = meta[peaks, , drop = FALSE],
                                 features.match = "GC.percent", n = N_BACKGROUND, verbose = FALSE)
  res <- Signac::FindMotifs(ec, features = peaks, background = bg, assay = "peaks", verbose = FALSE)
  cbind(ligand = ligand, res)
}

# Single-factor motifs only; "(var.n)" stripped; case-insensitive symbol match.
map_motif_genes <- function(pfm, genes) {
  nm <- vapply(pfm, TFBSTools::name, character(1))
  clean <- sub("\\(var\\.[0-9]+\\)$", "", nm)
  single <- !grepl("::", clean, fixed = TRUE)
  idx <- match(toupper(clean), toupper(genes))
  keep <- single & !is.na(idx)
  data.frame(motif_id = names(pfm)[keep], motif_name = nm[keep], gene = genes[idx[keep]], row.names = NULL)
}

# Per gene: motif hits in the locus peaks joined to the DE statistics.
tf_table <- function(motif_mat, peaks, mt, de) {
  hits <- motif_mat[peaks, mt$motif_id, drop = FALSE] > 0
  per_motif <- data.frame(mt, n_peaks = Matrix::colSums(hits),
                          peaks_with_motif = vapply(seq_len(ncol(hits)), function(j) {
                            paste(peaks[as.logical(hits[, j])], collapse = ";")
                          }, character(1)))
  present <- per_motif[per_motif$n_peaks > 0, ]
  by_gene <- split(present, present$gene)
  out <- do.call(rbind, lapply(names(by_gene), function(g) {
    d <- by_gene[[g]]
    data.frame(gene = g, motifs = paste(d$motif_name, collapse = ", "),
               peaks_with_motif = paste(unique(unlist(strsplit(d$peaks_with_motif, ";"))), collapse = ";"))
  }))
  out <- merge(out, de, by = "gene", all.x = TRUE)
  out$tested <- !is.na(out$p_val)
  out$up[is.na(out$up)] <- FALSE
  out[order(!out$up, out$p_val), ]
}

shared_table <- function(tfs) {
  s <- intersect(tfs$Wnt2$gene[tfs$Wnt2$up], tfs$Wnt9b$gene[tfs$Wnt9b$up])
  if (length(s) == 0) return(data.frame(gene = character(0)))
  cols <- c("gene", "avg_log2FC", "pct_pos", "pct_neg", "p_val", "p_val_adj_bonferroni")
  w2 <- tfs$Wnt2[tfs$Wnt2$gene %in% s, cols]
  w9 <- tfs$Wnt9b[tfs$Wnt9b$gene %in% s, cols]
  names(w2)[-1] <- paste0("Wnt2_", names(w2)[-1])
  names(w9)[-1] <- paste0("Wnt9b_", names(w9)[-1])
  merge(w2, w9, by = "gene")
}

# Focal TFs: statistics from the genome-wide tests (NA = not tested, < min.pct).
focal_table <- function(de) {
  do.call(rbind, lapply(LIGANDS, function(l) {
    d <- de[[l]]
    k <- match(FOCAL_TFS, d$gene)
    data.frame(ligand = l, gene = FOCAL_TFS, d[k, setdiff(names(d), c("ligand", "gene"))], row.names = NULL)
  }))
}

# TF sets and their intersection under alternative P-value rules (same rule
# for both ligands); "Unadjusted" is the primary analysis, adjusted rules are reference only.
P_RULES <- c(Unadjusted = "p_val", BH = "p_val_adj_BH", Bonferroni = "p_val_adj_bonferroni")

tf_rule_sensitivity <- function(tfs) {
  do.call(rbind, lapply(names(P_RULES), function(r) {
    s <- lapply(tfs, function(t) {
      p <- t[[P_RULES[[r]]]]
      t$gene[!is.na(p) & p < P_CUT & t$avg_log2FC > 0]
    })
    sh <- intersect(s$Wnt2, s$Wnt9b)
    data.frame(rule = r, threshold = paste0(P_RULES[[r]], " < ", P_CUT, " & avg_log2FC > 0"),
               n_Wnt2 = length(s$Wnt2), n_Wnt9b = length(s$Wnt9b), n_shared = length(sh),
               Wnt2_TFs = paste(s$Wnt2, collapse = ", "), Wnt9b_TFs = paste(s$Wnt9b, collapse = ", "),
               shared_TFs = paste(sh, collapse = ", "))
  }))
}

# Per-nucleus log-normalised values for plotting.
expression_long <- function(ec, genes) {
  x <- as.matrix(SeuratObject::LayerData(ec, assay = "RNA", layer = "data")[genes, , drop = FALSE])
  md <- ec@meta.data
  do.call(rbind, lapply(LIGANDS, function(l) {
    g <- md[[paste0(l, "_group")]]
    data.frame(ligand = l, gene = rep(genes, times = ncol(x)), barcode = rep(colnames(x), each = length(genes)),
               group = rep(ifelse(g == "pos", paste0(l, "+"), paste0(l, "-")), each = length(genes)),
               genotype = rep(as.character(md$Group), each = length(genes)), expr = as.vector(x))
  }))
}

# ---------------- run ----------------
run_wnt_tf_analysis <- function(ec) {
  ec <- add_ligand_groups(ec)
  de <- stats::setNames(lapply(LIGANDS, function(l) run_de(ec, l)), LIGANDS)
  pfm <- jaspar_pfm()
  ec <- add_motifs(ec, pfm)
  motif_mat <- Signac::GetMotifData(ec, assay = "peaks", slot = "data")
  mt <- map_motif_genes(pfm, rownames(ec[["RNA"]]))
  peaks <- stats::setNames(lapply(LIGANDS, function(l) locus_peaks(ec, l)), LIGANDS)
  enr <- stats::setNames(lapply(LIGANDS, function(l) motif_enrichment(ec, peaks[[l]]$peak, l)), LIGANDS)
  tfs <- stats::setNames(lapply(LIGANDS, function(l) tf_table(motif_mat, peaks[[l]]$peak, mt, de[[l]])), LIGANDS)
  shared <- shared_table(tfs)
  list(groups = group_summary(ec), de = de, peaks = peaks, enrichment = enr, tfs = tfs, shared = shared,
       rules = tf_rule_sensitivity(tfs), focal = focal_table(de), n_motifs = length(pfm), n_single_mapped = nrow(mt),
       long = expression_long(ec, unique(c(FOCAL_TFS, shared$gene))))
}
