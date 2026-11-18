# ============================================================
# Chromatin landscape of the multiome (helpers for atac_landscape_run.R).
#
# Cell types: figure labels with hepatocyte states pooled (R04
#   broad_celltype); types with < MIN_NUCLEI nuclei are not peak-called.
# Peaks per cell type: MACS3 3.1.0 via Signac::CallPeaks(group.by),
#   fragments from data/raw/atac_fragments.tsv.gz, effective genome size
#   1.87e9, Signac defaults (BED, --nomodel, extsize 200, shift -100,
#   q 0.05). Union = GenomicRanges::reduce of all cell-type peak sets.
# Peak annotation (peak centre, EnsDb v79 / object annotation, protein-
#   coding and lncRNA genes), priority Promoter > Exonic > Intronic > Distal:
#   Promoter = TSS -2,000/+100 bp (ArchR convention); Exonic = exon;
#   Intronic = gene body; Distal = otherwise.
# Linked UMAPs: stored RNA UMAP (harmony_rna) and the R04 ATAC UMAP; both
#   rescaled to the unit square; links drawn for a random subset of nuclei.
# Zonation heatmaps: control hepatocyte nuclei (R01 zonation coordinate,
#   0 = CV, 1 = PV) in N_ZONE_BINS equal-frequency bins; pseudobulk CPM per
#   bin (peaks: / nCount_peaks; RNA: / nCount_RNA); features with >=
#   MIN_FEATURE_COUNTS counts; Spearman rho of log2(CPM + 1) vs bin;
#   zonated = |rho| >= RHO_MIN and BH < 0.05. Rows ordered by the bin of the
#   loess-smoothed maximum. Reference columns: control Glul+ Hep and cKO
#   eGFP+ Hep pseudobulks on the same features.
# ============================================================

MIN_NUCLEI   <- 100
GENOME_SIZE  <- 1.87e9
PROMOTER_UP  <- 2000
PROMOTER_DN  <- 100
N_LINKS      <- 4000
N_ZONE_BINS  <- 40
MIN_FEATURE_COUNTS <- 50
RHO_MIN      <- 0.7
# No hepatocyte column: genome-wide, no peak is > ~1.5-fold more accessible in hepatocytes than in
# every other cell type (NPC barcodes carry hepatocyte chromatin background); hepatocytes stay as a row.
TRACK_GENES <- c(Endothelial = "Kdr", "Kupffer cell" = "Clec4f", "Stellate cell" = "Lrat",
                 "Immune cell" = "Ptprc", "Blood cell" = "Hbb-bs", Endothelial = "Wnt2", Endothelial = "Wnt9b")
TRACK_FLANK  <- 5000
TRACK_BINS   <- 250

peak_call_groups <- function(obj) {
  n <- table(obj$broad)
  keep <- names(n)[n >= MIN_NUCLEI]
  factor(ifelse(obj$broad %in% keep, obj$broad, NA), levels = keep)
}

call_celltype_peaks <- function(obj, macs_path, outdir) {
  if (!file.exists(macs_path)) stop("MACS3 not found at ", macs_path)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  obj$peak_group <- peak_call_groups(obj)
  sub <- subset(obj, cells = colnames(obj)[!is.na(obj$peak_group)])
  SeuratObject::DefaultAssay(sub) <- "peaks"
  pk <- Signac::CallPeaks(sub, assay = "peaks", group.by = "peak_group", macs2.path = macs_path,
                          outdir = outdir, effective.genome.size = GENOME_SIZE, combine.peaks = FALSE,
                          cleanup = TRUE, verbose = FALSE)
  names(pk) <- levels(droplevels(sub$peak_group))
  lapply(pk, function(g) GenomeInfoDb::keepStandardChromosomes(g, pruning.mode = "coarse"))
}

annotation_features <- function(ann) {
  ann <- ann[ann$gene_biotype %in% c("protein_coding", "lincRNA", "lncRNA")]
  d <- data.frame(tx = ann$tx_id, chr = as.character(GenomicRanges::seqnames(ann)), start = GenomicRanges::start(ann),
                  end = GenomicRanges::end(ann), strand = as.character(GenomicRanges::strand(ann)))
  tx <- do.call(rbind, lapply(split(d, d$tx), function(x) x[1, c("chr", "strand")]))
  tx$start <- tapply(d$start, d$tx, min)[rownames(tx)]
  tx$end <- tapply(d$end, d$tx, max)[rownames(tx)]
  txg <- GenomicRanges::GRanges(tx$chr, IRanges::IRanges(tx$start, tx$end), strand = tx$strand)
  tss <- GenomicRanges::resize(txg, width = 1, fix = "start")
  list(tx = txg, prom = GenomicRanges::promoters(tss, upstream = PROMOTER_UP, downstream = PROMOTER_DN),
       exon = ann[ann$type == "exon"])
}

annotate_peaks <- function(peaks, feats) {
  centre <- GenomicRanges::resize(peaks, width = 1, fix = "center")
  GenomeInfoDb::seqlevelsStyle(centre) <- GenomeInfoDb::seqlevelsStyle(feats$prom)[1]
  cls <- ifelse(IRanges::overlapsAny(centre, feats$prom, ignore.strand = TRUE), "Promoter",
         ifelse(IRanges::overlapsAny(centre, feats$exon, ignore.strand = TRUE), "Exonic",
         ifelse(IRanges::overlapsAny(centre, feats$tx, ignore.strand = TRUE), "Intronic", "Distal")))
  factor(cls, levels = c("Distal", "Exonic", "Intronic", "Promoter"))
}

peak_composition <- function(peak_sets, ann) {
  union <- GenomicRanges::reduce(unlist(GenomicRanges::GRangesList(lapply(peak_sets, GenomicRanges::granges))))
  sets <- c(peak_sets, list(Union = union))
  feats <- annotation_features(ann)
  do.call(rbind, lapply(names(sets), function(k) {
    tab <- table(annotate_peaks(sets[[k]], feats))
    data.frame(celltype = k, category = names(tab), n = as.integer(tab), total = length(sets[[k]]))
  }))
}

linked_umaps <- function(obj, atac_umap) {
  rna <- SeuratObject::Embeddings(obj, "umap")
  d <- data.frame(barcode = rownames(rna), rna_x = rna[, 1], rna_y = rna[, 2], celltype = obj$broad)
  d <- merge(d, atac_umap[, c("barcode", "UMAP1", "UMAP2")], by = "barcode")
  names(d)[names(d) == "UMAP1"] <- "atac_x"
  names(d)[names(d) == "UMAP2"] <- "atac_y"
  unit <- function(v) (v - min(v)) / (max(v) - min(v))
  transform(d, rna_x = unit(rna_x), rna_y = unit(rna_y), atac_x = unit(atac_x), atac_y = unit(atac_y))
}

# Tn5 insertions in TRACK_BINS bins over one window, per cell type.
window_coverage <- function(frag_path, region, groups, nfrag) {
  gr <- Signac::StringToGRanges(region)
  lines <- Rsamtools::scanTabix(Rsamtools::TabixFile(frag_path), param = gr)[[1]]
  f <- data.table::fread(text = paste(lines, collapse = "\n"), header = FALSE,
                         col.names = c("chr", "start", "end", "barcode", "count"))
  f$group <- groups[f$barcode]
  f <- f[!is.na(f$group), ]
  ins <- data.frame(pos = c(f$start + 1, f$end), group = c(f$group, f$group))
  s <- GenomicRanges::start(gr)
  e <- GenomicRanges::end(gr)
  breaks <- seq(s, e, length.out = TRACK_BINS + 1)
  tot <- tapply(nfrag, groups, sum)
  do.call(rbind, lapply(levels(groups), function(g) {
    h <- graphics::hist(ins$pos[ins$group == g & ins$pos >= s & ins$pos <= e], breaks = breaks, plot = FALSE)$counts
    data.frame(group = g, pos = (utils::head(breaks, -1) + utils::tail(breaks, -1)) / 2,
               signal = h / tot[[g]] * 1e6 / (diff(breaks)[1] / 1e3))
  }))
}

gene_window <- function(ann, gene) {
  g <- ann[ann$gene_name == gene]
  if (length(g) == 0) stop("Gene not in annotation: ", gene)
  paste0(as.character(GenomicRanges::seqnames(g))[1], "-", min(GenomicRanges::start(g)) - TRACK_FLANK, "-",
         max(GenomicRanges::end(g)) + TRACK_FLANK)
}

marker_tracks <- function(obj) {
  ann <- Signac::Annotation(obj[["peaks"]])
  frag <- Signac::GetFragmentData(Signac::Fragments(obj[["peaks"]])[[1]], "path")
  types <- c("Endothelial", "Hepatocyte", "Kupffer cell", "Stellate cell", "Immune cell", "Blood cell")
  grp <- factor(ifelse(obj$broad %in% types, obj$broad, NA), levels = types)
  ok <- !is.na(grp)
  groups <- stats::setNames(grp[ok], colnames(obj)[ok])
  nfrag <- obj$nCount_peaks[ok]
  cov <- do.call(rbind, lapply(unname(TRACK_GENES), function(gn) {
    reg <- gene_window(ann, gn)
    cbind(gene = gn, region = reg, window_coverage(frag, reg, groups, nfrag))
  }))
  models <- do.call(rbind, lapply(unname(TRACK_GENES), function(gn) {
    reg <- Signac::StringToGRanges(gene_window(ann, gn))
    a <- as.data.frame(ann[IRanges::overlapsAny(ann, reg) & ann$gene_name == gn & ann$type == "exon"])
    tx_len <- tapply(a$end, a$tx_id, max) - tapply(a$start, a$tx_id, min)
    a <- a[a$tx_id == names(which.max(tx_len)), ]
    data.frame(gene = gn, start = a$start, end = a$end, strand = as.character(a$strand))
  }))
  list(coverage = cov, models = models, n = as.data.frame(table(celltype = grp)))
}

# ---------------- zonation heatmaps ----------------
bin_pseudobulk <- function(counts, depth, bins) {
  m <- vapply(levels(bins), function(b) Matrix::rowSums(counts[, which(bins == b), drop = FALSE]), numeric(nrow(counts)))
  lib <- tapply(depth, bins, sum)
  rownames(m) <- rownames(counts)
  sweep(m, 2, lib[colnames(m)], "/") * 1e6
}

zonated_features <- function(cpm, raw_total) {
  keep <- raw_total >= MIN_FEATURE_COUNTS
  x <- log2(cpm[keep, , drop = FALSE] + 1)
  idx <- seq_len(ncol(x))
  rho <- apply(x, 1, function(v) suppressWarnings(stats::cor(v, idx, method = "spearman")))
  p <- apply(x, 1, function(v) suppressWarnings(stats::cor.test(v, idx, method = "spearman", exact = FALSE)$p.value))
  res <- data.frame(feature = rownames(x), rho = rho, p = p, p_BH = stats::p.adjust(p, "BH"))
  res <- res[!is.na(res$rho) & abs(res$rho) >= RHO_MIN & res$p_BH < 0.05, ]
  sm <- t(apply(x[res$feature, , drop = FALSE], 1, function(v) stats::predict(stats::loess(v ~ idx, span = 0.5))))
  res$peak_bin <- apply(sm, 1, which.max)
  res[order(res$peak_bin, -res$rho), ]
}

zonation_heatmaps <- function(obj, cells) {
  ctrl <- cells[cells$genotype == "Control", ]
  q <- stats::quantile(ctrl$zonation, seq(0, 1, length.out = N_ZONE_BINS + 1))
  # PV -> CV left to right: bin 1 = most periportal.
  ctrl$bin <- factor(N_ZONE_BINS + 1 - as.integer(cut(ctrl$zonation, unique(q), include.lowest = TRUE)),
                     levels = seq_len(N_ZONE_BINS))
  ref <- cells[cells$state %in% c("Control Glul+ Hep", "cKO eGFP+ Hep"), ]
  out <- list()
  for (mod in c("peaks", "RNA")) {
    cnt <- SeuratObject::LayerData(obj, assay = mod, layer = "counts")
    depth_col <- if (mod == "peaks") "nCount_peaks" else "nCount_RNA"
    cc <- cnt[, ctrl$barcode]
    cpm <- bin_pseudobulk(cc, obj@meta.data[ctrl$barcode, depth_col], ctrl$bin)
    feat <- zonated_features(cpm, Matrix::rowSums(cc))
    rcpm <- bin_pseudobulk(cnt[feat$feature, ref$barcode], obj@meta.data[ref$barcode, depth_col],
                           factor(ref$state, levels = c("Control Glul+ Hep", "cKO eGFP+ Hep")))
    out[[mod]] <- list(features = feat, cpm = cpm[feat$feature, , drop = FALSE], ref = rcpm)
  }
  out$n_bins <- N_ZONE_BINS
  out$n_ctrl <- nrow(ctrl)
  out
}

write_landscape_source_data <- function(res, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res$composition, file.path(dir, "peaks_per_celltype.csv"), row.names = FALSE)
  utils::write.csv(res$tracks$coverage, file.path(dir, "marker_track_coverage.csv"), row.names = FALSE)
  utils::write.csv(res$tracks$n, file.path(dir, "marker_track_group_sizes.csv"), row.names = FALSE)
  for (m in c("peaks", "RNA")) {
    z <- res$zonation[[m]]
    utils::write.csv(z$features, file.path(dir, paste0("zonated_", m, "_features.csv")), row.names = FALSE)
    utils::write.csv(data.frame(feature = rownames(z$cpm), z$cpm, z$ref, check.names = FALSE),
                     file.path(dir, paste0("zonated_", m, "_bin_cpm.csv")), row.names = FALSE)
  }
  for (k in names(res$peak_sets)) {
    g <- res$peak_sets[[k]]
    utils::write.table(data.frame(GenomicRanges::seqnames(g), GenomicRanges::start(g) - 1, GenomicRanges::end(g)),
                       file.path(dir, paste0("peaks_", gsub(" ", "_", k), ".bed")), sep = "\t",
                       quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  invisible(dir)
}
