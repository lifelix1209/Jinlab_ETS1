# ============================================================
# RNA-ATAC integration at the Wnt2 / Wnt9b loci (helpers for rna_atac_run.R).
#
# ATAC UMAP : RunUMAP on the stored harmony_lsi embedding, dims 2:30
#             (LSI dim 1 tracks depth), seed 42.
# Gene activity: Signac::GeneActivity (gene body + 2 kb upstream) for
#             GA_GENES; counts / nCount_peaks x 1e4, ln(x + 1).
# Groups    : figure cell types; at each locus, EC split into ligand+
#             (>= 1 RNA UMI of that ligand) and ligand- EC.
# Locus peaks: MACS2 'peaks' assay peaks inside LOCI windows (cover the
#             manuscript locus spans and the linked distal peaks).
#             Accessibility per group = peak counts / total peak counts of
#             the group x 1e6 (CPM); row z-score across groups for the heatmap.
# Coverage  : Tn5 insertion sites (fragment ends) from atac_fragments.tsv.gz
#             in BIN_BP bins, per group, scaled to insertions per 1e6
#             fragments-in-peaks of that group, smoothed (running mean of
#             SMOOTH_BINS bins).
# ETS motifs: JASPAR2020 CORE vertebrate ETS-class matrices (ETS1/2, ERG,
#             FLI1, FEV, ETV1-7, ELK1/3/4, GABPA, ELF1-5, EHF, ERF);
#             motifmatchr positions in locus peaks, p < 5e-5, mm10;
#             overlapping hits merged into sites.
# Peak-gene links: Signac::LinkPeaks (all nuclei; RNA LogNormalize; GC
#             RegionStats; distance 5e5; Signac default pvalue 0.05,
#             |score| 0.05).
# ============================================================

LOCI <- list(
  Wnt2  = list(region = "chr6-17915000-18045000", gene = "Wnt2"),
  Wnt9b = list(region = "chr11-103680000-103765000", gene = "Wnt9b")
)
GA_GENES  <- c("Wnt2", "Wnt9b", "Ets1", "Kdr", "Stab2", "Alb", "Clec4f", "Lrat")
DOT_GENES <- c("Wnt2", "Wnt9b", "Kdr", "Stab2", "Alb", "Clec4f", "Lrat")
TRACK_NON_EC <- c("Hepatocyte", "Kupffer cell", "Stellate cell")
HEP_ALL   <- c("Hepatocyte", "Glul+ Hep", "eGFP+ Hep", "PeriportalHep")
ETS_REGEX <- "^(ETS1|ETS2|ERG|FLI1|FEV|ETV[1-7]|ELK[134]|GABPA|ELF[1-5]|EHF|ERF)$"
LSI_DIMS  <- 2:30
BIN_BP    <- 100
SMOOTH_BINS <- 5
MOTIF_P   <- 5e-5
SEED      <- 42

# Broad groups: hepatocyte states pooled; EC kept; others as labelled.
broad_celltype <- function(obj) {
  ct <- as.character(figure_celltype(obj))
  ct[ct %in% HEP_ALL] <- "Hepatocyte"
  ct
}

load_multiome <- function(path) {
  if (!file.exists(path)) stop("Missing ", path)
  obj <- readRDS(path)
  SeuratObject::DefaultAssay(obj) <- "RNA"
  for (a in intersect(c("SCT", "ATAC"), SeuratObject::Assays(obj))) obj[[a]] <- NULL
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  rna <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
  obj$figure_celltype <- as.character(figure_celltype(obj))
  obj$broad <- broad_celltype(obj)
  obj$Wnt2_pos <- as.numeric(rna["Wnt2", ]) > 0
  obj$Wnt9b_pos <- as.numeric(rna["Wnt9b", ]) > 0
  obj
}

atac_umap <- function(obj) {
  emb <- SeuratObject::Embeddings(obj, "harmony_lsi")[, LSI_DIMS]
  set.seed(SEED)
  um <- uwot::umap(emb, n_neighbors = 30, min_dist = 0.3, metric = "cosine", seed = SEED, verbose = FALSE)
  data.frame(barcode = colnames(obj), UMAP1 = um[, 1], UMAP2 = um[, 2], celltype = obj$figure_celltype,
             broad = obj$broad, Wnt2_pos = obj$Wnt2_pos, Wnt9b_pos = obj$Wnt9b_pos, row.names = NULL)
}

gene_activity <- function(obj) {
  SeuratObject::DefaultAssay(obj) <- "peaks"
  ga <- Signac::GeneActivity(obj, features = GA_GENES, assay = "peaks", verbose = FALSE)
  ga <- ga[intersect(GA_GENES, rownames(ga)), , drop = FALSE]
  log1p(t(t(as.matrix(ga)) / obj$nCount_peaks * 1e4))
}

dot_table <- function(obj, ga) {
  rna <- as.matrix(SeuratObject::LayerData(obj, assay = "RNA", layer = "data")[DOT_GENES, ])
  grp <- obj$broad
  one <- function(m, modality) {
    do.call(rbind, lapply(split(seq_len(ncol(m)), grp), function(i) {
      data.frame(gene = rownames(m), modality = modality, celltype = grp[i[1]], n = length(i),
                 pct = 100 * rowMeans(m[, i, drop = FALSE] > 0), mean = rowMeans(m[, i, drop = FALSE]))
    }))
  }
  out <- rbind(one(rna, "RNA expression"), one(ga[intersect(DOT_GENES, rownames(ga)), , drop = FALSE], "Gene activity"))
  out$scaled <- stats::ave(out$mean, out$gene, out$modality, FUN = function(v) (v - min(v)) / max(max(v) - min(v), 1e-9))
  rownames(out) <- NULL
  out
}

locus_groups <- function(obj, ligand) {
  pos <- obj[[paste0(ligand, "_pos"), drop = TRUE]]
  g <- ifelse(obj$broad == "Endothelial", ifelse(pos, paste0(ligand, "+ EC"), paste0(ligand, "− EC")),
              ifelse(obj$broad %in% TRACK_NON_EC, obj$broad, NA))
  factor(g, levels = c(paste0(ligand, "+ EC"), paste0(ligand, "− EC"), TRACK_NON_EC))
}

locus_peaks <- function(obj, locus) {
  span <- Signac::StringToGRanges(locus$region)
  gr <- Signac::granges(obj[["peaks"]])
  pk <- gr[IRanges::overlapsAny(gr, span)]
  names(pk) <- Signac::GRangesToString(pk)
  pk
}

peak_group_cpm <- function(obj, peaks, groups) {
  cnt <- SeuratObject::LayerData(obj, assay = "peaks", layer = "counts")[names(peaks), , drop = FALSE]
  ok <- !is.na(groups)
  tot <- tapply(obj$nCount_peaks[ok], droplevels(groups[ok]), sum)
  m <- vapply(levels(droplevels(groups[ok])), function(g) Matrix::rowSums(cnt[, which(groups == g), drop = FALSE]),
              numeric(length(peaks)))
  m <- matrix(m, nrow = length(peaks), dimnames = list(names(peaks), levels(droplevels(groups[ok]))))
  sweep(m, 2, tot[colnames(m)], "/") * 1e6
}

ets_sites <- function(peaks) {
  pfm <- TFBSTools::getMatrixSet(JASPAR2020::JASPAR2020,
                                 opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE))
  nm <- toupper(vapply(pfm, TFBSTools::name, character(1)))
  ets <- pfm[grepl(ETS_REGEX, nm)]
  if (length(ets) == 0) stop("No ETS-class JASPAR matrices matched")
  genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  pos <- motifmatchr::matchMotifs(ets, peaks, genome = genome, out = "positions", p.cutoff = MOTIF_P)
  hits <- unlist(pos, use.names = FALSE)
  if (length(hits) == 0) return(data.frame(start = numeric(0), end = numeric(0), peak = character(0)))
  sites <- GenomicRanges::reduce(GenomicRanges::granges(hits), ignore.strand = TRUE)
  ov <- GenomicRanges::findOverlaps(sites, peaks)
  data.frame(start = GenomicRanges::start(sites)[S4Vectors::queryHits(ov)],
             end = GenomicRanges::end(sites)[S4Vectors::queryHits(ov)],
             peak = names(peaks)[S4Vectors::subjectHits(ov)], row.names = NULL)
}

peak_links <- function(obj) {
  genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
  SeuratObject::DefaultAssay(obj) <- "peaks"
  obj <- Signac::RegionStats(obj, genome = genome, assay = "peaks", verbose = FALSE)
  set.seed(SEED)
  obj <- Signac::LinkPeaks(obj, peak.assay = "peaks", expression.assay = "RNA",
                           genes.use = names(LOCI), distance = 5e5, verbose = FALSE)
  lk <- as.data.frame(Signac::Links(obj))
  lk[, c("seqnames", "start", "end", "gene", "peak", "score", "zscore", "pvalue")]
}

# Tn5 insertion coverage for each group in a locus window.
insertion_coverage <- function(frag_path, region, groups, nfrag_peaks) {
  gr <- Signac::StringToGRanges(region)
  tbx <- Rsamtools::TabixFile(frag_path)
  lines <- Rsamtools::scanTabix(tbx, param = gr)[[1]]
  if (length(lines) == 0) stop("No fragments in ", region)
  f <- data.table::fread(text = paste(lines, collapse = "\n"), header = FALSE,
                         col.names = c("chr", "start", "end", "barcode", "count"))
  f$group <- groups[f$barcode]
  f <- f[!is.na(f$group), ]
  ins <- data.frame(pos = c(f$start + 1, f$end), group = c(f$group, f$group))
  s <- GenomicRanges::start(gr)
  e <- GenomicRanges::end(gr)
  ins <- ins[ins$pos >= s & ins$pos <= e, ]
  breaks <- seq(s, e + BIN_BP, by = BIN_BP)
  tot <- tapply(nfrag_peaks, groups, sum)
  out <- do.call(rbind, lapply(levels(groups), function(g) {
    h <- graphics::hist(ins$pos[ins$group == g], breaks = breaks, plot = FALSE)$counts
    v <- stats::filter(h, rep(1 / SMOOTH_BINS, SMOOTH_BINS), sides = 2)
    v[is.na(v)] <- 0
    data.frame(group = g, pos = utils::head(breaks, -1) + BIN_BP / 2, signal = as.numeric(v) / tot[[g]] * 1e6)
  }))
  out$group <- factor(out$group, levels = levels(groups))
  out
}

gene_models <- function(obj, region) {
  ann <- Signac::Annotation(obj[["peaks"]])
  gr <- Signac::StringToGRanges(region)
  a <- ann[IRanges::overlapsAny(ann, gr) & ann$gene_biotype == "protein_coding"]
  if (length(a) == 0) return(NULL)
  d <- as.data.frame(a)
  tx_len <- tapply(d$end, d$tx_id, max) - tapply(d$start, d$tx_id, min)
  keep_tx <- vapply(split(d$tx_id, d$gene_name), function(t) names(which.max(tx_len[unique(t)])), character(1))
  d <- d[d$tx_id %in% keep_tx, ]
  d[, c("gene_name", "tx_id", "start", "end", "strand", "type")]
}

locus_analysis <- function(obj, ligand, links) {
  locus <- LOCI[[ligand]]
  groups <- locus_groups(obj, ligand)
  pk <- locus_peaks(obj, locus)
  cpm <- peak_group_cpm(obj, pk, groups)
  g_named <- stats::setNames(as.character(groups), colnames(obj))
  g_named <- factor(g_named[!is.na(g_named)], levels = levels(groups))
  nfp <- obj$nCount_peaks[!is.na(groups)]
  frag <- Signac::GetFragmentData(Signac::Fragments(obj[["peaks"]])[[1]], "path")
  span <- Signac::StringToGRanges(locus$region)
  lk <- links[links$gene == ligand, ]
  lk <- lk[lk$start >= GenomicRanges::start(span) & lk$end <= GenomicRanges::end(span), ]
  list(ligand = ligand, region = locus$region, peaks = data.frame(peak = names(pk), start = GenomicRanges::start(pk),
       end = GenomicRanges::end(pk), linked = names(pk) %in% lk$peak, row.names = NULL),
       cpm = cpm, ets = ets_sites(pk), links = lk,
       coverage = insertion_coverage(frag, locus$region, g_named, nfp),
       genes = gene_models(obj, locus$region),
       n = as.data.frame(table(group = groups)))
}

run_rna_atac <- function(obj) {
  links <- peak_links(obj)
  list(umap = atac_umap(obj), dot = dot_table(obj, gene_activity(obj)), links = links,
       loci = stats::setNames(lapply(names(LOCI), function(l) locus_analysis(obj, l, links)), names(LOCI)))
}

write_atac_source_data <- function(res, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(res$umap, file.path(dir, "atac_umap.csv"), row.names = FALSE)
  utils::write.csv(res$dot, file.path(dir, "rna_vs_gene_activity.csv"), row.names = FALSE)
  utils::write.csv(res$links, file.path(dir, "peak_gene_links.csv"), row.names = FALSE)
  for (l in names(res$loci)) {
    x <- res$loci[[l]]
    utils::write.csv(data.frame(peak = rownames(x$cpm), x$cpm, check.names = FALSE),
                     file.path(dir, paste0(l, "_peak_group_cpm.csv")), row.names = FALSE)
    utils::write.csv(x$peaks, file.path(dir, paste0(l, "_peaks.csv")), row.names = FALSE)
    utils::write.csv(x$ets, file.path(dir, paste0(l, "_ETS_sites.csv")), row.names = FALSE)
    utils::write.csv(x$coverage, file.path(dir, paste0(l, "_coverage.csv")), row.names = FALSE)
    utils::write.csv(x$n, file.path(dir, paste0(l, "_group_sizes.csv")), row.names = FALSE)
  }
  invisible(dir)
}
