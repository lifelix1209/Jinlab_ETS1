# ============================================================
# Hepatocyte -> endothelial NicheNet analysis (helpers for nichenet_run.R).
#
# Prior model: NicheNet v2 mouse (Zenodo 7074291):
#   ligand_target_matrix_nsga2r_final_mouse.rds, lr_network_mouse_21122021.rds,
#   weighted_networks_nsga2r_final_mouse.rds  (data/external/nichenet/).
# Sender   : cKO eGFP+ Hep; ligand expressed in >= 10% of sender nuclei.
# Receiver : EC; receptor detected in >= 5% of EC nuclei and EC-intrinsic
#            (ambient fraction < 0.5, R02).
# Background: EC-intrinsic genes detected in >= 5% of EC nuclei.
# Gene sets (from R02, ambient-corrected):
#   cKO response = up in corrected pseudobulk (P < 0.05, logFC > 0) AND up in
#                  nucleus-level test (BH < 0.05, log2FC > 0);
#   Wnt-EC program = Wnt2+/Wnt9b+ vs other EC, P < 0.01, log2FC > 0.5.
# Ligand activity: nichenetr::predict_ligand_activities (AUPR corrected).
# Sender change : cKO eGFP+ Hep vs control Glul+ Hep, log2 ratio of mean
#            normalised expression; % nuclei detecting per hepatocyte group.
# Niche prior: regulatory potential of each potential ligand for Ets1, Erg,
#            Fli1, Wnt2, Wnt9b, Rspo3, Hgf, as percentile among all potential
#            ligands (prior knowledge only; independent of the cKO data).
# ============================================================

NN_DIR       <- "data/external/nichenet"
SENDER_PCT   <- 0.10
RECEIVER_PCT <- 0.05
N_TOP        <- 10
N_TARGETS    <- 200
CANDIDATES_NN <- c("Igf1", "Nampt", "Vegfa", "Nrg4")
NICHE_GENES  <- c("Ets1", "Erg", "Fli1", "Wnt2", "Wnt9b", "Rspo3", "Hgf")
HEP_GROUPS   <- c("Control Glul+ Hep", "cKO eGFP+ Hep", "cKO other Hep")

load_priors <- function(dir = NN_DIR) {
  f <- file.path(dir, c("ligand_target_matrix_nsga2r_final_mouse.rds", "lr_network_mouse_21122021.rds",
                        "weighted_networks_nsga2r_final_mouse.rds"))
  miss <- f[!file.exists(f)]
  if (length(miss) > 0) stop("Missing NicheNet prior files: ", paste(miss, collapse = ", "))
  list(ltm = readRDS(f[1]), lr = dplyr::distinct(readRDS(f[2]), from, to), wn = readRDS(f[3]))
}

load_nn_counts <- function(obj_path) {
  obj <- readRDS(obj_path)
  ct <- as.character(figure_celltype(obj))
  geno <- genotype_label(obj$Group)
  grp <- ifelse(ct == "Endothelial", "EC",
         ifelse(geno == "Control" & ct == "Glul+ Hep", "Control Glul+ Hep",
         ifelse(geno == "β-cat cKO" & ct == "eGFP+ Hep", "cKO eGFP+ Hep",
         ifelse(geno == "β-cat cKO" & ct %in% c("Hepatocyte", "PeriportalHep", "Glul+ Hep"), "cKO other Hep", NA))))
  keep <- !is.na(grp)
  counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")[, keep]
  rm(obj); invisible(gc(verbose = FALSE))
  list(counts = counts, group = grp[keep])
}

pct_detect <- function(counts, cols) Matrix::rowMeans(counts[, cols, drop = FALSE] > 0)

mean_norm <- function(counts, cols) {
  m <- counts[, cols, drop = FALSE]
  Matrix::rowMeans(t(t(m) / Matrix::colSums(m)) * 1e4)
}

gene_sets <- function(ec) {
  pb <- ec$corrected
  nuc <- ec$nucleus
  resp <- intersect(pb$gene[pb$PValue < 0.05 & pb$logFC > 0], nuc$gene[nuc$p_BH < 0.05 & nuc$log2FC > 0])
  prog <- ec$program$gene[ec$program$p < 0.01 & ec$program$log2FC > 0.5]
  list("cKO response" = resp, "Wnt-EC program" = prog)
}

potential_ligands <- function(nc, ec, pri) {
  s_pct <- pct_detect(nc$counts, nc$group == "cKO eGFP+ Hep")
  e_pct <- pct_detect(nc$counts, nc$group == "EC")
  amb <- ec$ambient
  rec_ok <- names(e_pct)[e_pct >= RECEIVER_PCT]
  rec_ok <- rec_ok[rec_ok %in% amb$gene[amb$ambient_fraction < 0.5]]
  lig_ok <- names(s_pct)[s_pct >= SENDER_PCT]
  lr <- pri$lr[pri$lr$from %in% lig_ok & pri$lr$to %in% rec_ok, ]
  list(ligands = intersect(unique(lr$from), colnames(pri$ltm)), lr = lr, receivers = rec_ok)
}

ligand_activity <- function(sets, background, ligands, ltm) {
  do.call(rbind, lapply(names(sets), function(k) {
    gs <- intersect(sets[[k]], rownames(ltm))
    la <- nichenetr::predict_ligand_activities(geneset = gs, background_expressed_genes = background,
                                               ligand_target_matrix = ltm, potential_ligands = ligands)
    la <- la[order(-la$aupr_corrected), ]
    data.frame(gene_set = k, n_genes = length(gs), ligand = la$test_ligand, aupr_corrected = la$aupr_corrected,
               auroc = la$auroc, pearson = la$pearson, rank = seq_len(nrow(la)), row.names = NULL)
  }))
}

select_ligands <- function(act) {
  a <- act[act$gene_set == "cKO response", ]
  unique(c(utils::head(a$ligand, N_TOP), intersect(CANDIDATES_NN, a$ligand)))
}

ligand_targets <- function(ligands, targets, ltm) {
  rows <- lapply(ligands, function(l) {
    v <- ltm[, l]
    top <- names(utils::head(sort(v, decreasing = TRUE), N_TARGETS))
    hit <- intersect(top, targets)
    if (length(hit) == 0) return(NULL)
    data.frame(ligand = l, target = hit, weight = unname(v[hit]))
  })
  do.call(rbind, rows)
}

ligand_receptors <- function(ligands, lr, wn) {
  sig <- wn$lr_sig
  d <- merge(lr[lr$from %in% ligands, ], sig, by = c("from", "to"), all.x = TRUE)
  d$weight[is.na(d$weight)] <- 0
  names(d)[1:2] <- c("ligand", "receptor")
  d
}

sender_table <- function(nc, ligands) {
  cols <- lapply(stats::setNames(HEP_GROUPS, HEP_GROUPS), function(g) nc$group == g)
  pct <- vapply(cols, function(c) pct_detect(nc$counts[ligands, , drop = FALSE], c), numeric(length(ligands)))
  mn <- vapply(cols, function(c) mean_norm(nc$counts, c)[ligands], numeric(length(ligands)))
  dim(pct) <- dim(mn) <- c(length(ligands), length(HEP_GROUPS))
  long <- data.frame(ligand = rep(ligands, times = length(HEP_GROUPS)),
                     group = rep(HEP_GROUPS, each = length(ligands)),
                     pct = as.vector(pct), mean = as.vector(mn), row.names = NULL)
  lfc <- data.frame(ligand = ligands, log2FC = log2((mn[, 2] + 1e-3) / (mn[, 1] + 1e-3)))
  list(expression = long, lfc = lfc)
}

niche_prior <- function(ligands_all, ltm) {
  genes <- intersect(NICHE_GENES, rownames(ltm))
  m <- vapply(genes, function(g) {
    v <- ltm[g, ligands_all]
    100 * rank(v, ties.method = "average") / length(v)
  }, numeric(length(ligands_all)))
  rownames(m) <- ligands_all
  m
}

run_nichenet <- function(nc, ec, pri) {
  sets <- gene_sets(ec)
  pot <- potential_ligands(nc, ec, pri)
  e_pct <- pct_detect(nc$counts, nc$group == "EC")
  bg <- intersect(ec$intrinsic_genes[ec$intrinsic_genes %in% names(e_pct)[e_pct >= RECEIVER_PCT]], rownames(pri$ltm))
  act <- ligand_activity(sets, bg, pot$ligands, pri$ltm)
  show <- select_ligands(act)
  prior <- niche_prior(pot$ligands, pri$ltm)
  prior_long <- data.frame(ligand = rep(rownames(prior), times = ncol(prior)),
                           target = rep(colnames(prior), each = nrow(prior)),
                           percentile = as.vector(prior), row.names = NULL)
  snd <- sender_table(nc, show)
  list(gene_sets = data.frame(gene_set = rep(names(sets), lengths(sets)), gene = unlist(sets), row.names = NULL),
       n = data.frame(n_potential_ligands = length(pot$ligands), n_background = length(bg),
                      n_response = length(sets[[1]]), n_program = length(sets[[2]])),
       activity = act, shown = show,
       targets = ligand_targets(show, sets[["cKO response"]], pri$ltm),
       receptors = ligand_receptors(show, pot$lr, pri$wn),
       sender = snd$expression, sender_lfc = snd$lfc, prior = prior_long)
}

write_nn_source_data <- function(res, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (k in c("gene_sets", "n", "activity", "targets", "receptors", "sender", "sender_lfc", "prior")) {
    utils::write.csv(res[[k]], file.path(dir, paste0("nichenet_", k, ".csv")), row.names = FALSE)
  }
  invisible(dir)
}
