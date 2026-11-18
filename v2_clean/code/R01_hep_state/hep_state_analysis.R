# ============================================================
# Hepatocyte identity loss and signalling-state redistribution
# (analysis helpers for hep_state_run.R).
#
# Nuclei   : Hepatocyte, Glul+ Hep, eGFP+ Hep, PeriportalHep (v1 labels),
#            libraries 1, 2 (control) and 4, 5 (cKO), scDblFinder singlets.
# Groups   : Control Glul+ Hep (control, Glul+ Hep cluster) and cKO eGFP+ Hep
#            (cKO, eGFP+ Hep clusters), as in Fig. 3; cKO Glul+ Hep = the
#            WNT-responsive neighbouring layer in cKO livers; all other
#            hepatocytes = "Other".
# Scores   : Seurat AddModuleScore on LogNormalized RNA (scale 1e4) for
#            (i) pericentral WNT/beta-catenin-responsive landmarks (PC_GENES),
#            (ii) periportal landmarks (PP_GENES), (iii) candidate factors.
#            Candidate factors are excluded from both landmark sets.
# Zonation : d = PP score - PC score; coordinate = empirical CDF of d in
#            control hepatocytes (0 = most pericentral, 1 = most periportal);
#            cKO nuclei are placed with the same control CDF, so control
#            deciles hold 10% of control nuclei by construction.
# Zone-adjusted effect: lm(expr ~ ns(zonation, 4) + group) on control
#            nuclei plus one cKO group; beta = difference at matched zonation
#            (nuclei as observations; anti-conservative, 2 mice per genotype).
# ============================================================

PC_GENES <- c("Glul", "Axin2", "Lgr5", "Tbx3", "Cyp2e1", "Cyp1a2", "Oat", "Rhbg", "Lect2", "Slc1a2",
              "Cyp27a1", "Notum", "Rnf43", "Gulo", "Slc13a3", "Cyp7a1", "Lhpp", "Nkd1", "Sp5")
PP_GENES <- c("Cyp2f2", "Sds", "Hal", "Hsd17b13", "Pck1", "Gls2", "Aldh1b1", "Ass1", "Arg1", "Cps1",
              "Hamp", "Igfbp2", "Cdh1", "Sdsl", "Ctsc", "Mfsd2a", "Pigr")
FACTOR_SCORE_GENES <- c("Igf1", "Nampt", "Vegfa")
FACTOR_GENES <- c("Igf1", "Nampt", "Vegfa", "Nrg4")
HEATMAP_GENES <- list(
  "Pericentral (WNT-responsive)" = c("Glul", "Axin2", "Lgr5", "Tbx3", "Cyp2e1", "Oat", "Lect2", "Rhbg", "Cyp1a2"),
  "Periportal" = c("Cyp2f2", "Sds", "Hal", "Hsd17b13", "Pck1", "Arg1", "Igfbp2"),
  "Candidate factors" = FACTOR_GENES
)
HEP_TYPES  <- c("Hepatocyte", "Glul+ Hep", "eGFP+ Hep", "PeriportalHep")
STATE_LEVELS <- c("Control Glul+ Hep", "cKO eGFP+ Hep", "cKO Glul+ Hep", "Other")
N_BINS <- 10
SEED <- 42

if (length(intersect(FACTOR_GENES, c(PC_GENES, PP_GENES))) > 0) stop("Candidate factors overlap landmarks")

state_group <- function(genotype, celltype) {
  out <- ifelse(genotype == "Control" & celltype == "Glul+ Hep", "Control Glul+ Hep",
         ifelse(genotype == "β-cat cKO" & celltype == "eGFP+ Hep", "cKO eGFP+ Hep",
         ifelse(genotype == "β-cat cKO" & celltype == "Glul+ Hep", "cKO Glul+ Hep", "Other")))
  factor(out, levels = STATE_LEVELS)
}

load_hepatocytes <- function(path) {
  if (!file.exists(path)) stop("Missing ", path)
  obj <- tryCatch(readRDS(path), error = function(e) stop("Failed to read ", path, ": ", conditionMessage(e)))
  SeuratObject::DefaultAssay(obj) <- "RNA"
  keep <- colnames(obj)[as.character(obj$celltype) %in% HEP_TYPES]
  h <- subset(obj, cells = keep)
  rm(obj); invisible(gc(verbose = FALSE))
  for (a in intersect(c("SCT", "ATAC", "peaks"), SeuratObject::Assays(h))) h[[a]] <- NULL
  Seurat::NormalizeData(h, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
}

add_state_scores <- function(h) {
  missing <- setdiff(c(PC_GENES, PP_GENES, FACTOR_GENES), rownames(h))
  if (length(missing) > 0) stop("Genes absent from RNA assay: ", paste(missing, collapse = ", "))
  set.seed(SEED)
  Seurat::AddModuleScore(h, features = list(PC_GENES, PP_GENES, FACTOR_SCORE_GENES), name = "STATE",
                         seed = SEED)
}

# Per-nucleus table: scores, zonation coordinate, candidate expression.
state_table <- function(h) {
  md <- h@meta.data
  geno <- genotype_label(md$Group)
  ct <- as.character(md$celltype)
  expr <- as.matrix(SeuratObject::LayerData(h, assay = "RNA", layer = "data")[FACTOR_GENES, , drop = FALSE])
  d <- md$STATE2 - md$STATE1
  ecdf_ctrl <- stats::ecdf(d[geno == "Control"])
  z <- ecdf_ctrl(d)
  out <- data.frame(barcode = colnames(h), Library = paste0("L", md$Library), genotype = geno, celltype = ct,
                    state = state_group(geno, ct), pc_score = md$STATE1, pp_score = md$STATE2,
                    factor_score = md$STATE3, zonation = z,
                    zone_bin = cut(z, breaks = seq(0, 1, length.out = N_BINS + 1), include.lowest = TRUE,
                                   labels = seq_len(N_BINS)),
                    row.names = NULL)
  cbind(out, stats::setNames(as.data.frame(t(expr)), FACTOR_GENES))
}

# Library-level medians of the scores per state (n = 2 mice per genotype).
library_state_summary <- function(cells) {
  cells |>
    dplyr::filter(state != "Other") |>
    dplyr::group_by(state, genotype, Library) |>
    dplyr::summarise(n = dplyr::n(), pc_score = stats::median(pc_score), pp_score = stats::median(pp_score),
                     zonation = stats::median(zonation), factor_score = stats::median(factor_score),
                     .groups = "drop")
}

# Share of hepatocyte nuclei per control-defined zonation decile, per library.
zone_occupancy <- function(cells) {
  cells |>
    dplyr::count(genotype, Library, zone_bin, .drop = FALSE) |>
    dplyr::group_by(genotype, Library) |>
    dplyr::mutate(pct = 100 * n / sum(n)) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(genotype))
}

# Mean log-normalised expression of the heatmap genes per zonation bin and
# genotype, plus the two pericentral reference groups.
heatmap_matrix <- function(h, cells) {
  genes <- unlist(HEATMAP_GENES, use.names = FALSE)
  x <- SeuratObject::LayerData(h, assay = "RNA", layer = "data")[genes, cells$barcode, drop = FALSE]
  col_id <- paste(cells$genotype, cells$zone_bin, sep = "|")
  ref <- c("Control Glul+ Hep", "cKO eGFP+ Hep")
  ids <- c(split(seq_len(nrow(cells)), col_id),
           stats::setNames(lapply(ref, function(s) which(cells$state == s)), ref))
  m <- vapply(ids, function(i) Matrix::rowMeans(x[, i, drop = FALSE]), numeric(length(genes)))
  rownames(m) <- genes
  long <- data.frame(gene = rep(genes, times = ncol(m)), column = rep(colnames(m), each = length(genes)),
                     mean = as.vector(m), row.names = NULL)
  long$z <- stats::ave(long$mean, long$gene, FUN = function(v) (v - mean(v)) / max(stats::sd(v), 1e-9))
  long$module <- factor(rep(rep(names(HEATMAP_GENES), lengths(HEATMAP_GENES)), times = ncol(m)),
                        levels = names(HEATMAP_GENES))
  long$block <- ifelse(long$column %in% ref, "Pericentral groups", sub("\\|.*$", "", long$column))
  long$bin <- ifelse(long$column %in% ref, NA_character_, sub("^.*\\|", "", long$column))
  long
}

# Mean +/- s.e.m. of candidate expression per zonation bin and genotype.
factor_profiles <- function(cells) {
  long <- tidyr::pivot_longer(cells, dplyr::all_of(FACTOR_GENES), names_to = "gene", values_to = "expr")
  long |>
    dplyr::group_by(gene, genotype, zone_bin) |>
    dplyr::summarise(n = dplyr::n(), mean = mean(expr), sem = stats::sd(expr) / sqrt(dplyr::n()),
                     pct = 100 * mean(expr > 0), zonation = mean(zonation), .groups = "drop") |>
    dplyr::mutate(gene = factor(gene, levels = FACTOR_GENES))
}

# Pericentral reference groups placed on the profile plot.
factor_reference_points <- function(cells) {
  long <- tidyr::pivot_longer(dplyr::filter(cells, state %in% STATE_LEVELS[1:2]), dplyr::all_of(FACTOR_GENES),
                              names_to = "gene", values_to = "expr")
  long |>
    dplyr::group_by(gene, state) |>
    dplyr::summarise(n = dplyr::n(), mean = mean(expr), sem = stats::sd(expr) / sqrt(dplyr::n()),
                     zonation = stats::median(zonation), z_lo = stats::quantile(zonation, 0.25),
                     z_hi = stats::quantile(zonation, 0.75), .groups = "drop") |>
    dplyr::mutate(gene = factor(gene, levels = FACTOR_GENES))
}

# Zone-adjusted difference: per gene and cKO group, a linear model on the
# control nuclei plus that group's nuclei,
#   expr ~ ns(zonation, df = 4) + group,
# so beta = mean difference from control at the same zonation coordinate.
ZONE_DF <- 4

zone_adjusted_effects <- function(cells) {
  ctrl <- cells$genotype == "Control"
  groups <- list("cKO eGFP+ Hep" = cells$state == "cKO eGFP+ Hep",
                 "cKO Glul+ Hep" = cells$state == "cKO Glul+ Hep",
                 "cKO other Hep" = cells$genotype == "β-cat cKO" & cells$state == "Other")
  out <- do.call(rbind, lapply(FACTOR_GENES, function(g) {
    do.call(rbind, lapply(names(groups), function(k) {
      d <- cells[ctrl | groups[[k]], ]
      d$is_group <- as.integer(groups[[k]][ctrl | groups[[k]]])
      fit <- stats::lm(stats::as.formula(paste0("`", g, "` ~ splines::ns(zonation, df = ZONE_DF) + is_group")), data = d)
      co <- summary(fit)$coefficients["is_group", ]
      ci <- stats::confint(fit)["is_group", ]
      data.frame(gene = g, group = k, n_group = sum(d$is_group), n_control = sum(ctrl), beta = co[["Estimate"]],
                 ci_lo = ci[[1]], ci_hi = ci[[2]], p = co[["Pr(>|t|)"]])
    }))
  }))
  out$p_BH <- stats::p.adjust(out$p, "BH")
  out$gene <- factor(out$gene, levels = rev(FACTOR_GENES))
  out$group <- factor(out$group, levels = names(groups))
  out
}

# Two-sided Wilcoxon of the scores, cKO eGFP+ Hep vs control Glul+ Hep (per nucleus).
state_score_tests <- function(cells) {
  a <- cells[cells$state == "cKO eGFP+ Hep", ]
  b <- cells[cells$state == "Control Glul+ Hep", ]
  do.call(rbind, lapply(c("pc_score", "pp_score", "zonation", "factor_score"), function(v) {
    t <- stats::wilcox.test(a[[v]], b[[v]], exact = FALSE)
    data.frame(score = v, median_control_glul = stats::median(b[[v]]), median_cko_egfp = stats::median(a[[v]]),
               n_control = nrow(b), n_cko = nrow(a), p = t$p.value)
  }))
}

run_hep_state_analysis <- function(h) {
  h <- add_state_scores(h)
  cells <- state_table(h)
  list(cells = cells, library = library_state_summary(cells), occupancy = zone_occupancy(cells),
       heatmap = heatmap_matrix(h, cells), profiles = factor_profiles(cells),
       reference = factor_reference_points(cells), residuals = zone_adjusted_effects(cells),
       tests = state_score_tests(cells),
       gene_sets = data.frame(set = c(rep("pericentral_WNT", length(PC_GENES)), rep("periportal", length(PP_GENES)),
                                      rep("candidate_factor_score", length(FACTOR_SCORE_GENES))),
                              gene = c(PC_GENES, PP_GENES, FACTOR_SCORE_GENES)))
}

write_hep_source_data <- function(res, dir) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (k in c("cells", "library", "occupancy", "heatmap", "profiles", "reference", "residuals", "tests",
              "gene_sets")) {
    utils::write.csv(res[[k]], file.path(dir, paste0("hep_state_", k, ".csv")), row.names = FALSE)
  }
  invisible(dir)
}
