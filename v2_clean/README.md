# v2_clean — current analysis

Input object: `data/processed/seu_bko_lib1245_singlets.rds`, built by
`code/00_preprocessing/02_subset_clean_object.R` from the annotated multiome object:

| Library | Genotype | Nuclei | Doublets removed | Kept |
|---|---|---|---|---|
| 1 | control | 4,620 | 312 | 4,308 |
| 2 | control | 3,929 | 251 | 3,678 |
| 3 | control | 321 | — | dropped (PCA outlier, low depth) |
| 4 | β-cat cKO | 3,832 | 242 | 3,590 |
| 5 | β-cat cKO | 4,459 | 280 | 4,179 |
| 6 | β-cat cKO | 1,276 | — | dropped (low depth) |

Clustering, UMAP and cell-type labels are inherited unchanged (no re-clustering).

## Code

| Folder | Script | What it does | Output |
|---|---|---|---|
| `00_preprocessing` | `01_annotation_doublets.R` | scDblFinder per library; v2 label table | `data/processed/cell_annotation_v2.csv`, `results/00_preprocessing/annotation_doublets/` |
| | `02_subset_clean_object.R` | libraries 1/2/4/5, singlets; QC plot | clean `.rds`, `results/00_preprocessing/clean_object_qc/` |
| | `03_reassess_unknown_immune.R` | marker evidence for clusters 7 and 14 | `results/00_preprocessing/annotation_reassessment/` |
| `Figure3` | `fig3_run.R` (+ `fig3_analysis.R`, `fig3_main_panels.R`, `fig3_supp_panels.R`) | Fig. 3a-c, Fig. S3 | `results/Figure3/` |
| `Figure4` | `fig4_run.R` (+ `fig4_analysis.R`, `fig4_plots.R`) | Fig. 4B-C, Fig. S4A-B | `results/Figure4/` |
| `shared` | `nature_theme.R` | theme, palettes, `save_nature()` export | — |

## Figure 3 (`results/Figure3/`)

| File | Panel |
|---|---|
| `Fig3_abc.*` | composite a-c (183 mm) |
| `Fig3a.*`, `Fig3b.*` | workflow schematic; UMAP + genotype split |
| `Fig3c.*` | candidate factors in cKO eGFP⁺ Hep: mean UMI (y), % nuclei ≥ 1 UMI (size); genes sorted by % |
| `Fig3c_control_Glul.*` | same layout for control Glul⁺ Hep |
| `Fig3c_combined.*` | control and cKO side by side per gene (sorted by cKO %) |
| `FigS3a_candidates.*` | log-normalised expression, two-sided Wilcoxon |
| `FigS3_markers_v1.*` | marker dot plot for the cell-type annotation |
| `source_data/` | CSV behind every panel |

## Figure 4 (`results/Figure4/`)

Method (identical for Wnt2 and Wnt9b):

1. Endothelial nuclei; ligand⁺ = ≥ 1 UMI of *Wnt2* (or *Wnt9b*); control and cKO pooled.
2. LogNormalize (scale 10,000); Seurat `FindMarkers`, two-sided Wilcoxon, min.pct 0.01.
   Up-regulated = nominal *P* < 0.1 and log2FC > 0 (candidate screen, no multiple-testing
   correction; Bonferroni and BH values are reported in the tables).
3. Locus peaks: MACS2 `peaks` assay inside chr6:17,925,014-18,032,048 (*Wnt2*) and
   chr11:103,728,955-103,759,011 (*Wnt9b*). Motifs: JASPAR2020 CORE vertebrates,
   Signac `AddMotifs` (p < 5 × 10⁻⁵); single-factor motifs mapped to mouse genes.
4. TFs in ligand⁺ EC = motif in ≥ 1 locus peak ∩ up-regulated. Shared = intersection.

| File | Panel |
|---|---|
| `Fig4B_TF_venn.*` | TF overlap: Wnt2⁺ (37) vs Wnt9b⁺ (17); shared ETS1, TEAD1 |
| `Fig4B_TF_filter_matrix.*` | motif / up-regulation criteria for the candidate TFs |
| `Fig4C_Ets1_expression.*` | *Ets1*, ligand⁻ vs ligand⁺ EC (Wnt9b *P* = 0.088; Wnt2 *P* = 7.1 × 10⁻⁴) |
| `FigS4A_TF_venn.*` | motif TFs vs up-regulated genes for each ligand, and TF overlap |
| `FigS4B_Tead1_expression.*` | *Tead1* (Wnt9b *P* = 0.0097; Wnt2 *P* = 0.028) |
| `FigS4_volcano_Wnt2.*`, `FigS4_volcano_Wnt9b.*` | genome-wide ligand⁺ vs ligand⁻ (reference) |
| `FigS4_Nr2f2_expression.*`, `FigS4_other_TF_expression.*` | other v1 candidates (reference) |
| `Wnt_EC_TF_clean.xlsx`, `tables/` | DE tables, locus peaks, TF tables, motif enrichment |

### Why the v1 shared set (5 TFs) became 2

Re-running the v1 Wilcoxon test on the v1 endothelial object with and without the
removed nuclei (158 doublets, 330 nuclei from libraries 3/6):

| TF | Wnt9b⁺ vs Wnt9b⁻, v1 *P* | after removal | Wnt9b⁺ nuclei expressing (doublets / kept) |
|---|---|---|---|
| Klf2 | 0.047 | 0.50 | 3/6 / 4/31 |
| Nr2f2 | 0.033 | 0.33 | 4/6 / 12/31 |
| Prdm1 | 0.016 | 0.72 | 1/6 / 0/31 |
| Ets1 | 0.072 | 0.088 | 1/6 / 11/31 |
| Tead1 | 0.009 | 0.0097 | 1/6 / 9/31 |

The Wnt9b signal for KLF2, NR2F2 and PRDM1 came from predicted doublets; ETS1 and
TEAD1 are unchanged by the cleaning.

## Caveats

- Nuclei are treated as independent observations (2 mice per genotype); P values are
  anti-conservative and the Fig. 4 screen uses an unadjusted threshold.
- Ligand⁺ nuclei are deeper (median 2,106 vs 1,667 UMI for Wnt2); the Wilcoxon test
  does not adjust for depth.
- The ETS core motif (GGAA) is shared by many ETS-family factors; motif presence does
  not identify ETS1 specifically.

## Revision analyses R01–R04 (Nov 2026)

Run order (from repo root): `00_preprocessing/04_decontx_ambient.R` → `R01` → `R02` → `R03` (needs R02 cache) → `R04`.
Each `R0x_*_run.R` caches to `results/R0x_*/source_data/analysis_results.rds` (`FORCE_ANALYSIS=1` recomputes) and writes a
183-mm composite, single panels, source-data CSVs and a `*.collision-audit.json` (nature-figure QA).
Extra packages: celda (DecontX), nichenetr 2.2.1 + NicheNet v2 mouse priors in `data/external/nichenet/` (Zenodo 7074291),
msigdbr 26.1 (mouse-native MH, M2 CP:REACTOME, M5 GO:BP), fgsea, edgeR.

| Module | Question | Key result |
|---|---|---|
| `R01_hep_state` | Hepatocyte identity loss + signalling-state map (keeps Fig. 3c bubble plot) | cKO eGFP⁺ Hep lose the pericentral WNT score (median 0.83 → −0.06) and sit mid-lobular on a control-defined zonation axis (0.04 → 0.51). *Igf1* and *Nrg4* rise CV→PV in control; in cKO *Igf1* and *Nampt* are higher in every zone. Zone-adjusted (lm, ns(zonation,4)): *Nampt* ↑ in cKO eGFP⁺ (β = 0.36, P = 4 × 10⁻⁸); *Igf1* ↑ in cKO Glul⁺/other but not beyond zone in eGFP⁺ (P = 0.32); *Vegfa* not increased in snRNA (β = −0.27, P = 0.007). |
| `R02_ec_response` | Control vs cKO EC transcriptome, unbiased GSEA | EC nuclei carry ~70–77% hepatocyte ambient RNA (DecontX), lower in cKO. Naive pseudobulk "xenobiotic/bile-acid down, angiogenesis/VEGF/MAPK up" is an ambient/normalisation artefact. After removing ambient genes before normalisation (1,743 EC-intrinsic genes; edgeR QL, 2 vs 2) no gene or pathway passes FDR 0.05. Consistent nominal ↑: *Hgf* (P = 0.001), *Rspo3*, *Gab2*, *Nos3*, ETS family *Ets1/Erg/Fli1* (P ≈ 0.08); ↓ *Nrp1*. |
| `R03_nichenet` | Hepatocyte → EC ligands | Ligand activities are near zero for both gene sets (cKO response n = 23; Wnt-EC program n = 38; corrected AUPR ≤ 0.04) → supplementary only. Prior model: Vegfa has the top regulatory potential toward *Ets1* (100th percentile of 49 ligands), 98th for *Rspo3*, 94th for *Wnt9b*; Igf1/Nampt 88–100th. |
| `R04_rna_atac` | ATAC populations, Wnt2/Wnt9b gene activity, locus specificity, ETS motifs, peak–gene links | EC form a distinct ATAC cluster (948 nuclei). *Wnt2* gene activity is EC-specific; *Wnt9b* whole-gene activity is not (gene body has a ubiquitous peak), but its promoter peak (chr11:103,749 kb, 3 ETS sites) is highest in Wnt9b⁺ EC (50 vs 16 CPM in Wnt9b⁻ EC, ≤ 11 elsewhere). All 4 *Wnt2* linked peaks and both *Wnt9b* linked peaks are EC-enriched; 3/4 and 1/2 carry ETS-class motifs. |

Caveats: nuclei-level P values are anti-conservative (2 mice per genotype); Wnt9b⁺ EC n = 31; ETS motifs (GGAA) do not identify ETS1 specifically; NicheNet prior percentiles are knowledge-based, not data-derived.

### R05_atac_landscape (reference-figure panels on our data)

Run: `MACS3_PATH=<macs3> Rscript v2_clean/code/R05_atac_landscape/atac_landscape_run.R` (MACS3 3.1.0 in a Python venv; peak sets cached in `source_data/celltype_peaks.rds`, BED files per cell type).

| Panel | Reference | Result |
|---|---|---|
| a | Nature liver atlas Fig. b (transcriptome ↔ epigenome) | RNA (harmony_rna) and ATAC (harmony-LSI) UMAPs of the same 15,755 nuclei, 4,000 links |
| b, c | RPE Fig. B; brain Fig. C | MACS3 per cell type (≥ 100 nuclei): EC 59.2k, Hep 55.3k, Kupffer 66.0k, Stellate 21.0k, Immune 40.1k, Blood 96.3k; union 118.6k (12% promoter, 46% intronic, 37% distal) |
| d | RPE Fig. E | Marker-locus coverage per cell type (*Kdr, Clec4f, Lrat, Ptprc, Hbb-bs*) plus *Wnt2*/*Wnt9b*. No hepatocyte marker column: a genome-wide screen found no peak > ~1.5-fold more accessible in hepatocytes than in every other type; hepatocytes are kept as a control row. |
| e, f | Nature liver atlas Fig. e | Loess-smoothed (span 0.35), row z-scored, clipped ±1.5. Control hepatocytes, 40 PV→CV bins of the R01 coordinate: 1,465 zonated regions, 492 zonated genes (Spearman \|ρ\| ≥ 0.7, BH < 0.05). Reference columns: CV-high regions z = +1.33 in control Glul⁺ Hep vs −0.28 in cKO eGFP⁺ Hep (PV-high −1.11 vs −0.17) → identity loss is also visible in chromatin (axis defined from RNA only). |

Caveats: per-cell-type peak numbers scale with nuclei number and depth (Blood cell 96k from 215 nuclei is likely depth/quality-driven), so b/c are descriptive, not a biological comparison. Hepatocyte marker loci (*Alb, Hnf4a, Ttr, Fgb*) are almost as accessible in NPC barcodes as in hepatocytes (ratio 0.7–1.2): NPC nuclei carry hepatocyte chromatin background (ATAC counterpart of the RNA ambient load). EC-specific loci (*Kdr*, *Wnt2*) remain EC > hepatocyte.
