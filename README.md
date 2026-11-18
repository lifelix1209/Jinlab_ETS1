# Hepatocyte state disruption reinforces the endothelial Wnt niche through ETS1

Analysis code for the single-nucleus multiome (snRNA-seq + snATAC-seq) part of the
study: pericentral β-catenin deletion (`Axin2-rtTA; TetO-Cre; β-cat f/del`, "β-cat cKO")
versus control mouse liver.

Dataset analysed here: libraries 1, 2 (control) and 4, 5 (cKO), scDblFinder singlets
only (15,755 nuclei; 948 endothelial nuclei, of which Wnt2⁺ 75 and Wnt9b⁺ 31).
Fig. 4 TF screen: Wilcoxon *P* < 0.1, log2FC > 0, with the same threshold for Wnt2 and
Wnt9b. Shared TFs: **ETS1, TEAD1**.

## Repository layout

```
.
├── R/                          shared helpers (Arial fonts, annotation table path)
├── v2_clean/
│   ├── README.md               methods, key results and caveats per module
│   └── code/
│       ├── shared/             Nature-style theme and export helpers
│       ├── 00_preprocessing/   doublet calling, clean object, annotation check, DecontX
│       ├── Figure3/            Fig. 3a-c, Fig. S3
│       ├── Figure4/            Fig. 4B-C, Fig. S4A-B
│       ├── R01_hep_state/      hepatocyte identity loss / signalling-state map
│       ├── R02_ec_response/    EC response with ambient-RNA correction
│       ├── R03_nichenet/       hepatocyte -> EC NicheNet
│       ├── R04_rna_atac/       RNA-ATAC integration, Wnt2/Wnt9b loci
│       └── R05_atac_landscape/ per-cell-type ATAC landscape
└── data/                       not tracked; see data/README.md
```

Each analysis folder follows the same pattern: `*_run.R` is the entry point,
`*_analysis.R` holds the computations and `*_plots.R` (or `*_panels.R`) the figure
builders.

## Running the analysis

All scripts are run **from the repository root** (paths are relative to it). Inputs
are described in [`data/README.md`](data/README.md).

```bash
# 0. preprocessing
Rscript v2_clean/code/00_preprocessing/01_annotation_doublets.R    # scDblFinder -> data/processed/cell_annotation_v2.csv
Rscript v2_clean/code/00_preprocessing/02_subset_clean_object.R    # -> data/processed/seu_bko_lib1245_singlets.rds
Rscript v2_clean/code/00_preprocessing/03_reassess_unknown_immune.R
Rscript v2_clean/code/00_preprocessing/04_decontx_ambient.R        # -> data/processed/rna_decontx_*

# 1. main figures
Rscript v2_clean/code/Figure3/fig3_run.R
Rscript v2_clean/code/Figure4/fig4_run.R

# 2. revision analyses (R03 needs the R02 cache; R05 needs MACS3)
Rscript v2_clean/code/R01_hep_state/hep_state_run.R
Rscript v2_clean/code/R02_ec_response/ec_response_run.R
Rscript v2_clean/code/R03_nichenet/nichenet_run.R
Rscript v2_clean/code/R04_rna_atac/rna_atac_run.R
MACS3_PATH=/path/to/macs3 Rscript v2_clean/code/R05_atac_landscape/atac_landscape_run.R
```

Outputs are written to `v2_clean/results/<module>/`. Run scripts cache their analysis
in `results/<module>/.../analysis_results.rds`; set `FORCE_ANALYSIS=1` to recompute.

## Environment

R 4.4.3 (macOS, aarch64). Main packages: Seurat 5.4.0, SeuratObject 5.3.0,
Signac 1.16.0, scDblFinder 1.20.2, JASPAR2020 0.99.10, TFBSTools 1.44.0,
motifmatchr 1.28.0, BSgenome.Mmusculus.UCSC.mm10 1.4.3, EnsDb.Mmusculus.v79 2.99.0,
presto 1.0.0, ggplot2 4.0.1, patchwork 1.3.2, ggrepel 0.9.6, ggtext 0.1.2,
ggforce 0.5.0, ggrastr 1.0.2, svglite 2.2.2, ragg 1.5.0, openxlsx 4.2.8.1, dplyr 1.1.4.
Revision modules additionally use celda (DecontX), nichenetr 2.2.1 (NicheNet v2 mouse
priors, Zenodo 7074291), msigdbr 26.1, fgsea, edgeR and MACS3 3.1.0.
Full session information is written to `v2_clean/results/Figure4/sessionInfo.txt`
when Fig. 4 is run.

Figures use Arial; PDFs are written with the macOS Quartz device (editable text),
with `cairo_pdf` as a fallback on other systems.

## Data and outputs

Only code and documentation are version-controlled. Raw and processed data (`data/`),
generated figures and tables (`v2_clean/results/`) and manuscript files are excluded
by `.gitignore`.

## Author

Han Zhang (hanzhang)
