# data/ (not tracked by git)

Raw and processed data are too large for the repository. Place them here with this layout.

## raw/ (Cell Ranger ARC output)

| File | Description |
|---|---|
| `filtered_feature_bc_matrix.h5` | gene expression matrix, all libraries |
| `atac_fragments.tsv.gz` (+ `.tbi`) | ATAC fragments, tabix-indexed |

## processed/

| File | Made by | Used by |
|---|---|---|
| `BCatKO_multiome_clustered.rds` | upstream pipeline (QC, WNN clustering of all six libraries) | `00_preprocessing/01_annotation_doublets.R` |
| `seu_bko_annotated_with_ATAC_fragments.rds` | upstream pipeline (annotated object, fragments attached) | `00_preprocessing/02_subset_clean_object.R` |
| `cell_annotation_v2.csv` | `00_preprocessing/01_annotation_doublets.R` | doublet calls, v2 labels |
| `seu_bko_lib1245_singlets.rds` | `00_preprocessing/02_subset_clean_object.R` | all figure and revision scripts |
| `rna_decontx_counts.rds`, `rna_decontx_meta.csv` | `00_preprocessing/04_decontx_ambient.R` | `R02_ec_response` |

## external/

| Folder | Content | Used by |
|---|---|---|
| `nichenet/` | NicheNet v2 mouse priors (Zenodo 7074291): `ligand_target_matrix_nsga2r_final_mouse.rds`, `lr_network_mouse_21122021.rds`, `weighted_networks_nsga2r_final_mouse.rds` | `R03_nichenet` |
