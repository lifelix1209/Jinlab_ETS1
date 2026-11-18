# ============================================================
# Shared helpers for the v2 annotation (v2_clean/code/00_preprocessing/01_annotation_doublets.R)
# ============================================================

ANNOT_V2_CSV <- "data/processed/cell_annotation_v2.csv"

# Read the v2 annotation table; stops with a clear message if absent.
read_annotation_v2 <- function(path = ANNOT_V2_CSV) {
  if (!file.exists(path)) {
    stop("Annotation not found: ", path, ". Run v2_clean/code/00_preprocessing/01_annotation_doublets.R first.")
  }
  annot <- tryCatch(read.csv(path, stringsAsFactors = FALSE),
                    error = function(e) stop("Failed to read ", path, ": ", conditionMessage(e)))
  required <- c("barcode", "celltype_v2", "scDblFinder_class")
  missing <- setdiff(required, colnames(annot))
  if (length(missing) > 0) stop("Annotation lacks columns: ", paste(missing, collapse = ", "))
  annot
}

# Barcodes of nuclei called singlets by scDblFinder.
singlet_barcodes <- function(path = ANNOT_V2_CSV) {
  annot <- read_annotation_v2(path)
  annot$barcode[annot$scDblFinder_class == "singlet"]
}
