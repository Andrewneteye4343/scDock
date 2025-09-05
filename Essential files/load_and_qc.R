load_and_qc <- function(rds_files,
                        min_features = 200,
                        min_cells = 3,
                        names_delim = "_",
                        max_mito = 0.3,
                        load_and_qc_verbose = TRUE) {
  library(Seurat)

  seurat_list <- list()
  message("Loading your data...")
  for (i in seq_along(rds_files)) {
    counts <- readRDS(rds_files[i])

    seurat_obj <- CreateSeuratObject(
      counts = counts,
      min.features = min_features,
      min.cells = min_cells,
      names.delim = names_delim
    )

    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") # The "^MT" pattern infeaturesrepresents Use in human samples.
    if (all(is.na(seurat_obj$percent.mt)) || all(seurat_obj$percent.mt == 0)) {
      warning("⚠️ All percent.mt values are 0 or NA. Please check gene naming (e.g. Ensembl ID vs Symbol).")
    }

    raw_ncells <- ncol(seurat_obj)
    raw_ngenes <- nrow(seurat_obj)

    seurat_obj <- subset(
      seurat_obj,
      subset = nFeature_RNA > min_features & percent.mt < max_mito * 100
    )

    qc_ncells <- ncol(seurat_obj)
    qc_ngenes <- nrow(seurat_obj)

    if (load_and_qc_verbose) message("QC summary for dataset [", basename(rds_files[i]), "]")
    if (load_and_qc_verbose) message("Cells:  ", raw_ncells, " → ", qc_ncells)
    if (load_and_qc_verbose) message("Genes:  ", raw_ngenes, " → ", qc_ngenes)

    seurat_list[[i]] <- seurat_obj
  }

  if (load_and_qc_verbose) message("Merging all datasets ...")
  merged_obj <- seurat_list[[1]]
  if (length(seurat_list) > 1) {
    for (i in 2:length(seurat_list)) {
      merged_obj <- merge(merged_obj, y = seurat_list[[i]])
    }
  }

  if (load_and_qc_verbose) message("Merged Seurat object: ", ncol(merged_obj), " cells × ", nrow(merged_obj), " genes")
  message("✅ QC process complete.")
  return(merged_obj)
}