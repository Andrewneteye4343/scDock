# functions/Normalization_Scale.R
Normalization_Scale <- function(seurat_obj,
                                Normalization_Scale_use_assay = NULL,
                                Normalization_Scale_normalization_method = "LogNormalize",
                                Normalization_Scale_scale_factor = 1e4,
                                Normalization_Scale_CLR_margin = NULL,
                                Normalization_Scale_hvg_method = "vst",
                                Normalization_Scale_nVariableFeatures = 2000,
                                Normalization_Scale_split_by = NULL,
                                Normalization_Scale_model_use = "linear",
                                Normalization_Scale_scale_max = 10,
                                Normalization_Scale_scale_features = "variable",
                                Normalization_Scale_verbose = TRUE) {
  library(Seurat)
  
  # -------------------------
  # NormalizeData
  # -------------------------
  if (Normalization_Scale_verbose) message("Step 1: Normalizing data using method: ", Normalization_Scale_normalization_method)
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = Normalization_Scale_normalization_method,
    scale.factor = Normalization_Scale_scale_factor,
    margin = Normalization_Scale_CLR_margin,
    verbose = Normalization_Scale_verbose,
    assay = Normalization_Scale_use_assay
  )
  
  # -------------------------
  # FindVariableFeatures
  # -------------------------
  if (Normalization_Scale_verbose) message("Step 2: Finding top ", Normalization_Scale_nVariableFeatures, " variable features using method: ", Normalization_Scale_hvg_method)
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = Normalization_Scale_hvg_method,
    nfeatures = Normalization_Scale_nVariableFeatures,
    verbose = Normalization_Scale_verbose,
    assay = Normalization_Scale_use_assay
  )
  
  # -------------------------
  # ScaleData
  # -------------------------
  if (Normalization_Scale_verbose) message("Step 3: Scaling data...")
  
  # 決定要 scale 的基因
  if (is.character(Normalization_Scale_scale_features) && Normalization_Scale_scale_features == "variable") {
    features_to_scale <- VariableFeatures(seurat_obj)
    
  } else if (is.character(Normalization_Scale_scale_features) && Normalization_Scale_scale_features == "all") {
    features_to_scale <- rownames(seurat_obj)
    
  } else if (is.vector(Normalization_Scale_scale_features)) {
    missing_genes <- setdiff(Normalization_Scale_scale_features, rownames(seurat_obj))
    if (length(missing_genes) > 0) {
      warning("⚠️ Some genes in Normalization_Scale_scale_features were not found and will be ignored: ",
              paste(head(missing_genes), collapse = ", "),
              if (length(missing_genes) > 5) " ..." else "")
    }
    features_to_scale <- intersect(Normalization_Scale_scale_features, rownames(seurat_obj))
  }
  
  # 執行 scale
  seurat_obj <- ScaleData(
    seurat_obj,
    split.by = Normalization_Scale_split_by,
    model.use = Normalization_Scale_model_use,
    scale.max = Normalization_Scale_scale_max,
    features = features_to_scale,
    verbose = Normalization_Scale_verbose
  )
  return(seurat_obj)
}