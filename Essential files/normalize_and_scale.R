# functions/normalize_and_scale.R

normalize_and_scale <- function(seurat_obj,
                                use_assay = NULL,
                                normalization_method = "LogNormalize",
                                scale_factor = 1e4,
                                CLR_margin = NULL,
                                hvg_method = "vst",
                                nVariableFeatures = 2000,
                                split_by = NULL,
                                model_use = "linear",
                                scale_max = 10,
                                scale_features = "variable",
                                normalize_and_scale_verbose = TRUE) {
  library(Seurat)
  
  # -------------------------
  # NormalizeData
  # -------------------------
  message("Starting normalization and scaling...")
  if (normalize_and_scale_verbose) message("Step 1: Normalizing data using method: ", normalization_method)
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = normalization_method,
    scale.factor = scale_factor,
    margin = CLR_margin,
    verbose = normalize_and_scale_verbose,
    assay = use_assay
  )
  
  # -------------------------
  # FindVariableFeatures
  # -------------------------
  if (normalize_and_scale_verbose) message("Step 2: Finding top ", nVariableFeatures, " variable features using method: ", hvg_method)
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = hvg_method,
    nfeatures = nVariableFeatures,
    verbose = normalize_and_scale_verbose,
    assay = use_assay
  )
  
  # -------------------------
  # ScaleData
  # -------------------------
  if (normalize_and_scale_verbose) message("Step 3: Scaling data...")
  
  # 決定要 scale 的基因
  if (is.character(scale_features) && scale_features == "variable") {
    features_to_scale <- VariableFeatures(seurat_obj)
    
  } else if (is.character(scale_features) && scale_features == "all") {
    features_to_scale <- rownames(seurat_obj)
    
  } else if (is.vector(scale_features)) {
    missing_genes <- setdiff(scale_features, rownames(seurat_obj))
    if (length(missing_genes) > 0) {
      warning("⚠️ Some genes in scale_features were not found and will be ignored: ",
              paste(head(missing_genes), collapse = ", "),
              if (length(missing_genes) > 5) " ..." else "")
    }
    features_to_scale <- intersect(scale_features, rownames(seurat_obj))
  }
  
  # 執行 scale
  seurat_obj <- ScaleData(
    seurat_obj,
    split.by = split_by,
    model.use = model_use,
    scale.max = scale_max,
    features = features_to_scale,
    verbose = normalize_and_scale_verbose
  )
  
  message("✅ Normalization, HVG selection, and scaling complete.")
  return(seurat_obj)
}