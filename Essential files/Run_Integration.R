Run_Integration <- function(seurat_obj,
                            Run_Integration_group_by,
                            Run_Integration_method = "cca",   # 可選 "cca", "rpca", "harmony"
                            Run_Integration_nfeatures = 2000,
                            Run_Integration_normalization_method = "SCT",  # 可選 "SCT" 或 "LogNormalize"
                            Run_Integration_dims = 20,
                            Run_Integration_verbose = TRUE) {
  if (Run_Integration_verbose) message("========== Starting data integration ... ==========")

  if (is.null(seurat_obj@meta.data[[Run_Integration_group_by]])) {
    stop(paste0("❌ Metadata column '", Run_Integration_group_by, "' not found in Seurat object."))
  }
  # -------------------------
  # Harmony-based integration
  # -------------------------
  if (tolower(Run_Integration_method) == "harmony") {
    if (Run_Integration_verbose) message("▶ Using Harmony integration ...")

    # Step 1: Normalization & HVG
    if (Run_Integration_normalization_method == "SCT") {
      if (Run_Integration_verbose) message("Normalization: SCT")
      seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
    } else if (Run_Integration_normalization_method == "LogNormalize") {
      if (Run_Integration_verbose) message("Normalization: LogNormalize")
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = Run_Integration_nfeatures, verbose = FALSE)
    } else {
      stop("❌ Unsupported normalization method for Harmony: ", Run_Integration_normalization_method)
    }

    # Step 2: PCA
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = Run_Integration_dims, verbose = FALSE)

    # Step 3: Harmony integration
    library(harmony)
    seurat_obj <- RunHarmony(seurat_obj,
                             group.by.vars = Run_Integration_group_by,
                             reduction.use = "pca",
                             dims.use = 1:Run_Integration_dims,
                             verbose = Run_Integration_verbose)

    # Step 4: Downstream clustering
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:Run_Integration_dims)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:Run_Integration_dims, reduction = "harmony")

    if (Run_Integration_verbose) message("✅ Harmony integration completed successfully!")

  } else {
    # -------------------------
    # Seurat anchor-based integration (CCA / RPCA)
    # -------------------------
    if (Run_Integration_verbose) message("▶ Using Seurat integration (", Run_Integration_method, ") ...")

    # Step 1: Split object
    seurat_list <- SplitObject(seurat_obj, split.by = Run_Integration_group_by)

    if (Run_Integration_normalization_method == "SCT") {
      if (Run_Integration_verbose) message("Using SCT normalization ...")
      seurat_list <- lapply(seurat_list, function(x) SCTransform(x, verbose = FALSE))
      
      # Step 2: Feature selection
      features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = Run_Integration_nfeatures)
      seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)

      # Step 3: Find anchors
      anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                        normalization.method = "SCT",
                                        anchor.features = features,
                                        reduction = Run_Integration_method,
                                        dims = 1:Run_Integration_dims)

      # Step 4: Integrate
      seurat_obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

    } else if (Run_Integration_normalization_method == "LogNormalize") {
      if (Run_Integration_verbose) message("Using LogNormalize ...")
      seurat_list <- lapply(seurat_list, function(x) {
        x <- NormalizeData(x, verbose = FALSE)
        x <- FindVariableFeatures(x, nfeatures = Run_Integration_nfeatures, verbose = FALSE)
        x
      })

      # Step 2: Feature selection
      features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = Run_Integration_nfeatures)

      # Step 3: Find anchors
      anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                        normalization.method = "LogNormalize",
                                        anchor.features = features,
                                        reduction = Run_Integration_method,
                                        dims = 1:Run_Integration_dims)

      # Step 4: Integrate
      seurat_obj <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize")

    } else {
      stop("❌ Unsupported normalization method for Seurat integration: ", Run_Integration_normalization_method)
    }

    # Step 5: Downstream clustering
    DefaultAssay(seurat_obj) <- "integrated"
    seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    seurat_obj <- RunPCA(seurat_obj, npcs = Run_Integration_dims, verbose = FALSE)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:Run_Integration_dims)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.1)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:Run_Integration_dims, reduction = "pca")

    if (Run_Integration_verbose) message("✅ Seurat integration (", Run_Integration_method, ") completed successfully!")
  }

  return(seurat_obj)
}
