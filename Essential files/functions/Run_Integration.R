Run_Integration <- function(seurat_obj,
                            Run_Integration_group_by,
                            Run_Integration_method,
                            Run_Integration_nfeatures,
                            Run_Integration_normalization_method,
                            Run_Integration_dims,
                            Run_Integration_resolution,
                            Run_Integration_verbose) {

  suppressPackageStartupMessages({
    library(Seurat)
    library(harmony)
  })

  if (is.null(seurat_obj@meta.data[[Run_Integration_group_by]])) {
    stop(paste0("[Run_Integration] Metadata column '", Run_Integration_group_by, "' not found in Seurat object."))
  }

  # Branch 1: Harmony
  if (tolower(Run_Integration_method) == "harmony") {
    message("[Run_Integration] Using Harmony integration...")

    # Normalization
    DefaultAssay(seurat_obj) <- "RNA"
    if (toupper(Run_Integration_normalization_method) == "SCT") {
      if (Run_Integration_verbose) message("[Run_Integration] Global SCTransform (whole object)")
      seurat_obj <- SCTransform(seurat_obj, verbose = FALSE, assay = "RNA")
      assay_to_use <- "SCT"
    } else if (toupper(Run_Integration_normalization_method) == "LOGNORMALIZE") {
      if (Run_Integration_verbose) message("[Run_Integration] Global LogNormalize (whole object)")
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = Run_Integration_nfeatures, verbose = FALSE)
      seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
      assay_to_use <- "RNA"
    } else {
      stop("[Run_Integration] Unsupported normalization method: ", Run_Integration_normalization_method)
    }

    # PCA
    DefaultAssay(seurat_obj) <- assay_to_use
    if (Run_Integration_verbose) message("[Run_Integration] Running PCA for Harmony ...")
    seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)

    # Auto-detect optimal PC number (Elbow method)
    if (identical(Run_Integration_dims, "auto")) {
      if (Run_Integration_verbose) message("[Run_Integration] Automatically determining optimal PC number using geometric elbow method...")

      stdev <- seurat_obj[["pca"]]@stdev
      total_pcs <- length(stdev)
      x <- 1:total_pcs
      y <- stdev

      # Geometric method: Find the point farthest from the start-end line.
      p1 <- c(x[1], y[1])
      p2 <- c(x[total_pcs], y[total_pcs])
      distances <- abs((p2[2] - p1[2]) * x - (p2[1] - p1[1]) * y + p2[1]*p1[2] - p2[2]*p1[1]) /
                   sqrt((p2[2] - p1[2])^2 + (p2[1] - p1[1])^2)
      elbow_point <- which.max(distances)
      if (is.na(elbow_point) || elbow_point < 2) elbow_point <- min(30, total_pcs)

      Run_Integration_dims <- elbow_point

      # Save ElbowPlot to PDF
      elbow_plot <- ElbowPlot(seurat_obj)
      plot_filename <- "ElbowPlot_Integration.pdf"
      if (!is.null(seurat_obj@project.name))
        plot_filename <- paste0("ElbowPlot_Integration_", seurat_obj@project.name, ".pdf")
      pdf(plot_filename, width = 6, height = 4)
      print(elbow_plot)
      dev.off()

      if (Run_Integration_verbose) {
        message("[Run_Integration] Elbow plot saved as ", plot_filename)
        message("[Run_Integration] Auto-selected ", Run_Integration_dims, " PCs using geometric elbow method.")
      }
    }

    # Harmony integration
    if (Run_Integration_verbose) message("[Run_Integration] Running RunHarmony ...")
    seurat_obj <- RunHarmony(seurat_obj,
                             group.by.vars = Run_Integration_group_by,
                             reduction.use = "pca",
                             dims.use = 1:Run_Integration_dims,
                             theta = 2,
                             verbose = Run_Integration_verbose)

    # Downstream
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:Run_Integration_dims, verbose = FALSE)
    seurat_obj <- FindClusters(seurat_obj, resolution = Run_Integration_resolution, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:Run_Integration_dims, verbose = FALSE)

    if (Run_Integration_verbose) message("[Run_Integration] Harmony integration finished.")
    return(seurat_obj)
  }

  # Branch 2: Anchor-based integration (CCA / RPCA)
  message("[Run_Integration] Using Seurat anchor-based integration (", toupper(Run_Integration_method), ") ...")

  # Split data
  seurat_list <- SplitObject(seurat_obj, split.by = Run_Integration_group_by)
  if (Run_Integration_verbose) message("[Run_Integration] Split into ", length(seurat_list), " objects by ", Run_Integration_group_by)

  # Normalize
  if (toupper(Run_Integration_normalization_method) == "SCT") {
    if (Run_Integration_verbose) message("[Run_Integration] Performing SCTransform on each sample...")
    seurat_list <- lapply(names(seurat_list), function(nm) {
      x <- seurat_list[[nm]]
      DefaultAssay(x) <- "RNA"
      x <- SCTransform(x, verbose = FALSE, assay = "RNA")
      return(x)
    })
    norm_method_for_anchor <- "SCT"
  } else if (toupper(Run_Integration_normalization_method) == "LOGNORMALIZE") {
    if (Run_Integration_verbose) message("[Run_Integration] Performing LogNormalize + HVG on each sample...")
    seurat_list <- lapply(seurat_list, function(x) {
      DefaultAssay(x) <- "RNA"
      x <- NormalizeData(x, verbose = FALSE)
      x <- FindVariableFeatures(x, nfeatures = Run_Integration_nfeatures, verbose = FALSE)
      return(x)
    })
    norm_method_for_anchor <- "LogNormalize"
  } else {
    stop("[Run_Integration] Unsupported normalization method: ", Run_Integration_normalization_method)
  }

  # Integration features
  if (Run_Integration_verbose) message("[Run_Integration] Selecting integration features ...")
  features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = Run_Integration_nfeatures)

  # Prep SCT integration
  if (norm_method_for_anchor == "SCT") {
    if (Run_Integration_verbose) message("[Run_Integration] Preparing SCT integration ...")
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
  }

  # PCA per object (for anchor finding)
  if (Run_Integration_verbose) message("[Run_Integration] Running PCA on each object for anchor detection ...")
  seurat_list <- lapply(seurat_list, function(x) RunPCA(x, features = features, verbose = FALSE))

  # Auto PC detection for integration
  if (identical(Run_Integration_dims, "auto")) {
    if (Run_Integration_verbose) message("[Run_Integration] Automatically determining optimal PC number for anchor-based integration ...")

    stdev <- seurat_list[[1]][["pca"]]@stdev
    total_pcs <- length(stdev)
    x <- 1:total_pcs
    y <- stdev

    p1 <- c(x[1], y[1])
    p2 <- c(x[total_pcs], y[total_pcs])
    distances <- abs((p2[2] - p1[2]) * x - (p2[1] - p1[1]) * y + p2[1]*p1[2] - p2[2]*p1[1]) /
                 sqrt((p2[2] - p1[2])^2 + (p2[1] - p1[1])^2)
    elbow_point <- which.max(distances)
    if (is.na(elbow_point) || elbow_point < 2) elbow_point <- min(30, total_pcs)

    Run_Integration_dims <- elbow_point

    # Save ElbowPlot (from first object)
    elbow_plot <- ElbowPlot(seurat_list[[1]])
    plot_filename <- "ElbowPlot_Integration.pdf"
    if (!is.null(seurat_obj@project.name))
      plot_filename <- paste0("ElbowPlot_Integration_", seurat_obj@project.name, ".pdf")
    pdf(plot_filename, width = 6, height = 4)
    print(elbow_plot)
    dev.off()

    if (Run_Integration_verbose) {
      message("[Run_Integration] Elbow plot saved as ", plot_filename)
      message("[Run_Integration] Auto-selected ", Run_Integration_dims, " PCs using geometric elbow method.")
    }
  }

  # Find anchors
  if (Run_Integration_verbose) message("[Run_Integration] Finding integration anchors ...")
  anchors <- FindIntegrationAnchors(object.list = seurat_list,
                                    normalization.method = norm_method_for_anchor,
                                    anchor.features = features,
                                    reduction = tolower(Run_Integration_method),
                                    dims = 1:Run_Integration_dims,
                                    verbose = Run_Integration_verbose)

  # Integrate
  if (Run_Integration_verbose) message("[Run_Integration] Integrating data ...")
  seurat_obj <- IntegrateData(anchorset = anchors, normalization.method = norm_method_for_anchor, verbose = Run_Integration_verbose)

  # Downstream
  DefaultAssay(seurat_obj) <- "integrated"
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = Run_Integration_dims, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:Run_Integration_dims, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, resolution = Run_Integration_resolution, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:Run_Integration_dims, reduction = "pca", verbose = FALSE)

  if (Run_Integration_verbose) message("[Run_Integration] Anchor-based integration finished.")
  return(seurat_obj)
}
