# functions/DR_Cluster.R
DR_Cluster <- function(seurat_obj,
                       DR_Cluster_pca_features,
                       DR_Cluster_seed,
                       DR_Cluster_dims,
                       DR_Cluster_k_param,
                       DR_Cluster_n_trees,
                       DR_Cluster_resolution,
                       DR_Cluster_reduction_method,
                       DR_Cluster_reduction_assay,
                       DR_Cluster_clustering_algorithm,
                       DR_Cluster_verbose) {
  
  message("[DR_Cluster] Set seed with ", DR_Cluster_seed)
  set.seed(DR_Cluster_seed)

  # PCA
  if (DR_Cluster_pca_features == "variable") {
    if (DR_Cluster_verbose) message("[DR_Cluster] Using variable features for PCA...")
    seurat_obj <- RunPCA(seurat_obj,
                         features = VariableFeatures(seurat_obj),
                         verbose = DR_Cluster_verbose)
  } else if (DR_Cluster_pca_features == "all") {
    if (DR_Cluster_verbose) message("[DR_Cluster] Using all features for PCA...")
    seurat_obj <- RunPCA(seurat_obj,
                         features = rownames(seurat_obj),
                         verbose = DR_Cluster_verbose)
  } else if (is.vector(DR_Cluster_pca_features)) {
    if (DR_Cluster_verbose) message("[DR_Cluster] Using self-defined features for PCA...")
    seurat_obj <- RunPCA(seurat_obj,
                         features = DR_Cluster_pca_features,
                         verbose = DR_Cluster_verbose)
  }

  # Auto-detect optimal PC number (Elbow method)
  if (identical(DR_Cluster_dims, "auto")) {
    if (DR_Cluster_verbose) message("[DR_Cluster] Automatically determining optimal PC number using geometric elbow method...")

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

    # Prevent mistake
    if (is.na(elbow_point) || elbow_point < 2) {
      elbow_point <- min(30, total_pcs)
    }

    DR_Cluster_dims <- 1:elbow_point

    # Save ElbowPlot to PDF
    elbow_plot <- ElbowPlot(seurat_obj)
    plot_filename <- "ElbowPlot.pdf"
    if (!is.null(seurat_obj@project.name)) {
      plot_filename <- paste0("ElbowPlot_", seurat_obj@project.name, ".pdf")
    }
    pdf(plot_filename, width = 6, height = 4)
    print(elbow_plot)
    dev.off()

    if (DR_Cluster_verbose) {
      message("[DR_Cluster] Elbow plot saved as ", plot_filename)
      message("[DR_Cluster] Auto-selected ", length(DR_Cluster_dims), 
              " PCs (1:", elbow_point, ") using geometric elbow method.")
    }
  }

  # Neighbors + Clustering
  if (DR_Cluster_verbose) message("[DR_Cluster] Using ", 
                                  ifelse(is.numeric(DR_Cluster_dims), length(DR_Cluster_dims), DR_Cluster_dims),
                                  " dimensions for FindNeighbors()...")

  seurat_obj <- FindNeighbors(seurat_obj,
                              dims = DR_Cluster_dims,
                              k.param = DR_Cluster_k_param,
                              n.trees = DR_Cluster_n_trees,
                              verbose = DR_Cluster_verbose)

  if (DR_Cluster_verbose) message("[DR_Cluster] Using resolution ", DR_Cluster_resolution, " for FindClusters()...")
  seurat_obj <- FindClusters(seurat_obj,
                             resolution = DR_Cluster_resolution,
                             algorithm = DR_Cluster_clustering_algorithm,
                             verbose = DR_Cluster_verbose)

  # Dimension Reduction
  if (tolower(DR_Cluster_reduction_method) == "umap") {
    if (DR_Cluster_verbose) message("[DR_Cluster] Using umap for RunUMAP()...")
    seurat_obj <- RunUMAP(seurat_obj,
                          assay = DR_Cluster_reduction_assay,
                          dims = DR_Cluster_dims,
                          verbose = DR_Cluster_verbose)
  } else if (tolower(DR_Cluster_reduction_method) == "tsne") {
    if (DR_Cluster_verbose) message("[DR_Cluster] Using tsne for RunTSNE()...")
    seurat_obj <- RunTSNE(seurat_obj,
                          assay = DR_Cluster_reduction_assay,
                          dims = DR_Cluster_dims,
                          verbose = DR_Cluster_verbose)
  }

  return(seurat_obj)
}
