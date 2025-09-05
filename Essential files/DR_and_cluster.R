# functions/DR_and_cluster.R

DR_and_cluster <- function(seurat_obj,
                           pca_features = "variable",
                           seed = 42,
                           dims = 1:20,
                           k_param = 20,
                           n_trees = 50,
                           resolution = 0.1,
                           reduction_method = "umap",
                           reduction_assay = NULL,
                           clustering_algorithm = 1,
                           DR_and_cluster_verbose = TRUE) {
  
  message("Starting dimensional reduction and clustering...")
  message("Set seed with", seed)
  set.seed(seed)

  # -------------------------
  # PCA
  # -------------------------
  if (pca_features == "variable") {
    if (DR_and_cluster_verbose) message("Using variabel features for PCA...")
    seurat_obj <- RunPCA(seurat_obj,
                         features = VariableFeatures(seurat_obj),
                         verbose = DR_and_cluster_verbose)
  } else if (pca_features == "all") {
    if (DR_and_cluster_verbose) message("Using all features for PCA...")
    seurat_obj <- RunPCA(seurat_obj,
                         features = rownames(seurat_obj),
                         verbose = DR_and_cluster_verbose)
  } else if (is.vector(pca_features)) {
    if (DR_and_cluster_verbose) message("Using self-defined features for PCA...")
    seurat_obj <- RunPCA(seurat_obj,
                         features = pca_features,
                         verbose = DR_and_cluster_verbose)
  }

  # -------------------------
  # Neighbors + Clustering
  # -------------------------
  if (DR_and_cluster_verbose) message("Using ", dims, " dimensions for FindNeighbors()...")
  seurat_obj <- FindNeighbors(seurat_obj,
                              dims = dims,
                              k.param = k_param,
                              n.trees = n_trees,
                              verbose = DR_and_cluster_verbose)

  if (DR_and_cluster_verbose) message("Using resolution ", resolution, " for FindClusters()...")
  seurat_obj <- FindClusters(seurat_obj,
                             resolution = resolution,
                             algorithm = clustering_algorithm,
                             verbose = DR_and_cluster_verbose)

  # -------------------------
  # Dimension Reduction
  # -------------------------
  if (tolower(reduction_method) == "umap") {
    if (DR_and_cluster_verbose) message("Using umap for RunUMAP()...")
    seurat_obj <- RunUMAP(seurat_obj,
                          assay = reduction_assay,
                          dims = dims,
                          verbose = DR_and_cluster_verbose)
  } else if (tolower(reduction_method) == "tsne") {
    seurat_obj <- RunTSNE(seurat_obj,
                          assay = reduction_assay,
                          dims = dims,
                          verbose = DR_and_cluster_verbose)
  }

  message("✅ Dimensional reduction and clustering complete.")
  return(seurat_obj)
}