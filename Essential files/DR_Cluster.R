# functions/DR_Cluster.R
DR_Cluster <- function(seurat_obj,
                           DR_Cluster_pca_features = "variable",
                           DR_Cluster_seed = 42,
                           DR_Cluster_dims = 1:20,
                           DR_Cluster_k_param = 20,
                           DR_Cluster_n_trees = 50,
                           DR_Cluster_resolution = 0.1,
                           DR_Cluster_reduction_method = "umap",
                           DR_Cluster_reduction_assay = NULL,
                           DR_Cluster_clustering_algorithm = 1,
                           DR_Cluster_verbose = TRUE) {
  
  message("Set seed with ", DR_Cluster_seed)
  set.seed(DR_Cluster_seed)

  # -------------------------
  # PCA
  # -------------------------
  if (DR_Cluster_pca_features == "variable") {
    if (DR_Cluster_verbose) message("Using variabel features for PCA...")
    seurat_obj <- RunPCA(seurat_obj,
                         features = VariableFeatures(seurat_obj),
                         verbose = DR_Cluster_verbose)
  } else if (DR_Cluster_pca_features == "all") {
    if (DR_Cluster_verbose) message("Using all features for PCA...")
    seurat_obj <- RunPCA(seurat_obj,
                         features = rownames(seurat_obj),
                         verbose = DR_Cluster_verbose)
  } else if (is.vector(DR_Cluster_pca_features)) {
    if (DR_Cluster_verbose) message("Using self-defined features for PCA...")
    seurat_obj <- RunPCA(seurat_obj,
                         features = DR_Cluster_pca_features,
                         verbose = DR_Cluster_verbose)
  }

  # -------------------------
  # Neighbors + Clustering
  # -------------------------
  if (DR_Cluster_verbose) message("Using ", DR_Cluster_dims, " dimensions for FindNeighbors()...")
  seurat_obj <- FindNeighbors(seurat_obj,
                              dims = DR_Cluster_dims,
                              k.param = DR_Cluster_k_param,
                              n.trees = DR_Cluster_n_trees,
                              verbose = DR_Cluster_verbose)

  if (DR_Cluster_verbose) message("Using resolution ", DR_Cluster_resolution, " for FindClusters()...")
  seurat_obj <- FindClusters(seurat_obj,
                             resolution = DR_Cluster_resolution,
                             algorithm = DR_Cluster_clustering_algorithm,
                             verbose = DR_Cluster_verbose)

  # -------------------------
  # Dimension Reduction
  # -------------------------
  if (tolower(DR_Cluster_reduction_method) == "umap") {
    if (DR_Cluster_verbose) message("Using umap for RunUMAP()...")
    seurat_obj <- RunUMAP(seurat_obj,
                          assay = DR_Cluster_reduction_assay,
                          dims = DR_Cluster_dims,
                          verbose = DR_Cluster_verbose)
  } else if (tolower(DR_Cluster_reduction_method) == "tsne") {
    seurat_obj <- RunTSNE(seurat_obj,
                          assay = DR_Cluster_reduction_assay,
                          dims = DR_Cluster_dims,
                          verbose = DR_Cluster_verbose)
  }
  return(seurat_obj)
}