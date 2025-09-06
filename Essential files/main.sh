#!/bin/bash

set -e

CONFIG_FILE=$1
if [ -z "$CONFIG_FILE" ]; then
  echo "Please provide path to config.yaml"
  exit 1
fi

Rscript --vanilla - "$CONFIG_FILE" <<'EOF'

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]

options(
  warn = -1,          # Don't display warning
  cli.unicode = FALSE # Seurat info/warning
)

suppressPackageStartupMessages({
  library(yaml)
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(scMayoMap)
  library(ggplot2)
  library(CellChat)
  library(patchwork)
  library(igraph)
})

# 讀取 YAML
config <- yaml::read_yaml(config_file)

# -------------------------
# Step 0: Input files
# -------------------------
input_files <- config$input_paths
if (length(input_files) == 0) {
  stop("Error: No input paths specified in config.yaml")
}

# 匯入分析函數
source("functions/Load_QC.R")
source("functions/Normalization_Scale.R")
source("functions/DR_Cluster.R")
source("functions/Markers_Annotation.R")
source("functions/DR_Plot.R")
source("functions/Run_CellChat.R")
source("functions/Vina_Docking.R")

# -------------------------
# Step 1: Load and QC
# -------------------------
message("========== Starting loading data and QC ... ==========")
seurat_obj <- Load_QC(
  Load_QC_input_files   = config$Load_QC$Load_QC_input_files,
  Load_QC_input_type    = config$Load_QC$Load_QC_input_type,
  Load_QC_min_features  = config$Load_QC$Load_QC_min_features,
  Load_QC_min_cells     = config$Load_QC$Load_QC_min_cells,
  Load_QC_names_delim   = config$Load_QC$Load_QC_names_delim,
  Load_QC_max_mito      = config$Load_QC$Load_QC_max_mito,
  Load_QC_metadata_file = config$Load_QC$Load_QC_metadata_file,
  Load_QC_verbose       = config$Load_QC$Load_QC_verbose
)

base_input <- dirname(input_files[[1]])
message("Saving the QC-completed file to seurat_qc.rds")
saveRDS(seurat_obj, file = file.path(base_input, "seurat_qc.rds"))
message("========== QC process complete. ==========")
message("")


# -------------------------
# Step 2: Check for Normalization_Scale
# -------------------------
if (!inherits(seurat_obj, "Seurat")) {
  stop("Error: Input is not a Seurat object. Please run Load_QC() first.")
}

valid_hvg_methods <- c("vst", "mean.var.plot", "dispersion")
if (!(config$Normalization_Scale$Normalization_Scale_hvg_method %in% valid_hvg_methods)) {
  stop(paste0("Error: Invalid hvg.method (", config$Normalization_Scale$Normalization_Scale_hvg_method,
              "). Choose one of: ", paste(valid_hvg_methods, collapse = ", ")))
}

if (config$Normalization_Scale$Normalization_Scale_nVariableFeatures < 100 || config$Normalization_Scale$Normalization_Scale_nVariableFeatures > 10000) {
  warning("⚠️ nVariableFeatures = ", config$Normalization_Scale_nVariableFeatures,
          " is unusually small or large. Please double-check.")
}

if (!(config$Normalization_Scale$Normalization_Scale_scale_features %in% c("variable", "all")) && !is.vector(config$Normalization_Scale$Normalization_Scale_scale_features)) {
  stop("Error: Invalid value for scale_features. Must be 'variable', 'all', or a list of gene names.")
}

# -------------------------
# Step 3: Normalize, HVG, and Scale
# -------------------------
message("========== Starting normalization and scaling... ==========")
config$Normalization_Scale$Normalization_Scale_scale_factor <- as.numeric(config$Normalization_Scale$Normalization_Scale_scale_factor)
config$Normalization_Scale$Normalization_Scale_nVariableFeatures <- as.integer(config$Normalization_Scale$Normalization_Scale_nVariableFeatures)

seurat_obj <- Normalization_Scale(
  seurat_obj                                 = seurat_obj,
  Normalization_Scale_use_assay              = config$Normalization_Scale$Normalization_Scale_use_assay,
  Normalization_Scale_normalization_method   = config$Normalization_Scale$Normalization_Scale_normalization_method,
  Normalization_Scale_scale_factor           = config$Normalization_Scale$Normalization_Scale_scale_factor,
  Normalization_Scale_CLR_margin             = config$Normalization_Scale$Normalization_Scale_CLR_margin,
  Normalization_Scale_hvg_method             = config$Normalization_Scale$Normalization_Scale_hvg_method,
  Normalization_Scale_nVariableFeatures      = config$Normalization_Scale$Normalization_Scale_nVariableFeatures,
  Normalization_Scale_split_by               = config$Normalization_Scale$Normalization_Scale_split_by,
  Normalization_Scale_model_use              = config$Normalization_Scale$Normalization_Scale_model_use,
  Normalization_Scale_scale_max              = config$Normalization_Scale$Normalization_Scale_scale_max,
  Normalization_Scale_scale_features         = config$Normalization_Scale$Normalization_Scale_scale_features,
  Normalization_Scale_verbose                = config$Normalization_Scale$Normalization_Scale_verbose
)

message("Saving the normalized and scaled file to seurat_norm_scale.rds")
saveRDS(seurat_obj, file = file.path(base_input, "seurat_norm_scale.rds"))
message("========== Normalization and scaling complete ==========")
message("")

# -------------------------
# Step 4: Check for DR_Cluster
# -------------------------
message("========== Starting dimensional reduction and clustering... ==========")
if (!inherits(seurat_obj, "Seurat")) {
  stop("Error: Input is not a Seurat object. Please run previous steps first.")
}

scaled_ok <- FALSE
if ("scale.data" %in% slotNames(seurat_obj@assays$RNA)) {
  # Seurat v4
  scaled_ok <- length(seurat_obj@assays$RNA@scale.data) > 0
} else {
  # Seurat v5
  scaled_ok <- length(LayerData(seurat_obj[["RNA"]], "scale.data")) > 0
}
if (!scaled_ok) {
  stop("Error: Data has not been scaled. Please run Normalization_Scale() first.")
}

if (!(config$DR_Cluster$DR_Cluster_pca_features %in% c("variable", "all")) &&
    !is.vector(config$DR_Cluster$DR_Cluster_pca_features)) {
  stop("Error: Invalid DR_Cluster_pca_features. Must be 'variable', 'all', or a vector of genes.")
}

if (!(tolower(config$DR_Cluster$DR_Cluster_reduction_method) %in% c("umap", "tsne"))) {
  stop("Error: Invalid DR_Cluster_reduction_method. Must be 'umap' or 'tsne'.")
}

if (!(config$DR_Cluster$DR_Cluster_clustering_algorithm %in% 1:4)) {
  stop("Error: DR_Cluster_clustering_algorithm must be 1, 2, 3, or 4.")
}

if (!(is.numeric(config$DR_Cluster$DR_Cluster_resolution) && config$DR_Cluster$DR_Cluster_resolution > 0)) {
  stop("Error: resolution must be a positive number.")
}

# -------------------------
# Step 5: Run DR_Cluster
# -------------------------
config$DR_Cluster$DR_Cluster_dims <- seq_len(as.integer(config$DR_Cluster$DR_Cluster_dims))

seurat_obj <- DR_Cluster(
  seurat_obj                        = seurat_obj,
  DR_Cluster_pca_features           = config$DR_Cluster$DR_Cluster_pca_features,
  DR_Cluster_seed                   = config$DR_Cluster$DR_Cluster_seed,
  DR_Cluster_dims                   = config$DR_Cluster$DR_Cluster_dims,
  DR_Cluster_k_param                = config$DR_Cluster$DR_Cluster_k_param,
  DR_Cluster_n_trees                = config$DR_Cluster$DR_Cluster_n_trees,
  DR_Cluster_resolution             = config$DR_Cluster$DR_Cluster_resolution,
  DR_Cluster_reduction_method       = config$DR_Cluster$DR_Cluster_reduction_method,
  DR_Cluster_reduction_assay        = config$DR_Cluster$DR_Cluster_reduction_assay,
  DR_Cluster_clustering_algorithm   = config$DR_Cluster$DR_Cluster_clustering_algorithm,
  DR_Cluster_verbose                = config$DR_Cluster$DR_Cluster_verbose
)

message("Saving the clustered file to seurat_dr_cluster.rds")
saveRDS(seurat_obj, file = file.path(base_input, "seurat_dr_cluster.rds"))
message("========== Dimensional reduction and clustering complete. ==========")
message("")

# -------------------------
# Step 6: Find markers and annotate
# -------------------------
message("========== Starting cluster markers calculating... ==========")

# ---- Check ----
if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  stop("Error: 'seurat_clusters' not found. Please run DR_Cluster() first.")
}

if (!is.numeric(config$Markers_Annotation$Markers_Annotation_top_n) ||
    config$Markers_Annotation$Markers_Annotation_top_n < 0) {
  stop("Error: Markers_Annotation_top_n must be a non-negative number.")
}

# Check the output directory
for (out_file in c(config$Markers_Annotation$Markers_Annotation_output_marker_csv,
                   config$Markers_Annotation$Markers_Annotation_output_annotation_csv)) {
  if (!is.null(out_file)) {
    out_dir <- dirname(out_file)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  }
}

seurat_obj <- Markers_Annotation(
  seurat_obj                               = seurat_obj,
  Markers_Annotation_use_assay             = config$Markers_Annotation$Markers_Annotation_use_assay,
  Markers_Annotation_group_by              = config$Markers_Annotation$Markers_Annotation_group_by,
  Markers_Annotation_logfc_threshold       = config$Markers_Annotation$Markers_Annotation_logfc_threshold,
  Markers_Annotation_test_use              = config$Markers_Annotation$Markers_Annotation_test_use,
  Markers_Annotation_min_pct               = config$Markers_Annotation$Markers_Annotation_min_pct,
  Markers_Annotation_only_pos              = config$Markers_Annotation$Markers_Annotation_only_pos,
  Markers_Annotation_verbose               = config$Markers_Annotation$Markers_Annotation_verbose,
  Markers_Annotation_top_n                 = config$Markers_Annotation$Markers_Annotation_top_n,
  Markers_Annotation_output_marker_csv     = config$Markers_Annotation$Markers_Annotation_output_marker_csv,
  Markers_Annotation_output_annotation_csv = config$Markers_Annotation$Markers_Annotation_output_annotation_csv,
  Markers_Annotation_tissue_type           = config$Markers_Annotation$Markers_Annotation_tissue_type,
  Markers_Annotation_label_column          = ifelse(
    is.null(config$Markers_Annotation$Markers_Annotation_label_column),
    "celltype",
    config$Markers_Annotation$Markers_Annotation_label_column
  )
)

message("Saving the annotated file to seurat_markers_annotated.rds")
saveRDS(seurat_obj, file = file.path(base_input, "seurat_markers_annotated.rds"))
message("========== Find marker and annotation complete. ==========")
message("")

# -------------------------
# Step 7: Check + Dimensional Reduction Plot
# -------------------------
message("========== Starting visualization with Dimplot... ==========")

if (!is.null(config$DR_Plot)) {
  # Check reduction method
  available_reductions <- Reductions(seurat_obj)
  reduction_method <- config$DR_Plot$DR_Plot_reduction_method
  if (!(reduction_method %in% available_reductions)) {
    warning(paste0("⚠️ Reduction '", reduction_method, "' not found in Seurat object. Skipping plot."))
    reduction_method <- NULL
  }

  # Check group_by
  group_by <- config$DR_Plot$DR_Plot_group_by
  if (!is.null(reduction_method) && !(group_by %in% colnames(seurat_obj@meta.data))) {
    warning(paste0("⚠️ Metadata column '", group_by, "' not found. Skipping plot."))
    reduction_method <- NULL
  }

  if (!is.null(reduction_method)) {
    DR_Plot(
      seurat_obj               = seurat_obj,
      DR_Plot_output_path      = config$DR_Plot$DR_Plot_output_path,
      DR_Plot_reduction_method = reduction_method,
      DR_Plot_group_by         = group_by,
      DR_Plot_pt_size          = config$DR_Plot$DR_Plot_pt_size,
      DR_Plot_label            = config$DR_Plot$DR_Plot_label,
      DR_Plot_label_size       = config$DR_Plot$DR_Plot_label_size,
      DR_Plot_label_color      = config$DR_Plot$DR_Plot_label_color,
      DR_Plot_label_box        = config$DR_Plot$DR_Plot_label_box,
      DR_Plot_alpha            = config$DR_Plot$DR_Plot_alpha,
      DR_Plot_shuffle          = config$DR_Plot$DR_Plot_shuffle,
      DR_Plot_raster           = config$DR_Plot$DR_Plot_raster,
      DR_Plot_width            = config$DR_Plot$DR_Plot_width,
      DR_Plot_height           = config$DR_Plot$DR_Plot_height
    )
  }
}
message("========== Visualization complete. ==========")
message("")

# -------------------------
# Step 8: CellChat analysis
# -------------------------
message("========== Running CellChat on Seurat object ... ==========")

if (!is.null(config$Run_CellChat)) {
  seurat_obj <- readRDS(file.path(base_input, "seurat_markers_annotated.rds"))

  Run_CellChat(
    seurat_obj                   = seurat_obj,
    Run_CellChat_output_dir      = config$Run_CellChat$Run_CellChat_output_dir,
    Run_CellChat_species         = config$Run_CellChat$Run_CellChat_species,
    Run_CellChat_group_by        = config$Run_CellChat$Run_CellChat_group_by,
    Run_CellChat_source_celltype = config$Run_CellChat$Run_CellChat_source_celltype,
    Run_CellChat_target_celltype = config$Run_CellChat$Run_CellChat_target_celltype,
    Run_CellChat_plot_heatmap    = ifelse(is.null(config$Run_CellChat$Run_CellChat_plot_heatmap), TRUE, config$Run_CellChat$Run_CellChat_plot_heatmap),
    Run_CellChat_ntop_pathway    = config$Run_CellChat$Run_CellChat_ntop_pathway,
    Run_CellChat_ntop_signaling  = config$Run_CellChat$Run_CellChat_ntop_signaling,
    Run_CellChat_pdf_width       = ifelse(is.null(config$Run_CellChat$Run_CellChat_pdf_width), 8, config$Run_CellChat$Run_CellChat_pdf_width),
    Run_CellChat_pdf_height      = ifelse(is.null(config$Run_CellChat$Run_CellChat_pdf_height), 6, config$Run_CellChat$Run_CellChat_pdf_height)
  )
}
message("========== CellChat complete. ==========")
message("")

# -------------------------
# Step 9: AutoDock Vina docking
# -------------------------
message("========== Starting docking... ==========")

if (!is.null(config$Vina_Docking)) {
  result <- Vina_Docking(
    Vina_Docking_input_dir              = config$Vina_Docking$Vina_Docking_input_dir,
    Vina_Docking_ligand_ref_file        = config$Vina_Docking$Vina_Docking_ligand_ref_file,
    Vina_Docking_receptor_ref_file      = config$Vina_Docking$Vina_Docking_receptor_ref_file,
    Vina_Docking_output_dir             = config$Vina_Docking$Vina_Docking_output_dir,
    Vina_Docking_cas_txt_file           = config$Vina_Docking$Vina_Docking_cas_txt_file,
    Vina_Docking_docking_ligand_dir     = config$Vina_Docking$Vina_Docking_docking_ligand_dir,
    Vina_Docking_use_fda                = config$Vina_Docking$Vina_Docking_use_fda,
    Vina_Docking_fda_txt_path           = config$Vina_Docking$Vina_Docking_fda_txt_path,
    Vina_Docking_docking_receptor_dir   = config$Vina_Docking$Vina_Docking_docking_receptor_dir,
    Vina_Docking_vina_exhaustiveness    = ifelse(is.null(config$Vina_Docking$Vina_Docking_vina_exhaustiveness), 8, config$Vina_Docking$Vina_Docking_vina_exhaustiveness),
    Vina_Docking_vina_num_modes         = ifelse(is.null(config$Vina_Docking$Vina_Docking_vina_num_modes), 9, config$Vina_Docking$Vina_Docking_vina_num_modes),
    Vina_Docking_vina_seed              = ifelse(is.null(config$Vina_Docking$Vina_Docking_vina_seed), 42, config$Vina_Docking$Vina_Docking_vina_seed),
    Vina_Docking_vina_cpu               = ifelse(is.null(config$Vina_Docking$Vina_Docking_vina_cpu), 1, config$Vina_Docking$Vina_Docking_vina_cpu)
  )
}
message("========== Docking complete. ==========")
EOF