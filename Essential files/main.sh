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

# Read YAML
config <- yaml::read_yaml(config_file)

# Establish temp folder in user workspace
base_input <- config$input_paths[[1]]
tmp_dir <- file.path(dirname(base_input), "tmp_counts")
if (!dir.exists(tmp_dir)) {
  dir.create(tmp_dir, recursive = TRUE)
}

# Input files
rds_files <- c()
index <- 1
for (input_path in config$input_paths) {
  message("Check input: ", input_path)
  out_file <- file.path(tmp_dir, paste0("counts_", index, ".rds"))

  if (dir.exists(input_path)) {
    if (file.exists(file.path(input_path, "barcodes.tsv.gz")) &&
        file.exists(file.path(input_path, "features.tsv.gz")) &&
        file.exists(file.path(input_path, "matrix.mtx.gz"))) {
      message("Detected 10X folder → Read10X()")
      saveRDS(Read10X(data.dir = input_path), out_file)
    } else {
      stop("Invalid folder format: expected 10X files.")
    }

  } else if (file.exists(input_path)) {
    ext <- tools::file_ext(input_path)
    if (ext %in% c("h5", "hdf5")) {
      message("Detected HDF5 → Read10X_h5()")
      saveRDS(Read10X_h5(filename = input_path), out_file)
    } else if (ext == "csv") {
      message("Detected CSV → read.table()")
      raw_mat <- read.table(input_path, sep = ",", header = TRUE, row.names = 1, check.names = FALSE)
      saveRDS(as(as.matrix(raw_mat), "dgCMatrix"), out_file)
    } else if (ext %in% c("txt", "tsv")) {
      message("Detected TXT/TSV → read.table()")
      raw_mat <- read.table(input_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
      saveRDS(as(as.matrix(raw_mat), "dgCMatrix"), out_file)
    } else {
      stop("Unsupported file format: ", ext)
    }
  } else {
    stop("Path does not exist: ", input_path)
  }

  rds_files <- c(rds_files, out_file)
  index <- index + 1
}

# 匯入分析函數
source("functions/load_and_qc.R")
source("functions/normalize_and_scale.R")
source("functions/DR_and_cluster.R")
source("functions/find_markers_and_annotate.R")
source("functions/plot_dim_reduction.R")
source("functions/run_cellchat.R")
source("functions/vina_docking.R")

# -------------------------
# Step 1: Load and QC
# -------------------------
seurat_obj <- load_and_qc(
  rds_files     = rds_files,
  min_features  = config$load_and_qc$min_features,
  min_cells     = config$load_and_qc$min_cells,
  names_delim   = config$load_and_qc$names_delim,
  max_mito      = config$load_and_qc$max_mito,
  load_and_qc_verbose = config$load_and_qc$load_and_qc_verbose
)

message("Saving the QC-completed file to seurat_qc.rds")
saveRDS(seurat_obj, file = "/home/andrew/docking_work/Tasks/scDock/scData/seurat_qc.rds")

# -------------------------
# Step 2: Check for normalize_and_scale
# -------------------------
if (!inherits(seurat_obj, "Seurat")) {
  stop("Error: Input is not a Seurat object. Please run load_and_qc() first.")
}

valid_hvg_methods <- c("vst", "mean.var.plot", "dispersion")
if (!(config$normalize_and_scale$hvg_method %in% valid_hvg_methods)) {
  stop(paste0("Error: Invalid hvg.method (", config$normalize_and_scale$hvg_method,
              "). Choose one of: ", paste(valid_hvg_methods, collapse = ", ")))
}

if (config$normalize_and_scale$nVariableFeatures < 100 || config$normalize_and_scale$nVariableFeatures > 10000) {
  warning("⚠️ nVariableFeatures = ", config$nVariableFeatures,
          " is unusually small or large. Please double-check.")
}

if (!(config$normalize_and_scale$scale_features %in% c("variable", "all")) && !is.vector(config$normalize_and_scale$scale_features)) {
  stop("Error: Invalid value for scale_features. Must be 'variable', 'all', or a list of gene names.")
}

# -------------------------
# Step 3: Normalize, HVG, and Scale
# -------------------------
config$normalize_and_scale$scale_factor <- as.numeric(config$normalize_and_scale$scale_factor)
config$normalize_and_scale$nVariableFeatures <- as.integer(config$normalize_and_scale$nVariableFeatures)

seurat_obj <- normalize_and_scale(
  seurat_obj             = seurat_obj,
  use_assay              = config$normalize_and_scale$use_assay,
  normalization_method   = config$normalize_and_scale$normalization_method,
  scale_factor           = config$normalize_and_scale$scale_factor,
  CLR_margin             = config$normalize_and_scale$CLR_margin,
  hvg_method             = config$normalize_and_scale$hvg_method,
  nVariableFeatures      = config$normalize_and_scale$nVariableFeatures,
  split_by               = config$normalize_and_scale$split_by,
  model_use              = config$normalize_and_scale$model_use,
  scale_max              = config$normalize_and_scale$scale_max,
  scale_features         = config$normalize_and_scale$scale_features,
  normalize_and_scale_verbose = config$normalize_and_scale$normalize_and_scale_verbose
)

message("Saving the normalized and scaled file to seurat_norm_scale.rds")
saveRDS(seurat_obj, file = "/home/andrew/docking_work/Tasks/scDock/scData/seurat_norm_scale.rds")

# -------------------------
# Step 4: Check for DR_and_cluster
# -------------------------
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
  stop("Error: Data has not been scaled. Please run normalize_and_scale() first.")
}

if (!(config$DR_and_cluster$pca_features %in% c("variable", "all")) &&
    !is.vector(config$DR_and_cluster$pca_features)) {
  stop("Error: Invalid pca_features. Must be 'variable', 'all', or a vector of genes.")
}

if (!(tolower(config$DR_and_cluster$reduction_method) %in% c("umap", "tsne"))) {
  stop("Error: Invalid reduction_method. Must be 'umap' or 'tsne'.")
}

if (!(config$DR_and_cluster$clustering_algorithm %in% 1:4)) {
  stop("Error: clustering_algorithm must be 1, 2, 3, or 4.")
}

if (!(is.numeric(config$DR_and_cluster$resolution) && config$DR_and_cluster$resolution > 0)) {
  stop("Error: resolution must be a positive number.")
}

# -------------------------
# Step 5: Run DR_and_cluster
# -------------------------
config$DR_and_cluster$dims <- seq_len(as.integer(config$DR_and_cluster$dims))

seurat_obj <- DR_and_cluster(
  seurat_obj             = seurat_obj,
  pca_features           = config$DR_and_cluster$pca_features,
  seed                   = config$DR_and_cluster$seed,
  dims                   = config$DR_and_cluster$dims,
  k_param                = config$DR_and_cluster$k_param,
  n_trees                = config$DR_and_cluster$n_trees,
  resolution             = config$DR_and_cluster$resolution,
  reduction_method       = config$DR_and_cluster$reduction_method,
  clustering_algorithm   = config$DR_and_cluster$clustering_algorithm,
  DR_and_cluster_verbose = config$DR_and_cluster$DR_and_cluster_verbose
)

message("Saving the clustered file to seurat_dr_cluster.rds")
saveRDS(seurat_obj, file = "/home/andrew/docking_work/Tasks/scDock/scData/seurat_dr_cluster.rds")

# -------------------------
# Step 6: Find markers and annotate
# -------------------------

# ---- Check ----
if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  stop("Error: 'seurat_clusters' not found. Please run DR_and_cluster() first.")
}

if (!is.numeric(config$find_markers_and_annotate$top_n) ||
    config$find_markers_and_annotate$top_n < 0) {
  stop("Error: top_n must be a non-negative number.")
}

# Check the output directory
for (out_file in c(config$find_markers_and_annotate$output_marker_csv,
                   config$find_markers_and_annotate$output_annotation_csv)) {
  if (!is.null(out_file)) {
    out_dir <- dirname(out_file)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  }
}

seurat_obj <- find_markers_and_annotate(
  seurat_obj            = seurat_obj,
  use_assay             = config$find_markers_and_annotate$use_assay,
  group_by              = config$find_markers_and_annotate$group_by,
  logfc_threshold       = config$find_markers_and_annotate$logfc_threshold,
  test_use              = config$find_markers_and_annotate$test_use,
  min_pct               = config$find_markers_and_annotate$min_pct,
  only_pos              = config$find_markers_and_annotate$only_pos,
  find_markers_and_annotate_verbose = config$find_markers_and_annotate$find_markers_and_annotate_verbose,
  top_n                 = config$find_markers_and_annotate$top_n,
  output_marker_csv     = config$find_markers_and_annotate$output_marker_csv,
  output_annotation_csv = config$find_markers_and_annotate$output_annotation_csv,
  tissue_type           = config$find_markers_and_annotate$tissue_type,
  label_column          = ifelse(
    is.null(config$find_markers_and_annotate$label_column),
    "celltype",
    config$find_markers_and_annotate$label_column
  )
)

message("Saving the annotated file to seurat_markers_annotated.rds")
saveRDS(seurat_obj, file = "/home/andrew/docking_work/Tasks/scDock/scData/seurat_markers_annotated.rds")

# -------------------------
# Step 7: Check + Dimensional Reduction Plot
# -------------------------
if (!is.null(config$plot_dim_reduction)) {
  # Check reduction method
  available_reductions <- Reductions(seurat_obj)
  reduction_method <- config$plot_dim_reduction$reduction_method
  if (!(reduction_method %in% available_reductions)) {
    warning(paste0("⚠️ Reduction '", reduction_method, "' not found in Seurat object. Skipping plot."))
    reduction_method <- NULL
  }

  # Check group_by
  group_by <- config$plot_dim_reduction$group_by
  if (!is.null(reduction_method) && !(group_by %in% colnames(seurat_obj@meta.data))) {
    warning(paste0("⚠️ Metadata column '", group_by, "' not found. Skipping plot."))
    reduction_method <- NULL
  }

  if (!is.null(reduction_method)) {
    plot_dim_reduction(
      seurat_obj = seurat_obj,
      output_path = config$plot_dim_reduction$output_path,
      reduction_method = reduction_method,
      group_by = group_by,
      pt_size = config$plot_dim_reduction$pt_size,
      label = config$plot_dim_reduction$label,
      label_size = config$plot_dim_reduction$label_size,
      label_color = config$plot_dim_reduction$label_color,
      label_box = config$plot_dim_reduction$label_box,
      alpha = config$plot_dim_reduction$alpha,
      shuffle = config$plot_dim_reduction$shuffle,
      raster = config$plot_dim_reduction$raster,
      width = config$plot_dim_reduction$width,
      height = config$plot_dim_reduction$height
    )
  }
}

# -------------------------
# Step 8: CellChat analysis
# -------------------------
if (!is.null(config$run_cellchat)) {
  seurat_obj <- readRDS("/home/andrew/docking_work/Tasks/scDock/scData/seurat_markers_annotated.rds")

  run_cellchat(
    seurat_obj      = seurat_obj,
    output_dir      = config$run_cellchat$output_dir,
    species         = config$run_cellchat$species,
    group_by        = config$run_cellchat$group_by,
    source_celltype = config$run_cellchat$source_celltype,
    target_celltype = config$run_cellchat$target_celltype,
    plot_heatmap    = ifelse(is.null(config$run_cellchat$plot_heatmap), TRUE, config$run_cellchat$plot_heatmap),
    ntop_pathway    = config$run_cellchat$ntop_pathway,
    ntop_signaling  = config$run_cellchat$ntop_signaling,
    pdf_width       = ifelse(is.null(config$run_cellchat$pdf_width), 8, config$run_cellchat$pdf_width),
    pdf_height      = ifelse(is.null(config$run_cellchat$pdf_height), 6, config$run_cellchat$pdf_height)
  )
}

# -------------------------
# Step 9: AutoDock Vina docking
# -------------------------
if (!is.null(config$vina_docking)) {
  result <- vina_docking(
    input_dir         = config$vina_docking$input_dir,
    ligand_ref_file   = config$vina_docking$ligand_ref_file,
    receptor_ref_file = config$vina_docking$receptor_ref_file,
    output_dir        = config$vina_docking$output_dir,
    python_bin        = config$vina_docking$python_bin,
    cas_txt_file      = config$vina_docking$cas_txt_file,
    docking_ligand_dir     = config$vina_docking$docking_ligand_dir,
    use_fda                = config$vina_docking$use_fda,
    fda_txt_path           = config$vina_docking$fda_txt_path,
    docking_receptor_dir   = config$vina_docking$docking_receptor_dir,
    vina_exhaustiveness     = ifelse(is.null(config$vina_docking$vina_exhaustiveness), 8, config$vina_docking$vina_exhaustiveness),
    vina_num_modes          = ifelse(is.null(config$vina_docking$vina_num_modes), 9, config$vina_docking$vina_num_modes),
    vina_seed               = ifelse(is.null(config$vina_docking$vina_seed), 42, config$vina_docking$vina_seed),
    vina_cpu                = ifelse(is.null(config$vina_docking$vina_cpu), 1, config$vina_docking$vina_cpu)
  )
}

# 清理暫存
unlink(tmp_dir, recursive = TRUE)
EOF
