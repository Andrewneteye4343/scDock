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
  library(ComplexHeatmap)
})

# Read YAML
config <- yaml::read_yaml(config_file)

# -------------------------
# Validation function
# -------------------------
check_string_path <- function(value, arg_name = "path") {
  if (is.null(value)) {
    stop(paste0("Error: ", arg_name, " is missing in config.yaml"))
  }
  if (!is.character(value)) {
    stop(paste0("Error: ", arg_name, " must be a string or a vector of strings"))
  }
  bad <- value[!grepl("^/", value)]
  if (length(bad) > 0) {
    stop(paste0("Error: ", arg_name, " contains invalid path(s): ", paste(bad, collapse = ", ")))
  }
  return(value)
}

check_string_path_optional <- function(value, arg_name = "path") {
  if (is.null(value)) {
    return(value)
  }
  if (!is.character(value)) {
    stop(paste0("Error: ", arg_name, " must be a string or a vector of strings"))
  }
  bad <- value[!grepl("^/", value)]
  if (length(bad) > 0) {
    stop(paste0("Error: ", arg_name, " contains invalid path(s): ", paste(bad, collapse = ", ")))
  }
  return(value)
}

check_string <- function(value, arg_name = "value") { 
  if (!is.null(value)) {
    if (!is.character(value)) {
      stop(paste0("Error: ", arg_name, " must be NULL or a string"))
    }
    if (!grepl("^[A-Za-z0-9._-]+$", value)) {
      stop(paste0("Error: ", arg_name, " must be NULL or a single string (no slashes allowed). Got: ", value))
    }
    return(value)
  }
}

check_integer_string <- function(value, arg_name = "value", default = 10) {
  if (is.null(value)) {
    message(paste0("Warning: ", arg_name, " is missing in config.yaml"))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  bad <- value[!grepl("^-?[0-9]+$", value)]
  if (length(bad) > 0) {
    warning(paste0("Warning: ", arg_name, " must be an integer string. Got invalid: ", paste(bad, collapse = ", ")))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  return(value)
}

check_numeric_string <- function(value, arg_name = "value", default = 10.0) {
  if (is.null(value)) {
    warning(paste0("Warning: ", arg_name, " is missing"))
    message(paste0("Using default value ", default, " for ", arg_name))
    return(default)
  }
  pattern <- "^-?[0-9]+(\\.[0-9]+)?$"
  if (!grepl(pattern, value)) {
    warning(paste0("Warning: ", arg_name, " must be a numeric string (integer or float). Got: ", value))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  return(value)
}

check_boolean_string <- function(value, arg_name = "value", default = FALSE) {
  if (is.null(value)){
    warning(paste0("Warning: ", arg_name, " must not be NULL"))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  pattern <- "^(?i)(TRUE|FALSE)$"
  if (!grepl(pattern, value)) {
    warning(paste0("Warning: ", arg_name, " must be either TRUE or FALSE. Got: ", value))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  return(value)
}

check_boolean_string_optional <- function(value, arg_name = "value", default = FALSE) {
  if (is.null(value)){
    return(value)
  }
  pattern <- "^(?i)(TRUE|FALSE)$"
  if (!grepl(pattern, value)) {
    stop(paste0("Error: ", arg_name, " must be NULL, TRUE or FALSE. Got: ", value))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  return(value)
}

check_split_by <- function(value, arg_name = "split.by", default = NULL) {
  if (is.null(value)) {
    return(value)
  }
  if (!is.character(value)) {
    warning(paste0("Warning: ", arg_name, " must be a string or NULL"))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  if (!grepl("^[A-Za-z][A-Za-z0-9._]*$", value)) {
    warning(paste0("Warning: ", arg_name, " must be a valid metadata column name. Got: ", value))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  } 
  return(value)
}

check_group_by <- function(value, arg_name = "group.by", default = "seurat_clusters") {
  if (!is.character(value)) {
    warning(paste0("Warning: ", arg_name, " must be a string"))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  if (!grepl("^[A-Za-z][A-Za-z0-9._]*$", value)) {
    warning(paste0("Warning: ", arg_name, " must be a valid metadata column name. Got: ", value))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  } 
  return(value)
}

check_features <- function(value, arg_name = "features", default = "variable") {
  if (is.null(value)) {
    return(value)
  }
  if (!is.character(value)) {
    warning(paste0("Warning: ", arg_name, " must be a string or character vector"))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  if (length(value) == 1 && value %in% c("variable", "all")) {
    return(value)
  }
  invalid <- value[!grepl("^[A-Za-z0-9._:-]+$", value)]
  if (length(invalid) > 0) {
    warning(paste0("Warning: ", arg_name, " contains invalid gene names: ", paste(invalid, collapse = ", ")))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  return(value)
}

check_assay <- function(value, arg_name = "assay", default = "RNA") {
  if (is.null(value)) {
    return(value)
  }
  if (!is.character(value)) {
    warning(paste0("Warning: ", arg_name, " must be a string or NULL"))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  if (!grepl("^[A-Za-z][A-Za-z0-9._]*$", value)) {
    warning(paste0("Warning: ", arg_name, " must be a valid assay name (letters, numbers, ., _; cannot start with number). Got: ", value))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  return(value)
}

check_txt_path <- function(value, arg_name = "file path") {
  if (is.null(value)) {
    return(value)
  }
  if (!is.character(value)) {
    warning(paste0("Warning: ", arg_name, " must be a string or NULL"))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  if (!grepl("\\.txt$", value, ignore.case = FALSE)) {
    warning(paste0("Warning: ", arg_name, " must end with '.txt'. Got: ", value))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  return(value)
}

check_csv_path <- function(value, arg_name = "file path") {
  if (is.null(value)) {
    warning(paste0("Warning: ", arg_name, " is missing in config.yaml"))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }  
  if (!is.character(value)) {
    warning(paste0("Warning: ", arg_name, " must be a string"))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  if (!grepl("^.+\\.csv$", value, ignore.case = TRUE)) {
    warning(paste0("Warning: ", arg_name, " must be a valid path ending with '.csv'. Got: ", value))
    message(paste0("Using default value: ", default, " for ", arg_name))
    return(default)
  }
  return(value)
}

validate_config <- function(config) {
  # --- work_path ---
  arg_work_path <- check_string_path(config$work_path, arg_name = "work_path")
  
  # --- Seurat_output_dir ---
  if (is.null(config$Seurat_output_dir)) {
    stop(paste0("Error: Seurat_output_dir is missing in config.yaml"))
  }
  arg_Seurat_output_dir <- check_string(config$Seurat_output_dir, "Seurat_output_dir")

  arg_Seurat_output_path <- file.path(config$work_path, arg_Seurat_output_dir)
  if (!dir.exists(arg_Seurat_output_path)) {
    dir.create(arg_Seurat_output_path, recursive = TRUE)
  }

  # --- Load_QC ---
  qc <- config$Load_QC
  if (is.null(qc)) {
    stop(paste0("Error: Load_QC block is missing in config.yaml"))
  }

  # Load_QC_input_files
  arg_Load_QC_input_files <- check_string_path(qc$Load_QC_input_files, arg_name = "Load_QC_input_files")

  # Load_QC_input_type
  valid_types <- c("10x", "multi10x", "h5", "txt")
  if (!(qc$Load_QC_input_type %in% valid_types)) {
    stop(paste0("Error: Load_QC_input_type must be one of: ", paste(valid_types, collapse = ", ")))
  } else {
    arg_Load_QC_input_type <- qc$Load_QC_input_type
  }

  # Load_QC_species
  valid_types <- c("human","mouse")
  if (!(qc$Load_QC_species %in% valid_types)) {
    stop(paste0("Error: Load_QC_species must be one of: ", paste(valid_types, collapse = ", ")))
  } else {
    arg_Load_QC_species <- qc$Load_QC_species
  }

  # Load_QC_min_features
  arg_Load_QC_min_features <- check_integer_string(qc$Load_QC_min_features, "Load_QC_min_features", default = 200)

  # Load_QC_min_cells
  arg_Load_QC_min_cells <- check_integer_string(qc$Load_QC_min_cells, "Load_QC_min_cells", default = 3)

  # Load_QC_names_delim
  if (!is.character(qc$Load_QC_names_delim)) {
    warning(paste0("Warning: Load_QC_names_delim must be a string"))
    message(paste0("Using default value: _ for Load_QC_names_delim"))
    arg_Load_QC_names_delim <- "_"
  }
  if (grepl("^[._-]+$", qc$Load_QC_names_delim)) {
    arg_Load_QC_names_delim <- qc$Load_QC_names_delim
  } else {
    warning(paste0("Warning: Load_QC_names_delim must be a symbol string like . or _ or -"))
    message(paste0("Using default value: _ for Load_QC_names_delim"))
    arg_Load_QC_names_delim <- "_"
    }

  # Load_QC_max_mito
  arg_Load_QC_max_mito <- check_numeric_string(qc$Load_QC_max_mito, arg_name = "Load_QC_max_mito", default = 0.1)

  # Load_QC_verbose
  arg_Load_QC_verbose <- check_boolean_string(qc$Load_QC_verbose, arg_name = "Load_QC_verbose", default = TRUE)

  # Load_QC_metadata_file
  arg_Load_QC_metadata_file <- check_txt_path(qc$Load_QC_metadata_file, arg_name = "Load_QC_metadata_file")
  
  # --- Normalization_Scale ---
  ns <- config$Normalization_Scale
  if (is.null(ns)) {
    stop(paste0("Error: Normalization_Scale block is missing in config.yaml"))
  }

  # Normalization_Scale_use_assay
  arg_Normalization_Scale_use_assay <- check_assay(ns$Normalization_Scale_use_assay, arg_name = "Normalization_Scale_use_assay", default = "RNA")
  

  # Normalization_Scale_normalization_method
  valid_normalization_method <- c("LogNormalize","CLR","RC")
  if (!(ns$Normalization_Scale_normalization_method %in% valid_normalization_method)) {
    warning(paste0("Warning: Normalization_Scale_normalization_method must be one of: ", paste(valid_normalization_method, collapse = ", ")))
    message(paste0("Using default value: LogNormalize for Normalization_Scale_normalization_method"))
    arg_Normalization_Scale_normalization_method <- "LogNormalize"
  } else {
    arg_Normalization_Scale_normalization_method <- ns$Normalization_Scale_normalization_method
  }

  # Normalization_Scale_CLR_margin
  valid_CLR_margin <- c(1,2)
  if ((ns$Normalization_Scale_normalization_method == "CLR") && (!(ns$Normalization_Scale_CLR_margin %in% valid_CLR_margin))) {
    warning(paste0("Warning: Normalization_Scale_CLR_margin must be numeric"))
    message(paste0("Using default value: 1 [for features] for Normalization_Scale_CLR_margin"))
    arg_Normalization_Scale_CLR_margin <- 1
  } else {
    arg_Normalization_Scale_CLR_margin <- ns$Normalization_Scale_CLR_margin
  }

  # Normalization_Scale_scale_factor
  arg_Normalization_Scale_scale_factor <- check_numeric_string(ns$Normalization_Scale_scale_factor, arg_name = "Normalization_Scale_scale_factor", default = 10000)

  # Normalization_Scale_hvg_method
  valid_hvg_method <- c("vst", "mean.var.plot", "dispersion")
  if (!(ns$Normalization_Scale_hvg_method %in% valid_hvg_method)) {
    warning(paste0("Warning: Normalization_Scale_hvg_method must be one of: ", paste(valid_hvg_method, collapse = ", ")))
    message(paste0("Using default value: vst for Normalization_Scale_hvg_method"))
    arg_Normalization_Scale_hvg_method <- "vst"
  } else {
    arg_Normalization_Scale_hvg_method <- ns$Normalization_Scale_hvg_method
  }

  # Normalization_Scale_nVariableFeatures
  arg_Normalization_Scale_nVariableFeatures <- check_integer_string(ns$Normalization_Scale_nVariableFeatures, arg_name = "Normalization_Scale_nVariableFeatures", default = 2000)
  
  # Normalization_Scale_split_by
  arg_Normalization_Scale_split_by <- check_split_by(ns$Normalization_Scale_split_by, arg_name = "Normalization_Scale_split_by", default = NULL)

  # Normalization_Scale_model_use
  valid_model_use <- c("linear", "poisson", "negbinom")
  if (!(ns$Normalization_Scale_model_use %in% valid_model_use)) {
    warning(paste0("Warning: Normalization_Scale_model_use must be one of: ", paste(valid_model_use, collapse = ", ")))
    message(paste0("Using default value: linear for Normalization_Scale_model_use"))
    arg_Normalization_Scale_model_use <- "linear"
  } else {
    arg_Normalization_Scale_model_use <- ns$Normalization_Scale_model_use
  }
  
  # Normalization_Scale_scale_max
  arg_Normalization_Scale_scale_max <- check_numeric_string(ns$Normalization_Scale_scale_max, arg_name = "Normalization_Scale_scale_max", default = 10)

  # Normalization_Scale_scale_features
  arg_Normalization_Scale_scale_features <- check_features(ns$Normalization_Scale_scale_features, arg_name = "Normalization_Scale_scale_features", default = "variable")

  # Normalization_Scale_verbose
  arg_Normalization_Scale_verbose <- check_boolean_string(ns$Normalization_Scale_verbose, arg_name = "Normalization_Scale_verbose", default = TRUE)
  
  # --- DR_Cluster ---
  drc <- config$DR_Cluster
  if (is.null(drc)) {
    stop("Error: DR_Cluster block is missing in config.yaml")
  }

  # DR_Cluster_pca_features
  arg_DR_Cluster_pca_features <- check_features(drc$DR_Cluster_pca_features, arg_name = "DR_Cluster_pca_features", default = "variable")

  # DR_Cluster_seed
  arg_DR_Cluster_seed <- check_integer_string(drc$DR_Cluster_seed, arg_name = "DR_Cluster_seed", default = 42)

  # DR_Cluster_dims
  arg_DR_Cluster_dims <- check_integer_string(drc$DR_Cluster_dims, arg_name = "DR_Cluster_dims", default = "auto")

  # DR_Cluster_k_param
  arg_DR_Cluster_k_param <- check_integer_string(drc$DR_Cluster_k_param, arg_name = "DR_Cluster_k_param", default = 20)

  # DR_Cluster_n_trees
  arg_DR_Cluster_n_trees <- check_integer_string(drc$DR_Cluster_n_trees, arg_name = "DR_Cluster_n_trees", default = 50)

  # DR_Cluster_resolution
  arg_DR_Cluster_resolution <- check_numeric_string(drc$DR_Cluster_resolution, arg_name = "DR_Cluster_resolution", default = 0.5)

  # DR_Cluster_reduction_method
  valid_reduction_method <- c("umap","tsne")
  if (!(drc$DR_Cluster_reduction_method %in% valid_reduction_method)) {
    warning(paste0("Error: DR_Cluster_reduction_method must be one of: ", paste(valid_reduction_method, collapse = ", ")))
    message(paste0("Using default value: umap for DR_Cluster_reduction_method"))
    arg_DR_Cluster_reduction_method <- "umap"
  } else {
    arg_DR_Cluster_reduction_method <- drc$DR_Cluster_reduction_method
  }

  # DR_Cluster_reduction_assay
  arg_DR_Cluster_reduction_assay <- check_assay(drc$DR_Cluster_reduction_assay, arg_name = "DR_Cluster_reduction_assay", default = "RNA")

  # DR_Cluster_clustering_algorithm
  valid_algorithm <- c(1,2,3,4)
  if (!(drc$DR_Cluster_clustering_algorithm %in% valid_algorithm)) {
    warning(paste0("Warning: DR_Cluster_clustering_algorithm must be one of: ", paste(valid_algorithm, collapse = ", ")))
    message(paste0("Using default value: 1 [Louvain] for DR_Cluster_clustering_algorithm"))
    arg_DR_Cluster_clustering_algorithm <- 1
  } else {
    arg_DR_Cluster_clustering_algorithm <- drc$DR_Cluster_clustering_algorithm
  }

  # DR_Cluster_verbose
  arg_DR_Cluster_verbose <- check_boolean_string(drc$DR_Cluster_verbose, arg_name = "DR_Cluster_verbose")

  # --- Run_Integration ---
  ri <- config$Run_Integration
  if (is.null(ri)) stop("Error: Run_Integration block is missing in config.yaml")

  # Run_Integration_run_integration
  arg_Run_Integration_run_integration <- check_boolean_string(ri$Run_Integration_run_integration, arg_name = "Run_Integration_run_integration")

  # Run_Integration_group_by
  arg_Run_Integration_group_by <- check_group_by(ri$Run_Integration_group_by, arg_name = "Run_Integration_group_by", default = "orig.ident")

  # Run_Integration_method
  valid_integration_method <- c("cca","rpca","harmony")
  if (!(ri$Run_Integration_method %in% valid_integration_method)) {
    warning(paste0("Error: Run_Integration_method must be one of: ", paste(valid_integration_method, collapse = ", ")))
    message(paste0("Using default value: cca for Run_Integration_method"))
    arg_Run_Integration_method <- "cca"
  } else {
    arg_Run_Integration_method <- ri$Run_Integration_method
  }

  # Run_Integration_nfeatures
  arg_Run_Integration_nfeatures <- check_integer_string(ri$Run_Integration_nfeatures, "Run_Integration_nfeatures", default = 2000)

  # Run_Integration_normalization_method
  valid_integration_normalization_method <- c("SCT","LogNormalize")
  if (!(ri$Run_Integration_normalization_method %in% valid_integration_normalization_method)) {
    warning(paste0("Error: Run_Integration_normalization_method must be one of: ", paste(valid_integration_normalization_method, collapse = ", ")))
    message(paste0("Using default value: umap for Run_Integration_normalization_method"))
    arg_Run_Integration_normalization_method <- "SCT"
  } else {
    arg_Run_Integration_normalization_method <- ri$Run_Integration_normalization_method
  }

  # Run_Integration_dims
  arg_Run_Integration_dims <- check_integer_string(ri$Run_Integration_dims, "Run_Integration_dims", default = "auto")
  
  # Run_Integration_resolution
  arg_Run_Integration_resolution <- check_numeric_string(ri$Run_Integration_resolution, arg_name = "Run_Integration_resolution", default = 0.5)

  # Run_Integration_verbose
  arg_Run_Integration_verbose <- check_boolean_string(ri$Run_Integration_verbose, arg_name = "Run_Integration_verbose", default = TRUE)

  # --- Markers_Annotation ---
  ma <- config$Markers_Annotation
  if (is.null(ma)) stop("Error: Markers_Annotation block is missing in config.yaml")

  # Markers_Annotation_output_path
  arg_Markers_Annotation_output_path <- check_string_path(ma$Markers_Annotation_output_path, arg_name = "Markers_Annotation_output_path")

  # Markers_Annotation_use_assay
  arg_Markers_Annotation_use_assay <- check_assay(ma$Markers_Annotation_use_assay, arg_name = "Markers_Annotation_use_assay",default = "RNA")

  # Markers_Annotation_group_by
  arg_Markers_Annotation_group_by <- check_group_by(ma$Markers_Annotation_group_by, arg_name = "Markers_Annotation_group_by", default = "seurat_clusters")

  # Markers_Annotation_logfc_threshold
  arg_Markers_Annotation_logfc_threshold <- check_numeric_string(ma$Markers_Annotation_logfc_threshold, arg_name = "Markers_Annotation_logfc_threshold", default = 0.1)

  # Markers_Annotation_min_pct
  arg_Markers_Annotation_min_pct <- check_numeric_string(ma$Markers_Annotation_min_pct, arg_name = "Markers_Annotation_min_pct", default = 0.01)

  # Markers_Annotation_test_use
  valid_test_use <- c("wilcox","wilcox_limma","bimod","roc","t","negbinom","poisson","LR","MAST","DESeq2")
  if (!(ma$Markers_Annotation_test_use %in% valid_test_use)) {
    warning(paste0("Warning: Markers_Annotation_test_use must be one of: ", paste(valid_test_use, collapse = ", ")))
    message(paste0("Using default value: wilcox for Markers_Annotation_test_use"))
    arg_Markers_Annotation_test_use <- "wilcox"
  } else {
    arg_Markers_Annotation_test_use <- ma$Markers_Annotation_test_use
  }

  # Markers_Annotation_only_pos
  arg_Markers_Annotation_only_pos <- check_boolean_string(ma$Markers_Annotation_only_pos, arg_name = "Markers_Annotation_only_pos", default = TRUE)

  # Markers_Annotation_verbose
  arg_Markers_Annotation_verbose <- check_boolean_string(ma$Markers_Annotation_verbose, arg_name = "Markers_Annotation_verbose", default = TRUE)

  # Markers_Annotation_top_n
  arg_Markers_Annotation_top_n <- check_integer_string(ma$Markers_Annotation_top_n, arg_name = "Markers_Annotation_top_n", default = 150)

  # Markers_Annotation_tissue_type
  valid_tissue_type <- c("adipose tissue", "bladder", "blood", "bone", "bone marrow", "brain", "breast", "embryo", "eye", "gastrointestinal tract", "heart",
                         "kidney", "liver", "lung", "mammary gland", "muscle", "other", "ovary", "pancreas", "placenta", "prostate", "skin", "spleen", "stomach",
                         "testis", "thymus", "tooth", "uterus","neuroblastoma", "breast_cancer")
  if (!(ma$Markers_Annotation_tissue_type %in% valid_tissue_type)) {
    stop(paste0("Error: Markers_Annotation_tissue_type must be one of: ", paste(valid_tissue_type, collapse = ", ")))
  } else {
    arg_Markers_Annotation_tissue_type <- ma$Markers_Annotation_tissue_type
  }

  # Markers_Annotation_label_column
  arg_Markers_Annotation_label_column <- check_group_by(ma$Markers_Annotation_label_column, arg_name = "Markers_Annotation_label_column", default = "Celltype")

  # --- DR_Plot ---
  drp <- config$DR_Plot
  if (is.null(drp)) {
    stop("Error: DR_Plot block is missing in config.yaml")
  }

  # DR_Plot_output_path
  arg_DR_Plot_output_path <- check_string_path(drp$DR_Plot_output_path, arg_name = "DR_Plot_output_path")

  # DR_Plot_group_by
  arg_DR_Plot_group_by <- check_group_by(drp$DR_Plot_group_by, arg_name = "DR_Plot_group_by", default = "Celltype")

  # DR_Plot_pt_size
  arg_DR_Plot_pt_size <- check_numeric_string(drp$DR_Plot_pt_size, arg_name = "DR_Plot_pt_size", default = 0.5)

  # DR_Plot_label
  arg_DR_Plot_label <- check_boolean_string(drp$DR_Plot_label, arg_name = "DR_Plot_label", default = TRUE)

  # DR_Plot_label_size
  arg_DR_Plot_label_size <- check_numeric_string(drp$DR_Plot_label_size, arg_name = "DR_Plot_label_size", default = 3)

  # DR_Plot_label_color
  arg_DR_Plot_label_color <- check_string(drp$DR_Plot_label_color, arg_name = "DR_Plot_label_color")

  # DR_Plot_label_box
  arg_DR_Plot_label_box <- check_boolean_string(drp$DR_Plot_label_box, arg_name = "DR_Plot_label_box", default = FALSE)

  # DR_Plot_alpha
  arg_DR_Plot_alpha <- check_numeric_string(drp$DR_Plot_alpha, arg_name = "DR_Plot_alpha", default = 1)

  # DR_Plot_shuffle
  arg_DR_Plot_shuffle <- check_boolean_string(drp$DR_Plot_shuffle, arg_name = "DR_Plot_shuffle", default = FALSE)

  # DR_Plot_raster
  arg_DR_Plot_raster <- check_boolean_string_optional(drp$DR_Plot_raster, arg_name = "DR_Plot_raster", default = FALSE)

  # DR_Plot_width
  arg_DR_Plot_width <- check_numeric_string(drp$DR_Plot_width, arg_name = "DR_Plot_width", default = 6)

  # DR_Plot_height
  arg_DR_Plot_height <- check_numeric_string(drp$DR_Plot_height, arg_name = "DR_Plot_height", default = 6)

  # --- Run_CellChat ---
  rc <- config$Run_CellChat
  if (is.null(rc)) {
    stop("Error: Run_CellChat block is missing in config.yaml")
  }

  # Run_CellChat_output_path
  arg_Run_CellChat_output_path <- check_string_path(rc$Run_CellChat_output_path, arg_name = "Run_CellChat_output_path")

  # Run_CellChat_group_by
  arg_Run_CellChat_group_by <- check_group_by(rc$Run_CellChat_group_by, arg_name = "Run_CellChat_group_by", default = "Celltype")

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - adipose_tissue
  valid_adipose_tissue_celltype = c("Adipocyte","Adipose-derived stem cell","B cell","Basophil","Brown fat cell","Dendritic cell","Endothelial cell","Hematopoietic cell",
  "Luminal epithelial cell","Macrophage","Mammary epithelial cell","Mast cell","Monocyte","Natural killer cell","Neuron","Pericyte","Platelet","T cell","T memory cell")
  if (ma$Markers_Annotation_tissue_type == "adipose tissue") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_adipose_tissue_celltype)) {
      stop(paste0("Error: For adipose tissue, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_adipose_tissue_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "adipose tissue") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_adipose_tissue_celltype)) {
        stop(paste0("Error: For adipose tissue, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_adipose_tissue_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}
  
  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - bladder
  valid_bladder_celltype = c("Basal epithelial cell","Dendritic cell","Endothelial cell","Macrophage","Smooth muscle cell","Smooth muscle progenitor cell","Umbrella cell",
  "Ureteric epithelium cell","Urine-derived stem cell","Urothelial cell")
  if (ma$Markers_Annotation_tissue_type == "bladder") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_bladder_celltype)) {
      stop(paste0("Error: For bladder, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_bladder_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "bladder") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_bladder_celltype)) {
        stop(paste0("Error: For bladder, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_bladder_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - blood
  valid_blood_celltype = c("B cell","Basophil","CD14 Monocyte","CD16 Monocyte","CD4 Central Memory T cell","CD4 Cytotoxic T cell","CD4 Effector Memory T cell","CD4 Naive T cell",
  "CD4 Proliferating T cell","CD4 T cell","CD56-bright natural killer cell","CD56-dim natural killer cell","CD8 Central Memory T cell","CD8 Effector Memory T cell",
  "CD8 Naive T cell","CD8 Proliferating T cell","CD8 T cell","Dendritic cell","Double-negative T cell","Embryonic stem cell","Endothelial cell","Endothelial progenitor cell",
  "Endothelial stem cell","Eosinophil","Epithelial cell","Erythroblast","Erythroid cell","Erythroid precursor cell","Gamma delta T cell","Granulocyte","Hematopoietic stem cell",
  "Immature myeloid cell","Innate lymphoid cell","Intermediate B cell","Leukocyte","Lymphoid cell","Megakaryocyte","Memory B cell","Mesenchymal stem cell","Monocyte",
  "Mucosal-associated invariant T cell","Myeloid cell","Myeloid dendritic cell","Myeloid-derived suppressor cell","Naive B cell","Natural killer cell","Neutrophil",
  "Plasma cell","Plasmablast","Plasmacytoid dendritic cell","Platelet","Regulatory T cell","Reticulocyte","T cell","Thymic emigrant cell")
  if (ma$Markers_Annotation_tissue_type == "blood") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_blood_celltype)) {
      stop(paste0("Error: For blood, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_blood_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "blood") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_blood_celltype)) {
        stop(paste0("Error: For blood, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_blood_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}
  
  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - bone
  valid_bone_celltype = c("Chondrocyte","Dendritic cell","Eosinophil","Erythroid cell","Granulocyte","Hematopoietic stem cell","Meniscus-derived progenitor cell",
  "Mesenchymal progenitor cell","Mesenchymal stem cell","Monocyte","Myeloid progenitor cell","Neutrophil","Osteoblast","Osteoclast","Osteoclast precursor cell",
  "Osteocyte","Periosteum-derived progenitor cell")
  if (ma$Markers_Annotation_tissue_type == "bone") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_bone_celltype)) {
      stop(paste0("Error: For bone, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_bone_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "bone") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_bone_celltype)) {
        stop(paste0("Error: For bone, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_bone_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - bone marrow
  valid_bone_marrow_celltype = c("B cell","Basophil","Blastema cell","Common lymphoid progenitor cell","Common myeloid progenitor","Dendritic cell","Dendritic progenitor cell",
  "Endothelial cell","Endothelial progenitor cell","Eosinophil","Erythroid cell","Erythroid megakaryocyte progenitor cell","Erythroid precursor cell",
  "Granulocyte monocyte progenitor cell","Hematopoietic stem cell","Immature dendritic cell","Innate lymphoid cell","Langerhans cell","Lymphoid-primed multipotent progenitor cell",
  "Macrophage","Mast cell","Mature T cell","Megakaryocyte","Megakaryocyte progenitor cell","Mesenchymal stem cell","Monocyte","Monocyte derived dendritic cell",
  "Mucosal-associated invariant T cell","Myeloid cell","Myeloid dendritic cell","Myeloid progenitor cell","Natural killer cell","Neutrophil","Osteoclast precursor cell",
  "Plasma cell","Plasmacytoid dendritic cell","Platelet","Precursor plasmacytoid dendritic cell","T cell")
  if (ma$Markers_Annotation_tissue_type == "bone marrow") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_bone_marrow_celltype)) {
      stop(paste0("Error: For bone marrow, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_bone_marrow_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "bone marrow") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_bone_marrow_celltype)) {
        stop(paste0("Error: For bone marrow, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_bone_marrow_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - brain
  valid_brain_celltype = c("Adrenergic neurons","Anterior pituitary gland cell","Astrocyte","B cell","Basket cell","Bergmann glial cell","Cajal-Retzius cell",
  "CCK basket cell","Chandelier cell","Cholinergic neuron","Choroid plexus cell","Dopaminergic neuron","Endothelial cell","Ependymal cell","GABAergic neuron",
  "Ganglion cell","Glutamatergic neuron","Hypothalamic ependymal cell","Immature neuron","Interneuron","Interneuron-selective cell","Lepotomeningeal cell",
  "Macrophage","Martinotti cell","Meningeal cell","Microglia cell","Motor neuron","Mural cell","Neural progenitor cell","Neural stem cell","Neuroblast",
  "Neuroendocrine cell","Neuron","Noradrenergic neuron","Olfactory ensheathing glia","Olfactory sensory neuron","Oligodendrocyte","Oligodendrocyte precursor cell",
  "Pan-gabaergic","Pericyte","Pinealocyte","Purkinje neuron","Pyramidal cell","Radial glia cell","Rhombic lip cell","Satellite glial cell","Schwann cell",
  "Serotonergic neuron","Smooth muscle cell","T cell","Tanycyte","Trigeminal neuron","Type IA spiral ganglion neuron","Type IB spiral ganglion neuron",
  "Type IC spiral ganglion neuron","Type II spiral ganglion neuron")
  if (ma$Markers_Annotation_tissue_type == "brain") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_brain_celltype)) {
      stop(paste0("Error: For brain, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_brain_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "brain") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_brain_celltype)) {
        stop(paste0("Error: For brain, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_brain_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - breast
  valid_breast_celltype = c("B cell","Basal epithelial cell","Dendritic cell","Endothelial cell","Eosinophil","Epithelial cell","Fibroblast","Hematopoietic stem cell",
  "Luminal epithelial cell","Luminal progenitor cell","Macrophage","Mast cell","Mesenchymal stem cell","Myoepithelial cell","Natural killer cell","Neutrophil",
  "Pericyte","T cell")
  if (ma$Markers_Annotation_tissue_type == "breast") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_breast_celltype)) {
      stop(paste0("Error: For breast, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_breast_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "breast") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_breast_celltype)) {
        stop(paste0("Error: For breast, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_breast_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - embryo
  valid_embryo_celltype = c("Arterial cell","Blastomere","Cardiomyocyte","Dorsal otocyst","Ectoderm cell","Embryonic stem cell","Endocardial cell","Endocrine cell",
  "Endoderm cell","Endothelial cell","Epiblast cell","Germ cell","Hemangioblast","Hematopoietic stem cell","Macrophage","Mesenchymal stem cell","Mesoderm cell",
  "Myeloblast","Neural crest cell","Neural stem cell","Neural tube cell","Neuron","Neutrophil","Oocyte","Primitive endoderm cell","Skeletal muscle cell",
  "Trophectoderm cell","Trophoblast cell","Unrestricted somatic stem cell","Venous cell","Ventral otocyst")
  if (ma$Markers_Annotation_tissue_type == "embryo") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_embryo_celltype)) {
      stop(paste0("Error: For embryo, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_embryo_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "embryo") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_embryo_celltype)) {
        stop(paste0("Error: For embryo, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_embryo_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - eye
  valid_eye_celltype = c("Bipolar cell","Endothelial cell","Epithelial cell","Epithelial stem cell","Erythroid cell","Ganglion cell","Hematopoietic stem cell",
  "Lymphocyte","Macrophage","Mesenchymal cell","Mesenchymal stem cell","Muller cell","Myoepithelial cell","Photoreceptor cell","Progenitor cell")
  if (ma$Markers_Annotation_tissue_type == "eye") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_eye_celltype)) {
      stop(paste0("Error: For eye, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_eye_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "eye") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_eye_celltype)) {
        stop(paste0("Error: For eye, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_eye_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - gastrointestinal tract
  valid_gastrointestinal_tract_celltype = c("Ciliated epithelial cell","Colonic stem cell","Crypt cell","Dendritic cell","Enteric glia cell","Enteric neuron",
  "Enterochromaffin cell","Enterocyte","Enterocyte progenitor cell","Enteroendocrine cell","Enteroendocrine precursor cell","Foveolar cell","Gastric chief cell",
  "Goblet cell","Goblet progenitor cell","Intestinal stem cell","Macrophage","Mast cell","Microfold cell","Paneth cell","Parietal cell","S cell","Secretory progenitor cell",
  "T cell","Tuft cell","Tuft progenitor cell")
  if (ma$Markers_Annotation_tissue_type == "gastrointestinal tract") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_gastrointestinal_tract_celltype)) {
      stop(paste0("Error: For gastrointestinal tract, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_gastrointestinal_tract_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "gastrointestinal tract") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_gastrointestinal_tract_celltype)) {
        stop(paste0("Error: For gastrointestinal tract, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_gastrointestinal_tract_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - heart
  valid_heart_celltype = c("Aorta valve interstitial cell","Cardiocyte","Cardiomyocyte","Cardiovascular progenitor cell","Endocardial cell","Endothelial cell",
  "Erythroid cell","Fibroblast","Macrophage","Myofibroblast","Pericyte","Progenitor cell","Smooth muscle cell")
  if (ma$Markers_Annotation_tissue_type == "heart") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_heart_celltype)) {
      stop(paste0("Error: For heart, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_heart_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "heart") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_heart_celltype)) {
        stop(paste0("Error: For heart, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_heart_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - kidney
  valid_kidney_celltype = c("B cell","Basophil","Connecting tubule epithelial cell","Cortical cell","Dendritic cell","Descending vasa recta endothelial cell",
  "Distal convoluted tubule cell","Endothelial cell","Epithelial cell","Erythroblast","Erythroid precursor cell","Fibroblast","Glomerular capillary endothelial cell",
  "Granular cell","Infiltrated mononuclear cell","Inner medulla collecting duct epithelial cell","Intercalated cell","Juxtaglomerular cell","Kidney progenitor cell",
  "Loop of Henle cell","Loop of Henle cortical thick ascending limb epithelial cell","Loop of Henle medullary thick ascending limb epithelial cell",
  "Loop of Henle thin descending limb epithelial cell","Lymphatic endothelial cell","M2 Macrophage","Macrophage","Macula densa","Mast cell","Medullary fibroblast",
  "Mesangial cell","Monocyte","Mononuclear phagocyte","Natural killer cell","Nephron epithelial cell","Neutrophil","Papillary tip epithelial cell","Parietal epithelial cell",
  "Pericyte","Peritubular capilary endothelial cell","Plasma cell","Podocyte","Principal cell","Proximal tubule brush border cell","Proximal tubule cell",
  "Proximal tubule epithelial cell","Renal alpha-intercalated cell","Schwann cell","T cell")
  if (ma$Markers_Annotation_tissue_type == "kidney") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_kidney_celltype)) {
      stop(paste0("Error: For kidney, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_kidney_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "kidney") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_kidney_celltype)) {
        stop(paste0("Error: For kidney, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_kidney_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - liver
  valid_liver_celltype = c("B cell","Cholangiocyte","Dendritic cell","Endothelial cell","Epithelial cell","Hematopoietic cell","Hepatoblast","Hepatocyte","Kupffer cell",
  "Liver bud hepatic cell","Mesenchymal cell","Mesenchymal stem cell","Monocyte","Myofibroblast","Natural killer cell","Neutrophil","Progenitor cell","Stellate cell",
  "Stem cell","T cell")
  if (ma$Markers_Annotation_tissue_type == "liver") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_liver_celltype)) {
      stop(paste0("Error: For liver, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_liver_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "liver") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_liver_celltype)) {
        stop(paste0("Error: For liver, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_liver_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - lung
  valid_lung_celltype = c("Adventitial fibroblast","B cell","Basal cell","Basophil","Bronchial epithelial cell","Brush Cell","Capillary aerocyte","Ciliated cell",
  "Club cell","Dendritic cell","Endothelial cell","Eosinophil","Epithelial cell","Fibroblast","Fibromyocyte","Goblet cell","Granulocyte","Lipofibroblast","Lymphatic endothelial cell",
  "Lymphoid cell","M1 Macrophage","M2 Macrophage","Macrophage","Mast cell","Mesenchymal progenitor cell","Mesothelial cell","Monocyte","Mucous cell","Myeloid cell",
  "Myeloid leukocyte","Myofibroblast","Neuroendocrine cell","Neutrophil","Pericyte","Plasma cell","Pulmonary alveolar type I cell","Pulmonary alveolar type II cell",
  "Pulmonary Ionocyte cell","Secretory cell","Serous cell","Smooth muscle cell","Stem cell","T cell")
  if (ma$Markers_Annotation_tissue_type == "lung") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_lung_celltype)) {
      stop(paste0("Error: For lung, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_lung_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "lung") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_lung_celltype)) {
        stop(paste0("Error: For lung, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_lung_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - mammary gland
  valid_mammary_gland_celltype = c("B cell","Basal cell","Endothelial cell","Epithelial cell","Hormone sensing differentiated cell","Hormone sensing progenitor cell",
  "Luminal cell","Luminal epithelial cell","Luminal progenitor cell","Macrophage","Mast cell","Muscle cell","Myoepithelial cell","Progenitor cell","Stem cell","T cell")
  if (ma$Markers_Annotation_tissue_type == "mammary gland") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_mammary_gland_celltype)) {
      stop(paste0("Error: For mammary gland, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_mammary_gland_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "mammary gland") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_mammary_gland_celltype)) {
        stop(paste0("Error: For mammary gland, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_mammary_gland_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - muscle
  valid_muscle_celltype = c("Adventitial cell","B cell","Endothelial cell","Erythroblast","Fibroblast","Granulocyte monocyte progenitor cell","Inflammatory cell",
  "Macrophage","Mesenchymal progenitor cell","Mesenchymal stem cell","Muscle-derived cell","Myoblast","Myocyte","Myocyte progenitor cell","Myoepithelial cell",
  "Myofibroblast","Myogenic endothelial cell","Neutrophil","Pericyte","Progenitor cell","Satellite cell","Schwann cell","Smooth muscle cell","T cell","Tenocyte")
  if (ma$Markers_Annotation_tissue_type == "muscle") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_muscle_celltype)) {
      stop(paste0("Error: For muscle, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_muscle_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "muscle") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_muscle_celltype)) {
        stop(paste0("Error: For muscle, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_muscle_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - other
  valid_other_celltype = c("Chromaffin cell","Epithelial cell","Follicular cell","Glomus cell","Hair cell","Invasive spongiotrophoblast","Langerhans cell","Mesothelial cell",
  "Nucleus pulposus cell","Parathyroid chief cell","Salivary mucous cell","Thymocyte","Type III taste bud cell","Vomeronasal sensory neuron")
  if (ma$Markers_Annotation_tissue_type == "other") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_other_celltype)) {
      stop(paste0("Error: For other, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_other_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "other") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_other_celltype)) {
        stop(paste0("Error: For other, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_other_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - ovary
  valid_ovary_celltype = c("Cumulus cell","Endothelial cell","Endothelium cell","Epithelial cell","Germ cell","Granulosa cell","Luteal cell","Macrophage","Mesenchymal cell",
  "Oocyte","Pluripotent stem cell","Stem cell","Theca interna cell","Thecal cell")
  if (ma$Markers_Annotation_tissue_type == "ovary") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_ovary_celltype)) {
      stop(paste0("Error: For ovary, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_ovary_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "ovary") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_ovary_celltype)) {
        stop(paste0("Error: For ovary, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_ovary_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - pancreas
  valid_pancreas_celltype = c("Acinar cell","Alpha cell","B cell","Beta cell","Delta cell","Dendritic cell","Ductal cell","Ductal stem cell","Endocrine cell",
  "Endothelial cell","Epithelial cell","Epsilon cell","Gamma cell","Glial cell","Granulocyte","Macrophage","Mast cell","Pancreatic progenitor cell","Polypeptide cell",
  "Schwann cell","Stellate cell","T cell")
  if (ma$Markers_Annotation_tissue_type == "pancreas") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_pancreas_celltype)) {
      stop(paste0("Error: For pancreas, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_pancreas_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "pancreas") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_pancreas_celltype)) {
        stop(paste0("Error: For pancreas, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_pancreas_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - placenta
  valid_placenta_celltype = c("Basophil","Decidual stem cell","Dendritic cell","Endodermal cell","Endothelial cell","Erythroid cell","Hofbauer cell","Invasive spongiotrophoblast",
  "Labyrinthine trophoblast","Macrophage","Megakaryocyte progenitor cell","Mesenchymal stem cell","Monocyte","Natural killer cell","Pericyte","Spongiotrophoblast",
  "Stem cell","Trophoblast cell","Basal cell","Epithelial cell","Luminal cell","Progenitor cell","Prostate epithelial cell","Prostate stem cell","Stem cell")
  if (ma$Markers_Annotation_tissue_type == "placenta") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_placenta_celltype)) {
      stop(paste0("Error: For placenta, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_placenta_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "placenta") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_placenta_celltype)) {
        stop(paste0("Error: For placenta, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_placenta_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - skin
  valid_skin_celltype = c("B cell","Basal cell","Dendritic cell","Endothelial cell","Epidermal stem cell","Epithelial cell","Keratinocyte","Langerhans cell",
  "Macrophage","Melanocyte","Merkel cell","Mesenchymal stem cell","Neural crest stem cell","Sebocyte","Stem cell","T cell","Trichocyte")
  if (ma$Markers_Annotation_tissue_type == "skin") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_skin_celltype)) {
      stop(paste0("Error: For skin, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_skin_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "skin") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_skin_celltype)) {
        stop(paste0("Error: For skin, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_skin_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - spleen
  valid_spleen_celltype = c("B cell","Dendritic cell","Erythroblast","Granulocyte","Lymphocyte","Macrophage","Monocyte","Natural killer cell","Neutrophil",
  "Plasma cell","Plasmacytoid dendritic cell","T cell")
  if (ma$Markers_Annotation_tissue_type == "spleen") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_spleen_celltype)) {
      stop(paste0("Error: For spleen, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_spleen_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "spleen") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_spleen_celltype)) {
        stop(paste0("Error: For spleen, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_spleen_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - stomach
  valid_stomach_celltype = c("Dendritic cell","Gastric stem cell","Multipotent stem cell","Parietal cell","Parietal progenitor cell","Pit cell","Pit progenitor cell",
  "Smooth muscle cell","Tuft cell")
  if (ma$Markers_Annotation_tissue_type == "stomach") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_stomach_celltype)) {
      stop(paste0("Error: For stomach, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_stomach_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "stomach") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_stomach_celltype)) {
        stop(paste0("Error: For stomach, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_stomach_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - testis
  valid_testis_celltype = c("Leydig cell","Macrophage","Peritubular myoid cell","Pre sertoli cell","Sertoli cell","Spermatids","Spermatocyte","Spermatogonia",
  "Spermatogonial stem cell","Spermatogonium","Spermatozoa")
  if (ma$Markers_Annotation_tissue_type == "testis") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_testis_celltype)) {
      stop(paste0("Error: For testis, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_testis_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "testis") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_testis_celltype)) {
        stop(paste0("Error: For testis, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_testis_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - thymus
  valid_thymus_celltype = c("Cortical thymic epithelial cell","Dendritic cell","Epithelial cell","Medullary thymic epithelial cell","T cell","Thymocyte")
  if (ma$Markers_Annotation_tissue_type == "thymus") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_thymus_celltype)) {
      stop(paste0("Error: For thymus, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_thymus_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "thymus") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_thymus_celltype)) {
        stop(paste0("Error: For thymus, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_thymus_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - tooth
  valid_tooth_celltype = c("Alveolar osteocyte","Dental follicle cell","Dental pulp cell","Dental pulp stem cell","Endothelial cell","Epithelial cell","Glial cell",
  "Macrophage","Mesenchymal stem cell","Natural killer cell","Periodontal ligament stem cell","Perivascular cell","T cell")
  if (ma$Markers_Annotation_tissue_type == "tooth") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_tooth_celltype)) {
      stop(paste0("Error: For tooth, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_tooth_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "tooth") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_tooth_celltype)) {
        stop(paste0("Error: For tooth, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_tooth_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - uterus
  valid_uterus_celltype = c("B cell","Decidual cell","Endometrial stem cell","Keratinocyte","Macrophage","Monocyte","Natural killer cell","Neutrophil","Smooth muscle cell",
  "Stem cell")
  if (ma$Markers_Annotation_tissue_type == "uterus") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_uterus_celltype)) {
      stop(paste0("Error: For uterus, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_uterus_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "uterus") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_uterus_celltype)) {
        stop(paste0("Error: For uterus, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_uterus_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - neuroblastoma
  valid_neuroblastoma_celltype = c("B cell","Endothelial cell","Fibroblast","Neuroendocrine tumor cell","NK cell","Plasmacytoid dendritic cell","Conventional dendritic cell",
  "Monocyte","Macrophage","Plasma cell","Red blood cell","Schwann cell","Mesenchyme","T cell","Other stromal cell")
  if (ma$Markers_Annotation_tissue_type == "neuroblastoma") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_neuroblastoma_celltype)) {
      stop(paste0("Error: For neuroblastoma, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_neuroblastoma_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "neuroblastoma") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_neuroblastoma_celltype)) {
        stop(paste0("Error: For muscle, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_neuroblastoma_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}

  # Run_CellChat_source_celltype & Run_CellChat_target_celltype - breast_cancer
  valid_breast_cancer_celltype = c("Malignant epithelial cell","Basal-like tumor cell","Luminal tumor cell","HER2-enriched tumor cell","Endothelial cell","Fibroblast",
  "Myofibroblast","Pericyte","Adipocyte","T cell","B cell","Plasma cell","Macrophage","Monocyte","Dendritic cell","NK cell","Mast cell","Red blood cell","Other stromal cell")
  if (ma$Markers_Annotation_tissue_type == "breast_cancer") {
  if (!is.null(rc$Run_CellChat_source_celltype)) {
    if (!(rc$Run_CellChat_source_celltype %in% valid_breast_cancer_celltype)) {
      stop(paste0("Error: For breast cancer, Run_CellChat_source_celltype must be NULL or one of: ", paste(valid_breast_cancer_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_source_celltype <- rc$Run_CellChat_source_celltype
    }
  } else {
    arg_Run_CellChat_source_celltype <- NULL
  }
}
  if (ma$Markers_Annotation_tissue_type == "breast_cancer") {
  if (!is.null(rc$Run_CellChat_target_celltype)) {
    if (!(rc$Run_CellChat_target_celltype %in% valid_breast_cancer_celltype)) {
        stop(paste0("Error: For breast cancer, Run_CellChat_target_celltype must be NULL or one of: ", paste(valid_breast_cancer_celltype, collapse = ", ")))
    } else {
      arg_Run_CellChat_target_celltype <- rc$Run_CellChat_target_celltype
    }
  } else {
    arg_Run_CellChat_target_celltype <- NULL
  }
}
  # Run_CellChat_plot_heatmap
  arg_Run_CellChat_plot_heatmap <- check_boolean_string(rc$Run_CellChat_plot_heatmap, arg_name = "Run_CellChat_plot_heatmap", default = TRUE)

  # Run_CellChat_ntop_signaling
  arg_Run_CellChat_ntop_signaling <- check_integer_string(rc$Run_CellChat_ntop_signaling, arg_name = "Run_CellChat_ntop_signaling", default = 5)

  # Run_CellChat_MaxGroup
  arg_Run_CellChat_MaxGroup <- check_string(rc$Run_CellChat_MaxGroup, arg_name = "Run_CellChat_MaxGroup")

  # --- Vina_Docking ---
  vd <- config$Vina_Docking
  if (is.null(vd)) {
    stop("Error: Vina_Docking block is missing in config.yaml")
  }

  # Vina_Docking_ligand_ref_file
  arg_Vina_Docking_ligand_ref_file <- check_csv_path(vd$Vina_Docking_ligand_ref_file, arg_name = "Vina_Docking_ligand_ref_file")

  # Vina_Docking_receptor_ref_file
  arg_Vina_Docking_receptor_ref_file <- check_csv_path(vd$Vina_Docking_receptor_ref_file, arg_name = "Vina_Docking_receptor_ref_file")

  # Vina_Docking_output_path
  arg_Vina_Docking_output_path <- check_string_path(vd$Vina_Docking_output_path, arg_name = "Vina_Docking_output_path")

  # Vina_Docking_cas_txt_file (Optional)
  arg_Vina_Docking_cas_txt_file <- check_txt_path(vd$Vina_Docking_cas_txt_file, arg_name = "Vina_Docking_cas_txt_file")

  # Vina_Docking_docking_ligand_dir (Optional)
  arg_Vina_Docking_docking_ligand_dir <- check_string_path_optional(vd$Vina_Docking_docking_ligand_dir, arg_name = "Vina_Docking_docking_ligand_dir")

  # Vina_Docking_use_fda (Optional)
  arg_Vina_Docking_use_fda <- check_boolean_string_optional(vd$Vina_Docking_use_fda, arg_name = "Vina_Docking_use_fda", default = FALSE)

  # Vina_Docking_fda_txt (Optional)
  arg_Vina_Docking_fda_txt <- check_txt_path(vd$Vina_Docking_fda_txt, arg_name = "Vina_Docking_fda_txt")

  # Vina_Docking_docking_receptor_dir
  arg_Vina_Docking_docking_receptor_dir <- check_string_path_optional(vd$Vina_Docking_docking_receptor_dir, arg_name = "Vina_Docking_docking_receptor_dir")

  # Vina_Docking_vina_exhaustiveness
  arg_Vina_Docking_vina_exhaustiveness <- check_integer_string(vd$Vina_Docking_vina_exhaustiveness, arg_name = "Vina_Docking_vina_exhaustiveness", default = 8)

  # Vina_Docking_vina_num_modes
  arg_Vina_Docking_vina_num_modes <- check_integer_string(vd$Vina_Docking_vina_num_modes, arg_name = "Vina_Docking_vina_num_modes", default = 9)

  # Vina_Docking_vina_seed
  arg_Vina_Docking_vina_seed <- check_integer_string(vd$Vina_Docking_vina_seed, arg_name = "Vina_Docking_vina_seed", default = 42)

  # Vina_Docking_vina_cpu
  arg_Vina_Docking_vina_cpu <- check_integer_string(vd$Vina_Docking_vina_cpu, arg_name = "Vina_Docking_vina_cpu", default = 1)

  return(list(
  Seurat_output_path = arg_Seurat_output_path,
  Load_QC_input_files = arg_Load_QC_input_files,
  Load_QC_input_type= arg_Load_QC_input_type,
  Load_QC_species = arg_Load_QC_species,
  Load_QC_min_features= arg_Load_QC_min_features,
  Load_QC_min_cells = arg_Load_QC_min_cells,
  Load_QC_names_delim = arg_Load_QC_names_delim,
  Load_QC_max_mito = arg_Load_QC_max_mito,
  Load_QC_verbose = arg_Load_QC_verbose,
  Load_QC_metadata_file = arg_Load_QC_metadata_file,
  Normalization_Scale_use_assay = arg_Normalization_Scale_use_assay,
  Normalization_Scale_normalization_method = arg_Normalization_Scale_normalization_method,
  Normalization_Scale_scale_factor = arg_Normalization_Scale_scale_factor,
  Normalization_Scale_CLR_margin = arg_Normalization_Scale_CLR_margin,
  Normalization_Scale_hvg_method = arg_Normalization_Scale_hvg_method,
  Normalization_Scale_nVariableFeatures = arg_Normalization_Scale_nVariableFeatures,
  Normalization_Scale_split_by = arg_Normalization_Scale_split_by,
  Normalization_Scale_model_use = arg_Normalization_Scale_model_use,
  Normalization_Scale_scale_max = arg_Normalization_Scale_scale_max,
  Normalization_Scale_scale_features = arg_Normalization_Scale_scale_features,
  Normalization_Scale_verbose = arg_Normalization_Scale_verbose,
  DR_Cluster_pca_features = arg_DR_Cluster_pca_features,
  DR_Cluster_seed = arg_DR_Cluster_seed,
  DR_Cluster_dims = arg_DR_Cluster_dims,
  DR_Cluster_k_param = arg_DR_Cluster_k_param,
  DR_Cluster_n_trees = arg_DR_Cluster_n_trees,
  DR_Cluster_resolution = arg_DR_Cluster_resolution,
  DR_Cluster_reduction_method = arg_DR_Cluster_reduction_method,
  DR_Cluster_reduction_assay = arg_DR_Cluster_reduction_assay,
  DR_Cluster_clustering_algorithm = arg_DR_Cluster_clustering_algorithm,
  DR_Cluster_verbose = arg_DR_Cluster_verbose,
  Run_Integration_run_integration = arg_Run_Integration_run_integration,
  Run_Integration_group_by = arg_Run_Integration_group_by,
  Run_Integration_method = arg_Run_Integration_method,
  Run_Integration_nfeatures = arg_Run_Integration_nfeatures,
  Run_Integration_normalization_method = arg_Run_Integration_normalization_method,
  Run_Integration_dims = arg_Run_Integration_dims,
  Run_Integration_resolution = arg_Run_Integration_resolution,
  Run_Integration_verbose = arg_Run_Integration_verbose,
  Markers_Annotation_output_path = arg_Markers_Annotation_output_path,
  Markers_Annotation_use_assay = arg_Markers_Annotation_use_assay,
  Markers_Annotation_group_by = arg_Markers_Annotation_group_by,
  Markers_Annotation_logfc_threshold = arg_Markers_Annotation_logfc_threshold,
  Markers_Annotation_min_pct = arg_Markers_Annotation_min_pct,
  Markers_Annotation_test_use = arg_Markers_Annotation_test_use,
  Markers_Annotation_only_pos = arg_Markers_Annotation_only_pos,
  Markers_Annotation_verbose = arg_Markers_Annotation_verbose,
  Markers_Annotation_top_n = arg_Markers_Annotation_top_n,
  Markers_Annotation_tissue_type = arg_Markers_Annotation_tissue_type,
  Markers_Annotation_label_column = arg_Markers_Annotation_label_column,
  DR_Plot_output_path = arg_DR_Plot_output_path,
  DR_Plot_group_by = arg_DR_Plot_group_by,
  DR_Plot_pt_size = arg_DR_Plot_pt_size,
  DR_Plot_label = arg_DR_Plot_label,
  DR_Plot_label_size = arg_DR_Plot_label_size,
  DR_Plot_label_color = arg_DR_Plot_label_color,
  DR_Plot_label_box = arg_DR_Plot_label_box,
  DR_Plot_alpha = arg_DR_Plot_alpha,
  DR_Plot_shuffle = arg_DR_Plot_shuffle,
  DR_Plot_raster = arg_DR_Plot_raster,
  DR_Plot_width = arg_DR_Plot_width,
  DR_Plot_height = arg_DR_Plot_height,
  Run_CellChat_output_path = arg_Run_CellChat_output_path,
  Run_CellChat_group_by = arg_Run_CellChat_group_by,
  Run_CellChat_source_celltype = arg_Run_CellChat_source_celltype,
  Run_CellChat_target_celltype = arg_Run_CellChat_target_celltype,
  Run_CellChat_plot_heatmap = arg_Run_CellChat_plot_heatmap,
  Run_CellChat_ntop_signaling = arg_Run_CellChat_ntop_signaling,
  Run_CellChat_MaxGroup = arg_Run_CellChat_MaxGroup,
  Vina_Docking_output_path = arg_Vina_Docking_output_path,
  Vina_Docking_ligand_ref_file = arg_Vina_Docking_ligand_ref_file,
  Vina_Docking_receptor_ref_file = arg_Vina_Docking_receptor_ref_file,
  Vina_Docking_cas_txt_file = arg_Vina_Docking_cas_txt_file,
  Vina_Docking_use_fda = arg_Vina_Docking_use_fda,
  Vina_Docking_fda_txt = arg_Vina_Docking_fda_txt,
  Vina_Docking_docking_ligand_dir = arg_Vina_Docking_docking_ligand_dir,
  Vina_Docking_docking_receptor_dir = arg_Vina_Docking_docking_receptor_dir,
  Vina_Docking_vina_exhaustiveness = arg_Vina_Docking_vina_exhaustiveness,
  Vina_Docking_vina_num_modes = arg_Vina_Docking_vina_num_modes,
  Vina_Docking_vina_seed = arg_Vina_Docking_vina_seed,
  Vina_Docking_vina_cpu = arg_Vina_Docking_vina_cpu
  ))
}


# -------------------------
# Step 0: Run validation
# -------------------------

args <- validate_config(config)
seurat_output_dir <- args$Seurat_output_path
message("✅ Config validation passed.")

# Load R scripts
source("functions/Load_QC.R")
source("functions/Normalization_Scale.R")
source("functions/DR_Cluster.R")
source("functions/Run_Integration.R")
source("functions/Markers_Annotation.R")
source("functions/DR_Plot.R")
source("functions/Run_CellChat.R")
source("functions/Vina_Docking.R")

# -------------------------
# Step 1: Load and QC
# -------------------------
message("========== Start loading data and quality control ... ==========")
seurat_obj <- Load_QC(
  Load_QC_input_files   = args$Load_QC_input_files,
  Load_QC_input_type    = args$Load_QC_input_type,
  Load_QC_species       = args$Load_QC_species,
  Load_QC_min_features  = args$Load_QC_min_features,
  Load_QC_min_cells     = args$Load_QC_min_cells,
  Load_QC_names_delim   = args$Load_QC_names_delim,
  Load_QC_max_mito      = args$Load_QC_max_mito,
  Load_QC_metadata_file = args$Load_QC_metadata_file,
  Load_QC_verbose       = args$Load_QC_verbose
)

message("Saving the QC-completed file to seurat_qc.rds")
# saveRDS(seurat_obj, file = file.path(seurat_output_dir, "seurat_qc.rds"))
message("========== Loading data and QC process complete. ==========")
message("")

# -------------------------
# Step 2: Normalization and scaling
# -------------------------
message("========== Start normalization and scaling ... ==========")
args$Normalization_Scale_scale_factor <- as.numeric(args$Normalization_Scale_scale_factor)
args$Normalization_Scale_nVariableFeatures <- as.integer(args$Normalization_Scale_nVariableFeatures)

seurat_obj <- Normalization_Scale(
  seurat_obj                                 = seurat_obj,
  Normalization_Scale_use_assay              = args$Normalization_Scale_use_assay,
  Normalization_Scale_normalization_method   = args$Normalization_Scale_normalization_method,
  Normalization_Scale_scale_factor           = args$Normalization_Scale_scale_factor,
  Normalization_Scale_CLR_margin             = args$Normalization_Scale_CLR_margin,
  Normalization_Scale_hvg_method             = args$Normalization_Scale_hvg_method,
  Normalization_Scale_nVariableFeatures      = args$Normalization_Scale_nVariableFeatures,
  Normalization_Scale_split_by               = args$Normalization_Scale_split_by,
  Normalization_Scale_model_use              = args$Normalization_Scale_model_use,
  Normalization_Scale_scale_max              = args$Normalization_Scale_scale_max,
  Normalization_Scale_scale_features         = args$Normalization_Scale_scale_features,
  Normalization_Scale_verbose                = args$Normalization_Scale_verbose
)

message("Saving the normalized and scaled file to seurat_norm_scale.rds")
# saveRDS(seurat_obj, file = file.path(seurat_output_dir, "seurat_norm_scale.rds"))
message("========== Normalization and scaling complete. ==========")
message("")

# -------------------------
# Step 3.1: Run DR_Cluster
# -------------------------
message("========== Start dimensional reduction and clustering ... ==========")
if (!identical(args$DR_Cluster_dims, "auto")) {
  args$DR_Cluster_dims <- seq_len(as.integer(args$DR_Cluster_dims))
}

seurat_obj <- DR_Cluster(
  seurat_obj                        = seurat_obj,
  DR_Cluster_pca_features           = args$DR_Cluster_pca_features,
  DR_Cluster_seed                   = args$DR_Cluster_seed,
  DR_Cluster_dims                   = args$DR_Cluster_dims,
  DR_Cluster_k_param                = args$DR_Cluster_k_param,
  DR_Cluster_n_trees                = args$DR_Cluster_n_trees,
  DR_Cluster_resolution             = args$DR_Cluster_resolution,
  DR_Cluster_reduction_method       = args$DR_Cluster_reduction_method,
  DR_Cluster_reduction_assay        = args$DR_Cluster_reduction_assay,
  DR_Cluster_clustering_algorithm   = args$DR_Cluster_clustering_algorithm,
  DR_Cluster_verbose                = args$DR_Cluster_verbose
)

message("Saving the clustered file to seurat_dr_cluster.rds")
# saveRDS(seurat_obj, file = file.path(seurat_output_dir, "seurat_dr_cluster.rds"))
message("========== Dimensional reduction and clustering complete. ==========")
message("")

# -------------------------
# Step 3.2: Run Integration (Optional)
# -------------------------
if (!is.null(config$Run_Integration) && args$Run_Integration_run_integration == TRUE) {
  message("========== Start Seurat integration ... ==========")

  seurat_obj <- Run_Integration(
    seurat_obj                           = seurat_obj,
    Run_Integration_group_by             = args$Run_Integration_group_by,
    Run_Integration_method               = args$Run_Integration_method,
    Run_Integration_nfeatures            = args$Run_Integration_nfeatures,
    Run_Integration_normalization_method = args$Run_Integration_normalization_method,
    Run_Integration_dims                 = args$Run_Integration_dims,
    Run_Integration_resolution           = args$Run_Integration_resolution,
    Run_Integration_verbose              = args$Run_Integration_verbose
  )

  message("Saving the integrated file to seurat_integrated.rds")
  # saveRDS(seurat_obj, file = file.path(seurat_output_dir, "seurat_integrated.rds"))
  message("========== Integration complete. ==========")
  message("")
}

# -------------------------
# Step 4: Find markers and annotate
# -------------------------
message("========== Start cluster markers calculation and cell annotation... ==========")

seurat_obj <- Markers_Annotation(
  seurat_obj                               = seurat_obj,
  Markers_Annotation_use_assay             = args$Markers_Annotation_use_assay,
  Markers_Annotation_group_by              = args$Markers_Annotation_group_by,
  Markers_Annotation_logfc_threshold       = args$Markers_Annotation_logfc_threshold,
  Markers_Annotation_test_use              = args$Markers_Annotation_test_use,
  Markers_Annotation_min_pct               = args$Markers_Annotation_min_pct,
  Markers_Annotation_only_pos              = args$Markers_Annotation_only_pos,
  Markers_Annotation_verbose               = args$Markers_Annotation_verbose,
  Markers_Annotation_top_n                 = args$Markers_Annotation_top_n,
  Markers_Annotation_output_path           = args$Markers_Annotation_output_path,
  Markers_Annotation_tissue_type           = args$Markers_Annotation_tissue_type,
  Markers_Annotation_label_column          = args$Markers_Annotation_label_column
)

message("Saving the annotated file to seurat_markers_annotated.rds")
# saveRDS(seurat_obj, file = file.path(seurat_output_dir, "seurat_markers_annotated.rds"))
message("========== Find marker and annotation complete. ==========")
message("")

# -------------------------
# Step 5: Check and Dimensional Reduction Plot
# -------------------------
message("========== Start visualization with Dimplot... ==========")

DR_Plot(
    seurat_obj                  = seurat_obj,
    DR_Plot_output_path         = args$DR_Plot_output_path,
    DR_Cluster_reduction_method = args$DR_Cluster_reduction_method,
    DR_Plot_group_by            = args$DR_Plot_group_by,
    DR_Plot_pt_size             = args$DR_Plot_pt_size,
    DR_Plot_label               = args$DR_Plot_label,
    DR_Plot_label_size          = args$DR_Plot_label_size,
    DR_Plot_label_color         = args$DR_Plot_label_color,
    DR_Plot_label_box           = args$DR_Plot_label_box,
    DR_Plot_alpha               = args$DR_Plot_alpha,
    DR_Plot_shuffle             = args$DR_Plot_shuffle,
    DR_Plot_raster              = args$DR_Plot_raster,
    DR_Plot_width               = args$DR_Plot_width,
    DR_Plot_height              = args$DR_Plot_height
    )

message("========== Visualization complete. ==========")
message("")

# -------------------------
# Step 6: CellChat analysis
# -------------------------
message("========== Run CellChat on Seurat object ... ==========")

if (!is.null(config$Run_CellChat)) {
  # seurat_obj <- readRDS(file.path(seurat_output_dir, "seurat_markers_annotated.rds"))

  Run_CellChat(
    seurat_obj                   = seurat_obj,
    Run_CellChat_output_path     = args$Run_CellChat_output_path,
    Load_QC_species              = args$Load_QC_species,
    Run_CellChat_group_by        = args$Run_CellChat_group_by,
    Run_CellChat_source_celltype = args$Run_CellChat_source_celltype,
    Run_CellChat_target_celltype = args$Run_CellChat_target_celltype,
    Run_CellChat_plot_heatmap    = args$Run_CellChat_plot_heatmap,
    Run_CellChat_ntop_signaling  = args$Run_CellChat_ntop_signaling,
    Run_CellChat_MaxGroup        = args$Run_CellChat_MaxGroup
  )
}
message("========== CellChat complete. ==========")
message("")

# -------------------------
# Step 7: AutoDock Vina docking
# -------------------------
message("========== Start molecular docking ... ==========")

Vina_Docking(
    Run_CellChat_output_path            = args$Run_CellChat_output_path,
    Vina_Docking_output_path            = args$Vina_Docking_output_path,
    Vina_Docking_ligand_ref_file        = args$Vina_Docking_ligand_ref_file,
    Vina_Docking_receptor_ref_file      = args$Vina_Docking_receptor_ref_file,
    Vina_Docking_cas_txt_file           = args$Vina_Docking_cas_txt_file,
    Vina_Docking_use_fda                = args$Vina_Docking_use_fda,
    Vina_Docking_fda_txt                = args$Vina_Docking_fda_txt,
    Vina_Docking_docking_ligand_dir     = args$Vina_Docking_docking_ligand_dir,
    Vina_Docking_docking_receptor_dir   = args$Vina_Docking_docking_receptor_dir,
    Vina_Docking_vina_exhaustiveness    = args$Vina_Docking_vina_exhaustiveness,
    Vina_Docking_vina_num_modes         = args$Vina_Docking_vina_num_modes,
    Vina_Docking_vina_seed              = args$Vina_Docking_vina_seed,
    Vina_Docking_vina_cpu               = args$Vina_Docking_vina_cpu
  )
message("========== Molecular docking complete. ==========")
EOF
