Load_QC <- function(Load_QC_input_files,
                    Load_QC_input_type = c("10x", "multi10x", "rds", "h5", "txt"),
                    Load_QC_min_features = 200,
                    Load_QC_min_cells = 3,
                    Load_QC_names_delim = "_",
                    Load_QC_max_mito = 0.3,
                    Load_QC_verbose = TRUE,
                    Load_QC_metadata_file = NULL) {
  library(Seurat)
  library(Matrix)

  Load_QC_input_type <- match.arg(Load_QC_input_type)

  # ---- Read metadata if provided ----
  meta_info <- NULL
  if (!is.null(Load_QC_metadata_file)) {
    if (!file.exists(Load_QC_metadata_file)) stop("Metadata file not found: ", Load_QC_metadata_file)
    meta_info <- read.table(Load_QC_metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (!all(c("file_name", "group") %in% colnames(meta_info))) {
      stop("Metadata file must contain columns: file_name and group")
    }
  }

  seurat_list <- list()
  message("Loading your data...")

  for (i in seq_along(Load_QC_input_files)) {
    input_path <- Load_QC_input_files[i]
    counts <- NULL

    # ---- handle different input types ----
    if (Load_QC_input_type == "10x") {
      message("Input type = 10X folder → Read10X(): ", basename(input_path))
      counts <- Read10X(data.dir = input_path)
    } else if (Load_QC_input_type == "multi10x") {
      subdirs <- list.dirs(input_path, recursive = FALSE, full.names = TRUE)
      valid_subdirs <- subdirs[
        file.exists(file.path(subdirs, "barcodes.tsv.gz")) &
          file.exists(file.path(subdirs, "features.tsv.gz")) &
          file.exists(file.path(subdirs, "matrix.mtx.gz"))
      ]
      if (length(valid_subdirs) == 0) warning("No valid 10X subfolders found in: ", input_path)
      for (subdir in valid_subdirs) {
        message("Input type = 10X subfolder → Read10X(): ", basename(subdir))
        counts <- Read10X(data.dir = subdir)
        seurat_obj <- CreateSeuratObject(
          counts = counts,
          min.features = Load_QC_min_features,
          min.cells = Load_QC_min_cells,
          names.delim = Load_QC_names_delim,
          project = basename(subdir)
        )
        # Add group if metadata exists
        if (!is.null(meta_info)) {
          seurat_obj$sample_group <- meta_info$group[match(basename(subdir), meta_info$file_name)]
        }
        seurat_list[[length(seurat_list) + 1]] <- seurat_obj
      }
      next
    } else if (Load_QC_input_type == "rds") {
      rds_files <- if (dir.exists(input_path)) list.files(input_path, pattern = "\\.rds$", full.names = TRUE) else input_path
      for (rds_file in rds_files) {
        message("Input type = RDS → readRDS(): ", basename(rds_file))
        obj <- readRDS(rds_file)
        if (inherits(obj, "Seurat")) {
          message("  ✅ Detected Seurat object → using directly")
          seurat_obj <- obj
        } else if (is.matrix(obj) || inherits(obj, "dgCMatrix")) {
          message("  📦 Detected Matrix object → wrapping into Seurat")
          seurat_obj <- CreateSeuratObject(
            counts = obj,
            min.features = Load_QC_min_features,
            min.cells = Load_QC_min_cells,
            names.delim = Load_QC_names_delim,
            project = tools::file_path_sans_ext(basename(rds_file))
          )
        } else {
          stop("❌ Unsupported object type in RDS file: ", paste(class(obj), collapse = ", "))
        }
        if (!is.null(meta_info)) {
          seurat_obj$sample_group <- meta_info$group[match(basename(rds_file), meta_info$file_name)]
        }
        seurat_list[[length(seurat_list) + 1]] <- seurat_obj
      }
      next
    } else if (Load_QC_input_type == "h5") {
      h5_files <- if (dir.exists(input_path)) list.files(input_path, pattern = "\\.h5$", full.names = TRUE) else input_path
      for (h5_file in h5_files) {
        message("Input type = HDF5 → Read10X_h5(): ", basename(h5_file))
        counts <- Read10X_h5(filename = h5_file)
        if (is.list(counts)) {
          for (libname in names(counts)) {
            seurat_obj <- CreateSeuratObject(
              counts = counts[[libname]],
              min.features = Load_QC_min_features,
              min.cells = Load_QC_min_cells,
              names.delim = Load_QC_names_delim,
              project = paste0(tools::file_path_sans_ext(basename(h5_file)), "_", libname)
            )
            if (!is.null(meta_info)) {
              seurat_obj$sample_group <- meta_info$group[match(paste0(tools::file_path_sans_ext(basename(h5_file)), "_", libname),
                                                               meta_info$file_name)]
            }
            seurat_list[[length(seurat_list) + 1]] <- seurat_obj
          }
        } else {
          seurat_obj <- CreateSeuratObject(
            counts = counts,
            min.features = Load_QC_min_features,
            min.cells = Load_QC_min_cells,
            names.delim = Load_QC_names_delim,
            project = tools::file_path_sans_ext(basename(h5_file))
          )
          if (!is.null(meta_info)) {
            seurat_obj$sample_group <- meta_info$group[match(basename(h5_file), meta_info$file_name)]
          }
          seurat_list[[length(seurat_list) + 1]] <- seurat_obj
        }
      }
      next
    } else if (Load_QC_input_type == "txt") {
      tab_files <- if (dir.exists(input_path)) list.files(input_path, pattern = "\\.txt$", full.names = TRUE) else input_path
      for (tab_file in tab_files) {
        message("Input type = TXT → read.table(): ", basename(tab_file))
        raw_mat <- read.table(tab_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
        counts <- as(as.matrix(raw_mat), "dgCMatrix")
        seurat_obj <- CreateSeuratObject(
          counts = counts,
          min.features = Load_QC_min_features,
          min.cells = Load_QC_min_cells,
          names.delim = Load_QC_names_delim,
          project = tools::file_path_sans_ext(basename(tab_file))
        )
        if (!is.null(meta_info)) {
          seurat_obj$sample_group <- meta_info$group[match(basename(tab_file), meta_info$file_name)]
        }
        seurat_list[[length(seurat_list) + 1]] <- seurat_obj
      }
      next
    }

    # fallback
    if (!is.null(counts)) {
      seurat_obj <- CreateSeuratObject(
        counts = counts,
        min.features = Load_QC_min_features,
        min.cells = Load_QC_min_cells,
        names.delim = Load_QC_names_delim,
        project = tools::file_path_sans_ext(basename(input_path))
      )
      if (!is.null(meta_info)) {
        seurat_obj$sample_group <- meta_info$group[match(basename(input_path), meta_info$file_name)]
      }
      seurat_list[[length(seurat_list) + 1]] <- seurat_obj
    }
  }

  # ---- QC + Merge ----
  for (j in seq_along(seurat_list)) {
    seurat_obj <- seurat_list[[j]]
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

    raw_ncells <- ncol(seurat_obj)
    raw_ngenes <- nrow(seurat_obj)

    seurat_obj <- subset(
      seurat_obj,
      subset = nFeature_RNA > Load_QC_min_features & percent.mt < Load_QC_max_mito * 100
    )

    qc_ncells <- ncol(seurat_obj)
    qc_ngenes <- nrow(seurat_obj)

    if (Load_QC_verbose) {
      message("QC summary for dataset [", seurat_obj@project.name, "]")
      message("Cells:  ", raw_ncells, " → ", qc_ncells)
      message("Genes:  ", raw_ngenes, " → ", qc_ngenes)
    }

    seurat_list[[j]] <- seurat_obj
  }

  if (Load_QC_verbose) message("Merging all datasets ...")
  merged_obj <- seurat_list[[1]]
  if (length(seurat_list) > 1) {
    for (i in 2:length(seurat_list)) {
      merged_obj <- merge(merged_obj, y = seurat_list[[i]])
    }
  }

  if (Load_QC_verbose) {
    message("✅ Merged Seurat object: ", ncol(merged_obj), " cells × ", nrow(merged_obj), " genes")
  }

  merged_obj <- UpdateSeuratObject(merged_obj)
  return(merged_obj)
}