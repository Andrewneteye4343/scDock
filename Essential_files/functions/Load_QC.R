Load_QC <- function(Load_QC_input_files,
                    Load_QC_input_type = c("10x", "multi10x", "h5", "txt"),
                    Load_QC_min_features = 200,
                    Load_QC_min_cells = 3,
                    Load_QC_names_delim = "_",
                    Load_QC_max_mito = 0.3,
                    Load_QC_verbose = TRUE,
                    Load_QC_metadata_file = NULL,
                    Load_QC_species = c("human", "mouse")) {

  suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(org.Mm.eg.db)
  })

  Load_QC_input_type <- match.arg(Load_QC_input_type)
  species <- match.arg(Load_QC_species)

  # Load metadata if provided
  meta_info <- NULL
  if (!is.null(Load_QC_metadata_file)) {
    if (!file.exists(Load_QC_metadata_file)) stop("[Load_QC] Metadata file not found: ", Load_QC_metadata_file)
    meta_info <- read.table(Load_QC_metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (!all(c("file_name", "group") %in% colnames(meta_info))) {
      stop("[Load_QC] Metadata file must contain columns: file_name and group")
    }
  }

  # Function: Convert Ensembl IDs → Gene Symbols, and ensure legal, unique feature names
  convert_ensembl_if_needed <- function(seurat_obj, species) {
    gene_ids <- rownames(seurat_obj)
    is_human_ens <- all(grepl("^ENSG[0-9]+", gene_ids))
    is_mouse_ens <- all(grepl("^ENSMUSG[0-9]+", gene_ids))

    if ((species == "human" && is_human_ens) || (species == "mouse" && is_mouse_ens)) {
      if (Load_QC_verbose) message("[Load_QC] Detected Ensembl IDs → Converting to Gene Symbols ...")
      OrgDb <- if (species == "human") org.Hs.eg.db else org.Mm.eg.db

      symbol_map <- mapIds(
        OrgDb,
        keys = gene_ids,
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )

      new_names <- ifelse(!is.na(symbol_map) & symbol_map != "", symbol_map, gene_ids)
      dup_symbols <- new_names[duplicated(new_names)]
      if (length(dup_symbols) > 0) {
        for (dup in dup_symbols) {
          idx <- which(new_names == dup)
          new_names[idx] <- paste0(dup, "_", seq_along(idx))
        }
        if (Load_QC_verbose) message("[Load_QC] Duplicated gene symbols were made unique.")
      }

      rownames(seurat_obj) <- new_names
    }

    # Enforce Seurat-compatible gene names
    clean_names <- rownames(seurat_obj)
    clean_names <- gsub("_", "-", clean_names)
    clean_names <- gsub("[^A-Za-z0-9\\-\\.]", "-", clean_names)
    clean_names <- make.unique(clean_names)
    rownames(seurat_obj) <- clean_names

    return(seurat_obj)
  }

  # Load and process all input datasets
  seurat_list <- list()
  message("[Load_QC] Loading your data...")

  for (i in seq_along(Load_QC_input_files)) {
    input_path <- Load_QC_input_files[i]
    counts <- NULL

    if (Load_QC_input_type == "10x") {
      message("[Load_QC] Input type = 10X folder → Read10X(): ", basename(input_path))
      counts <- Read10X(data.dir = input_path)
    } else if (Load_QC_input_type == "multi10x") {
      subdirs <- list.dirs(input_path, recursive = FALSE, full.names = TRUE)
      valid_subdirs <- subdirs[
        file.exists(file.path(subdirs, "barcodes.tsv.gz")) &
          file.exists(file.path(subdirs, "features.tsv.gz")) &
          file.exists(file.path(subdirs, "matrix.mtx.gz"))
      ]
      for (subdir in valid_subdirs) {
        message("[Load_QC] Input type = 10X subfolder → Read10X(): ", basename(subdir))
        counts <- Read10X(data.dir = subdir)
        seurat_obj <- CreateSeuratObject(counts = counts,
                                         min.features = Load_QC_min_features,
                                         min.cells = Load_QC_min_cells,
                                         names.delim = Load_QC_names_delim,
                                         project = basename(subdir))
        seurat_obj <- convert_ensembl_if_needed(seurat_obj, species)
        seurat_obj$orig.ident <- seurat_obj@project.name
        if (!is.null(meta_info))
          seurat_obj$sample_group <- meta_info$group[match(basename(subdir), meta_info$file_name)]
        seurat_list[[length(seurat_list) + 1]] <- seurat_obj
      }
      next
    } else if (Load_QC_input_type == "h5") {
      h5_files <- if (dir.exists(input_path)) list.files(input_path, pattern = "\\.h5$", full.names = TRUE) else input_path
      for (h5_file in h5_files) {
        message("[Load_QC] Input type = HDF5 → Read10X_h5(): ", basename(h5_file))
        counts <- Read10X_h5(filename = h5_file)
        if (is.list(counts)) {
          for (libname in names(counts)) {
            seurat_obj <- CreateSeuratObject(counts = counts[[libname]],
                                             min.features = Load_QC_min_features,
                                             min.cells = Load_QC_min_cells,
                                             names.delim = Load_QC_names_delim,
                                             project = paste0(tools::file_path_sans_ext(basename(h5_file)), "_", libname))
            seurat_obj <- convert_ensembl_if_needed(seurat_obj, species)
            seurat_obj$orig.ident <- seurat_obj@project.name
            if (!is.null(meta_info))
              seurat_obj$sample_group <- meta_info$group[match(paste0(tools::file_path_sans_ext(basename(h5_file)), "_", libname),
                                                               meta_info$file_name)]
            seurat_list[[length(seurat_list) + 1]] <- seurat_obj
          }
        } else {
          seurat_obj <- CreateSeuratObject(counts = counts,
                                           min.features = Load_QC_min_features,
                                           min.cells = Load_QC_min_cells,
                                           names.delim = Load_QC_names_delim,
                                           project = tools::file_path_sans_ext(basename(h5_file)))
          seurat_obj <- convert_ensembl_if_needed(seurat_obj, species)
          seurat_obj$orig.ident <- seurat_obj@project.name
          if (!is.null(meta_info))
            seurat_obj$sample_group <- meta_info$group[match(basename(h5_file), meta_info$file_name)]
          seurat_list[[length(seurat_list) + 1]] <- seurat_obj
        }
      }
      next
    } else if (Load_QC_input_type == "txt") {
      tab_files <- if (dir.exists(input_path)) list.files(input_path, pattern = "\\.txt$", full.names = TRUE) else input_path
      for (tab_file in tab_files) {
        message("[Load_QC] Input type = TXT → read.table(): ", basename(tab_file))
        raw_mat <- read.table(tab_file, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE)
        counts <- as(as.matrix(raw_mat), "dgCMatrix")
        seurat_obj <- CreateSeuratObject(counts = counts,
                                         min.features = Load_QC_min_features,
                                         min.cells = Load_QC_min_cells,
                                         names.delim = Load_QC_names_delim,
                                         project = tools::file_path_sans_ext(basename(tab_file)))
        seurat_obj <- convert_ensembl_if_needed(seurat_obj, species)
        seurat_obj$orig.ident <- seurat_obj@project.name
        if (!is.null(meta_info))
          seurat_obj$sample_group <- meta_info$group[match(basename(tab_file), meta_info$file_name)]
        seurat_list[[length(seurat_list) + 1]] <- seurat_obj
      }
      next
    }

    if (!is.null(counts)) {
      seurat_obj <- CreateSeuratObject(counts = counts,
                                       min.features = Load_QC_min_features,
                                       min.cells = Load_QC_min_cells,
                                       names.delim = Load_QC_names_delim,
                                       project = tools::file_path_sans_ext(basename(input_path)))
      seurat_obj <- convert_ensembl_if_needed(seurat_obj, species)
      seurat_obj$orig.ident <- seurat_obj@project.name
      if (!is.null(meta_info))
        seurat_obj$sample_group <- meta_info$group[match(basename(input_path), meta_info$file_name)]
      seurat_list[[length(seurat_list) + 1]] <- seurat_obj
    }
  }

  # Quality Control & Merge
  for (j in seq_along(seurat_list)) {
    seurat_obj <- seurat_list[[j]]
    # Human vs mouse mitochondrial gene prefix
    mito_pattern <- if (species == "mouse") "^mt-" else "^MT-"
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mito_pattern)

    raw_ncells <- ncol(seurat_obj)
    raw_ngenes <- nrow(seurat_obj)

    seurat_obj <- subset(seurat_obj,
                         subset = nFeature_RNA > Load_QC_min_features & percent.mt < Load_QC_max_mito * 100)

    if (Load_QC_verbose) {
      message("[Load_QC] QC summary for dataset [", seurat_obj@project.name, "]")
      message("[Load_QC] Cells: ", raw_ncells, " → ", ncol(seurat_obj))
      message("[Load_QC] Genes: ", raw_ngenes, " → ", nrow(seurat_obj))
    }

    seurat_list[[j]] <- seurat_obj
  }

  if (Load_QC_verbose) message("[Load_QC] Merging all datasets ...")
  merged_obj <- seurat_list[[1]]
  if (length(seurat_list) > 1) {
    for (i in 2:length(seurat_list)) merged_obj <- merge(merged_obj, y = seurat_list[[i]])
  }

  if (Load_QC_verbose)
    message("[Load_QC] Merged Seurat object: ", ncol(merged_obj), " cells × ", nrow(merged_obj), " genes")

  merged_obj <- UpdateSeuratObject(merged_obj)
  return(merged_obj)
}
