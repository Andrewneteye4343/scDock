# functions/Markers_Annotation.R
Markers_Annotation <- function(seurat_obj,
                               Markers_Annotation_use_assay = NULL,
                               Markers_Annotation_group_by = "seurat_clusters",
                               Markers_Annotation_logfc_threshold = 0.25,
                               Markers_Annotation_test_use = "wilcox",
                               Markers_Annotation_min_pct = 0.01,
                               Markers_Annotation_only_pos = TRUE,
                               Markers_Annotation_verbose = TRUE,
                               Markers_Annotation_top_n = 10,
                               Markers_Annotation_output_path = NULL,
                               Markers_Annotation_tissue_type = NULL,
                               Markers_Annotation_label_column = "CellType") {
  library(Seurat)
  library(dplyr)
  library(scMayoMap)

  # ----------- Step 0: Assay selection -----------
  assay_name <- if (!is.null(Markers_Annotation_use_assay)) {
    Markers_Annotation_use_assay
  } else {
    Seurat::DefaultAssay(seurat_obj)
  }

  assay_obj <- seurat_obj[[assay_name]]

  # 判斷 assay 類型
  assay_class <- class(assay_obj)[1]
  if (Markers_Annotation_verbose) {
    message("Using assay [", assay_name, "] of class: ", assay_class)
  }

  # ----------- Step 1: JoinLayers only if Assay5 -----------
  if (inherits(assay_obj, "Assay5")) {
    if (Markers_Annotation_verbose) {
      message("Detected Assay5 with layers → applying JoinLayers()")
    }
    if (exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
      seurat_obj[[assay_name]] <- SeuratObject::JoinLayers(assay_obj)
    } else {
      message("⚠️ JoinLayers() not found in SeuratObject namespace → skipping")
    }
  } else {
    if (Markers_Annotation_verbose) {
      message("Assay class [", assay_class, "] does not support JoinLayers() → skipping")
    }
  }

  # ----------- Step 2: Find markers -----------
  if (Markers_Annotation_verbose) message("Finding markers using FindAllMarkers() ...")

  markers <- FindAllMarkers(seurat_obj,
                            assay = assay_name,
                            group.by = Markers_Annotation_group_by,
                            logfc.threshold = Markers_Annotation_logfc_threshold,
                            test.use = Markers_Annotation_test_use,
                            min.pct = Markers_Annotation_min_pct,
                            only.pos = Markers_Annotation_only_pos)

  if (Markers_Annotation_verbose) message("Marker finding complete. Total markers: ", nrow(markers))

  if (!is.null(Markers_Annotation_output_path)) {
    Marker_output <- file.path(Markers_Annotation_output_path, "Markers.csv")
    write.csv(markers, Marker_output, row.names = FALSE)
    message("Markers.csv saved to: ", Markers_Annotation_output_path)
  }

  # ----------- Step 3: Select top markers per cluster -----------
  if (!is.null(Markers_Annotation_top_n) && Markers_Annotation_top_n > 0) {
    top_markers <- markers %>%
      filter(avg_log2FC > 0) %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = Markers_Annotation_top_n, with_ties = FALSE) %>%
      ungroup()
  } else {
    top_markers <- markers %>% filter(avg_log2FC > 0)
  }

  if (Markers_Annotation_verbose) message("Running scMayoMap for cell annotation...")
  annotation <- scMayoMap::scMayoMap(data = top_markers, tissue = Markers_Annotation_tissue_type)

  # ----------- Step 4: Save annotation result -----------
  if (!is.null(Markers_Annotation_output_path) && "res" %in% names(annotation)) {
    Annotation_output <- file.path(Markers_Annotation_output_path, "Annotations.csv")
    write.csv(annotation$res, Annotation_output, row.names = FALSE)
    message("Annotation.csv saved to: ", Markers_Annotation_output_path)
  }

  # ----------- Step 5: Add annotation to Seurat object -----------
  if ("markers" %in% names(annotation) && is.data.frame(annotation$markers)) {
    cluster_map_df <- annotation$markers %>%
      mutate(cluster = as.character(cluster),
             score   = suppressWarnings(as.numeric(score))) %>%
      group_by(cluster) %>%
      slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(cluster = as.character(cluster),
                celltype = as.character(celltype))
  } else if ("res" %in% names(annotation) && is.data.frame(annotation$res)) {
    tmp <- annotation$res
    names(tmp) <- tolower(names(tmp))
    if (all(c("cluster", "celltype") %in% names(tmp))) {
      cluster_map_df <- tmp %>%
        mutate(cluster = as.character(cluster),
               celltype = as.character(celltype)) %>%
        distinct(cluster, .keep_all = TRUE) %>%
        select(cluster, celltype)
    } else {
      stop("scMayoMap returned an unexpected structure: cannot find 'markers' or ('cluster','celltype') in 'res'.")
    }
  } else {
    stop("scMayoMap returned an unexpected structure without 'markers' or 'res'.")
  }

  # ----------- Step 6: Map annotation to Seurat metadata -----------
  celltype_map <- setNames(cluster_map_df$celltype, cluster_map_df$cluster)
  clusters_in_obj <- as.character(seurat_obj$seurat_clusters)
  celltype_per_cell <- unname(celltype_map[clusters_in_obj])
  celltype_per_cell[is.na(celltype_per_cell)] <- "Unknown"

  meta_df <- data.frame(
    tmp = celltype_per_cell,
    row.names = colnames(seurat_obj),
    stringsAsFactors = FALSE
  )
  colnames(meta_df) <- Markers_Annotation_label_column

  seurat_obj <- SeuratObject::AddMetaData(
    object   = seurat_obj,
    metadata = meta_df
  )

  return(seurat_obj)
}
