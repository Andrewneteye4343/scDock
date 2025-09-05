# functions/find_markers_and_annotate.R

find_markers_and_annotate <- function(seurat_obj,
                                      use_assay = NULL,
                                      group_by = "seurat_clusters",
                                      logfc_threshold = 0.25,
                                      test_use = "wilcox",
                                      min_pct = 0.01,
                                      min_diff_pct = -Inf,
                                      only_pos = TRUE,
                                      find_markers_and_annotate_verbose = TRUE,
                                      top_n = 10,
                                      output_marker_csv = NULL,
                                      tissue_type = NULL,
                                      output_annotation_csv = NULL,
                                      label_column = "CellType") {
  library(Seurat)
  library(dplyr)
  library(scMayoMap)
  
  message("Starting cluster markers calculating...")
  # ----------- Step 1: Find markers -----------
  if (find_markers_and_annotate_verbose) message("Finding markers using FindAllMarkers() ...")
  markers <- FindAllMarkers(seurat_obj,
                            assay = use_assay,
                            group.by = group_by,
                            logfc.threshold = logfc_threshold,
                            test.use = test_use,
                            min.pct = min_pct,
                            min.diff.pct =min_diff_pct,
                            only.pos = only_pos)
  
  if (find_markers_and_annotate_verbose) message("Marker finding complete. Total markers: ", nrow(markers))
  
  if (!is.null(output_marker_csv)) {
    write.csv(markers, output_marker_csv, row.names = FALSE)
    message("Markers saved to: ", output_marker_csv)
  }
  
  # ----------- Step 2: Select top markers per cluster -----------
  if (!is.null(top_n) && top_n > 0) {
    top_markers <- markers %>%
      filter(avg_log2FC > 0) %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = top_n, with_ties = FALSE) %>%
      ungroup()
  } else {
    top_markers <- markers %>%
      filter(avg_log2FC > 0)
  }
  
  if (find_markers_and_annotate_verbose) message("Running scMayoMap for cell annotation...")
  annotation <- scMayoMap::scMayoMap(data = top_markers, tissue = tissue_type)
  # ----------- Step 3: Save annotation result -----------
  if (!is.null(output_annotation_csv) && "res" %in% names(annotation)) {
    write.csv(annotation$res, output_annotation_csv, row.names = FALSE)
    message("Annotation saved to: ", output_annotation_csv)
  }
  
  # ----------- Step 4: Add annotation to Seurat object -----------
  # For each cluster, it is annotated as the celltype with highest score, which is obtained from annotation$markers.
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
    # Or just use the results from annotation$res.
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

  # Mapping the annotation back to seura_obj.
  celltype_map <- setNames(cluster_map_df$celltype, cluster_map_df$cluster)
  clusters_in_obj <- as.character(seurat_obj$seurat_clusters)

  # Each seurat_clusters will get their annotated celltype name from scMayoMap. If not, the cluster(s) is annotated as "Unknown".
  celltype_per_cell <- unname(celltype_map[clusters_in_obj])
  celltype_per_cell[is.na(celltype_per_cell)] <- "Unknown"

  # Create metadata
  meta_df <- data.frame(
    tmp = celltype_per_cell,
    row.names = colnames(seurat_obj),
    stringsAsFactors = FALSE
  )
  colnames(meta_df) <- label_column

  # ----------- Step 5: Add to Seurat object -----------
  seurat_obj <- SeuratObject::AddMetaData(
    object   = seurat_obj,
    metadata = meta_df
  )
  
  message("✅ Find marker and annotation complete.")
  return(seurat_obj)
}
