# functions/Markers_Annotation.R
Markers_Annotation <- function(seurat_obj,
                               Markers_Annotation_use_assay,
                               Markers_Annotation_group_by,
                               Markers_Annotation_logfc_threshold,
                               Markers_Annotation_test_use,
                               Markers_Annotation_min_pct,
                               Markers_Annotation_only_pos,
                               Markers_Annotation_verbose,
                               Markers_Annotation_top_n,
                               Markers_Annotation_output_path,
                               Markers_Annotation_tissue_type,
                               Markers_Annotation_label_column) {
  suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(scMayoMap)
  library(homologene)
})
  # Step 1: Assay selection
  assay_name <- if (!is.null(Markers_Annotation_use_assay)) {
    Markers_Annotation_use_assay
  } else {
    Seurat::DefaultAssay(seurat_obj)
  }

  assay_obj <- seurat_obj[[assay_name]]

  # Check assay type
  assay_class <- class(assay_obj)[1]
  if (Markers_Annotation_verbose) {
    message("[Markers_Annotation] Using assay [", assay_name, "] of class: ", assay_class)
  }

  # Step 2: JoinLayers only if Assay5
  if (inherits(assay_obj, "Assay5")) {
    if (Markers_Annotation_verbose) {
      message("[Markers_Annotation] Detected Assay5 with layers → applying JoinLayers()")
    }
    if (exists("JoinLayers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
      seurat_obj[[assay_name]] <- SeuratObject::JoinLayers(assay_obj)
    } else {
      message("[Markers_Annotation] JoinLayers() not found in SeuratObject namespace → skipping")
    }
  } else {
    if (Markers_Annotation_verbose) {
      message("[Markers_Annotation] Assay class [", assay_class, "] does not support JoinLayers() → skipping")
    }
  }

  # Step 3: Find markers
  if (Markers_Annotation_verbose) message("[Markers_Annotation] Finding markers using FindAllMarkers() ...")

  markers <- FindAllMarkers(seurat_obj,
                            assay = assay_name,
                            group.by = Markers_Annotation_group_by,
                            logfc.threshold = Markers_Annotation_logfc_threshold,
                            test.use = Markers_Annotation_test_use,
                            min.pct = Markers_Annotation_min_pct,
                            only.pos = Markers_Annotation_only_pos)

  if (Markers_Annotation_verbose) message("[Markers_Annotation] Marker finding complete. Total markers: ", nrow(markers))

  if (!is.null(Markers_Annotation_output_path)) {
    Marker_output <- file.path(Markers_Annotation_output_path, "Cluster_markers.csv")
    write.csv(markers, Marker_output, row.names = FALSE)
    message("[Markers_Annotation] Cluster_markers.csv saved to: ", Markers_Annotation_output_path)
  }

  # Step 4: Select top markers per cluster
  if (!is.null(Markers_Annotation_top_n) && Markers_Annotation_top_n > 0) {
    top_markers <- markers %>%
      filter(avg_log2FC > 0) %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = Markers_Annotation_top_n, with_ties = FALSE) %>%
      ungroup()
  } else {
    top_markers <- markers %>% filter(avg_log2FC > 0)
  }

  # Step 4.5: Gene name conversion (mouse → human) 
  convert_FindAllMarkers_mouse_to_human_offline <- function(marker_df, gene_col = "gene", keep_unique = TRUE) {
    if (!gene_col %in% colnames(marker_df))
      stop(paste("❌ Can't find:", gene_col))

    mouse_genes <- unique(marker_df[[gene_col]])
    # Mouse NCBI Taxonomy ID: 10090; Human NCBI Taxonomy ID: 9606
    gene_map <- homologene(mouse_genes, inTax = 10090, outTax = 9606)
    colnames(gene_map) <- c("MouseGene", "HumanGene")
    map_dict <- setNames(gene_map$HumanGene, gene_map$MouseGene)

    marker_df[[gene_col]] <- ifelse(marker_df[[gene_col]] %in% names(map_dict),
                                    map_dict[marker_df[[gene_col]]],
                                    marker_df[[gene_col]])
    if (keep_unique)
      marker_df <- marker_df[!duplicated(marker_df[[gene_col]]), ]
    rownames(marker_df) <- NULL
    return(marker_df)
  }

  gene_names <- unique(top_markers$gene)
  proportion_lower <- mean(grepl("[a-z]", gene_names))
  if (proportion_lower > 0.5) {
    message("[Markers_Annotation] Detected mouse-like gene symbols (", round(proportion_lower*100), "% lowercase) → converting to human ...")
    top_markers <- convert_FindAllMarkers_mouse_to_human_offline(top_markers)
    message("[Markers_Annotation] Gene names converted to human symbols.")
  } else {
    message("[Markers_Annotation] Detected human symbols.")
  }

# Extra step: Extend annotation database of scMayoMap
# === Neuroblastoma start ===
original_db <- scMayoMapDatabase

# Set Neuroblastoma cell types and marker
cell_types_new <- c(
  "neuroblastoma:B cell",
  "neuroblastoma:Endothelial cell",
  "neuroblastoma:Fibroblast",
  "neuroblastoma:Neuroendocrine tumor cell",
  "neuroblastoma:NK cell",
  "neuroblastoma:Plasmacytoid dendritic cell",
  "neuroblastoma:Conventional dendritic cell",
  "neuroblastoma:Monocyte",
  "neuroblastoma:Macrophage",
  "neuroblastoma:Plasma cell",
  "neuroblastoma:Red blood cell",
  "neuroblastoma:Schwann cell",
  "neuroblastoma:Mesenchyme",
  "neuroblastoma:T cell",
  "neuroblastoma:Other stromal cell"
)

marker_list <- list(
  "neuroblastoma:B cell" = c(
    "MS4A1","CD19","CD20" ,"CD79A","CD79B","PAX5","CD37","BANK1","TCL1A",
    "IGKC","IGHM","IGHD","CD24","FCER2","CR2"
  ),
  "neuroblastoma:Endothelial cell" = c(
    "PECAM1","VWF","ENG","KDR","CD34","TEK","CLDN5","CDH5","TIE1","ESM1",
    "ACKR1","FLT1","SELP"
  ),
  "neuroblastoma:Fibroblast" = c(
    "COL1A1","COL1A2","FAP","PDGFRA","DCN","ACTA2","PDGFRB","THY1","S100A4",
    "MMP11","MMP2","ELN","FBLN1","COL3A1","POSTN"
  ),
  "neuroblastoma:Neuroendocrine tumor cell" = c(
    "PHOX2B","PHOX2A","ASCL1","LMO1","LIN28B","TH","DBH","DDC","CHGA","CHGB",
    "SYP","ENO2","DCX","ISL1","GATA3","MYCN","ALK","SLC6A2","SLC18A1",
    "DLK1","NTRK1","NTRK2","SOX11","MKI67","TOP2A","PCNA","B4GALNT1","ST8SIA1",
    "UCHL1","GAP43","NEFM","NEFL","TUBB3","MAP2"
  ),
  "neuroblastoma:Neuron" = c(
    "RBFOX3","MAP2","NEFL","NEFM","TUBB3","SYN1","SYN2","SNAP25","SYT1",
    "SLC17A7","SLC17A6","GRIN1","GRIA1","GABRA1","GAD1","GAD2","SCN1A","SCN2A",
    "KCNQ2","NLGN1"
  ),
  "neuroblastoma:NK cell" = c(
    "NCAM1","KLRD1","NKG7","GNLY","PRF1","KLRF1","KLRC1","KLRC2","FCGR3A",
    "NCR1","NCR3","GZMA","GZMB","KIR2DL1"
  ),
  "neuroblastoma:Plasmacytoid dendritic cell" = c(
    "CLEC4C","IL3RA","TCF4","LILRA4","IRF7","GZMB","LILRB4","SERPINF1"
  ),
  "neuroblastoma:Conventional dendritic cell" = c(
    "ITGAX","CLEC9A","XCR1","CD1C","FCER1A","ZBTB46","BATF3","FLT3","CCR7"
  ),
  "neuroblastoma:Monocyte" = c(
    "CD14","LYZ","FCN1","S100A8","S100A9","VCAN","FCGR3A","CST3","SELPLG"
  ),
  "neuroblastoma:Macrophage" = c(
    "CD68","CD163","MRC1","CSF1R","CD14","CD86","CD80","MSR1","MARCO","LYVE1"
  ),
  "neuroblastoma:Plasma cell" = c(
    "MZB1","XBP1","PRDM1","SDC1","IGHG1","IGKC","IGLL5","JCHAIN","TXNDC5"
  ),
  "neuroblastoma:Red blood cell" = c(
    "HBB","HBA1","HBA2","GYPA","ALAS2","EPB41","SLC4A1"
  ),
  "neuroblastoma:Schwann cell" = c(
    "SOX10","S100B","MPZ","PLP1","PMP22","ERBB3","NGFR","FOXD3","GJB1"
  ),
  "neuroblastoma:Mesenchyme" = c(
    "PDGFRA","COL3A1","VIM","FN1","LUM","POSTN","PDGFRB","THY1","TAGLN"
  ),
  "neuroblastoma:T cell" = c(
    "CD3D","CD3E","CD2","CD4","IL7R","CD8A","CD8B","GZMB","PRF1",
    "TRAC","TRBC1","CCR7","SELL","FOXP3"
  ),
  "neuroblastoma:Other stromal cell" = c(
    "COL6A1","PDPN","LUM","MMP2","SPARC","RGS5","CSPG4","MCAM","CD146"
  )
)

# Establish Neuroblastoma binary matrix
all_genes_new <- unique(unlist(marker_list))

new_db <- matrix(0, nrow = length(all_genes_new), ncol = length(cell_types_new),
                 dimnames = list(all_genes_new, cell_types_new))

for(ct in cell_types_new){
  genes <- marker_list[[ct]]
  new_db[genes, ct] <- 1
}

new_db_df <- as.data.frame(new_db)
new_db_df$gene <- rownames(new_db_df)
new_db_df$tissue <- "neuroblastoma"
new_db_df <- new_db_df[, c("tissue","gene", cell_types_new)]

# Align all cell type column
all_cell_types <- union(colnames(original_db)[3:ncol(original_db)], cell_types_new)

expand_cols <- function(db, all_cell_types){
  missing_cols <- setdiff(all_cell_types, colnames(db))
  for(mc in missing_cols){
    db[[mc]] <- 0
  }
  db <- db[, c("tissue","gene", all_cell_types)]
  return(db)
}

original_db_exp <- expand_cols(original_db, all_cell_types)
new_db_exp <- expand_cols(new_db_df, all_cell_types)

# Combination
combined_db <- rbind(original_db_exp, new_db_exp)
rownames(combined_db) <- NULL
# === Neuroblastoma end ===

# === Breast cancer start ===
# Set Breast cancer cell types and marker
cell_types_new <- c(
  "breast_cancer:Malignant epithelial cell",
  "breast_cancer:Basal-like tumor cell",
  "breast_cancer:Luminal tumor cell",
  "breast_cancer:HER2-enriched tumor cell",
  "breast_cancer:Endothelial cell",
  "breast_cancer:Fibroblast",
  "breast_cancer:Myofibroblast",
  "breast_cancer:Pericyte",
  "breast_cancer:Adipocyte",
  "breast_cancer:T cell",
  "breast_cancer:B cell",
  "breast_cancer:Plasma cell",
  "breast_cancer:Macrophage",
  "breast_cancer:Monocyte",
  "breast_cancer:Dendritic cell",
  "breast_cancer:NK cell",
  "breast_cancer:Mast cell",
  "breast_cancer:Red blood cell",
  "breast_cancer:Other stromal cell"
)

marker_list <- list(
  "breast_cancer:Malignant epithelial cell" = c(
    "EPCAM","KRT8","KRT18","KRT19","CDH1","MUC1","CLDN4","AGR2","TFF1",
    "TP63","SOX9","PROM1","CEACAM6","S100A6","CSTB","MKI67","PCNA"
  ),
  "breast_cancer:Basal-like tumor cell" = c(
    "KRT5","KRT6A","KRT14","KRT17","TP63","EGFR","CD44","ITGA6",
    "VIM","LAMC2","CLDN1","FABP5","S100A2","MIA","LGALS7"
  ),
  "breast_cancer:Luminal tumor cell" = c(
    "ESR1","PGR","GATA3","FOXA1","KRT8","KRT18","XBP1","BCL2",
    "SCGB2A1","TFF1","AGR2","EPCAM","CDH1","MUC1","CLDN3"
  ),
  "breast_cancer:HER2-enriched tumor cell" = c(
    "ERBB2","GRB7","PGAP3","STARD3","CDK12","TCAP","CD24","CDH1",
    "KRT8","KRT18","MUC1","CLDN4","CYP24A1","MKI67"
  ),
  "breast_cancer:Endothelial cell" = c(
    "PECAM1","VWF","ENG","KDR","CD34","TEK","CLDN5","CDH5","TIE1","ESM1",
    "ACKR1","FLT1","SELP"
  ),
  "breast_cancer:Fibroblast" = c(
    "COL1A1","COL1A2","FAP","PDGFRA","DCN","ACTA2","PDGFRB","THY1","S100A4",
    "MMP11","MMP2","ELN","FBLN1","COL3A1","POSTN"
  ),
  "breast_cancer:Myofibroblast" = c(
    "ACTA2","TAGLN","MYH11","CNN1","TPM2","PDGFRB","DES","CALD1","RGS5","VIM"
  ),
  "breast_cancer:Pericyte" = c(
    "RGS5","PDGFRB","MCAM","CSPG4","CD146","ACTA2","TAGLN","NOTCH3","KCNJ8"
  ),
  "breast_cancer:Adipocyte" = c(
    "PLIN1","PLIN4","ADIPOQ","FABP4","LPL","LIPE","CIDEC","LEP","PPARG","CEBPA"
  ),
  "breast_cancer:T cell" = c(
    "CD3D","CD3E","CD2","CD4","IL7R","CD8A","CD8B","GZMB","PRF1",
    "TRAC","TRBC1","CCR7","SELL","FOXP3"
  ),
  "breast_cancer:B cell" = c(
    "MS4A1","CD19","CD20" ,"CD79A","CD79B","PAX5","CD37","BANK1","TCL1A",
    "IGKC","IGHM","IGHD","CD24","FCER2","CR2"
  ),
  "breast_cancer:Plasma cell" = c(
    "MZB1","XBP1","PRDM1","SDC1","IGHG1","IGKC","IGLL5","JCHAIN","TXNDC5"
  ),
  "breast_cancer:Macrophage" = c(
    "CD68","CD163","MRC1","CSF1R","CD14","CD86","CD80","MSR1","MARCO","LYVE1"
  ),
  "breast_cancer:Monocyte" = c(
    "CD14","LYZ","S100A8","S100A9","VCAN","FCN1","FCGR3A","CST3","SELL"
  ),
  "breast_cancer:Conventional dendritic cell" = c(
    "ITGAX","CLEC9A","XCR1","CD1C","FCER1A","ZBTB46","BATF3","FLT3","CCR7"
  ),
  "breast_cancer:Plasmacytoid dendritic cell" = c(
    "CLEC4C","IL3RA","TCF4","LILRA4","IRF7","GZMB","LILRB4","SERPINF1"
  ),
  "breast_cancer:NK cell" = c(
    "NCAM1","KLRD1","NKG7","GNLY","PRF1","KLRF1","KLRC1","KLRC2","FCGR3A",
    "NCR1","NCR3","GZMA","GZMB","KIR2DL1"
  ),
  "breast_cancer:Mast cell" = c(
    "TPSAB1","TPSB2","CPA3","KIT","MS4A2","FCER1A","HDC","HPGDS","CD63"
  ),
  "breast_cancer:Red blood cell" = c(
    "HBB","HBA1","HBA2","GYPA","ALAS2","EPB41","SLC4A1"
  ),
  "breast_cancer:Other stromal cell" = c(
    "COL6A1","PDPN","SPARC","FN1","MMP2","CSPG4","RGS5","CD146","VIM"
  )
)
# Establish breast cancer binary matrix
all_genes_new <- unique(unlist(marker_list))

new_db <- matrix(0, nrow = length(all_genes_new), ncol = length(cell_types_new),
                 dimnames = list(all_genes_new, cell_types_new))

for(ct in cell_types_new){
  genes <- marker_list[[ct]]
  new_db[genes, ct] <- 1
}

new_db_df <- as.data.frame(new_db)
new_db_df$gene <- rownames(new_db_df)
new_db_df$tissue <- "breast_cancer"
new_db_df <- new_db_df[, c("tissue","gene", cell_types_new)]

# Align all cell type column
all_cell_types <- union(colnames(combined_db)[3:ncol(combined_db)], cell_types_new)
combined_db <- expand_cols(combined_db, all_cell_types)
new_db_exp <- expand_cols(new_db_df, all_cell_types)

# Combination
combined_db <- rbind(combined_db, new_db_exp)
rownames(combined_db) <- NULL
# === Breast cancer end ===

if (Markers_Annotation_verbose) message("[Markers_Annotation] Running scMayoMap for cell annotation...")
annotation <- scMayoMap::scMayoMap(data = top_markers, tissue = Markers_Annotation_tissue_type, database = combined_db)

# Step 5: Save annotation result
if (!is.null(Markers_Annotation_output_path) && "res" %in% names(annotation)) {
  Annotation_output <- file.path(Markers_Annotation_output_path, "Annotation_results.csv")
  write.csv(annotation$res, Annotation_output, row.names = FALSE)
  message("[Markers_Annotation] Annotation_results.csv saved to: ", Markers_Annotation_output_path)
}

# Step 6: Add annotation to Seurat object
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
    stop("[Markers_Annotation] scMayoMap returned an unexpected structure: cannot find 'markers' or ('cluster','celltype') in 'res'.")
  }
} else {
  stop("[Markers_Annotation] scMayoMap returned an unexpected structure without 'markers' or 'res'.")
}

# Step 7: Map annotation to Seurat metadata
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
