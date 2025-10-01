# functions/Run_CellChat.R
Run_CellChat <- function(seurat_obj,
                         Run_CellChat_output_path = "path/to/dir/",
                         Run_CellChat_species = "human",
                         Run_CellChat_group_by = "Celltype",
                         Run_CellChat_source_celltype = NULL,
                         Run_CellChat_target_celltype = NULL,
                         Run_CellChat_plot_heatmap = TRUE,
                         Run_CellChat_ntop_signaling = 5,
                         Run_CellChat_pdf_width = 6,
                         Run_CellChat_pdf_height = 6,
                         Run_CellChat_MaxGroup = NULL) {  # 新增參數

  library(CellChat)
  library(patchwork)
  library(scales)
  library(ComplexHeatmap)

  if (Run_CellChat_species == "human") {
    CellChatDB.use <- CellChatDB.human
  } else if (Run_CellChat_species == "mouse") {
    CellChatDB.use <- CellChatDB.mouse
  } else {
    stop("❌ Unsupported species: ", Run_CellChat_species)
  }

  dir.create(Run_CellChat_output_path, showWarnings = FALSE, recursive = TRUE)

  # ================================
  # 檢查 assay → 若有 SCT 用 SCT，否則用 RNA
  # ================================
  if ("SCT" %in% names(seurat_obj@assays)) {
    DefaultAssay(seurat_obj) <- "SCT"
    message("📌 Using SCT assay for CellChat input")
  } else if ("RNA" %in% names(seurat_obj@assays)) {
    DefaultAssay(seurat_obj) <- "RNA"
    message("📌 Using RNA assay for CellChat input")
  } else {
    stop("❌ Neither SCT nor RNA assay found in Seurat object. Cannot run CellChat.")
  }

  # ================================
  # 判斷單組別或多組別
  # ================================
  is_multi_group <- FALSE
  groups <- NULL
  cellchat <- NULL
  cellchat_list <- NULL

  if (Run_CellChat_group_by %in% colnames(seurat_obj@meta.data) &&
      length(unique(seurat_obj@meta.data[[Run_CellChat_group_by]])) > 1 &&
      Run_CellChat_group_by == "sample_group") {

    is_multi_group <- TRUE
    groups <- unique(seurat_obj@meta.data[[Run_CellChat_group_by]])
    message("🔄 Multi-group CellChat analysis based on: ", Run_CellChat_group_by, "; groups: ", paste(groups, collapse = ", "))

    # 檢查 Run_CellChat_MaxGroup 是否有效
    if (!is.null(Run_CellChat_MaxGroup)) {
      if (!all(Run_CellChat_MaxGroup %in% groups)) {
        stop("❌ Run_CellChat_MaxGroup contains invalid group names: ", paste(setdiff(Run_CellChat_MaxGroup, groups), collapse = ", "))
      }
    } else {
      Run_CellChat_MaxGroup <- groups  # 若未指定，預設使用全部組別
    }

    # ================================
    # 建立各組 CellChat 對象
    # ================================
    cellchat_list <- list()
    for (grp in groups) {
      message("▶ Running CellChat for group: ", grp)
      sub_obj <- subset(seurat_obj,
                        cells = rownames(seurat_obj@meta.data[seurat_obj@meta.data[[Run_CellChat_group_by]] == grp, ]))
      sub_cellchat <- createCellChat(object = GetAssayData(sub_obj, slot = "data"),
                                     meta = sub_obj@meta.data,
                                     group.by = "Celltype")
      sub_cellchat@DB <- CellChatDB.use
      sub_cellchat <- subsetData(sub_cellchat)
      sub_cellchat <- identifyOverExpressedGenes(sub_cellchat)
      sub_cellchat <- identifyOverExpressedInteractions(sub_cellchat)
      sub_cellchat <- computeCommunProb(sub_cellchat)
      sub_cellchat <- filterCommunication(sub_cellchat, min.cells = 10)
      sub_cellchat <- computeCommunProbPathway(sub_cellchat)
      sub_cellchat <- aggregateNet(sub_cellchat)
      sub_cellchat <- netAnalysis_computeCentrality(sub_cellchat, slot.name = "netP")
      cellchat_list[[grp]] <- sub_cellchat
    }

    # merge multi-group CellChat object
    cellchat <- mergeCellChat(cellchat_list, add.names = groups)
    cellchat@DB <- CellChatDB.use
    cellchat <- updateCellChat(cellchat)
    saveRDS(cellchat, file = file.path(Run_CellChat_output_path, "cellchat_merged.rds"))
    message("✅ Merged multi-group cellchat object saved to: ", Run_CellChat_output_path)

  } else {
    # ================================
    # 單組別
    # ================================
    message("Single-group CellChat analysis, grouping by: ", Run_CellChat_group_by)
    cellchat <- createCellChat(object = GetAssayData(seurat_obj, slot = "data"),
                               meta = seurat_obj@meta.data,
                               group.by = Run_CellChat_group_by)
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    saveRDS(cellchat, file = file.path(Run_CellChat_output_path, "cellchat.rds"))
    message("✅ Single-group cellchat object saved to: ", Run_CellChat_output_path)
  }

  # ================================
  # signaling role heatmap
  # ================================
  if (Run_CellChat_plot_heatmap) {
    message("Drawing netAnalysis_signalingRole_heatmap ...")

    if (Run_CellChat_plot_heatmap) {
    message("Drawing netAnalysis_signalingRole_heatmap ...")

    if (is_multi_group) {
      message("Multi-group mode: drawing heatmaps for each group with unified pathway set")

      # 取得所有 pathways 的 union
      pathway.union <- unique(unlist(lapply(cellchat_list, function(x) x@netP$pathways)))

      # --- outgoing ---
      ht_list_out <- list()
      for (i in seq_along(cellchat_list)) {
        grp_name <- names(cellchat_list)[i]
        ht_list_out[[i]] <- netAnalysis_signalingRole_heatmap(
          cellchat_list[[i]],
          pattern = "outgoing",
          signaling = pathway.union,
          title = grp_name,
          width = 6, height = 6
        )
      }

      pdf(file = file.path(Run_CellChat_output_path, "signalingRole_heatmap_outgoing.pdf"),
          width = Run_CellChat_pdf_width * length(cellchat_list), 
          height = Run_CellChat_pdf_height * 1.5)
      draw(Reduce(`+`, ht_list_out), ht_gap = unit(0.5, "cm"))
      dev.off()

      # --- incoming ---
      ht_list_in <- list()
      for (i in seq_along(cellchat_list)) {
        grp_name <- names(cellchat_list)[i]
        ht_list_in[[i]] <- netAnalysis_signalingRole_heatmap(
          cellchat_list[[i]],
          pattern = "incoming",
          signaling = pathway.union,
          title = grp_name,
          width = 6, height = 6,
          color.heatmap = "GnBu"
        )
      }

      pdf(file = file.path(Run_CellChat_output_path, "signalingRole_heatmap_incoming.pdf"),
          width = Run_CellChat_pdf_width * length(cellchat_list), 
          height = Run_CellChat_pdf_height * 1.5)
      draw(Reduce(`+`, ht_list_in), ht_gap = unit(0.5, "cm"))
      dev.off()

      message("✅ signalingRole_heatmap (multi-group, union pathways) saved: incoming & outgoing")

    } else {
      # --- 單組別 ---
      ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
      ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.heatmap = "GnBu")
      pdf(file = file.path(Run_CellChat_output_path, "signalingRole_heatmap.pdf"),
          width = Run_CellChat_pdf_width * 1.5, height = Run_CellChat_pdf_height * 1.5)
      print(ht1 + ht2)
      dev.off()
      message("✅ signalingRole_heatmap (single-group) saved")
      }
    }
  }
  # ================================
  # L-R pair selection
  # ================================
  if (!is_multi_group) {
    # -------- 單組別 --------
    comm <- subsetCommunication(cellchat, slot.name = "net")
    pcol <- if ("pval" %in% colnames(comm)) "pval" else if ("p.value" %in% colnames(comm)) "p.value" else NULL
    if (!is.null(pcol)) comm <- comm[comm[[pcol]] < 0.05, ]
    if ("prob" %in% colnames(comm)) comm <- comm[order(-comm$prob), ]
    if (nrow(comm) > 0) {
      top_comm <- head(comm, Run_CellChat_ntop_signaling)
      csv_file <- file.path(Run_CellChat_output_path, paste0("top", Run_CellChat_ntop_signaling, "_LR_singleGroup.csv"))
      write.csv(top_comm, file = csv_file, row.names = FALSE)
      message("✅ Saved top L-R pairs to: ", csv_file)

      pdf(file = file.path(Run_CellChat_output_path, "netVisual_bubble_singleGroup.pdf"),
          width = Run_CellChat_pdf_width,
          height = Run_CellChat_pdf_height)
      print(netVisual_bubble(cellchat,
                             signaling = unique(as.character(top_comm$pathway_name)),
                             sources.use = Run_CellChat_source_celltype,
                             targets.use = Run_CellChat_target_celltype,
                             remove.isolate = TRUE))
      dev.off()
      message("✅ Single-group bubble plot saved")
    }
  } else {
    # -------- 多組別 --------
    message("Performing multi-group bubble plots based on Run_CellChat_MaxGroup ...")
    for (max_grp in Run_CellChat_MaxGroup) {
      max_idx <- which(groups == max_grp)
      comparison_vec <- c(1:length(groups))
      n_groups <- length(comparison_vec)
      color_text <- RColorBrewer::brewer.pal(min(n_groups, 8), "Set2")
      if (n_groups > 8) {
        color_text <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_groups)
      }
      
      pdf(file = file.path(Run_CellChat_output_path, paste0("netVisual_bubble_maxGroup_", max_grp, ".pdf")),
          width = Run_CellChat_pdf_width * 1.5, height = Run_CellChat_pdf_height * 1.5)
      print(netVisual_bubble(cellchat,
                             comparison = comparison_vec,
                             max.dataset = max_idx,
                             sources.use = Run_CellChat_source_celltype,
                             targets.use = Run_CellChat_target_celltype,
                             remove.isolate = TRUE,
                             title.name = paste("Max signaling in", max_grp),
                             color.text = color_text))
      dev.off()
      message("✅ Bubble plot saved for max group: ", max_grp)
    }

    # -------- 計算多組別 L-R 差異 --------
    comm_list <- lapply(cellchat_list, function(cc) {
      df <- subsetCommunication(cc, slot.name = "net")
      if ("p.value" %in% colnames(df) && !"pval" %in% colnames(df)) df$pval <- df[["p.value"]]
      if (!("interaction_name" %in% colnames(df)) && "interaction" %in% colnames(df)) df$interaction_name <- df$interaction
      if (!("source" %in% colnames(df)) && "from" %in% colnames(df)) df$source <- df$from
      if (!("target" %in% colnames(df)) && "to" %in% colnames(df)) df$target <- df$to
      if ("pval" %in% colnames(df)) df <- df[df$pval < 0.05, , drop = FALSE]
      if (nrow(df) > 0) df$key <- paste(df$interaction_name, df$source, df$target, sep = "|")
      df
    })
    names(comm_list) <- groups

    sig_details <- list()
    for (i in seq_along(groups)) {
      for (j in seq_along(groups)) {
        if (i >= j) next
        g1 <- groups[i]; g2 <- groups[j]
        comm1 <- comm_list[[g1]]; comm2 <- comm_list[[g2]]
        keys_to_check <- unique(c(comm1$key, comm2$key))
        if (length(keys_to_check) == 0) next
        for (k in keys_to_check) {
          row1 <- comm1[comm1$key == k, , drop = FALSE]
          row2 <- comm2[comm2$key == k, , drop = FALSE]
          prob1 <- if (nrow(row1) > 0) as.numeric(row1$prob[1]) else 0
          prob2 <- if (nrow(row2) > 0) as.numeric(row2$prob[1]) else 0
          pval1 <- if (nrow(row1) > 0) as.numeric(row1$pval[1]) else NA
          pval2 <- if (nrow(row2) > 0) as.numeric(row2$pval[1]) else NA
          parts <- strsplit(k, "\\|")[[1]]
          sig_details[[length(sig_details) + 1]] <- data.frame(
            interaction = parts[1],
            source = parts[2],
            target = parts[3],
            key = k,
            group1 = g1,
            group2 = g2,
            prob1 = prob1,
            prob2 = prob2,
            pval1 = pval1,
            pval2 = pval2,
            prob_diff = prob1 - prob2,
            direction = paste(g1, "vs", g2),
            higher_group = ifelse(prob1 > prob2, g1, g2),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (length(sig_details) > 0) {
      diff_df <- do.call(rbind, sig_details)
  
      # 只保留 higher_group 在 Run_CellChat_MaxGroup
      diff_df <- diff_df[diff_df$higher_group %in% Run_CellChat_MaxGroup, , drop = FALSE]
  
      # 如果指定 source_celltype / target_celltype，過濾
      if (!is.null(Run_CellChat_source_celltype)) {
        diff_df <- diff_df[diff_df$source %in% Run_CellChat_source_celltype, , drop = FALSE]
      }
      if (!is.null(Run_CellChat_target_celltype)) {
        diff_df <- diff_df[diff_df$target %in% Run_CellChat_target_celltype, , drop = FALSE]
      }
  
      # 按 prob_diff 絕對值排序
      diff_df <- diff_df[order(-abs(diff_df$prob_diff)), ]
  
      # 取前 Run_CellChat_ntop_signaling
      if (nrow(diff_df) > Run_CellChat_ntop_signaling) {
        diff_df <- diff_df[1:Run_CellChat_ntop_signaling, ]
      }

      # === 新增：從 interaction 拆出 ligand / receptor ===
      if ("interaction" %in% colnames(diff_df)) {
        parts <- strsplit(diff_df$interaction, "_")
        diff_df$ligand <- sapply(parts, `[`, 1)
        diff_df$receptor <- sapply(parts, `[`, 2)
      }
  
      # 輸出 CSV
      out_csv <- file.path(Run_CellChat_output_path, "multiGroup_significant_LR_by_prob_diff.csv")
      write.csv(diff_df, file = out_csv, row.names = FALSE)
      message("✅ Multi-group significant L-R pairs saved to: ", out_csv)
    } else {
      message("⚠️ No significant L-R pairs found across group comparisons.")
    }
}
  return(cellchat)
}
