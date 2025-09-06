# functions/Run_CellChat.R
Run_CellChat <- function(seurat_obj,
                         Run_CellChat_output_dir = "path/to/dir/",
                         Run_CellChat_species = "human",
                         Run_CellChat_group_by = "Celltype",
                         Run_CellChat_source_celltype = NULL,
                         Run_CellChat_target_celltype = NULL,
                         Run_CellChat_plot_heatmap = TRUE,
                         Run_CellChat_ntop_pathway = 5,
                         Run_CellChat_ntop_signaling = 1,
                         Run_CellChat_pdf_width = 6,
                         Run_CellChat_pdf_height = 6) {
  library(CellChat)
  library(patchwork)

  if (Run_CellChat_species == "human") {
    CellChatDB.use <- CellChatDB.human
  } else if (Run_CellChat_species == "mouse") {
    CellChatDB.use <- CellChatDB.mouse
  } else {
    stop("❌ Unsupported species: ", Run_CellChat_species)
  }

  dir.create(Run_CellChat_output_dir, showWarnings = FALSE, recursive = TRUE)

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

  saveRDS(cellchat, file = file.path(Run_CellChat_output_dir, "cellchat.rds"))
  message("Output cellchat object saved to: ", Run_CellChat_output_dir)

  # -----------------------------
  # signaling role heatmap (outgoing & incoming)
  # -----------------------------
  if (Run_CellChat_plot_heatmap) {
    message("Drawing netAnalysis_signalingRole_heatmap ...")
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

    pdf(file = file.path(Run_CellChat_output_dir, "signalingRole_heatmap.pdf"),
        width = Run_CellChat_pdf_width, height = Run_CellChat_pdf_height)
    print(ht1 + ht2)
    dev.off()
    message("Outgoing and incoming signalingRole_heatmap plot saved: ",
            file.path(Run_CellChat_output_dir, "signalingRole_heatmap.pdf"))
  }

  # -----------------------------
  # Pathway ranking (netP) + bubble
  # -----------------------------
  top_paths <- character(0)

  if (!is.null(Run_CellChat_ntop_pathway)) {
    message("Ranking signaling pathways using rankNet ...")

    # Case 1: If provide only source:
    if (!is.null(Run_CellChat_source_celltype) && is.null(Run_CellChat_target_celltype)) {
      rn <- CellChat::rankNet(
        object = cellchat,
        slot.name = "netP",
        measure = "weight",
        mode = "single",
        sources.use = Run_CellChat_source_celltype,
        return.data = TRUE
      )
      top_paths <- head(rn$signaling.contribution$name, Run_CellChat_ntop_pathway)
      message("Top signaling pathways from ", Run_CellChat_source_celltype, ": ", paste(top_paths, collapse = ", "))

      if (length(top_paths) > 0) {
        pdf(file = file.path(Run_CellChat_output_dir, paste0("netVisual_bubble_top", Run_CellChat_ntop_pathway, "_from_", Run_CellChat_source_celltype, ".pdf")),
            width = Run_CellChat_pdf_width,
            height = Run_CellChat_pdf_height)
        print(CellChat::netVisual_bubble(cellchat,
                                         signaling = top_paths,
                                         sources.use = Run_CellChat_source_celltype,
                                         targets.use = NULL,
                                         remove.isolate = TRUE))
        dev.off()
      }

    # Case 2: If provide only target:
    } else if (is.null(Run_CellChat_source_celltype) && !is.null(Run_CellChat_target_celltype)) {
      rn <- CellChat::rankNet(
        object = cellchat,
        slot.name = "netP",
        measure = "weight",
        mode = "single",
        targets.use = Run_CellChat_target_celltype,
        return.data = TRUE
      )
      top_paths <- head(rn$signaling.contribution$name, Run_CellChat_ntop_pathway)
      message("Top signaling pathways to ", Run_CellChat_target_celltype, ": ", paste(top_paths, collapse = ", "))

      if (length(top_paths) > 0) {
        pdf(file = file.path(Run_CellChat_output_dir, paste0("netVisual_bubble_top", Run_CellChat_ntop_pathway, "_to_", Run_CellChat_target_celltype, ".pdf")),
            width = Run_CellChat_pdf_width,
            height = Run_CellChat_pdf_height)
        print(CellChat::netVisual_bubble(cellchat,
                                         signaling = top_paths,
                                         sources.use = NULL,
                                         targets.use = Run_CellChat_target_celltype,
                                         remove.isolate = TRUE))
        dev.off()
      }

    # Case 3: If provide both source and target:
    } else if (!is.null(Run_CellChat_source_celltype) && !is.null(Run_CellChat_target_celltype)) {

      # Check valid communication
      comm_subset <- subsetCommunication(cellchat,
                                         sources.use = Run_CellChat_source_celltype,
                                         targets.use = Run_CellChat_target_celltype)
      if (nrow(comm_subset) == 0) {
        warning("⚠️ No inferred communications between ", Run_CellChat_source_celltype, " -> ", Run_CellChat_target_celltype, ". Skipping rankNet and netVisual_bubble.")
        top_paths <- character(0)
      } else {
        rn <- CellChat::rankNet(
          object = cellchat,
          slot.name = "netP",
          measure = "weight",
          mode = "single",
          sources.use = Run_CellChat_source_celltype,
          targets.use = Run_CellChat_target_celltype,
          return.data = TRUE
        )
        top_paths <- head(rn$signaling.contribution$name, Run_CellChat_ntop_pathway)
        message("Top signaling pathways from ", Run_CellChat_source_celltype, " to ", Run_CellChat_target_celltype, ": ", paste(top_paths, collapse = ", "))

        if (length(top_paths) > 0) {
          pdf(file = file.path(Run_CellChat_output_dir, paste0("netVisual_bubble_top", Run_CellChat_ntop_pathway, "_", Run_CellChat_source_celltype, "_to_", Run_CellChat_target_celltype, ".pdf")),
              width = Run_CellChat_pdf_width,
              height = Run_CellChat_pdf_height)
          print(CellChat::netVisual_bubble(cellchat,
                                           signaling = top_paths,
                                           sources.use = Run_CellChat_source_celltype,
                                           targets.use = Run_CellChat_target_celltype,
                                           remove.isolate = TRUE))
          dev.off()
        }
      }

    # Default: If provide no source or target: Global top pathways
    } else {
      rn <- CellChat::rankNet(
        object = cellchat,
        slot.name = "netP",
        measure = "weight",
        mode = "single",
        return.data = TRUE
      )
      top_paths <- head(rn$signaling.contribution$name, Run_CellChat_ntop_pathway)
      message("Overall top signaling pathways: ", paste(top_paths, collapse = ", "))

      if (length(top_paths) > 0) {
        pdf(file = file.path(Run_CellChat_output_dir, paste0("netVisual_bubble_top", Run_CellChat_ntop_pathway, "_signaling.pdf")),
            width = Run_CellChat_pdf_width,
            height = Run_CellChat_pdf_height)
        print(CellChat::netVisual_bubble(cellchat,
                                         signaling = top_paths,
                                         remove.isolate = TRUE))
        dev.off()
      }
    }

    message("netVisual_bubble plots saved (top ", Run_CellChat_ntop_pathway, " signaling)")
  }

  # -----------------------------
  # Extract top L-R pairs (net) for each top pathway
  # Only works when user provides both source and target celltype
  # -----------------------------
  if (!is.null(Run_CellChat_source_celltype) && !is.null(Run_CellChat_target_celltype) && length(top_paths) > 0) {
    for (path in top_paths) {
    message("Extracting top L-R pairs for pathway: ", path)
      
    # Extract L-R pair 資料
    comm_df <- subsetCommunication(
      cellchat,
      signaling = path,
      sources.use = Run_CellChat_source_celltype,
      targets.use = Run_CellChat_target_celltype
    )
      
    if (nrow(comm_df) == 0) {
      message("⚠️ No communication found for pathway ", path)
      next
    }
      
    comm_df <- comm_df[order(comm_df$prob, decreasing = TRUE), ]
    top_lr <- head(comm_df, Run_CellChat_ntop_signaling)
      
    # save CSV file
    csv_file <- file.path(Run_CellChat_output_dir, paste0("top", Run_CellChat_ntop_signaling, "_LR_", path, ".csv"))
    write.csv(top_lr, file = csv_file, row.names = FALSE)
    message("Saved top ", Run_CellChat_ntop_signaling," L-R pairs for ", path, " to: ", csv_file)
      
    # draw bubble plot
    pdf(file = file.path(Run_CellChat_output_dir, paste0("netVisual_bubble_top", Run_CellChat_ntop_signaling, "_", path, ".pdf")),
        width = Run_CellChat_pdf_width,
        height = Run_CellChat_pdf_height)
    print(CellChat::netVisual_bubble(
      cellchat,
      signaling = path,
      sources.use = Run_CellChat_source_celltype,
      targets.use = Run_CellChat_target_celltype,
      remove.isolate = TRUE
    ))
    dev.off()
    message("Bubble plot saved for top ", Run_CellChat_ntop_signaling," L-R pairs of ", path)
    }
  }
  return(cellchat)
}
