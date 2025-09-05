# functions/run_cellchat.R

run_cellchat <- function(seurat_obj,
                         output_dir = "output/cellchat",
                         species = "human",
                         group_by = "Celltype",
                         source_celltype = NULL,
                         target_celltype = NULL,
                         plot_heatmap = TRUE,
                         ntop_pathway = 5,
                         ntop_signaling = 1,
                         pdf_width = 6,
                         pdf_height = 6) {
  library(CellChat)
  library(patchwork)

  if (species == "human") {
    CellChatDB.use <- CellChatDB.human
  } else if (species == "mouse") {
    CellChatDB.use <- CellChatDB.mouse
  } else {
    stop("❌ Unsupported species: ", species)
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  message("Running CellChat on Seurat object ...")
  cellchat <- createCellChat(object = GetAssayData(seurat_obj, slot = "data"),
                             meta = seurat_obj@meta.data,
                             group.by = group_by)
  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

  saveRDS(cellchat, file = file.path(output_dir, "cellchat.rds"))
  message("Output cellchat object saved to: ", output_dir)

  # -----------------------------
  # signaling role heatmap (outgoing & incoming)
  # -----------------------------
  if (plot_heatmap) {
    message("Drawing netAnalysis_signalingRole_heatmap ...")
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

    pdf(file = file.path(output_dir, "signalingRole_heatmap.pdf"),
        width = pdf_width, height = pdf_height)
    print(ht1 + ht2)
    dev.off()
    message("Outgoing and incoming signalingRole_heatmap plot saved: ",
            file.path(output_dir, "signalingRole_heatmap.pdf"))
  }

  # -----------------------------
  # Pathway ranking (netP) + bubble
  # -----------------------------
  top_paths <- character(0)

  if (!is.null(ntop_pathway)) {
    message("Ranking signaling pathways using rankNet ...")

    # Case 1: If provide only source:
    if (!is.null(source_celltype) && is.null(target_celltype)) {
      rn <- CellChat::rankNet(
        object = cellchat,
        slot.name = "netP",
        measure = "weight",
        mode = "single",
        sources.use = source_celltype,
        return.data = TRUE
      )
      top_paths <- head(rn$signaling.contribution$name, ntop_pathway)
      message("Top signaling pathways from ", source_celltype, ": ", paste(top_paths, collapse = ", "))

      if (length(top_paths) > 0) {
        pdf(file = file.path(output_dir, paste0("netVisual_bubble_top", ntop_pathway, "_from_", source_celltype, ".pdf")),
            width = pdf_width,
            height = pdf_height)
        print(CellChat::netVisual_bubble(cellchat,
                                         signaling = top_paths,
                                         sources.use = source_celltype,
                                         targets.use = NULL,
                                         remove.isolate = TRUE))
        dev.off()
      }

    # Case 2: If provide only target:
    } else if (is.null(source_celltype) && !is.null(target_celltype)) {
      rn <- CellChat::rankNet(
        object = cellchat,
        slot.name = "netP",
        measure = "weight",
        mode = "single",
        targets.use = target_celltype,
        return.data = TRUE
      )
      top_paths <- head(rn$signaling.contribution$name, ntop_pathway)
      message("Top signaling pathways to ", target_celltype, ": ", paste(top_paths, collapse = ", "))

      if (length(top_paths) > 0) {
        pdf(file = file.path(output_dir, paste0("netVisual_bubble_top", ntop_pathway, "_to_", target_celltype, ".pdf")),
            width = pdf_width,
            height = pdf_height)
        print(CellChat::netVisual_bubble(cellchat,
                                         signaling = top_paths,
                                         sources.use = NULL,
                                         targets.use = target_celltype,
                                         remove.isolate = TRUE))
        dev.off()
      }

    # Case 3: If provide both source and target:
    } else if (!is.null(source_celltype) && !is.null(target_celltype)) {

      # Check valid communication
      comm_subset <- subsetCommunication(cellchat,
                                         sources.use = source_celltype,
                                         targets.use = target_celltype)
      if (nrow(comm_subset) == 0) {
        warning("⚠️ No inferred communications between ", source_celltype, " -> ", target_celltype, ". Skipping rankNet and netVisual_bubble.")
        top_paths <- character(0)
      } else {
        rn <- CellChat::rankNet(
          object = cellchat,
          slot.name = "netP",
          measure = "weight",
          mode = "single",
          sources.use = source_celltype,
          targets.use = target_celltype,
          return.data = TRUE
        )
        top_paths <- head(rn$signaling.contribution$name, ntop_pathway)
        message("Top signaling pathways from ", source_celltype, " to ", target_celltype, ": ", paste(top_paths, collapse = ", "))

        if (length(top_paths) > 0) {
          pdf(file = file.path(output_dir, paste0("netVisual_bubble_top", ntop_pathway, "_", source_celltype, "_to_", target_celltype, ".pdf")),
              width = pdf_width,
              height = pdf_height)
          print(CellChat::netVisual_bubble(cellchat,
                                           signaling = top_paths,
                                           sources.use = source_celltype,
                                           targets.use = target_celltype,
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
      top_paths <- head(rn$signaling.contribution$name, ntop_pathway)
      message("Overall top signaling pathways: ", paste(top_paths, collapse = ", "))

      if (length(top_paths) > 0) {
        pdf(file = file.path(output_dir, paste0("netVisual_bubble_top", ntop_pathway, "_signaling.pdf")),
            width = pdf_width,
            height = pdf_height)
        print(CellChat::netVisual_bubble(cellchat,
                                         signaling = top_paths,
                                         remove.isolate = TRUE))
        dev.off()
      }
    }

    message("netVisual_bubble plots saved (top ", ntop_pathway, " signaling)")
  }

  # -----------------------------
  # Extract top L-R pairs (net) for each top pathway
  # Only works when user provides both source and target celltype
  # -----------------------------
  if (!is.null(source_celltype) && !is.null(target_celltype) && length(top_paths) > 0) {
    for (path in top_paths) {
    message("Extracting top L-R pairs for pathway: ", path)
      
    # Extract L-R pair 資料
    comm_df <- subsetCommunication(
      cellchat,
      signaling = path,
      sources.use = source_celltype,
      targets.use = target_celltype
    )
      
    if (nrow(comm_df) == 0) {
      message("⚠️ No communication found for pathway ", path)
      next
    }
      
    comm_df <- comm_df[order(comm_df$prob, decreasing = TRUE), ]
    top_lr <- head(comm_df, ntop_signaling)
      
    # save CSV file
    csv_file <- file.path(output_dir, paste0("top", ntop_signaling, "_LR_", path, ".csv"))
    write.csv(top_lr, file = csv_file, row.names = FALSE)
    message("Saved top ", ntop_signaling," L-R pairs for ", path, " to: ", csv_file)
      
    # draw bubble plot
    pdf(file = file.path(output_dir, paste0("netVisual_bubble_top", ntop_signaling, "_", path, ".pdf")),
        width = pdf_width,
        height = pdf_height)
    print(CellChat::netVisual_bubble(
      cellchat,
      signaling = path,
      sources.use = source_celltype,
      targets.use = target_celltype,
      remove.isolate = TRUE
    ))
    dev.off()
    message("Bubble plot saved for top ", ntop_signaling," L-R pairs of ", path)
    }
  }

  message("✅ CellChat complete.")
  return(cellchat)
}
