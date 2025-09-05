# functions/plot_dim_reduction.R

plot_dim_reduction <- function(seurat_obj,
                               output_path,
                               reduction_method = "umap",
                               group_by = "Celltype",
                               pt_size = 0.5, 
                               label = TRUE,
                               label_size = 4,
                               label_color = "black",
                               label_box = FALSE,
                               alpha = alpha,
                               shuffle = FALSE,
                               raster = NULL,
                               width = 6,
                               height = 6) {
  library(Seurat)
  library(ggplot2)
  
  message("Starting visualization with Dimplot...")
  dim_plot <- DimPlot(seurat_obj,
                      reduction = reduction_method,
                      group.by = group_by,
                      pt.size = pt_size,
                      label = label,
                      label.size = label_size,
                      label.color = label_color,
                      label.box = label_box,
                      alpha = alpha,
                      shuffle = shuffle,
                      raster = raster
                      ) +
    ggtitle(paste0(toupper(reduction_method), " by ", group_by))
  
  ggsave(output_path, plot = dim_plot, width = width, height = height)
  message(toupper(reduction_method), " plot saved to: ", output_path)
  message("✅ Visualization complete.")
}
