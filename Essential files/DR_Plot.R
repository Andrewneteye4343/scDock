# functions/DR_Plot.R
DR_Plot <- function(seurat_obj,
                               DR_Plot_output_path = "/path/to/dir/",
                               DR_Plot_reduction_method = "umap",
                               DR_Plot_group_by = "Celltype",
                               DR_Plot_pt_size = 0.5, 
                               DR_Plot_label = TRUE,
                               DR_Plot_label_size = 4,
                               DR_Plot_label_color = "black",
                               DR_Plot_label_box = FALSE,
                               DR_Plot_alpha = 1,
                               DR_Plot_shuffle = FALSE,
                               DR_Plot_raster = NULL,
                               DR_Plot_width = 6,
                               DR_Plot_height = 6) {
  library(Seurat)
  library(ggplot2)
  
  dim_plot <- DimPlot(seurat_obj,
                      reduction = DR_Plot_reduction_method,
                      group.by = DR_Plot_group_by,
                      pt.size = DR_Plot_pt_size,
                      label = DR_Plot_label,
                      label.size = DR_Plot_label_size,
                      label.color = DR_Plot_label_color,
                      label.box = DR_Plot_label_box,
                      alpha = DR_Plot_alpha,
                      shuffle = DR_Plot_shuffle,
                      raster = DR_Plot_raster
                      ) +
    ggtitle(paste0(toupper(DR_Plot_reduction_method), " by ", DR_Plot_group_by))
  
  ggsave(DR_Plot_output_path, plot = dim_plot, width = DR_Plot_width, height = DR_Plot_height)
  message(toupper(DR_Plot_reduction_method), " plot saved to: ", DR_Plot_output_path)
}
