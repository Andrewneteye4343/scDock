# functions/DR_Plot.R
DR_Plot <- function(seurat_obj,
                               DR_Plot_output_path,
                               DR_Cluster_reduction_method,
                               DR_Plot_group_by,
                               DR_Plot_pt_size, 
                               DR_Plot_label,
                               DR_Plot_label_size,
                               DR_Plot_label_color,
                               DR_Plot_label_box,
                               DR_Plot_alpha,
                               DR_Plot_shuffle,
                               DR_Plot_raster,
                               DR_Plot_width,
                               DR_Plot_height) {
  suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})
  dim_plot <- DimPlot(seurat_obj,
                      reduction = DR_Cluster_reduction_method,
                      group.by = DR_Plot_group_by,
                      pt.size = DR_Plot_pt_size,
                      label = DR_Plot_label,
                      label.size = DR_Plot_label_size,
                      label.color = DR_Plot_label_color,
                      label.box = DR_Plot_label_box,
                      alpha = DR_Plot_alpha,
                      shuffle = DR_Plot_shuffle,
                      raster = DR_Plot_raster
                      )
  if (DR_Cluster_reduction_method == "umap"){
      dim_plot <- dim_plot + labs(title = "", x = "UMAP_1", y = "UMAP_2")
  } else {
    dim_plot <- dim_plot + labs(title = "", x = "TSNE_1", y = "TSNE_2")
  }
  
  DRplot = file.path(DR_Plot_output_path, "DimPlot.png")
  ggsave(DRplot, plot = dim_plot, width = DR_Plot_width, height = DR_Plot_height)
  message("[DR_Plot] ", toupper(DR_Cluster_reduction_method), " plot saved to: ", DRplot)
}
