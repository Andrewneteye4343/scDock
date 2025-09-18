In this tutorial, we will demonstrate how to utilize scDock with GSE149878 dataset from NCBI GEO website (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149878). This scRNA dataset contains 4 lung tissue samples from a COVID-19 severe patient.  

● Data preparation:  
Please manually download those data to your computer and should get:  
`GSM4516279_C166_filtered_feature_bc_matrix.h5`  
`GSM4516280_C168_filtered_feature_bc_matrix.h5`  
`GSM4516281_C170_filtered_feature_bc_matrix.h5`  
`GSM4516282_C172_filtered_feature_bc_matrix.h5`  

● Modify config.yaml:  
In the tutorial, the scDock workpath is set as `/home/andrew/scDock`. 
Users should adjust the paths according to your requirement.  Except the path arguments, users can stay default settings if you are not sure how to adjust them.     
`work_path: /home/andrew/scDock`  
`Seurat_output_dir: seurat_output`  
`Load_QC_input_files: /home/andrew/scDock/datasets`  
`Load_QC_input_type: h5`  
`Load_QC_min_features: 200`  
`Load_QC_min_cells: 3`  
`Load_QC_names_delim: "_"`  
`Load_QC_max_mito: 0.1`  
`Load_QC_metadata_file: NULL`  
`Load_QC_verbose: TRUE`  
`Normalization_Scale_use_assay: RNA`  
`Normalization_Scale_normalization_method: LogNormalize`  
`Normalization_Scale_CLR_margin: 2`  
`Normalization_Scale_scale_factor: 10000`  
`Normalization_Scale_hvg_method: vst`  
`Normalization_Scale_nVariableFeatures: 2000`  
`Normalization_Scale_split_by: NULL`  
`Normalization_Scale_model_use: linear`  
`Normalization_Scale_scale_max: 10`  
`Normalization_Scale_scale_features: variable`  
`Normalization_Scale_verbose: TRUE`  
`DR_Cluster_pca_features: variable`  
`DR_Cluster_seed: 42`  
`DR_Cluster_dims: 20`  
`DR_Cluster_k_param: 20`  
`DR_Cluster_n_trees: 50`  
`DR_Cluster_resolution: 0.1`  
`DR_Cluster_reduction_method: umap`  
`DR_Cluster_reduction_assay: NULL`  
`DR_Cluster_clustering_algorithm: 1`  
`DR_Cluster_verbose: TRUE`  
`Markers_Annotation_output_path: /home/andrew/scDock`  
`Markers_Annotation_use_assay: RNA`  
`Markers_Annotation_group_by: seurat_clusters`  
`Markers_Annotation_logfc_threshold: 0.1`  
`Markers_Annotation_min_pct:  0.01`  
`Markers_Annotation_test_use: wilcox`  
`Markers_Annotation_only_pos: TRUE`  
`Markers_Annotation_verbose: TRUE`  
`Markers_Annotation_top_n: 100`  
`Markers_Annotation_tissue_type: lung`  
`Markers_Annotation_label_column: Celltype`  
`DR_Plot_output_path: /home/andrew/scDock`  
`DR_Plot_reduction_method: umap`  
`DR_Plot_group_by: Celltype`  
`DR_Plot_pt_size: 0.5`  
`DR_Plot_label: TRUE`  
`DR_Plot_label_size: 4`  
`DR_Plot_label_color: black`  
``
