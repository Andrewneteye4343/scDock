In this tutorial, we will demonstrate how to utilize scDock with GSE149878 dataset from NCBI GEO website (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149878). This scRNA dataset contains 4 lung tissue samples from a COVID-19 severe patient.  

● Data preparation:  
Please manually download those data to your computer and should get:  
`GSM4516279_C166_filtered_feature_bc_matrix.h5`  
`GSM4516280_C168_filtered_feature_bc_matrix.h5`  
`GSM4516281_C170_filtered_feature_bc_matrix.h5`  
`GSM4516282_C172_filtered_feature_bc_matrix.h5`  

● Modify config.yaml (read in the above GSE149878_config.yaml file):  
In the tutorial, the scDock workpath is set as `/home/andrew/scDock`. 
Users should adjust the paths according to your requirement.  Except the path arguments, users can stay default settings if you are not sure how to adjust them.  
Here, we mention some important arguments that users might want to adjust:  
1. `Load_QC_max_mito` This argument defines the persent of mitochondrial counts cut off in subset(). The valid range locates from 0 to 1. Default is 0.1. Lower cut off might remove too much more cells in the further analysis.
2. `Load_QC_metadata_file` This path defines where to find metadata for grouping your samples. Should not be the same path with Load_QC_input_files when you input txt as Load_QC_input_type. Only use it if you want to group your samples.
3. `DR_Cluster_resolution` This argument defines the resoltion in Seurat::FindClusters(). Default is 0.1. The higher resolution might bring more clusters in results.
