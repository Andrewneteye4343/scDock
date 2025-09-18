In this tutorial, we will demonstrate how to utilize scDock with GSE149878 dataset from NCBI GEO website (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149878). This scRNA dataset contains 4 lung tissue samples from a COVID-19 severe patient.  

● Data preparation:  
Please manually download those data to your computer and should get:  
`GSM4516279_C166_filtered_feature_bc_matrix.h5`  
`GSM4516280_C168_filtered_feature_bc_matrix.h5`  
`GSM4516281_C170_filtered_feature_bc_matrix.h5`  
`GSM4516282_C172_filtered_feature_bc_matrix.h5`  

● Modify config.yaml:  
In the tutorial, the scDock workpath is set as `/home/andrew/scDock`. 
Users should adjust the paths according to your requirement.    
Here are the demonstratation for config.yaml:  
`work_path: /home/andrew/scDock` # This is the workpath where you conduct `bash main.sh config.yaml`.  
`Seurat_output_dir: seurat_output` # This will be the output directory under the work_path.  
`Load_QC_input_files: /home/andrew/scDock/datasets` # Users should place your scRNA in this directory.  
`Load_QC_input_type: h5` # Users should set the correct file format type.  
`Load_QC_min_features: 200` # Set the min.features for Seurat::CreateSeuratObject().  
`Load_QC_min_cells: 3` # Set the min.cells for Seurat::CreateSeuratObject().  
`Load_QC_names_delim: "_" # Set the names.delim for Seurat::CreateSeuratObject().  
`Load_QC_max_mito: 0.1` # Set the QC of mitochondrial counts percent.  
`Load_QC_metadata_file: NULL` # Whether to use meta data. Should be set path/to/metadata.txt to read the meta data (Please refer README.md in MetaData).
