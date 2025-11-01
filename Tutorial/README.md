In this tutorial, we will demonstrate how to utilize scDock with GSE218563 dataset from NCBI GEO website (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218563). This scRNA dataset contains 16 kidney tissue samples from control and diabetic nephropathy (DN) mice.  

● Data preparation:  
Please manually download those data to your computer. We recommend users rename their files as shown in the example below for 10X Cell Ranger data. For example, each sample directory should contain the following files:  
`barcodes.tsv.gz`  
`features.tsv.gz`  
`matrix.mtx.gz`  

● Meta data preparation (Optional):  
If you are interested in the difference of cellular interaction between control and DN mice, you can provide metadata.txt as we have prepared.  

● Modify config.yaml (read in the above GSE218563_config.yaml file):  
In the tutorial, the scDock workpath is set as `/home/andrew/scDock`. 
Users can modify the paths according to your requirement.  Except the path arguments, users can stay default settings if you are not sure how to modify them.  

Here, we mention some important arguments that users might want to modify:  
1. `Load_QC_max_mito` This argument defines the persent of mitochondrial counts cut off in subset(). The valid range locates from 0 to 1. Default is 0.1. Lower cut off will remove more cells in the further analysis.
2. `Run_Integration_resolution` This argument defines the resoltion in Seurat::FindClusters(). This clustering is processed after data integration. Default is 0.5. Higher resolution might bring more clusters in results.
3. `Run_Integration_dims` This argument defines how many PCs are used in the dimensional reduction analysis. Default is "auto". The default setting will use an proper PC number based on geometric Elbow method.  

