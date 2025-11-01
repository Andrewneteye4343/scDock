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
2. `Run_Integration_method` This argument defines the method in Seurat::IntegrateLayers(). The valid value are "cca", "harmony" and "rpca". Default is "cca".  
3. `Run_Integration_dims` This argument defines how many PCs will be used in the dimensional reduction analysis. Default is "auto". The default setting will use a proper PC number based on geometric Elbow method. Should be modified only if you intend to run integration.  
4. `Run_Integration_resolution` This argument defines the resoltion in Seurat::FindClusters(). Default is 0.5. Higher resolution might bring more clusters in results. This argument is optional. Should be modified only if you intend to run integration.  
5. `Markers_Annotation_tissue_type` This argument defines the tissue type of your scRNA/snRNA data. Please remember to set this correctly.
6. `Run_CellChat_group_by` This argument defines the group.by in CellChat::createCellChat(). Default is "Celltype". If you are interested in groups from metadata.txt, you can set this argument with "sample_group".  
7. `Run_CellChat_ntop_signaling` Choose how many top signalings (detailed L-R pairs in each pathway) you want to explore. More signalings will increase the time cost for further docking. Defualt is 5.
8. `Vina_Docking_vina_exhaustiveness` This argument defines the exhaustiveness in AudoDock Vina. Default is 8.
9. `Vina_Docking_vina_cpu` This argument defines the cpu in AutoDock Vina. Default is 1.  

● Start tutorial:  
1. If you havn't installed scDock, you can conduct `bash install.sh`. The installation might be finished in minutes.
2. First, you should place "main.sh", "config.yaml", and "functions" directory in your work place (Figure 1). (images/Figure 1.jpg)
 
