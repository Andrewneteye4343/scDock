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
2. `Load_QC_metadata_file` This path defines where to find metadata for grouping your samples (Optional).  
3. `Run_Integration_method` This argument defines the method in Seurat::IntegrateLayers(). The valid value are "cca", "harmony" and "rpca". Default is "cca" (Optional).  
4. `Run_Integration_dims` This argument defines how many PCs will be used in the dimensional reduction analysis. Default is "auto". The default setting will use a proper PC number based on geometric Elbow method. Should be modified only if you intend to run integration (Optional).  
5. `Run_Integration_resolution` This argument defines the resoltion in Seurat::FindClusters(). Default is 0.5. Higher resolution might bring more clusters in results. This argument is optional. Should be modified only if you intend to run integration (Optional).  
6. `Markers_Annotation_tissue_type` This argument defines the tissue type of your scRNA/snRNA data. Please remember to set this correctly.
7. `Run_CellChat_group_by` This argument defines the group.by in CellChat::createCellChat(). Default is "Celltype". If you are interested in groups from metadata.txt, you can set this argument with "sample_group".  
8. `Run_CellChat_ntop_signaling` Choose how many top signalings (detailed L-R pairs in each pathway) you want to explore. More signalings will increase the time cost for further docking. Defualt is 5.
9. `Vina_Docking_vina_exhaustiveness` This argument defines the exhaustiveness in AudoDock Vina. Default is 8.
10. `Vina_Docking_vina_cpu` This argument defines the cpu in AutoDock Vina. Default is 1.  

● Start tutorial:  
1. If you havn't installed scDock, you can conduct `bash install.sh`. The installation might be finished in minutes.
2. Before analysis, you should place "main.sh", "config.yaml", and "functions" directory in your workpath.  
<img width="119" height="91" alt="Figure 1" src="https://github.com/user-attachments/assets/6596c326-d52b-4651-aed3-173fef6a198f" />  

3. Condut `bash main.sh config.yaml` to start scDock analysis.  First, scDock will check the validity of arguments in config.yaml. We provide default setting for some conditions if the setting goes wrong. After the configuration, scDock will start load your data.  
<img width="885" height="203" alt="Figure 2" src="https://github.com/user-attachments/assets/168130be-c6c7-49c7-84f6-64a4b86e322e" />

4. During QC process, the cell and gene number before and after qulaity control will be displayed.  
<img width="463" height="303" alt="Figure 3" src="https://github.com/user-attachments/assets/5d2a64e7-aa93-4dc7-a7dd-00e76e319f7f" />

5. During dimensional reduction, if users set "auto" to find a proper PCs number. You will get the information below and output a Elbow plot.pdf.
  <img width="891" height="81" alt="Figure 4" src="https://github.com/user-attachments/assets/f881a352-f7ea-4e94-aa34-1249a129cfe3" />  

6. If users set Run_Integration_run_integration = true in config.yaml, you will see the information below.
<img width="532" height="65" alt="Figure 5" src="https://github.com/user-attachments/assets/a7718b3a-5158-423e-843b-39d148cab5ec" />  
<img width="931" height="416" alt="Figure 6" src="https://github.com/user-attachments/assets/8f37f278-d12c-416c-8e9c-2be97190e48f" />  

7. 
