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
3. `Run_Integration_method` This argument defines the method in Seurat::IntegrateLayers(). The valid values are "cca", "harmony" and "rpca". Default is "cca". Should be modified only if you intend to run integration (Optional).  
4. `Run_Integration_dims` This argument defines how many PCs will be used in the dimensional reduction analysis. Default is "auto". The default setting will use a proper PC number based on geometric Elbow method. Should be modified only if you intend to run integration (Optional).  
5. `Run_Integration_resolution` This argument defines the resoltion in Seurat::FindClusters(). Default is 0.5. Higher resolution might bring more clusters in results. Should be modified only if you intend to run integration (Optional).  
6. `Markers_Annotation_tissue_type` This argument defines the tissue type of your scRNA/snRNA data. Please remember to set this correctly.
7. `Run_CellChat_group_by` This argument defines the group.by in CellChat::createCellChat(). Default is "Celltype". If you are interested in groups from metadata.txt, you can set this argument with "sample_group".  
8. `Run_CellChat_ntop_signaling` Choose how many top signalings (detailed L-R pairs in each pathway) you want to explore. More signalings will increase the time cost for further molecular docking. Defualt is 5.
9. `Vina_Docking_vina_exhaustiveness` This argument defines the exhaustiveness in AudoDock Vina. The higher value will increase the exhaustiveness of the global search, but also increase the time cost. Default is 8.
10. `Vina_Docking_vina_cpu` This argument defines the how many cpu(s) will be used in AutoDock Vina. Default is 1.  

● Start tutorial:  
1. If you haven’t installed scDock, run `bash install.sh`. The installation should complete within a few minutes.  
2. Before starting the analysis, make sure to place main.sh, config.yaml, and the functions directory in your working directory.  <img width="119" height="91" alt="Figure 1" src="https://github.com/user-attachments/assets/6596c326-d52b-4651-aed3-173fef6a198f" />  

3. Please make sure you are working in the scDock virtual environment. If not, activate it by running `conda activate scDock` first. Then, execute bash main.sh config.yaml to start the scDock analysis.  <img width="885" height="203" alt="Figure 2" src="https://github.com/user-attachments/assets/168130be-c6c7-49c7-84f6-64a4b86e322e" />

4. During the QC process, the numbers of cells and genes before and after quality control will be displayed.  <img width="463" height="303" alt="Figure 3" src="https://github.com/user-attachments/assets/5d2a64e7-aa93-4dc7-a7dd-00e76e319f7f" />

5. During dimensional reduction, if the user sets "auto" to automatically determine the optimal number of PCs, scDock will display the following information and generate an Elbow plot (ElbowPlot_Integration_SeuratProject.pdf).  <img width="891" height="81" alt="Figure 4" src="https://github.com/user-attachments/assets/f881a352-f7ea-4e94-aa34-1249a129cfe3" />  

6. If users set Run_Integration_run_integration = true in config.yaml, the following information may be displayed.  <img width="532" height="65" alt="Figure 5" src="https://github.com/user-attachments/assets/a7718b3a-5158-423e-843b-39d148cab5ec" />  <img width="931" height="416" alt="Figure 6" src="https://github.com/user-attachments/assets/8f37f278-d12c-416c-8e9c-2be97190e48f" />  

7. After marker identification and cell annotation are completed, the corresponding results will be saved. If the data originates from mouse samples, gene names will be converted to their human homologs before cell annotation.  <img width="957" height="164" alt="Figure 7" src="https://github.com/user-attachments/assets/7e9b6996-43ed-4294-973f-a3944adc7cc5" />

8. If users set Run_CellChat_group_by = "sample_group", the interaction strength and probability between groups will be compared during the CellChat analysis. The results — including incoming and outgoing heatmaps, bubble plots (with different groups as the max group), and detailed probability differences in a CSV file — will be saved as:
netVisual_bubble_maxGroup_Control.pdf, netVisual_bubble_maxGroup_DN.pdf, signalingRole_heatmap_incoming.pdf, signalingRole_heatmap_outgoing.pdf, and multiGroup_significant_LR_by_prob_diff.csv.  <img width="1891" height="144" alt="Figure 8" src="https://github.com/user-attachments/assets/583ca4a4-40fb-4239-a6f2-fa7aa6008985" />  <img width="1896" height="142" alt="Figure 9" src="https://github.com/user-attachments/assets/2d1230d4-9607-48ef-a2bc-85b3d648e6e3" />  <img width="1224" height="346" alt="Figure 10" src="https://github.com/user-attachments/assets/f16737cc-afbb-4a43-8309-cda7298ef2d5" />  

9. Before molecular docking, scDock retrieves the PDB IDs of ligands and receptors involved in the top cellular signaling interactions between groups. It then downloads the corresponding protein structures from the SWISS-MODEL template library and performs preprocessing for docking preparation. Afterward, scDock downloads and preprocesses the compounds specified by the user. The final prepared structures are saved as .pdbqt files.  <img width="1172" height="391" alt="Figure 11" src="https://github.com/user-attachments/assets/62b1c0a9-1564-4f4f-9a2b-3bee1dc5604e" />

10. During molecular docking, scDock will display the current docking participators. When the docking is finished, scDock will automatically rank the binding affinity of compounds in the single receptor results and save it (AutoDockVina_score.csv).  <img width="1169" height="862" alt="Figure 12" src="https://github.com/user-attachments/assets/30d4c2d6-4bf9-4c1d-be15-7c9a9ae6a4f6" />  

11. In the end, you should obtain results similar to those shown in the figure.
<img width="374" height="393" alt="Figure 13" src="https://github.com/user-attachments/assets/3d43956f-f1da-4d8a-b2e6-ef6d5e454f89" />

● Result explanation:  
1. ElbowPlot_SeuratProject.pdf – Visualizes the ranking of principal components based on the percentage of variance explained in the Elbow plot.
The optimal number of PCs should be chosen at the point where an “elbow” is observed.
2. ElbowPlot_Integration_SeuratProject.pdf – Visualizes the ranking of principal components based on the percentage of variance explained in the Elbow plot during data integration. The optimal number of PCs should be chosen at the point where an “elbow” is observed.
3. DimPlot.png – Visualizes the distribution of cells with annotated labels.  <img width="3000" height="2400" alt="DimPlot" src="https://github.com/user-attachments/assets/475cffb8-b98e-4c35-95a8-04a38f1331d3" />
4. Markers.csv – Records the detailed information of differentially expressed genes (DEGs), including p-value, adjusted p-value and average log2 fold change, and additional related details in each cluster. 
5. Annotation.csv – Records detailed annotation results calculated by scMayoMap (https://github.com/chloelulu/scMayoMap).  
6. signalingRole_heatmap_incoming.pdf & signalingRole_heatmap_outgoing.pdf – Visualize incoming and outgoing signaling activities among groups.
7. netVisual_bubble_maxGroup_DN.pdf – Visualizes signaling pathways that are more active in a specific group (e.g., DN) compared to others.
8. multiGroup_significant_LR_by_prob_diff.csv – Records detailed information on top signaling pathways across all groups or within a specific group, including: interaction name, source (signaling-providing cell type), target (signaling-receiving cell type), group name, prob1 (probability in group 1), prob2 (probability in group 2), prob_diff (probability difference), higher_group (group with higher probability), and additional related details.  
9. ligands_with_PDB.csv & receptors_with_PDB.csv – Record ligand and receptor proteins involved in top signaling pathways across or within specific groups. Information includes protein name and their PDB model matched from Uniprot-PDB-mapper (https://github.com/iriziotis/Uniprot-PDB-mapper) (accessed August 8, 2025).
10. ligand_from_PDB_LR_pairs/ – Contains the protein structures and detailed docking results (LigandName_AutoDockVina_result_structure.pdbqt & AutoDockVina_score.csv) of ligand proteins listed in ligands_with_PDB.csv. In AutoDockVina_score.csv, scDock ranks the compounds by their binding affinity from lowest to highest.
11. receptor_from_PDB_LR_pairs/ – Contains the protein structures and detailed docking results (LigandName_AutoDockVina_result_structure.pdbqt & AutoDockVina_score.csv) of receptor proteins listed in receptors_with_PDB.csv.
12. ligand_structures_from_CAStxt_for_AutoDockVina/ – Stores the prepared compound structures for molecular docking provided in cas.txt.
[netVisual_bubble_maxGroup_DN.pdf](https://github.com/user-attachments/files/23565965/netVisual_bubble_maxGroup_DN.pdf)
