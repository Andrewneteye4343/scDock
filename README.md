# scDock
### Overview  
This repository provides an integrated analysis pipeline covering:  
● scRNA/snRNA-seq analysis  
● Cell–cell communication inference  
● Small-molecule molecular docking  

scDock streamlines the complex workflow from single-cell/nucleus RNA sequencing analysis to drug discovery, reducing manual preparation and making the process more efficient. With a user-friendly design, researchers only need to modify a single configuration file to run the entire workflow—without switching programming languages or platforms—greatly simplifying operation and enhancing reproducibility.  

### Input:
Available scRNA/snRNA format:  
● 10X Cell Ranger output (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz)  
● .h5  
● .txt (Column: Cell; Row: Feature)  

### Output:  
● Core scRNA-seq analysis (Dimensional reduction plot with cell annotation, marker genes, detailed annotation results, and Seurat_object.rds files of intermeidated Seurat analysis steps)  
● Intercellular communication inference (Incoming and outgoing signaling, global communication / group-associated communication, and CellChat.rds)  
● Molecular docking results (binding affinity and position) based on the top intercellular signaling pathways across cell types and the most differentially activated pathways between groups (you can also select the cell types and groups of interest), protein structure(s), and compound structure(s).  

### How to use:  
Enter scDock virtual environment and conduct the following command  
`conda activate scDock`  
`bash main.sh config.yaml` 

### Preparation:
1. [Recommended] Please follow the README.md in our Essential environment.
2. [Recommended] Please follow the README.md in our Essential files.
3. [Optional] If you want to use FDA compounds as compound library in molecular docking, you can follow the README.md in our FDA parent compounds_202509.
4. [Optional] If you want to group your input samples, you can follow the README.md in our Metadata.

### Note:  
We recommend clearing the files in argument "Vina_Docking_output_path" (in config.yaml) before each run, as scDock does not automatically remove the ligand and receptor directories. Otherwise, unnecessary structural files may accumulate.
