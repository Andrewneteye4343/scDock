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
● Core scRNA-seq analysis outputs include dimensional reduction plots with cell-type annotations, marker gene information, detailed annotation results, and Seurat_object.rds files from intermediate Seurat analysis steps.  
● Intercellular communication inference outputs include incoming and outgoing signaling, global and group-specific communication networks, and the CellChat.rds file.  
● Molecular docking results include the binding affinities and binding poses between compounds and proteins involved in the top intercellular signaling pathways across cell types, as well as in the most differentially activated pathways between groups. Users can also specify the cell types and groups of interest. The output includes the corresponding protein structure(s) and compound structure(s) used for docking.

### How to use:  
Enter scDock virtual environment and conduct the following command  
`conda activate scDock`  
`bash main.sh config.yaml` 

### Preparation:
1. [Recommended] Please follow the README.md in our [Essential environment](./Essential%20environment).
2. [Recommended] Please follow the README.md in our Essential files [Essential files](./Essential%20files).
3. [Optional] If you want to use FDA compounds as compound library in molecular docking, you can follow the README.md in our FDA parent compounds_202509.
4. [Optional] If you want to group your input samples, you can follow the README.md in our Metadata.

### Tutorial:
We provide a [Tutorial](./Tutorial) that demonstrates how to use scDock. By following the instructions, you can perform a comparative analysis using kidney samples from control and diabetic nephropathy (DN) mice.  

### Note:  
We recommend clearing the files in argument "Vina_Docking_output_path" (in config.yaml) before each run, as scDock does not automatically remove the ligand and receptor directories. Otherwise, unnecessary structural files may accumulate.
