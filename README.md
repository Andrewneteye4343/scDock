# scDock
### Overview  
This repository provides an integrated analysis pipeline covering:  
● scRNA/snRNA-seq analysis  
● Cell–cell communication inference  
● Small-molecule docking  

scDock streamlines the complex workflow from single-cell/nucleus RNA sequencing analysis to drug discovery, reducing manual preparation and making the process more efficient.

### Input:
Available scRNA/snRNA format:  
● barcodes.tsv.gz, features.tsv.gz, matrix.mtx  
● .h5  
● .txt  

### Output:  
● Seurat_object.rds files from Seurat analysis steps.  
● CellChat.rds file from CellChat analysis.  
● Protein structure(s) from CellChat analysis.  
● Compound strucutre(s) from cas.txt or fda.txt.  
● Molecular docking results (binding affinity and binding position) based on the strongest intercellular signaling pathways across cell types (you can select the cell types of interest).

### How to use:  
Enter scDock virtual environment and conduct the following command  
`conda activate scDock`  
`bash main.sh config.yaml 2>&1 | tee log.txt` 

### Preparation:
1. [Recommended] Please follow the README.md in our Essential environment.
2. [Recommended] Please follow the README.md in our Essential files.
3. [Optional] If you want to use FDA compounds as compound library in molecular docking, you can follow the README.md in our FDA parent compounds_202509.
4. [Optional] If you want to group your input samples, you can follow the README.md in our Metadata.

### Note:  
We recommend clearing the files in "Vina_Docking_output_path" (defined in config.yaml) before each run, as scDock does not automatically remove the ligand and receptor directories. Otherwise, unnecessary structural files may accumulate.
