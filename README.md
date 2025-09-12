# scDock
#### Overview  
This repository provides an integrated analysis pipeline covering:  
⚪ scRNA/snRNA-seq analysis  
⚪ Cell–cell communication inference  
⚪ Small-molecule docking  
scDock streamlines the complex workflow from single-cell/nucleus RNA sequencing analysis to drug discovery, reducing manual preparation and making the process more efficient.

#### Input:
scRNA/snRNA data (CellRanger output files [barcodes.tsv.gz, features.tsv.gz, matrix.mtx], .h5 and .txt are available in scDock)  
#### Output:
Small molecular docking results based on the top strongest cellular signalings across cell types (you can choose your interested cell type)  

#### How to use:  
`bash main.sh config.yaml 2>&1 | tee log.txt` 

#### Preparation:
1. [Mandatory] First, please follow the README in our Essential environment.
2. [Mandatory] Then, please follow the README in our Essential files.
3. [Optional] If you want to use FDA compounds as drug library, we can follow the README in our FDA parent compounds_202509.
4. [Optional] If you want to group your input samples, you can follow the README in our Metadata.

#### Note:  
We recommend clearing the files in "Vina_Docking_output_path" (defined in config.yaml) before each run, as scDock does not automatically remove the ligand and receptor directories. Otherwise, unnecessary structural files may accumulate.
