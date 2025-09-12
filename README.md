# scDock
Here, we establish a series analysis pipeline, including scRNA/snRNA sequencing analysis, cell-cell communication and small molecular docking.  
scDock simplify the complicated and preliminary process from single cell / nucleus RNA sequencing analysis to drug discovery.  
Input: scRNA/snRNA data (CellRanger output files [barcodes.tsv.gz, features.tsv.gz, matrix.mtx], .h5 and .txt are available in scDock)  
Output: Small molecular docking results based on the top strongest cellular signalings across cell types (you can choose your interested cell type)  

How to use:
conduct <bash main.sh config.yaml 2>&1 | tee log.txt>

Preparation:
1. [Mandatory] First, please follow the README in our Essential environment/.
2. [Mandatory] Then, please follow the README in our Essential files/.
3. [Optional] If you want to use FDA compounds as drug library, we can follow the README in our FDA parent compounds_202509/.
4. [Optional] If you want to group your input samples, you can follow the README in our Metadata/.
