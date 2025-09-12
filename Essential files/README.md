### To utilize scDock, users need to download the essential files from Github. Please follow the instructions:
#### You will need to adjust the path arguments (where belong to your computer) in config.yaml later.
1. Please create and move to your work directory for scDock analysis. We will use `/home/path/to/scDock` as an example.
2. Please place "main.sh" and "config.yaml" at `/home/path/to/scDock`
3. Please place the ligand_reference.csv and receptor_reference.csv in your computer. For example, `/home/path/to/scDock`
4. Please create a directory named "functions" at `/home/path/to/scDock`. You should get a directory `/home/path/to/scDock/functions`
5. Please place all the below scripts at `/home/to/your/path/scDock/functions`
##### R scripts
● Load_QC.R  
● Normalization_Scale.R  
● DR_Cluster.R  
● Markers_Annotation.R  
● DR_Plot.R  
● Run_CellChat.R  
● Vina_Docking.R  

##### Python scripts
● download_alphafold.py  
● download_swissmodel.py  
● download_cas_pubchem.py  
● prepare_receptor.py  

#### config.yaml:
Users only need to adjust the arguments in config.yaml according to their requirements.
If you are unsure about the meaning of certain argument, we recommend using the default settings.

Each argument in config.yaml is documented with its usage and original function. For more details, you can also refer to the documentation provided on the respective official websites.

#### Decide the compound library for molecular docking
We provides three options for ligand usage: FDA compounds, CAS number and your own structures.  
● For FDA compounds, you need to set Vina_Docking_use_fda = true and Vina_Docking_fda_txt = path/to/fda.txt. In this option, scDock will download and pre-process the parent compounds of FDA-approved compounds (version: September 2025). The fda.txt can be download from Github.  
● For CAS number, you need to set Vina_Docking_cas_txt_file = path/to/cas.txt. In this option, you can provide the CAS number of your interested compounds in the cas.txt file. The CAS number can be serached at PubChem webste (https://pubchem.ncbi.nlm.nih.gov/).  
● For your own structures, you need to set Vina_Docking_docking_ligand_dir = path/to/YourStucture. You can provide your compound structure.pdbqt file by using this argument.




#### For Markers_Annotation_tissue_type argument, you must input the below tissue name:
adipose tissue, bladder, blood, bone, bone marrow, brain, breast, embryo, eye, gastrointestinal tract, heart, kidney, liver, lung, mammary gland, muscle, other, ovary, pancreas, placenta, prostate, skin, spleen, stomach, testis, thymus, tooth, uterus

#### For Run_CellChat_source_celltype and Run_CellChat_target_celltype arguments, you must input NULL or the below cell name:
#### adipose tissue
Adipocyte, Adipose-derived stem cell, B cell, Basophil, Brown fat cell, Dendritic cell, Endothelial cell, Hematopoietic cell, Luminal epithelial cell, Macrophage, Mammary epithelial cell, Mast cell, Monocyte, Natural killer cell, Neuron, Pericyte, Platelet, T cell, T memory cell

#### bladder
Basal epithelial cell, Dendritic cell, Endothelial cell, Macrophage, Smooth muscle cell, Smooth muscle progenitor cell, Umbrella cell, Ureteric epithelium cell, Urine-derived stem cell, Urothelial cell

#### blood 
B cell, Basophil, CD14 Monocyte, CD16 Monocyte, CD4 Central Memory T cell, CD4 Cytotoxic T cell, CD4 Effector Memory T cell, CD4 Naive T cell, CD4 Proliferating T cell, CD4 T cell, CD56-bright natural killer cell, CD56-dim natural killer cell, CD8 Central Memory T cell, CD8 Effector Memory T cell, CD8 Naive T cell, CD8 Proliferating T cell, CD8 T cell, Dendritic cell, Double-negative T cell, Embryonic stem cell, Endothelial cell, Endothelial progenitor cell, Endothelial stem cell, Eosinophil, Epithelial cell, Erythroblast, Erythroid cell, Erythroid precursor cell, Gamma delta T cell, Granulocyte, Hematopoietic stem cell, Immature myeloid cell, Innate lymphoid cell, Intermediate B cell, Leukocyte, Lymphoid cell, Megakaryocyte, Memory B cell, Mesenchymal stem cell, Monocyte, Mucosal-associated invariant T cell, Myeloid cell, Myeloid dendritic cell, Myeloid-derived suppressor cell, Naive B cell, Natural killer cell, Neutrophil, Plasma cell, Plasmablast, Plasmacytoid dendritic cell, Platelet, Regulatory T cell, Reticulocyte, T cell, Thymic emigrant cell

#### bone
Chondrocyte, Dendritic cell, Eosinophil, Erythroid cell, Granulocyte, Hematopoietic stem cell, Meniscus-derived progenitor cell, Mesenchymal progenitor cell, Mesenchymal stem cell, Monocyte, Myeloid progenitor cell, Neutrophil, Osteoblast, Osteoclast, Osteoclast precursor cell, Osteocyte, Periosteum-derived progenitor cell

#### bone marrow 
B cell, Basophil, Blastema cell, Common lymphoid progenitor cell, Common myeloid progenitor, Dendritic cell, Dendritic progenitor cell, Endothelial cell, Endothelial progenitor cell, Eosinophil, Erythroid cell, Erythroid megakaryocyte progenitor cell, Erythroid precursor cell, Granulocyte monocyte progenitor cell, Hematopoietic stem cell, Immature dendritic cell, Innate lymphoid cell, Langerhans cell, Lymphoid-primed multipotent progenitor cell, Macrophage, Mast cell, Mature T cell, Megakaryocyte, Megakaryocyte progenitor cell, Mesenchymal stem cell, Monocyte, Monocyte derived dendritic cell, Mucosal-associated invariant T cell, Myeloid cell, Myeloid dendritic cell, Myeloid progenitor cell, Natural killer cell, Neutrophil, Osteoclast precursor cell, Plasma cell, Plasmacytoid dendritic cell, Platelet, Precursor plasmacytoid dendritic cell, T cell

#### brain 
Adrenergic neurons, Anterior pituitary gland cell, Astrocyte, B cell, Basket cell, Bergmann glial cell, Cajal-Retzius cell, CCK basket cell, Chandelier cell, Cholinergic neuron, Choroid plexus cell, Dopaminergic neuron, Endothelial cell, Ependymal cell, GABAergic neuron, Ganglion cell, Glutamatergic neuron, Hypothalamic ependymal cell, Immature neuron, Interneuron, Interneuron-selective cell, Lepotomeningeal cell, Macrophage, Martinotti cell, Meningeal cell, Microglia cell, Motor neuron, Mural cell, Neural progenitor cell, Neural stem cell, Neuroblast, Neuroendocrine cell, Neuron, Noradrenergic neuron, Olfactory ensheathing glia, Olfactory sensory neuron, Oligodendrocyte, Oligodendrocyte precursor cell, Pan-gabaergic, Pericyte, Pinealocyte, Purkinje neuron, Pyramidal cell, Radial glia cell, Rhombic lip cell, Satellite glial cell, Schwann cell, Serotonergic neuron, Smooth muscle cell, T cell, Tanycyte, Trigeminal neuron, Type IA spiral ganglion neuron, Type IB spiral ganglion neuron, Type IC spiral ganglion neuron, Type II spiral ganglion neuron

#### breast
B cell, Basal epithelial cell, Dendritic cell, Endothelial cell, Eosinophil, Epithelial cell, Fibroblast, Hematopoietic stem cell, Luminal epithelial cell, Luminal progenitor cell, Macrophage, Mast cell, Mesenchymal stem cell, Myoepithelial cell, Natural killer cell, Neutrophil, Pericyte, T cell

#### embryo
Arterial cell,Blastomere,Cardiomyocyte,Dorsal otocyst,Ectoderm cell,Embryonic stem cell,Endocardial cell,Endocrine cell,Endoderm cell,Endothelial cell,Epiblast cell,Germ cell,Hemangioblast,
Hematopoietic stem cell,Macrophage,Mesenchymal stem cell,Mesoderm cell,Myeloblast,Neural crest cell,Neural stem cell,Neural tube cell,Neuron,Neutrophil,Oocyte,Primitive endoderm cell,
Skeletal muscle cell,Trophectoderm cell,Trophoblast cell,Unrestricted somatic stem cell,Venous cell,Ventral otocyst

#### eye
Bipolar cell,Endothelial cell,Epithelial cell,Epithelial stem cell,Erythroid cell,Ganglion cell,Hematopoietic stem cell,Lymphocyte,Macrophage,Mesenchymal cell,Mesenchymal stem cell,Muller cell,
Myoepithelial cell,Photoreceptor cell,Progenitor cell

#### gastrointestinal tract
Ciliated epithelial cell,Colonic stem cell,Crypt cell,Dendritic cell,Enteric glia cell,Enteric neuron,Enterochromaffin cell,Enterocyte,Enterocyte progenitor cell,Enteroendocrine cell,
Enteroendocrine precursor cell,Foveolar cell,Gastric chief cell,Goblet cell,Goblet progenitor cell,Intestinal stem cell,Macrophage,Mast cell,Microfold cell,Paneth cell,Parietal cell,
S cell,Secretory progenitor cell,T cell,Tuft cell,Tuft progenitor cell

#### heart
Aorta valve interstitial cell,Cardiocyte,Cardiomyocyte,Cardiovascular progenitor cell,Endocardial cell,Endothelial cell,Erythroid cell,Fibroblast,Macrophage,Myofibroblast,Pericyte,
Progenitor cell,Smooth muscle cell

#### kidney
B cell,Basophil,Connecting tubule epithelial cell,Cortical cell,Dendritic cell,Descending vasa recta endothelial cell,Distal convoluted tubule cell,Endothelial cell,Epithelial cell,
Erythroblast,Erythroid precursor cell,Fibroblast,Glomerular capillary endothelial cell,Granular cell,Infiltrated mononuclear cell,Inner medulla collecting duct epithelial cell,
Intercalated cell,Juxtaglomerular cell,Kidney progenitor cell,Loop of Henle cell,Loop of Henle cortical thick ascending limb epithelial cell,
Loop of Henle medullary thick ascending limb epithelial cell,Loop of Henle thin descending limb epithelial cell,Lymphatic endothelial cell,M2 Macrophage,Macrophage,Macula densa,
Mast cell,Medullary fibroblast,Mesangial cell,Monocyte,Mononuclear phagocyte,Natural killer cell,Nephron epithelial cell,Neutrophil,Papillary tip epithelial cell,
Parietal epithelial cell,Pericyte,Peritubular capilary endothelial cell,Plasma cell,Podocyte,Principal cell,Proximal tubule brush border cell,Proximal tubule cell,
Proximal tubule epithelial cell,Renal alpha-intercalated cell,Schwann cell,T cell

#### liver
B cell,Cholangiocyte,Dendritic cell,Endothelial cell,Epithelial cell,Hematopoietic cell,Hepatoblast,Hepatocyte,Kupffer cell,Liver bud hepatic cell,Mesenchymal cell,
Mesenchymal stem cell,Monocyte,Myofibroblast,Natural killer cell,Neutrophil,Progenitor cell,Stellate cell,Stem cell,T cell

#### lung
Adventitial fibroblast,B cell,Basal cell,Basophil,Bronchial epithelial cell,Brush Cell,Capillary aerocyte,Ciliated cell,Club cell,Dendritic cell,Endothelial cell,Eosinophil,
Epithelial cell,Fibroblast,Fibromyocyte,Goblet cell,Granulocyte,Lipofibroblast,Lymphatic endothelial cell,Lymphoid cell,M1 Macrophage,M2 Macrophage,Macrophage,Mast cell,
Mesenchymal progenitor cell,Mesothelial cell,Monocyte,Mucous cell,Myeloid cell,Myeloid leukocyte,Myofibroblast,Neuroendocrine cell,Neutrophil,Pericyte,Plasma cell,
Pulmonary alveolar type I cell,Pulmonary alveolar type II cell,Pulmonary Ionocyte cell,Secretory cell,Serous cell,Smooth muscle cell,Stem cell,T cell

#### mammary gland
B cell,Basal cell,Endothelial cell,Epithelial cell,Hormone sensing differentiated cell,Hormone sensing progenitor cell,Luminal cell,Luminal epithelial cell,Luminal progenitor cell,
Macrophage,Mast cell,Muscle cell,Myoepithelial cell,Progenitor cell,Stem cell,T cell

#### muscle
Adventitial cell,B cell,Endothelial cell,Erythroblast,Fibroblast,Granulocyte monocyte progenitor cell,Inflammatory cell,Macrophage,Mesenchymal progenitor cell,Mesenchymal stem cell,
Muscle-derived cell,Myoblast,Myocyte,Myocyte progenitor cell,Myoepithelial cell,Myofibroblast,Myogenic endothelial cell,Neutrophil,Pericyte,Progenitor cell,Satellite cell,Schwann cell,
Smooth muscle cell,T cell,Tenocyte

#### other
Chromaffin cell,Epithelial cell,Follicular cell,Glomus cell,Hair cell,Invasive spongiotrophoblast,Langerhans cell,Mesothelial cell,Nucleus pulposus cell,Parathyroid chief cell,
Salivary mucous cell,Thymocyte,Type III taste bud cell,Vomeronasal sensory neuron

#### ovary
Cumulus cell,Endothelial cell,Endothelium cell,Epithelial cell,Germ cell,Granulosa cell,Luteal cell,Macrophage,Mesenchymal cell,Oocyte,Pluripotent stem cell,Stem cell,
Theca interna cell,Thecal cell

#### pancreas
Acinar cell,Alpha cell,B cell,Beta cell,Delta cell,Dendritic cell,Ductal cell,Ductal stem cell,Endocrine cell,Endothelial cell,Epithelial cell,Epsilon cell,Gamma cell,Glial cell,
Granulocyte,Macrophage,Mast cell,Pancreatic progenitor cell,Polypeptide cell,Schwann cell,Stellate cell,T cell

#### placenta
Basophil,Decidual stem cell,Dendritic cell,Endodermal cell,Endothelial cell,Erythroid cell,Hofbauer cell,Invasive spongiotrophoblast,Labyrinthine trophoblast,Macrophage,
Megakaryocyte progenitor cell,Mesenchymal stem cell,Monocyte,Natural killer cell,Pericyte,Spongiotrophoblast,Stem cell,Trophoblast cell,Basal cell,Epithelial cell,
Luminal cell,Progenitor cell,Prostate epithelial cell,Prostate stem cell,Stem cell

#### skin
B cell,Basal cell,Dendritic cell,Endothelial cell,Epidermal stem cell,Epithelial cell,Keratinocyte,Langerhans cell,Macrophage,Melanocyte,Merkel cell,Mesenchymal stem cell,
Neural crest stem cell,Sebocyte,Stem cell,T cell,Trichocyte

#### spleen
B cell,Dendritic cell,Erythroblast,Granulocyte,Lymphocyte,Macrophage,Monocyte,Natural killer cell,Neutrophil,Plasma cell,Plasmacytoid dendritic cell,T cell

#### stomach
Dendritic cell,Gastric stem cell,Multipotent stem cell,Parietal cell,Parietal progenitor cell,Pit cell,Pit progenitor cell,Smooth muscle cell,Tuft cell

#### testis
Leydig cell,Macrophage,Peritubular myoid cell,Pre sertoli cell,Sertoli cell,Spermatids,Spermatocyte,Spermatogonia,Spermatogonial stem cell,Spermatogonium,Spermatozoa

#### thymus
Cortical thymic epithelial cell,Dendritic cell,Epithelial cell,Medullary thymic epithelial cell,T cell,Thymocyte

#### tooth
Alveolar osteocyte,Dental follicle cell,Dental pulp cell,Dental pulp stem cell,Endothelial cell,Epithelial cell,Glial cell,Macrophage,Mesenchymal stem cell,
Natural killer cell,Periodontal ligament stem cell,Perivascular cell,T cell

#### uterus
B cell,Decidual cell,Endometrial stem cell,Keratinocyte,Macrophage,Monocyte,Natural killer cell,Neutrophil,Smooth muscle cell,Stem cell
