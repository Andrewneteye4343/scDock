We provide a list of compounds (FDA_202509_compounds.xlsx) obtained from the FDA Drugs Database
 (version: September 2, 2025; accessed: September 4, 2025).

The .xlsx file was generated using the script "Get CAS from compound name.py", which retrieves the PubChem CID and CAS number for each compound name in the FDA list, along with the corresponding parent compound’s name, PubChem CID, and CAS number.

The spreadsheet contains six columns:  
● Input Name  
● Input CID  
● Input CAS  
● Parent Name  
● Parent CID  
● Parent CAS  

Since molecular docking requires the major active component of a drug, only the parent compound is used. The CAS numbers of the parent compounds are recorded in "fda.txt".

Using the FDA compound list as a drug library
To enable this feature:
1. Place the provided fda.txt file in your computer.  
2. Set the argument Vina_Docking_use_fda = true and Vina_Docking_fda_txt to /path/to/fda.txt.  
Once configured, scDock will automatically download and preprocess the corresponding FDA parent compounds for subsequent molecular docking.
