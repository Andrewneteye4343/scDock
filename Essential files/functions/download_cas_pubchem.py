# functions/download_cas_pubchem.py
#!/usr/bin/env python3
# coding: utf-8

import sys
import os
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess

def download_sdf(cas_number, out_sdf):
    """Download 3D SDF from PubChem for given CAS number"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/SDF?record_type=3d"
    r = requests.get(url)
    if r.status_code == 200 and len(r.content) > 0:
        with open(out_sdf, "wb") as f:
            f.write(r.content)
        return True
    else:
        print(f"[WARN] CAS {cas_number} not found in PubChem.")
        return False

def sdf_to_pdbqt(sdf_file, out_dir):
    """Convert SDF to PDBQT using RDKit + OpenBabel"""
    mol_name = os.path.splitext(os.path.basename(sdf_file))[0]
    pdb_file = os.path.join(out_dir, f"{mol_name}.pdb")
    pdbqt_file = os.path.join(out_dir, f"{mol_name}.pdbqt")

    # Generate 3D conformer using RDKit
    mol = Chem.MolFromMolFile(sdf_file, removeHs=False)
    if mol is None:
        raise ValueError(f"RDKit failed to read {sdf_file}")

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol)
    Chem.MolToPDBFile(mol, pdb_file)

    # Convert into .pdbqt using OpenBabel
    subprocess.run(["obabel", pdb_file, "-O", pdbqt_file, "-xh"], check=True)

    return pdbqt_file

def main():
    if len(sys.argv) != 3:
        print("Usage: download_cas_pubchem.py <cas_txt_file> <output_dir>")
        sys.exit(1)
    
    cas_txt_file = sys.argv[1]
    output_dir = sys.argv[2]
    os.makedirs(output_dir, exist_ok=True)

    with open(cas_txt_file) as f:
        cas_numbers = [line.strip() for line in f if line.strip()]

    for cas in cas_numbers:
        sdf_file = os.path.join(output_dir, f"{cas}.sdf")
        if not os.path.exists(sdf_file):
            success = download_sdf(cas, sdf_file)
            if not success:
                continue
        try:
            pdbqt_file = sdf_to_pdbqt(sdf_file, output_dir)
            print(f"[INFO] Generated PDBQT: {pdbqt_file}")
        except subprocess.CalledProcessError:
            print(f"[ERROR] Failed to convert {cas} to PDBQT.")

if __name__ == "__main__":
    main()

