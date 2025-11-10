#!/usr/bin/env python3
# coding: utf-8
# functions/prepare_receptor.py

import sys, os, yaml
from MolKit import Read
import AutoDockTools.Utilities24.prepare_receptor4 as pr4
import numpy as np

# Function to measure the maximum grid box
def compute_grid_box(coords, padding=5.0):
    min_coords = np.min(coords, axis=0) - padding
    max_coords = np.max(coords, axis=0) + padding
    center = (min_coords + max_coords) / 2
    size = max_coords - min_coords
    return center, size

# Function to prepare receptor structure and generate grid box settings
def process_receptor(in_pdb, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    receptor = Read(in_pdb)[0]
    out_pdbqt = os.path.join(outdir, os.path.basename(in_pdb).replace(".pdb", "_prepared.pdbqt"))

    # Prepare receptor structure with AD4ReceptorPreparation from AutoDockTools
    prep = pr4.AD4ReceptorPreparation(
        receptor,
        repairs='checkhydrogens',
        charges_to_add='gasteiger',
        cleanup='nphs_lps_waters_nonstdres',
        outputfilename=out_pdbqt
    )
    print(f"[prepare_receptor] Generated: {out_pdbqt}")

    # Measuring grid box
    coords = receptor.allAtoms.coords
    center, size = compute_grid_box(coords)
    grid_file = out_pdbqt.replace("_prepared.pdbqt", "_prepared_grid.txt")
    with open(grid_file, "w") as f:
        f.write(f"center_x = {center[0]:.3f}\n")
        f.write(f"center_y = {center[1]:.3f}\n")
        f.write(f"center_z = {center[2]:.3f}\n")
        f.write(f"size_x = {size[0]:.3f}\n")
        f.write(f"size_y = {size[1]:.3f}\n")
        f.write(f"size_z = {size[2]:.3f}\n")
    print(f"[prepare_receptor] Grid info saved: {grid_file}")

# Generate their own grid box for all structures in the directory
def evaluate_prepared_receptors(receptor_dir):
    files = [f for f in os.listdir(receptor_dir) if f.endswith("_prepared.pdbqt")]
    if not files:
        print(f"[prepare_receptor] No *_prepared.pdbqt found in {receptor_dir}")
        return

    for f in files:
        pdbqt_path = os.path.join(receptor_dir, f)
        print(f"[prepare_receptor] Processing grid for: {f}")

        # Read coordinates of receptor structure
        receptor = Read(pdbqt_path)[0]
        coords = receptor.allAtoms.coords
        center, size = compute_grid_box(coords)

        grid_file = pdbqt_path.replace("_prepared.pdbqt", "_prepared_grid.txt")
        with open(grid_file, "w") as f_out:
            f_out.write(f"center_x = {center[0]:.3f}\n")
            f_out.write(f"center_y = {center[1]:.3f}\n")
            f_out.write(f"center_z = {center[2]:.3f}\n")
            f_out.write(f"size_x = {size[0]:.3f}\n")
            f_out.write(f"size_y = {size[1]:.3f}\n")
            f_out.write(f"size_z = {size[2]:.3f}\n")

        print(f"[prepare_receptor] Grid info saved: {grid_file}")

if __name__ == "__main__":
    # Check whether config.yaml exists
    config_path = "config.yaml"
    receptor_dir = None

    if os.path.exists(config_path):
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)

        # Read the Vina_Docking_docking_receptor_dir argument where stores customized receptor structure(s) 
        vina_cfg = config.get("Vina_Docking", {})
        receptor_dir = vina_cfg.get("Vina_Docking_docking_receptor_dir", None)

    # If Vina_Docking_docking_receptor_dir is provided, scDock will measure grid box for the customized receptor structure(s)
    if receptor_dir and str(receptor_dir).lower() not in ["null", "none", ""]:
        receptor_dir = os.path.expanduser(receptor_dir)
        print(f"[prepare_receptor] Evaluating receptors in directory: {receptor_dir}")
        evaluate_prepared_receptors(receptor_dir)
        sys.exit(0)

    # If not, use the downloaded receptor structure(s)
    if len(sys.argv) < 3:
        print("Usage: prepare_receptor.py input.pdb output_dir")
        sys.exit(1)

    in_pdb, outdir = sys.argv[1], sys.argv[2]
    process_receptor(in_pdb, outdir)