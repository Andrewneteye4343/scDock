# functions/prepare_receptor.py
#!/usr/bin/env python3
# coding: utf-8

import sys, os
from MolKit import Read
import AutoDockTools.Utilities24.prepare_receptor4 as pr4
import numpy as np

def compute_grid_box(coords, padding=5.0):
    """計算最大 grid box"""
    min_coords = np.min(coords, axis=0) - padding
    max_coords = np.max(coords, axis=0) + padding
    center = (min_coords + max_coords) / 2
    size = max_coords - min_coords
    return center, size

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: prepare_receptor.py input.pdb output_dir")
        sys.exit(1)

    in_pdb, outdir = sys.argv[1], sys.argv[2]
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    receptor = Read(in_pdb)[0]
    out_pdbqt = os.path.join(outdir, os.path.basename(in_pdb).replace(".pdb", "_prepared.pdbqt"))

    # Using AD4ReceptorPreparation from AutoDockTools
    prep = pr4.AD4ReceptorPreparation(
    receptor,
    repairs='checkhydrogens',
    charges_to_add='gasteiger',
    cleanup='nphs_lps_waters_nonstdres',
    outputfilename=out_pdbqt
)


    print(f"[prepare_receptor] Generated: {out_pdbqt}")

    # Get the maximum grid box
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
