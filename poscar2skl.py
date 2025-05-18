#!/usr/bin/env python3

import argparse
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# Parse argument for POSCAR file path
parser = argparse.ArgumentParser(description="Convert POSCAR to olcao.skl.new format")
parser.add_argument("-f", "--file", type=str, required=True, help="Path to POSCAR file")
args = parser.parse_args()

# Read structure from POSCAR
structure = Poscar.from_file(args.file).structure

# Get lattice parameters
lattice = structure.lattice
a, b, c = lattice.abc
alpha, beta, gamma = lattice.angles

# Get atomic data
cart_coords = structure.cart_coords
elements = structure.species
num_atoms = len(cart_coords)

# Get space group symbol
spacegroup = SpacegroupAnalyzer(structure, symprec=1e-3).get_space_group_symbol()

# Write to file
with open("olcao.skl.new", "w") as f:
    f.write("title\n")
    f.write("POSCAR\n")
    f.write("end\n")
    f.write("cell\n")
    f.write(f"{a:.6f} {b:.6f} {c:.6f} {alpha:.2f} {beta:.2f} {gamma:.2f}\n")
    f.write(f"cart {num_atoms}\n")
    for elem, coord in zip(elements, cart_coords):
        f.write(f"{elem} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")
    f.write(f"space 1_a\n")
    f.write("supercell 1 1 1\n")
    f.write("full\n")

print("Wrote olcao.skl.new")
