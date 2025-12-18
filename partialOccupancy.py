#!/usr/bin/env python3

import numpy as np
import random
import argparse
from collections import defaultdict
from pymatgen.core import Structure, Element
from pymatgen.io.cif import CifParser
from pymatgen.io.vasp import Poscar

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Convert CIF to POSCAR with probabilistic atom selection.")
parser.add_argument('-f', '--file', type=str, required=True, help='Input CIF file')
args = parser.parse_args()

# Load CIF file with occupancy tolerance
cif_file = args.file
parser = CifParser(cif_file, occupancy_tolerance=1.05)
structures = parser.parse_structures(primitive=False)
if not structures:
    raise ValueError(f"No valid structure parsed from {cif_file}")
structure_raw = structures[0]

# Extract raw data: label, fractional coords, occupancy, and element
raw_data = []
for site in structure_raw.sites:
    coords = site.frac_coords
    occu_dict = site.species
    for elem, occ in occu_dict.items():
        raw_data.append((site.species_string, coords, occ, elem.symbol))

# Group atoms by rounded coordinates
coord_dict = defaultdict(list)
for label, coords, occ, elem in raw_data:
    key = tuple(np.round(coords, decimals=5))
    coord_dict[key].append((elem, occ))

for i in range(10):
    # Probabilistic atom selection
    coords = []
    species = []
    for pos, occupants in coord_dict.items():
        total_occ = sum(o[1] for o in occupants)
        if total_occ == 0:
            continue
        normed_probs = [o[1]/total_occ for o in occupants]
        selected = random.choices([o[0] for o in occupants], weights=normed_probs, k=1)[0]
        species.append(Element(selected))
        coords.append(pos)

    # Build structure
    structure = Structure(structure_raw.lattice, species, coords)

    # Write to POSCAR
    poscar = Poscar(structure)
    poscar.write_file(f"POSCAR_PartialOcc_{i}.vasp")
    print(f"Wrote POSCAR to POSCAR_PartialOcc{i}.vasp")

    # Determine unique species in POSCAR order
    unique_species = poscar.site_symbols

    # Write matching POTCAR
    potcar_filename = f"POTCAR_PartialOcc_{i}.vasp"
    with open(potcar_filename, 'w') as f:
        for el in unique_species:
            f.write(f"TITEL = DUMMY DUMMY {el}_pv\n")
    print(f"Wrote POTCAR to {potcar_filename}")
