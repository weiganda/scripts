#!/usr/bin/env python3

import os

# Get the current directory
current_directory = os.path.basename(os.getcwd())


# Define the content to be written to the file
content1 = """

CalculationMode = gs
UnitsOutput = eV_Angstrom

XYZCoordinates = 'structure.xyz'

BoxShape = sphere
Radius = 5.0*angstrom
Spacing = 0.18*angstrom
"""

# Write the content to the file named 'inp'
with open("inp", "w") as file:
    file.write(content1)

print("File 'inp' has been created successfully.")



# Define the content to be written to the file
content2 = f"""

#!/bin/bash
#SBATCH -p general
#SBATCH -J {current_directory} 
#SBATCH -o {current_directory}.o%J
#SBATCH -e {current_directory}.e%J
#SBATCH -N 1
#SBATCH -t 03:00:00
#SBATCH --mem=10G
#

mpirun octopus

"""

# Write the content to the file named 'inp'
with open("slurm", "w") as file:
    file.write(content2)

print("File 'slurm' has been created successfully.")


