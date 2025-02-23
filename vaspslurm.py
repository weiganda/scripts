#!/usr/bin/env python3

import os

def create_slurm_file():
    # Get the current working directory
    current_directory = os.getcwd()

    # Extract the directory name to use as the job name
    directory_name = os.path.basename(current_directory)

    slurm_content = f"""#!/bin/bash
#SBATCH -p general
#SBATCH -J "{directory_name}"
#SBATCH -o "{directory_name}.o%J"
#SBATCH -e "{directory_name}.e%J"
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -t 47:00:00
#SBATCH --mem=30G
#

cd {current_directory}

export MV2_ENABLE_AFFINITY=0
mpirun -n 6 vasp_gam
"""

    # Write the content to a SLURM file
    slurm_filename = f"slurm"
    with open(slurm_filename, 'w') as file:
        file.write(slurm_content)

    print(f"SLURM file created successfully in {current_directory}!")

# Run the function
create_slurm_file()
