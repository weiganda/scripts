#!/bin/env python

# PROGRAM:  seriesData
# PURPOSE:  This program will take a time series of skeleton files
#           (currently Dipole specific) run olcao -dimo, copy the 
#           dipole moment to a csv file, move the .skl, scfV, and
#           iter files to a sub folder and run the next .skl file.
#           this will automate the calculation of the GHz frequency
#           dielectric calculation as well as rewrite over old
#           files to hopefully save data space.
#
# AUTHOR:   Alysse Weigand
# LAST MODIFIED: April 16, 2024
# USAGE: seriesData

import shutil
import os
import subprocess
import csv

def copy_file(src, dst):
    try:
        shutil.copy(src, dst)
        print(f"File '{src}' copied successfully to '{dst}'")
    except FileNotFoundError:
        print(f"Error: '{src}' not found.")
    except PermissionError:
        print(f"Error: Permission denied to copy '{src}'.")

def execute_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Command '{command}' executed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command '{command}': {e}")

def main():
    destination_file = "olcao.skl"
    makeinput_command = "makeinput -kp 1 1 1"
    olcao_command = "olcao -scfdimo"
    remove_files = "rm olcao.skl *.dat *.out"
    remove_folders = "rm -r i*"

    # Initial file counter for olcao.skl files
    file_counter = 1

        # If olcao.skl folder doesnt exit, exit
        #if not os.path.exists(source_file):
            #break

    if os.path.exists(destination_file):

            #copy_file(source_file, destination_file)

        # Execute makeinput
        execute_command(makeinput_command)

        # Execute olcao -scfdimo
        execute_command(olcao_command)

        # Read the cell volume of the original .skl file and save to file
        if file_counter == 1:
                
            with open("gs_main-fb.out", 'r') as file:
                lines = file.readlines()

            # Initialize a variable to hold the real cell volume
            real_cell_volume = None

            # Loop through each line to find the real cell volume
            for line in lines:
                if 'The real  cell volume is:' in line:
                    # Split the line to extract the numerical value
                    parts = line.split(':')
                    if len(parts) > 1:
                        real_cell_volume = parts[1].strip()

            # Check if the real cell volume was found
            if real_cell_volume is not None:
                # Open the output file in write mode and write the value
                with open('vol', 'w') as output_file:
                    output_file.write(real_cell_volume)
            else:
                print("The real cell volume was not found in the fort.20 file.")

           
        # Read data from gs_dimo.plot and write to dipole_series
        with open("gs_dimo-fb.t.plot", "r") as gs_dimo_file:
            lines = gs_dimo_file.readlines()
            # Extract dipole moment values
            dipole_values = []
            for line in lines:
                if "Dipole Moment" in line:
                    parts = line.split(":")
                    value = float(parts[1].strip())
                    dipole_values.append(value)
 
            with open("dipole_series.csv", "a", newline='') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow(dipole_values)
                    
        # Remove previous data files
        execute_command(remove_files)
        execute_command(remove_folders)
        # Increment file counter
        file_counter += 1



if __name__ == "__main__":
    main()

