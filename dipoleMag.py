#!/usr/bin/env python3

"""
Created By: Alysse Weigand
Last Updated: 12/20/25

Dipole Magnitude Calculator

Purpose:
This script calculates the magnitude of dipole moment vectors provided in a CSV file 
and converts the magnitudes into multiple units: atomic units (ea0), Debye, and Coulomb-meters.
It also writes the average dipole magnitude to an output CSV file for further analysis.

Student-Friendly Notes:
- Dipole moment vectors have components (x, y, z). The magnitude is computed as:
      |μ| = sqrt(x^2 + y^2 + z^2)
- Conversion factors:
      1 atomic unit (ea0) → 2.541746473 Debye
      1 atomic unit (ea0) → 8.478353619813e-30 Coulomb-meters
- Input format:
      CSV file with each line as x,y,z (no header required)
- Output:
      dipoleMag.csv containing average magnitude in ea0, Debye, and Coulomb-meters
- This script is intended for post-processing OLCAO electronic structure calculations
  where dipole moment vectors are generated at each step.

Usage Example:
    python dipole_magnitude.py -f dipole_vectors.csv
"""

import argparse
import math
import matplotlib.pyplot as plt
import csv

# Conversion factors from atomic units to Debye and Coulomb-meters
CONVERSION_FACTOR_Debye = 2.541746473
CONVERSION_FACTOR_Coulomb_Meter = 8.478353619813e-30

# -----------------------------
# Function to calculate vector magnitude
# -----------------------------
def calculate_magnitude(x, y, z):
    """Compute |μ| = sqrt(x^2 + y^2 + z^2)"""
    return math.sqrt(x**2 + y**2 + z**2)

# -----------------------------
# Function to read CSV vectors and compute magnitudes
# -----------------------------
def read_vectors_and_compute_magnitudes(filename):
    """
    Reads a CSV file containing x,y,z vectors (one vector per line)
    and returns a list of their magnitudes in atomic units (ea0).
    """
    magnitudes = []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip():  # skip empty lines
                x_str, y_str, z_str = line.strip().split(',')
                x, y, z = float(x_str), float(y_str), float(z_str)
                magnitude = calculate_magnitude(x, y, z)
                magnitudes.append(magnitude)
    return magnitudes

# -----------------------------
# Function to write magnitudes to CSV
# -----------------------------
def write_magnitudes_to_csv(magnitudes, output_filename="dipoleMag.csv"):
    """
    Writes the average dipole magnitude in ea0, Debye, and Coulomb-meters to a CSV file.
    """
    avg_mag = sum(magnitudes) / len(magnitudes)
    avg_converted_debye = avg_mag * CONVERSION_FACTOR_Debye
    avg_converted_coul = avg_mag * CONVERSION_FACTOR_Coulomb_Meter

    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Average (ea0)", "Debye", "Coulomb-Meter"])
        writer.writerow([f"{avg_mag:.5f}", f"{avg_converted_debye:.5f}", f"{avg_converted_coul:.5e}"])

# -----------------------------
# Main program
# -----------------------------
def main():
    parser = argparse.ArgumentParser(description="Calculate, convert, and save dipole vector magnitudes.")
    parser.add_argument("-f", "--file", required=True, help="Path to the CSV file containing dipole vectors (x,y,z)")
    args = parser.parse_args()

    # Compute magnitudes from input vectors
    magnitudes = read_vectors_and_compute_magnitudes(args.file)

    # Save average magnitude in multiple units
    write_magnitudes_to_csv(magnitudes, "dipoleMag.csv")

if __name__ == "__main__":
    main()

