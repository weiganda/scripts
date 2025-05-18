#!/usr/bin/env python3

import argparse
import math
import matplotlib.pyplot as plt
import csv

CONVERSION_FACTOR_Debye = 2.541746473
CONVERSION_FACTOR_Coulomb_Meter = 8.478353619813e-30

def calculate_magnitude(x, y, z):
    return math.sqrt(x**2 + y**2 + z**2)

def read_vectors_and_compute_magnitudes(filename):
    magnitudes = []
    with open(filename, 'r') as f:
        for line in f:
            if line.strip():  # skip empty lines
                x_str, y_str, z_str = line.strip().split(',')
                x, y, z = float(x_str), float(y_str), float(z_str)
                magnitude = calculate_magnitude(x, y, z)
                magnitudes.append(magnitude)
    return magnitudes

def write_magnitudes_to_csv(magnitudes, output_filename="dipoleMag.csv"):
    avg_mag = sum(magnitudes) / len(magnitudes)
    avg_converted_debye = avg_mag * CONVERSION_FACTOR_Debye
    avg_converted_coul = avg_mag * CONVERSION_FACTOR_Coulomb_Meter
    with open(output_filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write average
        writer.writerow(["Average", "Debye", "CoulombMeter"])
        writer.writerow([f"{avg_mag:.5f}", f"{avg_converted_debye:.5f}", f"{avg_converted_coul}"])

def main():
    parser = argparse.ArgumentParser(description="Calculate, convert, save, and plot vector magnitudes.")
    parser.add_argument("-f", "--file", required=True, help="Path to the file containing vectors")
    args = parser.parse_args()

    magnitudes = read_vectors_and_compute_magnitudes(args.file)

    write_magnitudes_to_csv(magnitudes, "dipoleMag.csv")

if __name__ == "__main__":
    main()

