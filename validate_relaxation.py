#!/usr/bin/env python3

"""
Created by: Alysse Weigand
Last Updated: 12/20/2025

Development Notes:
- Conceptual design, purpose, and validation by Alysse Weigand.
- All scientific reasoning, method choices, and interpretation of results by Alysse Weigand.
- Code implementation and structure assisted by ChatGPT.

Purpose:
This script is used to analyze VASP relaxation runs and determine if the system
has converged in energy, forces, and volume. It reads the OUTCAR file, extracts
total energies, maximum forces, RMS forces, and volume, and provides:

- Plots of energy and force evolution
- Summary of convergence status

Usage:
python validate_relaxation.py
"""

import re
import matplotlib.pyplot as plt
import math

filename = "OUTCAR"

total_energies = []
max_forces = []
rms_forces = []
volumes = []

# --- Parse OUTCAR ---
with open(filename, "r") as f:
    for line in f:

        # Total energy
        if "free energy" in line and "TOTEN" in line:
            match = re.search(r"=\s*([-\d\.Ee]+)", line)
            if match:
                total_energies.append(float(match.group(1)))

        # Max force from 'd Force' lines
        if line.strip().startswith("d Force"):
            # Example line: d Force = 0.0000000E+00[ 0.000E+00, 0.000E+00]
            force_match = re.findall(r"([-+]?\d+\.\d+E[-+]\d+)", line)
            if force_match:
                # Convert string to float
                fx, fy, fz = [float(f) for f in force_match[:3]]
                magnitude = math.sqrt(fx**2 + fy**2 + fz**2)
                max_forces.append(magnitude)
                rms_forces.append(math.sqrt((fx**2 + fy**2 + fz**2)/3))

        # Volume
        if "volume of cell" in line.lower():
            match = re.search(r"volume of cell\s*[:=]\s*([-\d\.Ee]+)", line, re.I)
            if match:
                volumes.append(float(match.group(1)))

# --- Steps ---
steps_energy = list(range(1, len(total_energies) + 1))
steps_forces = list(range(1, len(max_forces) + 1))
steps_volume = list(range(1, len(volumes) + 1))

# --- Convergence thresholds ---
energy_tol = 1e-5      # eV, EDIFF from INCAR
force_tol = 1e-3       # eV/Å, EDIFFG from INCAR
volume_tol = 1e-4      # relative change

# --- Compute final values ---
final_energy = total_energies[-1] if total_energies else None
delta_energy = abs(total_energies[-1] - total_energies[-2]) if len(total_energies) > 1 else None
final_max_force = max(max_forces) if max_forces else None
final_rms_force = rms_forces[-1] if rms_forces else None
final_volume = volumes[-1] if volumes else None
rel_delta_volume = (volumes[-1] - volumes[0])/volumes[0] if len(volumes) > 1 else None

# --- Convergence checks ---
energy_converged = delta_energy is not None and delta_energy < energy_tol
force_converged = final_max_force is not None and final_max_force < force_tol
volume_converged = rel_delta_volume is not None and abs(rel_delta_volume) < volume_tol

# --- Print summary ---
print("\n--- Relaxation Summary ---")
print(f"Final energy = {final_energy} eV")
print(f"Final ΔE = {delta_energy}")
print(f"Energy convergence: {'Yes' if energy_converged else 'No'}")

print(f"Final max force = {final_max_force} eV/Å")
print(f"Final RMS force = {final_rms_force} eV/Å")
print(f"Force convergence: {'Yes' if force_converged else 'No'}")

print(f"Final volume = {final_volume} Å³")
print(f"Final relative ΔV = {rel_delta_volume}")
print(f"Volume convergence: {'Yes' if volume_converged else 'No'}")

# --- Student-friendly guidance ---
if energy_converged and force_converged and volume_converged:
    print("\nRelaxation complete: energy, forces, and volume converged")
else:
    print("\nRelaxation not fully converged. Consider continuing relaxation or checking parameters.")

# --- Plotting ---
fig, axs = plt.subplots(3,1, figsize=(10,8), sharex=False)

# Energy
if total_energies:
    axs[0].plot(steps_energy, total_energies, marker='s', color='b')
    axs[0].set_ylabel("Total Energy (eV)")
    axs[0].set_title("Energy vs MD Step")
    axs[0].grid(True)

# Max Force
if max_forces:
    axs[1].plot(steps_forces, max_forces, marker='o', color='r', label="Max Force")
    axs[1].axhline(force_tol, color='k', linestyle='--', label="Force threshold")
    axs[1].set_ylabel("Max Force (eV/Å)")
    axs[1].set_title("Max Force vs MD Step")
    axs[1].grid(True)
    axs[1].legend()

# Volume
if volumes:
    axs[2].plot(steps_volume, volumes, marker='^', color='g')
    axs[2].set_ylabel("Volume (Å³)")
    axs[2].set_title("Volume vs MD Step")
    axs[2].grid(True)

plt.tight_layout()
plt.show()
