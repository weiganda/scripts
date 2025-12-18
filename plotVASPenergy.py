#!/usr/bin/env python3

import re
import matplotlib.pyplot as plt

filename = "OUTCAR"

temperatures = []
total_energies = []
volumes = []

with open(filename, "r") as f:
    for line in f:

        # Total energy
        if "free energy" in line and "TOTEN" in line:
            match = re.search(r"=\s*([-\d\.Ee]+)", line)
            if match:
                total_energies.append(float(match.group(1)))

# Steps based on available data
steps_energy = list(range(1, len(total_energies)+1))

# --- PLOTTING ---
fig, ax = plt.subplots(figsize=(8, 5))  # single Axes object

if total_energies:
    ax.plot(steps_energy, total_energies, marker='s', color='b')
    ax.set_ylabel("Total Energy (eV)")
    ax.set_xlabel("Step")
    ax.grid(True)
else:
    ax.text(0.5, 0.5, "No energy data found", ha='center', va='center')

plt.tight_layout()
plt.show()

