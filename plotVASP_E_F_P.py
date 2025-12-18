#!/usr/bin/env python3

import re
import math
import matplotlib.pyplot as plt

filename = "OUTCAR"

total_energies = []
max_forces = []
pressures = []

with open(filename, "r") as f:
    lines = f.readlines()

# Temporary storage for per-step force data
force_block = []
collect_force = False

for line in lines:

    # --- Total energy ---
    if "free energy" in line and "TOTEN" in line:
        match = re.search(r"=\s*([-\d\.Ee]+)", line)
        if match:
            total_energies.append(float(match.group(1)))

    # --- Forces ---
    if "TOTAL-FORCE" in line:
        force_block = []
        collect_force = True
        continue

    if collect_force:
        if line.strip() == "":
            collect_force = False
            if force_block:
                forces = []
                for l in force_block:
                    try:
                        # Compute vector magnitude of the force
                        fx, fy, fz = map(float, l.split()[-3:])
                        forces.append(math.sqrt(fx**2 + fy**2 + fz**2))
                    except ValueError:
                        continue
                if forces:
                    max_forces.append(max(forces))
        else:
            force_block.append(line)

    # --- Pressure ---
    if "external pressure" in line:
        match = re.search(r"external pressure\s*=\s*([-\d\.Ee]+)", line)
        if match:
            pressures.append(float(match.group(1)))

# Steps
energy_steps = list(range(1, len(total_energies)+1))
force_steps = list(range(1, len(max_forces)+1))
pressure_steps = list(range(1, len(pressures)+1))

# Print max force and pressure values
print("Step\tMax Force (eV/Å)\tPressure (kB)")
for i in range(max(len(max_forces), len(pressures))):
    f_val = max_forces[i] if i < len(max_forces) else ""
    p_val = pressures[i] if i < len(pressures) else ""
    print(f"{i+1}\t{f_val}\t{p_val}")

def get_ylims(data, padding=0.1):
    """Return min/max y-limits with padding."""
    min_val = min(data)
    max_val = max(data)
    delta = max_val - min_val
    if delta == 0:
        delta = abs(min_val) * 0.1 if min_val != 0 else 1
    return (min_val - padding*delta, max_val + padding*delta)

# --- PLOTTING ---
fig, axs = plt.subplots(3, 1, figsize=(10, 12), sharex=False)

# Total energy subplot
axs[0].plot(energy_steps, total_energies, 'b-', marker='s')
axs[0].set_ylabel("Total Energy (eV)", color='b')
axs[0].tick_params(axis='y', labelcolor='b')
axs[0].grid(True)
axs[0].set_title("VASP Relaxation Data")
axs[0].set_ylim(get_ylims(total_energies))

# Max force subplot
if max_forces:
    axs[1].plot(force_steps, max_forces, 'r--', marker='o')
    axs[1].set_ylim(get_ylims(max_forces))
axs[1].set_ylabel("Max Force (eV/Å)", color='r')
axs[1].tick_params(axis='y', labelcolor='r')
axs[1].grid(True)

# Pressure subplot
if pressures:
    axs[2].plot(pressure_steps, pressures, 'g-.', marker='^')
    axs[2].set_ylim(get_ylims(pressures))
axs[2].set_xlabel("Step")
axs[2].set_ylabel("Pressure (kB)", color='g')
axs[2].tick_params(axis='y', labelcolor='g')
axs[2].grid(True)

plt.tight_layout()
plt.show()

