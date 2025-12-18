#!/usr/bin/env python3

import re
import matplotlib.pyplot as plt
import numpy as np

# -------------------------
# Read N_atoms from OUTCAR
# -------------------------
N_atoms = None
with open("OUTCAR") as f:
    for line in f:
        if "NIONS" in line:
            match = re.search(r"NIONS\s*=\s*(\d+)", line)
            if match:
                N_atoms = int(match.group(1))
                break
if N_atoms is None:
    raise ValueError("Could not find NIONS in OUTCAR")
print(f"N_atoms = {N_atoms}")

# -------------------------
# Read POTIM from INCAR or OUTCAR
# -------------------------
POTIM_fs = None
# Try INCAR first
try:
    with open("INCAR") as f:
        for line in f:
            if "POTIM" in line:
                POTIM_fs = float(line.split("=")[1].split()[0])
                break
except FileNotFoundError:
    pass

# If not found in INCAR, try OUTCAR
if POTIM_fs is None:
    with open("OUTCAR") as f:
        for line in f:
            if "POTIM" in line:
                POTIM_fs = float(line.split("=")[1].split()[0])
                break

if POTIM_fs is None:
    raise ValueError("Could not find POTIM in INCAR or OUTCAR")
print(f"POTIM = {POTIM_fs} fs")

# -------------------------
# Read energies, temperatures, and volumes from OUTCAR
# -------------------------
temperatures = []
total_energies = []
volumes = []

with open("OUTCAR") as f:
    for line in f:
        # Temperature
        if "ekin_lat" in line.lower():
            match = re.search(r"\(temperature\s*([\d\.]+)\s*K\)", line, re.I)
            if match:
                temperatures.append(float(match.group(1)))

        # Total energy
        if "etotal" in line.lower():
            match = re.search(r"etotal\s*=\s*([-\d\.Ee]+)", line, re.I)
            if match:
                total_energies.append(float(match.group(1)))

        # Volume
        if "volume" in line.lower():
            match = re.search(r"volume\s*(?:of cell)?\s*[:=]\s*([-\d\.Ee]+)", line, re.I)
            if match:
                volumes.append(float(match.group(1)))

# -------------------------
# Steps
# -------------------------
steps_temp = list(range(1, len(temperatures) + 1))
steps_energy = list(range(1, len(total_energies) + 1))
steps_volume = list(range(1, len(volumes) + 1))

# -------------------------
# Drift calculation for NVE
# -------------------------
if total_energies:
    deltaE_basic = total_energies[-1] - total_energies[0]
    coeffs = np.polyfit(steps_energy, total_energies, 1)  # slope, intercept
    slope = coeffs[0]  # eV per step
    deltaE_slope = slope * (len(steps_energy)-1)

    N_steps = len(total_energies)
    t_total_ps = POTIM_fs * N_steps * 1e-3  # fs -> ps

    drift_basic = deltaE_basic / (N_atoms * t_total_ps)
    drift_slope = deltaE_slope / (N_atoms * t_total_ps)

    print(f"\n--- NVE Drift ---")
    print(f"ΔE (start to end) = {deltaE_basic:.6f} eV, drift = {drift_basic:.2e} eV/atom/ps")
    print(f"Slope-based ΔE = {deltaE_slope:.6f} eV, drift = {drift_slope:.2e} eV/atom/ps")

# -------------------------
# Temperature analysis (SMASS)
# -------------------------
if temperatures:
    avg_temp = np.mean(temperatures)
    std_temp = np.std(temperatures, ddof=1)
    coeffs_temp = np.polyfit(steps_temp, temperatures, 1)
    slope_temp = coeffs_temp[0] * 1e3 / POTIM_fs  # Convert slope to K/ps
    expected_sigma = avg_temp / np.sqrt(3*N_atoms)

    print(f"\n--- Temperature Analysis (SMASS) ---")
    print(f"Avg Temp = {avg_temp:.2f} K, Std = {std_temp:.2f} K")
    print(f"Expected fluctuation (σ_T) = {expected_sigma:.2f} K")
    print(f"Slope (systematic drift) = {slope_temp:.4f} K/ps")

# -------------------------
# Volume analysis (PMASS)
# -------------------------
if volumes:
    avg_vol = np.mean(volumes)
    std_vol = np.std(volumes, ddof=1)
    rel_fluct = std_vol / avg_vol  # relative fluctuation
    coeffs_vol = np.polyfit(steps_volume, volumes, 1)
    slope_vol = coeffs_vol[0] * 1e3 / POTIM_fs  # Convert slope to Å^3/ps
    rel_drift = abs(slope_vol) / avg_vol  # relative drift per ps

    print(f"\n--- Volume Analysis (PMASS) ---")
    print(f"Avg Volume = {avg_vol:.2f} Å³")
    print(f"Std Volume = {std_vol:.2f} Å³")
    print(f"Relative fluctuation = {rel_fluct:.2%}")
    print(f"Slope (systematic drift) = {slope_vol:.4f} Å³/ps")
    print(f"Relative volume drift = {rel_drift:.2%}/ps")  # new line added

    
# -------------------------
# Plotting
# -------------------------
fig, axs = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

# Temperature
if temperatures:
    axs[0].plot(steps_temp, temperatures, marker='o', color='r', label="Temperature")
    axs[0].axhline(avg_temp, color='k', linestyle='--', label=f"Avg = {avg_temp:.2f} K")
    axs[0].set_ylabel("Temperature (K)")
    axs[0].grid(True)
    axs[0].legend()
else:
    axs[0].text(0.5, 0.5, "No temperature data found", ha='center', va='center')

# Total energy
if total_energies:
    axs[1].plot(steps_energy, total_energies, marker='s', color='b', label="Total Energy")
    axs[1].axhline(np.mean(total_energies), color='k', linestyle='--', label=f"Avg = {np.mean(total_energies):.4f} eV")
    axs[1].set_ylabel("Total Energy (eV)")
    axs[1].grid(True)
    axs[1].legend()
else:
    axs[1].text(0.5, 0.5, "No energy data found", ha='center', va='center')

# Volume
if volumes:
    axs[2].plot(steps_volume, volumes, marker='^', color='g', label="Volume")
    axs[2].axhline(avg_vol, color='k', linestyle='--', label=f"Avg = {avg_vol:.2f} Å³")
    axs[2].set_ylabel("Cell Volume (Å³)")
    axs[2].set_xlabel("Step")
    axs[2].grid(True)
    axs[2].legend()
else:
    axs[2].text(0.5, 0.5, "No volume data found", ha='center', va='center')
    axs[2].set_xlabel("Step")

plt.tight_layout()
plt.show()

