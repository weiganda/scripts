#!/usr/bin/env python3
"""
Created by: Alysse Weigand
Last Updated: 12/20/2025

Development Notes:
- Conceptual design, purpose, and validation by Alysse Weigand.
- All scientific reasoning, method choices, and interpretation of results by Alysse Weigand.
- Code implementation and structure assisted by ChatGPT.

================================================================================
SMASS Optimization from Multiple NPT Runs (Robust Version)
================================================================================

This script analyzes multiple NPT MD simulation runs at different SMASS values
to determine an "ideal" thermostat mass (SMASS) that results in critically
damped temperature fluctuations.

--------------------------------------------------------------------------------
INPUT FILE FORMAT (CSV)
--------------------------------------------------------------------------------
The input CSV should contain the following columns (one row per simulation run):

SMASS   : Thermostat mass used in the simulation (float)
Std_T   : Observed temperature fluctuation during steady-state (K)
Sigma_T : Expected temperature fluctuation for N atoms (K)

This data can be collected from the validate_SMASS.py script. 

Example:

SMASS,Std_T,Sigma_T
0.05,40.1,12.6
0.10,38.2,12.6
0.15,25.0,12.6
0.20,15.5,12.6
0.30,11.0,12.6

--------------------------------------------------------------------------------
WHAT THE SCRIPT DOES
--------------------------------------------------------------------------------
1. Reads the CSV file.
2. Computes normalized overshoot for each run:
       overshoot = (Std_T - Sigma_T) / Sigma_T
   - Positive overshoot → underdamped (thermostat too light)
   - Negative overshoot → overdamped (thermostat too heavy)
3. Uses linear regression to estimate SMASS where overshoot ≈ 0,
   corresponding to critically damped behavior.
   - Works even if all points are underdamped or overdamped.
4. Plots overshoot vs SMASS and optionally shows predicted critically damped SMASS.
5. Prints a recommended SMASS for future simulations.

================================================================================

Note: This is a rough estimate to guide your SMASS choices. As you collect more
data points, the estimation will improve. 
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# User input: CSV filename
# -------------------------
input_file = "smass_runs.csv"

# -------------------------
# Read the CSV file
# -------------------------
df = pd.read_csv(input_file)

# -------------------------
# Compute overshoot: normalized deviation from expected fluctuation
# Positive → underdamped, Negative → overdamped
# -------------------------
df["overshoot"] = (df["Std_T"] - df["Sigma_T"]) / df["Sigma_T"]

# -------------------------
# Linear regression: SMASS vs overshoot
# -------------------------
# We fit overshoot = a * SMASS + b
a, b = np.polyfit(df["SMASS"], df["overshoot"], 1)

# Critically damped SMASS: overshoot = 0
SMASS_crit = -b / a

# -------------------------
# Print results
# -------------------------
print("\n================ Recommended SMASS =================")
print(f"Estimated critically damped SMASS = {SMASS_crit:.3f}")
print("====================================================\n")

# -------------------------
# Plotting
# -------------------------
fig, ax = plt.subplots(figsize=(8, 6))

ax.scatter(df["SMASS"], df["overshoot"], color="b", s=80, label="Simulation data")
ax.plot(df["SMASS"], a*df["SMASS"] + b, color="g", linestyle="--", label="Linear fit")
ax.axhline(0, color="k", linestyle="--", label="Critical damping (overshoot=0)")
ax.axvline(SMASS_crit, color="r", linestyle="--", label=f"Estimated SMASS = {SMASS_crit:.3f}")

ax.set_xlabel("SMASS")
ax.set_ylabel("Overshoot (Std_T / Sigma_T - 1)")
ax.set_title("Temperature Fluctuation Overshoot vs SMASS")
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.show()
