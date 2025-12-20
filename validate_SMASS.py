#!/usr/bin/env python3

"""
Created by: Alysse Weigand
Last Updated: 12/20/2025

Development Notes:
- Conceptual design, purpose, and validation by Alysse Weigand.
- All scientific reasoning, method choices, and interpretation of results by Alysse Weigand.
- Code implementation and structure assisted by ChatGPT.

SMASS Validation Script

This script is used to determine and validate the thermostat mass (SMASS) for NVT 
molecular dynamics simulations. It monitors the temperature fluctuations of the 
system and compares them to the expected statistical values to ensure the thermostat 
is neither too weak nor too strong. 

This script is intended to be used after the POTIM value has been validated using
an NVE ensemble and and validate_POTIM.py script. 

Usage:
- Run this script before production NVT simulations to choose an appropriate SMASS.
- Analyze temperature fluctuations to select a SMASS that maintains stable and accurate 
  temperature control.

Intuition for Students:

- SMASS is the “thermostat mass” that controls how quickly the system responds to temperature fluctuations.
- Too small → thermostat reacts too fast → temperature oscillates (underdamped).
- Too large → thermostat reacts too slowly → temperature barely fluctuates (overdamped).
- Ideal → critically damped: temperature fluctuates around the mean naturally, consistent with canonical ensemble.
- This script analyzes NVT temperature data to suggest whether SMASS is appropriate.

"""

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

# -------------------------
# Read SMASS from INCAR if available
# -------------------------
SMASS = None

try:
    with open("INCAR") as f:
        for line in f:
            if "SMASS" in line:
                SMASS = float(line.split("=")[1].split()[0])

except FileNotFoundError:
    pass

# Use defaults or indicate missing values
if SMASS is None:
    SMASS = 0.0
    print("SMASS not found in INCAR")


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


# -----------------------------
# Simulation summary
# -----------------------------
print("\n---------- Simulation Parameters -----------")
print(f"Number of atoms (N_atoms) = {N_atoms}")
print(f"POTIM (timestep) = {POTIM_fs:.3f} fs")
print(f"SMASS (thermostat mass) = {SMASS}")
print("------------------------------------------\n")



# -------------------------
# Autocorrelation and Effective Sample Size
# -------------------------
def autocorr(x):
    x = np.array(x)
    n = len(x)
    x_mean = np.mean(x)
    var = np.var(x, ddof=1)
    c = np.correlate(x - x_mean, x - x_mean, mode='full') / (var * n)
    return c[n-1:]

if temperatures:
    # Steady-state temperatures (last 50%)
    #   Using the last half avoids initial equilibrium effects
    n_temp = len(temperatures)
    start_idx = n_temp // 2
    temperatures_steady = np.array(temperatures[start_idx:])
    steps_temp_steady = np.array(steps_temp[start_idx:])
    
    c_temp = autocorr(temperatures_steady)
    # Find correlation time (first lag where autocorr < 0.1)

    print("\n---------------- Sample Size -----------------------")
    try:
        tau_corr = np.where(c_temp < 0.1)[0][0]
        N_eff = len(temperatures_steady) / tau_corr
        print(f"Estimated correlation time = {tau_corr} steps, N_eff = {N_eff:.1f} independent samples")
        if N_eff < 30:
            print("ALERT: Not enough independent samples for reliable statistics. Consider longer simulation.")
    except IndexError:
        print("WARNING: Autocorrelation did not drop below threshold — unable to estimate correlation time")
    print("----------------------------------------------------\n")
    




# -----------------------------
# Suggested starting parameters
# -----------------------------
# POTIM (fs): timestep; must resolve fastest vibrations
# SMASS: thermostat mass; adjust to control temperature fluctuations
#
# Examples for solid-state materials:
#   Ionic solids (NaCl, MgO):           POTIM=1.5,      SMASS=0.2,      PMASS=200
#   Covalent solids (Si, C diamond):    POTIM=1.0,      SMASS=0.2,      PMASS=200
#   Metals (Al, Cu, Au):                POTIM=2.0,      SMASS=0.5,      PMASS=500
#   Molecular solids (ice, NH3):        POTIM=0.5-1.0,  SMASS=0.1-0.5,  PMASS=200
#
# Users should always validate POTIM with NVE energy drift and then adjust SMASS/PMASS.
# 
#   ALYSSE GRAB CITIATIONS, RULIS WILL WANT THOSE

# -------------------------
# Temperature analysis (SMASS)
# -------------------------
# The temperature analysis gives insight into how well the thermostat is regulating the system.
# SMASS controls the "mass" of the thermostat and affects how quickly the system responds to temperature fluctuations.
# This test should be run on an NPT.
#
# Two types of temperature metrics are calculated here:
#   1. Average and standard deviation of temperature:
#       
#           Avg Temp = mean(T)
#           Std Temp = std(T)
#
#       The expected fluctuation for N atoms is:
#           
#           σ_T = Avg Temp / sqrt(3 * N)
#
#       3 degrees of freedom
#
#       If std_temp >> σ_T, the thermostat is too light (SMASS too small). The temperature fluctuates excessively.
#       If std_temp << σ_T, the thermostat is too heavy (SMASS too large). The temperature is too sluggish.
#
#   2. Temperature drift (slope_temp):
#       
#           T(step) ≈ slope_temp * step + intercept
#
#       Linear regression of temperature vs. step gives the systematic drift of temperature.
#       Converted to K/ps via:
#           
#           slope_temp = slope_step * 1e3 / POTIM_fs
#
#       Ideally, slope_temp should be close to zero. A significant slope indicates that SMASS may need adjustment.
#
# Interpretation guidelines:
#   - std_temp >> σ_T, Thermostat too light, consider increasing SMASS
#   - std_temp << σ_T, Thermostat too heavy, consider decreasing SMASS
#   - |slope_temp| large → System is drifting in temperature, may indicate SMASS or timestep issues
#
# Note: SMASS should be chosen to balance realistic fluctuations with smooth temperature control,
#       generally resulting in ~1–2% relative temperature fluctuations.

if temperatures:
    # --- Steady-state temperature statistics ---
    avg_temp = np.mean(temperatures_steady)
    std_temp = np.std(temperatures_steady, ddof=1)

    coeffs_temp = np.polyfit(steps_temp_steady, temperatures_steady, 1)
    slope_temp = coeffs_temp[0] * 1e3 / POTIM_fs  # K/ps

    # Expected statistical fluctuation in a canonical ensemble (3 Degrees of Freedom per atom)
    expected_sigma = avg_temp / np.sqrt(3 * N_atoms)

    print(f"\n--- Temperature Analysis (SMASS) [last 50% only] ---")
    print(f"Avg Temp = {avg_temp:.2f} K")
    print(f"Std = {std_temp:.2f} K")
    print(f"Expected fluctuation (σ_T) = {expected_sigma:.2f} K")
    print(f"Slope (systematic drift) = {slope_temp:.4f} K/ps")

    # --- Evaluate thermostat behavior ---
    if std_temp > 2 * expected_sigma:
        print("ALERT: Temperature fluctuates too much! Increase SMASS!!!")
    elif std_temp < 0.5 * expected_sigma:
        print("ALERT: Temperature fluctuations too small! Decrease SMASS!!!")
    else:
        print("Temperature fluctuations within expected range")

    # --- Evaluate temperature drift ---
    if abs(slope_temp) < 0.01:
        print("Temperature drift negligible. SMASS is fine.")
    elif abs(slope_temp) <= 0.05:
        print("Temperature drift moderate. Monitor simulation")
    else:
        print("ALERT: Significant temperature drift!!! Adjust SMASS.")

    print(f"----------------------------------------------------\n")


# -----------------------------
# IMPORTANT NOTE ON ADJUSTMENTS
# -----------------------------
# Although the ALERTS suggest adjusting the SMASS, PMASS, or timestep, the POTIM, once validated
# through an NVE simulation, should generally remain fixed.
# 
# - POTIM validation: Use NVE energy drift to ensure the timestep is small enough to capture
#   the fastest vibrations without causing excessive drift.
# - Once POTIM is validated, any subsequent tuning for temperature or volume stability should
#   focus on SMASS (thermostat) or PMASS (barostat), not the timestep.
#
# ------------------------------------------------------------------
# Suggestions Explanation
# ------------------------------------------------------------------
#
# This block provides guidance for adjusting SMASS based on observed
# temperature fluctuations compared to the expected theoretical fluctuations.
#
# 1. Compute expected fluctuation:
#      σ_T = Avg Temp / sqrt(3 * N_atoms)
#    - This is the statistical fluctuation for N atoms in the canonical ensemble.
#
# 2. Compare observed standard deviation (std_temp) to σ_T:
#    - If std_temp > 2 * σ_T:
#         - Temperature fluctuates more than expected
#         - Thermostat is too "light" (SMASS too small), reacts too quickly
#         - Suggest increasing SMASS proportionally to the severity of the overshoot
#
#    - If std_temp < 0.5 * σ_T:
#         - Temperature fluctuates less than expected
#         - Thermostat is too "heavy" (SMASS too large), reacts too slowly
#         - Suggest decreasing SMASS proportionally to the severity of the undershoot
#
#    - If 0.5 * σ_T <= std_temp <= 2 * σ_T:
#         - Temperature fluctuations are within expected range
#         - SMASS is considered appropriate
#
# 3. Severity factor:
#    - Quantifies how far the observed fluctuation is from the expected value
#    - Used to scale suggested SMASS up or down
#
# 4. Notes:
#    - These suggestions assume the simulation is already at steady state
#      (i.e., the last half of the temperature trajectory is analyzed)
#    - POTIM should already be validated using NVE drift; do not adjust timestep here
#    - SMASS tuning is aimed at achieving a critically damped thermostat response
#      where temperature autocorrelation decays smoothly and N_eff is sufficient
#
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# Severity Factor
# ------------------------------------------------------------------
#
# The severity factor quantifies **how far the observed temperature fluctuations
# deviate from the expected theoretical fluctuations**. It is used to scale
# the suggested SMASS adjustment proportionally. The severity factor is the
# percentage of deviation from the expected fluctuations. 
#
# Definition:
#   - If temperature fluctuates too much (std_temp > expected_sigma):
#         severity = (std_temp - expected_sigma) / expected_sigma
#         - Measures the relative excess in fluctuations
#         - Larger severity → greater overshoot → increase SMASS more
#
#   - If temperature fluctuates too little (std_temp < expected_sigma):
#         severity = (expected_sigma - std_temp) / expected_sigma
#         - Measures the relative deficit in fluctuations
#         - Larger severity → thermostat too heavy → decrease SMASS more
#
# Interpretation:
#   - severity ≈ 0 → temperature fluctuations are close to expected → small or no adjustment
#   - severity ≈ 1 → fluctuations are twice the expected deviation → substantial adjustment
#
# Application:
#   - For excessive fluctuations:
#         suggested_SMASS = SMASS * (1 + severity)
#   - For too stable fluctuations:
#         suggested_SMASS = SMASS / (1 + severity)
#
# Notes:
#   - This is not exact, it is a best guess, exact scaling may need fine-tuning
#     depending on material and system size
#   - Works best when applied to the steady-state portion of the trajectory
#   - Always cross-check with autocorrelation and effective sample size (N_eff)
#     to ensure the adjustment will improve thermostat performance
#
# ------------------------------------------------------------------

print("\n------------- Suggestions (Improved) -------------")

# Severity factor quantifies how far the observed fluctuations deviate from theory; used to scale
#   SMASS adjustments. Again, this is just a best guess to help guide future SMASS choices.

if temperatures_steady is not None:
    # -------------------------------
    # Base severity from std deviation
    # -------------------------------
    if std_temp > expected_sigma:  # fluctuations too high
        severity = (std_temp - expected_sigma) / expected_sigma
        direction = "increase"
    elif std_temp < expected_sigma:  # fluctuations too low
        severity = (expected_sigma - std_temp) / expected_sigma
        direction = "decrease"
    else:
        severity = 0
        direction = None

    # -------------------------------
    # Incorporate autocorrelation
    # -------------------------------
    # Weight based on effective correlation time (tau_corr)
    # Longer tau_corr → slower thermostat response → reduce adjustment
    weight = min(1.0, tau_corr / 10)  # cap weight at 1.0

    # Adjust for underdamping: if autocorr oscillates (first negative crossing), boost severity
    first_neg = np.where(c_temp < 0)[0]
    if len(first_neg) > 0 and first_neg[0] < len(c_temp)//2:
        damping_factor = 1.5  # underdamped → increase severity
        damping_note = "Oscillatory autocorrelation detected (underdamped)"
    else:
        damping_factor = 1.0  # critically or overdamped
        damping_note = "Autocorrelation smooth (critically/overdamped)"

    severity_adj = severity * weight * damping_factor

    # -------------------------------
    # Cap adjustment to avoid overshoot
    # -------------------------------
    max_change = 0.5  # maximum 50% change at a time
    severity_adj = min(severity_adj, max_change)

    # -------------------------------
    # Compute suggested SMASS
    # -------------------------------
    if direction == "increase":
        suggested_SMASS = SMASS * (1 + severity_adj)
        print(f"Temperature fluctuates too much (std={std_temp:.2f} K, expected={expected_sigma:.2f} K).")
        print(f"{damping_note}")
        print(f"Consider increasing SMASS to ~{suggested_SMASS:.2f} (incremental adjustment).")
    elif direction == "decrease":
        suggested_SMASS = SMASS / (1 + severity_adj)
        print(f"Temperature too stable (std={std_temp:.2f} K, expected={expected_sigma:.2f} K).")
        print(f"{damping_note}")
        print(f"Consider decreasing SMASS to ~{suggested_SMASS:.2f} (incremental adjustment).")
    else:
        print("Temperature fluctuations within expected range. SMASS is appropriate.")
        print(f"{damping_note}")

    # -------------------------------
    # Check if enough independent samples
    # -------------------------------
    if N_eff < 30:
        print("WARNING: Not enough independent samples (N_eff < 30). Consider running a longer simulation.")

print(f"---------------------------------------\n")
print("\nNote: Once POTIM is validated via NVE drift, only SMASS and PMASS should be adjusted.")

print("\n---------------- Tips for Students ----------------")
print("1. POTIM should already be validated using NVE drift; do NOT change timestep here.")
print("2. SMASS controls thermostat damping; aim for critically damped response.")
print("3. Check temperature fluctuation vs expected σ_T: too high → increase SMASS; too low → decrease SMASS.")
print("4. Use temperature autocorrelation plot to visually inspect damping behavior.")
print("5. Effective sample size (N_eff) should be > 30 for reliable statistics; otherwise run longer.")
print("6. Adjust SMASS incrementally; avoid large jumps.")
print("--------------------------------------------------\n")

    
# ------------------------------------------------
# Plotting, Desciptions, and Interpretation Guide
# ------------------------------------------------
fig, axs = plt.subplots(3, 1, figsize=(12, 6), sharex=False)

# 1. Total Energy vs MD Step
#    - Shows the evolution of the total energy throughout the simulation.
#    - The black dashed line indicates the average energy.
#
if total_energies:
    axs[0].plot(steps_energy, total_energies, marker='s', color='b', label="Total Energy")
    axs[0].set_ylabel("Total Energy (eV)")
    axs[0].set_title("Energy vs MD Step")
    axs[0].grid(True)
    axs[0].legend()
else:
    axs[0].text(0.5, 0.5, "No energy data found", ha='center', va='center')

# 2. Temperature vs MD Step
#    - Displays the instantaneous temperature over the trajectory.
#    - The black dashed line represents the average steady-state temperature.
#    - The gray shaded region represents the expected fluctuation range (±σ_T) for N atoms:
#          σ_T = Avg Temp / sqrt(3 * N_atoms)
#    - Observing temperature mostly within the gray band indicates the thermostat is correctly controlling T.
#    - Large deviations may indicate SMASS is too small (thermostat reacts too fast) or too large (too sluggish).
#    - Slope of the temperature (linear fit) shows systematic drift:
#          - |slope| < 0.01 K/ps → negligible drift
#          - 0.01 ≤ |slope| ≤ 0.05 K/ps → moderate drift
#          - |slope| > 0.05 K/ps → significant drift; consider adjusting SMASS
if temperatures:
    axs[1].plot(steps_temp, temperatures, marker='o', color='r', label="Temperature")
    axs[1].axhline(avg_temp, color='k', linestyle='--', label=f"Avg = {avg_temp:.2f} K")
    axs[1].fill_between(steps_temp_steady, avg_temp - expected_sigma, avg_temp + expected_sigma, color='gray', alpha=0.2, label='Expected ±σ_T')
    axs[1].set_ylabel("Temperature (K)")
    axs[1].set_title("Temperature vs MD Step")
    axs[1].grid(True)
    axs[1].legend()
else:
    axs[1].text(0.5, 0.5, "No temperature data found", ha='center', va='center')

# 3. Temperature Autocorrelation Function
#    - Quantifies how temperature values are correlated over time.
#    - C(0) = 1 and C(τ) → 0 as lag τ increases.
#    - The lag where C(τ) < 0.1 gives the correlation time (τ_corr).
#    - Effective number of independent samples: N_eff = N_steps / τ_corr
#        - Ensures statistical reliability of mean and std temperature.
#        - N_eff < 30 → insufficient independent data; consider longer simulation.
#    - Short correlation time → thermostat too light (fast fluctuations)
#    - Long correlation time → thermostat too heavy (slow fluctuations)
# NOTE: You do NOT want oscillation on this plot. You want to see a smooth decay
#   to the dashed line. If you see oscillation your system is underdamped. If you
#   see very very slow decay, your system is over damped. You should think about
#   this plot like a damped oscillator. You want critical dampening. 
if temperatures:
    c_temp = autocorr(temperatures_steady)
    axs[2].plot(c_temp, marker='o', color='purple')
    axs[2].axhline(0.1, color='r', linestyle='--', label="Corr < 0.1 threshold")
    axs[2].set_ylabel("Autocorrelation")
    axs[2].set_xlabel("Lag")
    axs[2].set_title("Temperature Autocorrelation Function")
    axs[2].grid(True)
    axs[2].legend()
else:
    axs[2].text(0.5, 0.5, "No temperature data found", ha='center', va='center')

plt.tight_layout()


# 4. Steady-State Temperature Distribution Histogram
#    - Compares the observed temperature distribution (last half of simulation) to the expected Gaussian (σ_T).
#    - Histogram shows actual fluctuations; dashed line shows ideal thermal distribution.
#    - Good agreement indicates SMASS yields physically realistic temperature sampling.
#    - Large deviations suggest the thermostat is not properly sampling the canonical ensemble.
#
if temperatures:
    plt.figure(figsize=(6,4))
    plt.hist(temperatures_steady, bins=30, density=True, alpha=0.6, color='red', label='Observed T')
    x = np.linspace(min(temperatures_steady), max(temperatures_steady), 200)
    from scipy.stats import norm
    y = norm.pdf(x, loc=avg_temp, scale=expected_sigma)
    plt.plot(x, y, 'k--', lw=2, label='Expected Gaussian (σ_T)')
    plt.xlabel("Temperature (K)")
    plt.ylabel("Probability density")
    plt.title("Steady-State Temperature Distribution")
    plt.legend()
    plt.grid(True)
    plt.show()



