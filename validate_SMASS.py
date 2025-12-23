#!/usr/bin/env python3

"""
Created by: Alysse Weigand
Last Updated: 12/22/2025

Development Notes:
    Conceptual design, purpose, and validation by Alysse Weigand.
    All scientific reasoning, method choices, and interpretation of results by Alysse Weigand.
    Code implementation and structure assisted by ChatGPT.

SMASS Validation Script

This script is used to determine and validate the thermostat mass (SMASS) for NVT 
molecular dynamics simulations. It monitors the temperature fluctuations of the 
system and compares them to the expected statistical values to ensure the thermostat 
is neither too weak nor too strong. 

This script is intended to be used after the POTIM value has been validated using
an NVE ensemble and and validate_POTIM.py script. 

Usage:
    Run this script before production NVT simulations to choose an appropriate SMASS.
    Analyze temperature fluctuations to select a SMASS that maintains stable and accurate 
        temperature control.

Intuition for Students:
    SMASS sets the characteristic timescale of the thermostat relative to the system’s intrinsic 
        vibrational frequencies. If SMASS is too small, the thermostat reacts too quickly. 
    Instead of gently correcting the temperature, it overcorrects, causing the temperature 
        to swing back and forth like an oscillator.
    If SMASS is too large, the thermostat reacts too slowly. The temperature then drifts around 
        for a long time before being corrected, so the thermostat doesn’t do its job very well.

    A good SMASS value makes the temperature fluctuate naturally around the target value on 
        realistic timescales, but only if the thermostat is behaving properly overall.

Important:
    For solid materials using a single Nose–Hoover thermostat, the temperature can keep oscillating
        no matter what SMASS you choose. When this happens, the system is not sampling temperature 
        correctly, and you cannot tell whether SMASS is good or bad just by looking at temperature 
        fluctuations.

    Because of this, the script first checks whether the temperature actually forgets its past 
        (using the temperature autocorrelation). Only if the temperature decorrelates quickly 
        and does not oscillate does the script try to decide whether SMASS is appropriate.

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
# There are two purposes for the autocorrelation function and the interpretation
#   can be a bit confusing and needs to be discussed here. 
#
# The first purpose of the autocorrelation function is to determine whether or
#   not any of the statistical information below can be considered valid.
#
#   The autocorrelation function is going to determine if there is enough
#       independent temperature data for meaningful statistics.
#
#   If the correlation time is too long (>20 steps) then the SMASS is likely too
#       large and due to the nature of the autocorrelation function (see below).
#       The data can be misinterpreted leading into a black hole (personal experience).
#       So, to save everyone time, a variety of alerts will be triggered.
#
#   It is important to note that this threshold may be too large and some logic
#       or intuition will need to be used. This threshold is a "best guess".
#
#   Ok, so interpreting the autocorrelation (ACF) data and plots. What we are looking for 
#       here is oscillation. But watch out, oscillation happens at both extremes
#
#       1. SMASS too small - High-frequency and large amplitude oscillations present
#           in the ACF plot. Underdamped. 
#
#       2. SMASS perfect - Smooth decay, no oscillations. Critically damped.
#
#       3. SMASS too large - Slow decay, stuck. Overdamped.
#
#       4. SMASS way too large - Low-frequency and low-amplitude oscillations. 
#
# The second purpose of the ACF is to determine effective sample size. Do we
#   have enough independet samples for statistical analysis. I set the value at 
#   30. The ACF tau * the effective sample size is going to be the number of steps
#   that the simulation has completed.
#
#   Additionally, it is important to note that we are only using the last half of
#       the temperature data that is being collected. This is to ensure that we
#       are ignoring any of the initial equilibration information and only looking
#       for the steady-state data.
#
def autocorr(x):
    x = np.array(x)
    n = len(x)
    x_mean = np.mean(x)
    var = np.var(x, ddof=1)
    c = np.correlate(x - x_mean, x - x_mean, mode='full') / (var * n)
    return c[n-1:]

if temperatures:
    # Steady-state temperatures (last half)
    #   Using the last half avoids initial equilibrium effects
    n_temp = len(temperatures)
    start_idx = n_temp // 2
    temperatures_steady = np.array(temperatures[start_idx:])
    steps_temp_steady = np.array(steps_temp[start_idx:])
   
    # Find the correlation time
    c_temp = autocorr(temperatures_steady)

    print("\n---------------- Sample Size -----------------------")
    try:
        tau_corr = np.where(c_temp < 0.1)[0][0]
        N_eff = len(temperatures_steady) / tau_corr
        print(f"Estimated correlation time = {tau_corr} steps, N_eff = {N_eff:.1f} independent samples")

        # If the correlation time is too long, all statistical information is invalid. The SMASS is too
        #   large. A tau of 20 was chosen as a testing value. 
        if tau_corr > 20:
            print("ALERT: The correlation time may too long.")
            print("ALERT: If you see low-frequency oscillations, DECREASE SMASS!!!")
            print("Ignore all following statistical data")   
        if N_eff < 30:
            print("ALERT: Not enough independent samples for reliable statistics. Consider longer simulation.")
    except IndexError:
        print("WARNING: Autocorrelation did not drop below threshold - unable to estimate correlation time")
    print("----------------------------------------------------\n")

# -------------------------
# Energy drift analysis (steady state)
# -------------------------
# This analysis checks for energy drift in a simulation, typically performed after the system has
#   reached steady state (e.g., second half of the run). It helps determine whether energy is being
#   artificially drained or injected. This can be used for secondary analysis.
#
# Calculated:
#   1. Energy fluctuations (std_E)
#
#       std_E = std(E_total_steady)
#
#       This measures the random fluctuations of energy in the steady-state portion of the simulation.
#
#   2. Energy drift (slope_E)
#
#       E(t) ~ slope_E * t + intercept_E
#
#       Linear regression of total energy vs. simulation time (ps) gives the systematic energy drift.
#       Converted to eV/ps using the timestep (POTIM_fs):
#
#           slope_E = slope of linear fit over steady-state portion
#
#   3. Total drift across window (total_drift)
#
#       total_drift = slope_E * (t_end - t_start)
#
#       Compares the magnitude of systematic drift to typical fluctuations (std_E).
#
# Interpretation guidelines:
#   drift_ratio = |total_drift| / std_E
#       > 1.0, Significant drift, energy being drained or injected
#       0.3-1.0, Possible weak drift, consider longer run
#       < 0.3, Energy appears stationary within statistical noise
#
# Notes:
#   Only the second half of the simulation is analyzed to focus on steady-state behavior.
#   Proper energy conservation ensures realistic dynamics and reliable fluctuation statistics.

if total_energies:
    n_energy = len(total_energies)
    start_idx_E = n_energy // 2

    energies_steady = np.array(total_energies[start_idx_E:])
    steps_energy_steady = np.array(steps_energy[start_idx_E:])

    # Convert steps to time (ps)
    time_ps = steps_energy_steady * POTIM_fs * 1e-3

    # Linear fit: E(t) = a t + b
    slope_E, intercept_E = np.polyfit(time_ps, energies_steady, 1)

    # Energy fluctuations
    std_E = np.std(energies_steady, ddof=1)

    # Total drift across window
    total_drift = slope_E * (time_ps[-1] - time_ps[0])
    print("\n------------- Energy Drift Analysis -------------")
    print(f"Energy drift slope = {slope_E:.3e} eV/ps")
    print(f"Energy std (steady) = {std_E:.3e} eV")
    print(f"Total drift over window = {total_drift:.3e} eV")

    drift_ratio = abs(total_drift) / std_E if std_E > 0 else np.inf
    print(f"Drift-to-noise ratio = {drift_ratio:.2f}")

    if drift_ratio > 1.0:
        print("ALERT: Systematic energy drift detected (thermostat draining or injecting energy).")
    elif drift_ratio > 0.3:
        print("WARNING: Possible weak energy drift, consider longer run.")
    else:
        print("Energy appears stationary within statistical noise.")
    print("-------------------------------------------------\n")



# -------------------------
# Temperature analysis (SMASS)
# -------------------------
# The temperature analysis gives insight into how well the thermostat is regulating the system.
# SMASS controls the "mass" of the thermostat and affects how quickly the system responds to temperature fluctuations.
# This test should be run on an NVT.
#
# Two types of temperature metrics are calculated here:
#   1. Average and standard deviation of temperature:
#       
#           Avg Temp = mean(T) - Are we staying at the temperature we want?
#           Std Temp = std(T) - How big are the fluctuations around this temperature?
#
#       The expected fluctuation for N atoms is:
#           
#           sigma_T = sqrt(2 * Avg Temp^2) / (3 * N))
#
#       See:
#           Hickman J, Mishin Y. Temperature fluctuations in canonical systems: 
#               Insights from molecular dynamics simulations. Phys Rev B. 2016 Nov 29;94(18):184311. 
#           
#       If std_temp >> sigma_T, the thermostat is too light (SMASS too small). The temperature fluctuates excessively.
#       If std_temp << sigma_T, the thermostat is too heavy (SMASS too large). The temperature is too sluggish.
#
#       Remember though, this is only valid after we have assessed the ACF data.
#
#   2. Temperature drift (slope_temp):
#       
#           T(step) ~ slope_temp * step + intercept
#
#       Linear regression of temperature vs. step gives the systematic drift of temperature.
#       Converted to K/ps via:
#           
#           slope_temp = slope_step * 1e3 / POTIM_fs
#
#       Ideally, slope_temp should be close to zero. A significant slope indicates that SMASS may need adjustment.
#
# Interpretation guidelines:
#       std_temp >> sigma_T, Thermostat too light, consider increasing SMASS
#       std_temp << sigma_T, Thermostat too heavy, consider decreasing SMASS
#       |slope_temp| large, System is drifting in temperature, may indicate SMASS or timestep issues
#
# Note: SMASS should be chosen to balance realistic fluctuations with smooth temperature control,
#       generally resulting in ~1-2% relative temperature fluctuations.

if temperatures:
    # --- Steady-state temperature statistics ---
    avg_temp = np.mean(temperatures_steady)
    std_temp = np.std(temperatures_steady, ddof=1)

    coeffs_temp = np.polyfit(steps_temp_steady, temperatures_steady, 1)
    slope_temp = coeffs_temp[0] * 1e3 / POTIM_fs  # K/ps


    # Expected statistical fluctuation in a canonical ensemble
    expected_sigma = np.sqrt( (2 * avg_temp**2) /(3 * N_atoms))
 
    print(f"\n--- Temperature Analysis (SMASS) [last 50% only] ---")
    print(f"Avg Temp = {avg_temp:.2f} K")
    print(f"Std = {std_temp:.2f} K")
    print(f"Expected fluctuation (sigma_T) = {expected_sigma:.2f} K")
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


# ------------------------------------------------------------------
# Suggestions Explanation
# ------------------------------------------------------------------
#
# This chunck provides guidance for adjusting SMASS based on observed
# temperature fluctuations compared to the expected theoretical fluctuations.
#
# 1. Compute expected fluctuation:
#      sigma_T = sqrt (( 2 * Avg Temp^2 ) / (3 * N_atoms))
#     This is the statistical fluctuation for N atoms in the canonical ensemble.
#
# 2. Compare observed standard deviation (std_temp) to sigma_T:
#       If std_temp > 2 * sigma_T:
#           Temperature fluctuates more than expected
#           Thermostat is too "light" (SMASS too small), reacts too quickly
#           Suggest increasing SMASS proportionally to the severity of the overshoot
#
#       If std_temp < 0.5 * sigma_T:
#           Temperature fluctuates less than expected
#           Thermostat is too "heavy" (SMASS too large), reacts too slowly
#           Suggest decreasing SMASS proportionally to the severity of the undershoot
#
#       If 0.5 * sigma_T <= std_temp <= 2 * sigma_T:
#           Temperature fluctuations are within expected range
#           SMASS is considered appropriate
#
# 3. Severity factor:
#       Quantifies how far the observed fluctuation is from the expected value
#       Used to scale suggested SMASS up or down
#
# 4. Notes:
#       These suggestions assume the simulation is already at steady state
#           (i.e., the last half of the temperature trajectory is analyzed)
#       These suggestions also assume that you have validated the ACF data.
#       SMASS tuning is aimed at achieving a critically damped thermostat response
#           where temperature autocorrelation decays smoothly and N_eff is sufficient
#
# ------------------------------------------------------------------
# Severity Factor
# ------------------------------------------------------------------
#
# The severity factor quantifies "how far the observed temperature fluctuations
#   deviate from the expected theoretical fluctuations". It is used to scale
#   the suggested SMASS adjustment proportionally. The severity factor is the
#   percentage of deviation from the expected fluctuations. 
#
# Definition:
#       If temperature fluctuates too much (std_temp > expected_sigma):
#           severity = (std_temp - expected_sigma) / expected_sigma
#           Larger severity, greater overshoot, increase SMASS more
#
#       If temperature fluctuates too little (std_temp < expected_sigma):
#           severity = abs(std_temp - expected_sigma) / expected_sigma
#           Larger severity, thermostat too heavy, decrease SMASS more
#
# Interpretation:
#       severity = 0, temperature fluctuations are close to expected, small or no adjustment
#       severity = 1, fluctuations are twice the expected deviation, substantial adjustment
#
# Application:
#       For excessive fluctuations:
#         suggested_SMASS = SMASS * (1 + severity)
#       For too stable fluctuations:
#         suggested_SMASS = SMASS / (1 + severity)
#
# Notes:
#   This is not exact, it is a best guess.
# ------------------------------------------------------------------

print("\n------------- Suggestions -------------")
if temperatures_steady is not None:
    # -------------------------------
    # Base severity from std deviation
    # -------------------------------
    if std_temp > expected_sigma:  # fluctuations too high
        severity = (std_temp - expected_sigma) / expected_sigma
        direction = "increase"
    elif std_temp < expected_sigma:  # fluctuations too low
        severity = np.abs(std_temp - expected_sigma) / expected_sigma
        direction = "decrease"
    else:
        severity = 0
        direction = None


    # -------------------------------
    # Compute suggested SMASS
    # -------------------------------
    if direction == "increase":
        suggested_SMASS = SMASS * severity
        print(f"Temperature fluctuates too much (std={std_temp:.2f} K, expected={expected_sigma:.2f} K).")
        print(f"{damping_note}")
        print(f"Consider increasing SMASS to ~{suggested_SMASS:.2f}")
    elif direction == "decrease":
        suggested_SMASS = SMASS / severity
        print(f"Temperature too stable (std={std_temp:.2f} K, expected={expected_sigma:.2f} K).")
        print(f"{damping_note}")
        print(f"Consider decreasing SMASS to ~{suggested_SMASS:.2f}")
    else:
        print("Temperature fluctuations within expected range. SMASS is appropriate.")
        print(f"{damping_note}")
    
# ------------------------------------------------
# Plotting, Desciptions, and Interpretation Guide
# ------------------------------------------------
fig, axs = plt.subplots(3, 1, figsize=(12, 6), sharex=False)

# 1. Temperature Autocorrelation Function
#    Quantifies how temperature values are correlated over time.
#    C(0) = 1 and C(tau), 0 as lag tau increases.
#    The lag where C(tau) < 0.1 gives the correlation time (tau_corr).
#    Effective number of independent samples: N_eff = N_steps / tau_corr
#        Ensures statistical reliability of mean and std temperature.
#        N_eff < 30, insufficient independent data; consider longer simulation.
#    Short correlation time, thermostat too light (fast fluctuations)
#    Long correlation time, thermostat too heavy (slow fluctuations)
# NOTE: You do NOT want oscillation on this plot. You want to see a smooth decay
#   to the dashed line. If you see oscillation your system is underdamped or SMASS
#   is extremly large. If you see very very slow decay, your system is over damped.
#   You should think about this plot like a damped oscillator. 
#   You want critical dampening. 
if temperatures:
    c_temp = autocorr(temperatures_steady)
    axs[0].plot(c_temp, marker='o', color='purple')
    axs[0].axhline(0.1, color='r', linestyle='--', label="Corr < 0.1 threshold")
    axs[0].set_ylabel("Autocorrelation")
    axs[0].set_xlabel("Lag")
    axs[0].set_title("Temperature Autocorrelation Function")
    axs[0].grid(True)
    axs[0].legend()
else:
    axs[0].text(0.5, 0.5, "No temperature data found", ha='center', va='center')

# 2. Temperature vs MD Step
#    Displays the instantaneous temperature over the trajectory.
#    The black dashed line represents the average steady-state temperature.
#    The gray shaded region represents the expected fluctuation range (+-sigma_T) for N atoms:
#    Observing temperature mostly within the gray band indicates the thermostat is correctly controlling T.
#    Large deviations may indicate SMASS is too small (thermostat reacts too fast) or too large (too sluggish).
#    Slope of the temperature (linear fit) shows systematic drift:
#          |slope| < 0.01 K/ps, negligible drift
#          0.01 <= |slope| <= 0.05 K/ps, moderate drift
#          |slope| > 0.05 K/ps, significant drift; consider adjusting SMASS
if temperatures:
    axs[1].plot(steps_temp, temperatures, marker='o', color='r', label="Temperature")
    axs[1].axhline(avg_temp, color='k', linestyle='--', label=f"Avg = {avg_temp:.2f} K")
    axs[1].fill_between(steps_temp_steady, avg_temp - expected_sigma, avg_temp + expected_sigma, color='gray', alpha=0.2, label='Expected sig_T')
    axs[1].set_ylabel("Temperature (K)")
    axs[1].set_title("Temperature vs MD Step")
    axs[1].grid(True)
    axs[1].legend()
else:
    axs[1].text(0.5, 0.5, "No temperature data found", ha='center', va='center')

# 3. Total Energy vs MD Step
#    Shows the evolution of the total energy throughout the simulation.
#    The black dashed line indicates the average energy.
#
if total_energies:
    axs[2].plot(steps_energy, total_energies, marker='s', color='b', label="Total Energy")
    axs[2].set_ylabel("Total Energy (eV)")
    axs[2].set_title("Energy vs MD Step")
    axs[2].grid(True)
    axs[2].legend()
else:
    axs[2].text(0.5, 0.5, "No energy data found", ha='center', va='center')


plt.tight_layout()


# 4. Steady-State Temperature Distribution Histogram
#    Compares the observed temperature distribution (last half of simulation) to the expected Gaussian (sigma_T).
#    Histogram shows actual fluctuations; dashed line shows ideal thermal distribution.
#    Good agreement indicates SMASS yields physically realistic temperature sampling.
#    Large deviations suggest the thermostat is not properly sampling the canonical ensemble.
#
if temperatures:
    plt.figure(figsize=(6,4))
    plt.hist(temperatures_steady, bins=30, density=True, alpha=0.6, color='red', label='Observed T')
    x = np.linspace(min(temperatures_steady), max(temperatures_steady), 200)
    from scipy.stats import norm
    y = norm.pdf(x, loc=avg_temp, scale=expected_sigma)
    plt.plot(x, y, 'k--', lw=2, label='Expected Gaussian (sig_T)')
    plt.xlabel("Temperature (K)")
    plt.ylabel("Probability density")
    plt.title("Steady-State Temperature Distribution")
    plt.legend()
    plt.grid(True)
    plt.show()
