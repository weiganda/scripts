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
# Read energies, temperatures, and volumes from OUTCAR
# -------------------------
temperatures = []
total_energies = []

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

# -------------------------
# Steps
# -------------------------
steps_temp = list(range(1, len(temperatures) + 1))
steps_energy = list(range(1, len(total_energies) + 1))

# -----------------------------
# Simulation summary
# -----------------------------
print("\n---------- Simulation Parameters ----------")
print(f"Number of atoms (N_atoms) = {N_atoms}")
print(f"POTIM (timestep) = {POTIM_fs:.3f} fs")
print("-------------------------------------------\n")

# -----------------------------
# Nyquist Frequency
# -----------------------------
dt_ps = POTIM_fs * 1e-3
nyquist_Hz = 1e12 / (2 * dt_ps)
print("\n------------ Nyquist Frequency ------------")
print(f"Nyquist frequency = {nyquist_Hz:.2e} Hz")
print("-------------------------------------------\n")



# ------------------------------------------
# Autocorrelation and Effective Sample Size
# ------------------------------------------
# This calculation estimates how many samples are present in the MD trajectory.
# Consecutive MD steps are correlated, meaning that not every temperature or energy value 
# represents a truly independent observation. Using all steps without accounting for correlation 
# would overestimate the statistical reliability of averages and fluctuations.
#
# The procedure is as follows:
#
# 1. Compute the autocorrelation function of the property of interest (e.g., temperature):
#       C(τ) = < (x_t - <x>) * (x_{t+τ} - <x>) > / < (x_t - <x>)^2 >
#    where τ is the lag in steps, and <> denotes the mean over time.
#
# 2. Determine the correlation time, tau_corr:
#       tau_corr = first lag τ where C(τ) < 1/e
#    This represents the typical number of steps over which the system retains memory of its previous state.
#
# 3. Compute the effective number of independent samples:
#       N_eff = N_total / tau_corr
#    This gives a more accurate count of independent data points for statistical analysis.
#
# 4. Interpretation:
#       - N_eff >> 30, sufficient data for reliable statistics
#       - N_eff < 30, insufficient independent samples; simulation may need to be longer
#       - N_eff >> 30 is considered very reliable statistics, below this is unreliable. 
#
# This analysis ensures that computed averages, standard deviations, or drift calculations
# are based on statistically meaningful data rather than raw step counts.
#
# Estimate correlation time for temperature
def autocorr(x):
    x = np.array(x)
    n = len(x)
    x_mean = np.mean(x)
    var = np.var(x)
    c = np.correlate(x - x_mean, x - x_mean, mode='full') / (var * n)
    return c[n-1:]

c_temp = autocorr(temperatures)
tau_corr = np.where(c_temp < 0.1)[0][0]  # first index where autocorr < 0.1
N_eff = len(temperatures) / tau_corr
print("\n--------------- Sample Size ---------------")
print(f"Estimated correlation time = {tau_corr} steps, N_eff = {N_eff:.1f} independent samples")
if N_eff < 30:
    print("WARNING: Not enough independent samples for reliable statistics.")
    print("This can be ignored if the plot shows steady-state energy.")
print("-------------------------------------------\n")




# -----------------------------
# Drift calculation for Energy
# -----------------------------
# The drift calculation will tell us how much the energy is changing per atom, ev, picosecond.
# There are two different drift calculations that are happening here. 
#   1. The first, drift_basic is determining the drift from the final and initial step energy values. 
#      
#           Drift (eV/atom/ps) = ΔE / (N * POTIM * 1e-3)   
#       
#       The demonimatior is POTIM (fs) * 1e-3 to get it in terms of ps.
#
# The drift_slope calculation gives an alternative measure of the energy drift per atom, eV, picosecond. 
# Instead of just using the first and last energy values, it fits a straight line to the last half of 
#   the energy data. The last half is used to identify the steady-state of the energy data. 
# 
#   2. drift_slope is determined by performing a linear regression on total energy vs. step number:
#       
#           E_tot(step) ≈ slope * step + intercept
#       
#       The slope from this fit represents the average energy change per MD step.
#
#       The total energy change over the simulation is then:
#           
#           ΔE = slope * (N_steps - 1)
#
#       Finally, the drift per atom per ps is calculated as:
#           
#           Drift (eV/atom/ps) = ΔE / (N * POTIM * 1e-3)
#
#
# The drift calculation (when run on NVE) will help us determine whether or not the POTIM value that we have choosen is
#       appropriate for the material. A general rule of thumb is as follows:
#           1. Great: Drift < 10e-5
#           2. OK:    Drift 10e-5 - 10e-4
#           3. Nope:  Drift > 10e-4
#
# If the drift is in the "Nope" range, than the POTIM is too large for the material. Ideally, your POTIM value should be 
#       able to capture the highest-frequency vibration of the material.

#   drift_slope is 
if total_energies:
    N_steps = len(total_energies)
    half_index = N_steps // 2  # start of steady-state region

    # Basic drift (still uses first and last of full trajectory)
    deltaE_basic = total_energies[-1] - total_energies[0]

    # Slope-based drift using only the last half
    ss_steps = steps_energy[half_index:]
    ss_energies = total_energies[half_index:]
    coeffs = np.polyfit(ss_steps, ss_energies, 1)  # slope, intercept
    slope = coeffs[0]  # eV per step
    deltaE_slope = slope * (len(ss_steps) - 1)

    # Total simulation time (full trajectory)
    t_total_ps = POTIM_fs * N_steps * 1e-3  # fs -> ps

    # Drift calculations
    drift_basic = abs(deltaE_basic / (N_atoms * t_total_ps))
    drift_slope = abs(deltaE_slope / (N_atoms * t_total_ps))
    
    print(f"\n--------------- Drift Energy --------------")
    print(f"ΔE (start to end) = {deltaE_basic:.6f} eV, drift = {drift_basic:.2e} eV/atom/ps")
    print(f"Slope-based ΔE = {deltaE_slope:.6f} eV, drift = {drift_slope:.2e} eV/atom/ps")


    # Evaluate drift quality
    if drift_slope < 1e-5:
        print("POTIM is Great")
    elif 1e-5 <= drift_slope <= 1e-4:
        print("POTIM is Good")
    else:  # drift_slope > 1e-4
        print("ALERT: POTIM is too large!!!")

    print(f"------------------------------------------\n")


# -----------------------------
# Energy Drift Linearity (R²)
# -----------------------------
# This section evaluates how well a straight line fits the steady-state (last half) 
#   of the energy trajectory.
# While the slope of the energy vs. step gives a measure of drift (eV/atom/ps),
#   the R² value quantifies how linear the energy evolution is over time. 
# A linear energy trend indicates a consistent, predictable drift.
# Steps:
# 1. Fit a straight line to the total energy vs. MD step number:
#       E_tot(step) ≈ slope * step + intercept
#    where slope gives the average energy change per step.
#
# 2. Compute the predicted energies from the linear fit:
#       E_fit(step) = slope * step + intercept
#
# 3. Compute the residual sum of squares (SS_res) and total sum of squares (SS_tot):
#       SS_res = Σ (E_actual - E_fit)^2
#       SS_tot = Σ (E_actual - mean(E_actual))^2
#
# 4. Compute R²:
#       R² = 1 - SS_res / SS_tot
#    Interpretation:
#       R² ≈ 1 → energy drift is highly linear; slope is reliable
#       R² < ~0.8 → energy exhibits significant nonlinearity; slope may be misleading
#       R² < 0 → fit is worse than using the mean energy; indicates instability

if total_energies:
    N_steps = len(total_energies)
    half_index = N_steps // 2  # start of steady-state region

    # Use only last half for slope and R²
    ss_steps = steps_energy[half_index:]
    ss_energies = total_energies[half_index:]

    # Linear fit
    coeffs = np.polyfit(ss_steps, ss_energies, 1)  # slope, intercept
    slope = coeffs[0]
    intercept = coeffs[1]

    # Predicted energies from linear fit
    fit_energy = np.polyval(coeffs, ss_steps)

    # Residual and total sum of squares
    ss_res = np.sum((ss_energies - fit_energy)**2)
    ss_tot = np.sum((ss_energies - np.mean(ss_energies))**2)

    # Coefficient of determination
    r2_energy = 1 - ss_res / ss_tot

    # Compute slope-based drift
    deltaE_slope = slope * (len(ss_steps) - 1)
    t_total_ps = POTIM_fs * N_steps * 1e-3  # use full simulation time
    drift_slope = abs(deltaE_slope / (N_atoms * t_total_ps))

    # Print results
    print(f"\n-------- Energy Drift and Linearity (steady-state) -------")
    print(f"Slope-based ΔE = {deltaE_slope:.6f} eV, drift = {drift_slope:.2e} eV/atom/ps")
    print(f"Energy drift linearity (R²) = {r2_energy:.4f}")

    # Interpret R²
    if r2_energy < 0.8:
        print("WARNING: Energy drift is not linear — slope may be unreliable.")
    elif r2_energy < 0.95:
        print("Energy drift reasonably linear, slope is acceptable.")
    else:
        print("Energy drift highly linear; slope is reliable.")
    print(f"-----------------------------------------------------------\n")


# -----------------------------
# Suggestions (still in testing mode)
# -----------------------------
print("\n--------------- Suggestions ---------------")

# 1. POTIM / energy drift
if drift_slope < 1e-5:
    print("POTIM is Great — no adjustment needed.")
elif 1e-5 <= drift_slope <= 1e-4:
    print("POTIM is Good — no adjustment needed.")
else:  # drift_slope > 1e-4
    print("ALERT: POTIM is too large!!! Needs a smaller timestep.")

print(f"-------------------------------------------\n")

   
# -------------------------
# Plotting
# -------------------------
fig, ax = plt.subplots(1, 1, figsize=(8, 4), sharex=True)

# Total energy
if total_energies:
    ax.plot(steps_energy, total_energies, marker='s', color='b', label="Total Energy")
    ax.axhline(np.mean(total_energies), color='k', linestyle='--', label=f"Avg = {np.mean(total_energies):.4f} eV")
    ax.set_ylabel("Total Energy (eV)")
    ax.grid(True)
    ax.legend()
else:
    ax.text(0.5, 0.5, "No energy data found", ha='center', va='center')
plt.tight_layout()
plt.show()

