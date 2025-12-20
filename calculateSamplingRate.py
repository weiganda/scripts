#!/usr/bin/env python3

"""
Created By: Alysse Weigand
Last Updated: 12/20/2025

AIMD Sampling Rate and Timestep Helper

Purpose:
This script helps determine the appropriate AIMD timestep (POTIM) and sampling rate 
for analyzing vibrational frequencies in molecular dynamics simulations. 

By providing a target frequency, the script calculates how often you should sample 
your trajectory (in MD steps) to accurately capture the oscillations without 
violating the Nyquist criterion.

Nyquist Criterion:

- To accurately capture a signal of frequency f, you must sample it at least twice per cycle.
- In MD, this means that if your system has vibrations at frequency f, your sampling interval 
  (MD steps) must be short enough that you take at least 2 samples per period.
- Sampling less frequently than this leads to "aliasing," where high-frequency oscillations 
  appear as slower, incorrect frequencies in your analysis (FFT).
- Formula: Minimum sampling frequency = 2 * f
- Consequence: The AIMD timestep and trajectory sampling rate should be chosen with this in mind.

Notes:
- Nyquist theorem: To properly capture a frequency f, you must sample at least 2*f.
- AIMD timestep (POTIM) is the time between MD steps. Smaller timestep → more accurate, 
  but more computationally expensive.
- Sampling rate tells you how many MD steps correspond to one “observation” at the 
  target frequency.
- This tool is intended for post-processing NVT trajectories to guide FFT analysis.
- Always ensure your timestep is small enough to resolve the fastest vibrations in your system.
"""

import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Calculate Sampling Rate and AIMD Timestep from target frequency."
    )
    parser.add_argument('-f', '--frequency', type=float, required=True,
                        help='Target vibrational frequency in Hz (e.g., 3.33e15)')
    parser.add_argument('-a', '--aimd-timestep', type=float, default=0.5e-15,
                        help='AIMD timestep in seconds (default: 0.5e-15 s)')
    args = parser.parse_args()

    # Constants from user input
    AIMD_TIMESTEP = args.aimd_timestep  # seconds per MD step
    frequency = args.frequency          # target vibrational frequency in Hz

    # -----------------------------------
    # Sampling Rate Calculation
    # -----------------------------------
    # Nyquist criterion: sample at least twice per period of the highest frequency
    nyquist_correction = 2 * frequency          # minimum sampling frequency (Hz)
    freq_to_s = 1 / nyquist_correction          # period in seconds between samples
    sampling_rate = freq_to_s / AIMD_TIMESTEP   # number of MD steps per sample

    # -----------------------------------
    # Output Results
    # -----------------------------------
    print(f"AIMD Timestep: {AIMD_TIMESTEP:.2e} s")
    print(f"Target Frequency: {frequency:.2e} Hz")
    print(f"Suggested Sampling Rate: 1 sample every {sampling_rate:.2f} MD steps")

if __name__ == "__main__":
    main()
