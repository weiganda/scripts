#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(description="Calculate Sampling Rate and AIMD Timestep from frequency.")
    parser.add_argument('-f', '--frequency', type=float, required=True,
                        help='Desired frequency in Hz (e.g., 3.33e15)')
    parser.add_argument('-a', '--aimd-timestep', type=float, default=0.5e-15,
                        help='AIMD timestep in seconds (default: 0.5e-15)')
    args = parser.parse_args()

    # Constants
    AIMD_TIMESTEP = args.aimd_timestep
    frequency = args.frequency

    # Calculations
    nyquist_correction = 2 * frequency
    freq_to_s = (1 / nyquist_correction)
    sampling_rate = freq_to_s / AIMD_TIMESTEP

    # Output
    print(f"AIMD Timestep: {AIMD_TIMESTEP:.2e} s")
    print(f"Sampling Rate: {sampling_rate:.2f}")

if __name__ == "__main__":
    main()

