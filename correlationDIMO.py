#!/usr/bin/env python3

# PROGRAM:  correlationDIMO
# PURPOSE:  This program will calculate the time dependent correlation function
#           of the dipole moment data produced by VASP and olcao -dimo. It will 
#           then apply a numerical derivative and fast fourier transform of the
#           correlation function to produce the real and imaginary electric 
#           susceptibility and convert that data into the real and imaginary 
#           parts of the permittivity/dielectric function for low frequency
#           calculations.
#
# AUTHOR:   Alysse Weigand
# LAST MODIFIED: April 1, 2025
# USAGE:    Low frequency dielectric function calculation


##################################################
# Import necessary modules
##################################################

import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.constants as const
from scipy.signal import correlate
from scipy.signal import savgol_filter
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.tsa.stattools import acf
import scipy as sp
import csv

##################################################
# Define routines
##################################################


#---------Tested and confirmed to work *Do NOT CHANGE 3/27/25
def normalize_vectors(vector1, vectors):

    # Normalization may not be the correct term for what this part
    #   of the script does. What this is doing is setting the first
    #   vector of the dipole moment as 0,0,0. Since we anticipate 
    #   using this on bulk materials that do not have an intrinsic
    #   dipole moment, this logically seems appropriate. Although
    #   it does not appear to make a difference. If we were
    #   to apply this to materials that may naturally have a dipole 
    #   moment, this would need to be changed to accommodate that. 
    #   Another idea would be to get rid of this completely, but we 
    #   will wait and see if that is necessary. avw 3/27/35

    # Translation of the first vector to the origin (0,0,0)
    vector1_float = np.array(vector1, dtype=float)    
    vectors_float = np.array(vectors, dtype=float)

    # Set the first vector to the origin and adjust the others accordingly
    vector1_norm = vector1_float #- vector1_float # I do not think we need this
    vectors_norm = vectors_float #- vector1_float
    
    # Combine vector1_norm with vectors_norm to include it as the first line
    norm_vectors = [vector1_norm] + list(vectors_norm)

    # Write to CSV
    with open("dipole_norm.csv", "a", newline='') as csvfile:
        writer = csv.writer(csvfile)
        for vec in norm_vectors:
            writer.writerow([float(x) for x in vec])  # Ensure each element is written as float

    return norm_vectors



#--------- Tested ---- 3/28/25
def dipole_correlation_function(vectors):
    
    # The unit of the dipole are a_0e, bohr radii * charge. 

    # Initialize arrays
    acf_xyz = []
    dipole_correlation = []

    # Ensure that the vectors are in array format
    vectors = np.array(vectors)
    
    # Extract data for x, y, and z dimensions
    x_values = vectors[:, 0]
    y_values = vectors[:, 1]
    z_values = vectors[:, 2]

    length = len(x_values) - 1 

    acf_x = acf(x_values, nlags=length)
    acf_y = acf(y_values, nlags=length)  
    acf_z = acf(z_values, nlags=length) 

    # Calculate the leading term ⟨Mi(0) Mj(0)⟩ for each component
    leading_term_x = x_values[0]**2
    leading_term_y = y_values[0]**2
    leading_term_z = z_values[0]**2
    leading_term_avg = np.mean([leading_term_x, leading_term_y, leading_term_z])

    print(f"⟨Mx(0) Mx(0)⟩: {leading_term_x}")
    print(f"⟨My(0) My(0)⟩: {leading_term_y}")
    print(f"⟨Mz(0) Mz(0)⟩: {leading_term_z}")
    print(f"⟨M(0) M(0)⟩: {leading_term_avg}")


    # Combine ACF values into a single array
    
    for i in range(length + 1):
        acf_vector = [acf_x[i], acf_y[i], acf_z[i]]
        acf_xyz.append(acf_vector)
        acf_avgs = np.mean(acf_vector) 
        dipole_correlation.append(acf_avgs) 

    # Save the average ACF values to a CSV file
    avg_acf_df = pd.DataFrame(dipole_correlation, columns=["ACF_tot"])
    avg_acf_df.to_csv("dipole_correlation_avg.csv", index=False)

    # Save xyz correlation data to CSV
    correlation_df = pd.DataFrame(acf_xyz, columns=["ACF_x", "ACF_y", "ACF_z"])
    output_csv = "dipole_correlation_xyz.csv"
    correlation_df.to_csv(output_csv, index=False)

    # Plotting X component (input and correlation)
    plt.figure(figsize=(12, 8))

    # X Component Plot
    plt.subplot(2, 2, 1)
    ax1 = plt.gca()
    ax1.plot(x_values, label="Input X", color="red", alpha=0.7)
    ax1.set_xlabel("Index/Lag")
    ax1.set_ylabel("Input Value", color="red")
    ax1.tick_params(axis='y', labelcolor="red")
    ax1.set_title("X Component - Input and Dipole Correlation")
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(acf_x, label="Dipole Correlation X", color="darkred", linestyle="--")
    ax2.set_ylabel("Correlation", color="darkred")
    ax2.tick_params(axis='y', labelcolor="darkred")
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")

    # Y Component Plot
    plt.subplot(2, 2, 2)
    ax1 = plt.gca()
    ax1.plot(y_values, label="Input Y", color="green", alpha=0.7)
    ax1.set_xlabel("Index/Lag")
    ax1.set_ylabel("Input Value", color="green")
    ax1.tick_params(axis='y', labelcolor="green")
    ax1.set_title("Y Component - Input and Dipole Correlation")
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(acf_y, label="Dipole Correlation Y", color="darkgreen", linestyle="--")
    ax2.set_ylabel("Correlation", color="darkgreen")
    ax2.tick_params(axis='y', labelcolor="darkgreen")
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")

    # Z Component Plot
    plt.subplot(2, 2, 3)
    ax1 = plt.gca()
    ax1.plot(z_values, label="Input Z", color="blue", alpha=0.7)
    ax1.set_xlabel("Index/Lag")
    ax1.set_ylabel("Input Value", color="blue")
    ax1.tick_params(axis='y', labelcolor="blue")
    ax1.set_title("Z Component - Input and Dipole Correlation")
    ax1.grid(True)

    ax2 = ax1.twinx()
    ax2.plot(acf_z, label="Dipole Correlation Z", color="darkblue", linestyle="--")
    ax2.set_ylabel("Correlation", color="darkblue")
    ax2.tick_params(axis='y', labelcolor="darkblue")
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")

    # Total Dipole Correlation Plot
    plt.subplot(2, 2, 4)
    plt.plot(dipole_correlation, label="Total Dipole Correlation", color="purple")
    plt.xlabel("Lag")
    plt.ylabel("Correlation")
    plt.title("Total Dipole Correlation")
    plt.legend()
    plt.grid(True)

    # Adjust layout and save the figure
    plt.tight_layout()
    plt.savefig("dipole_correlation_plot.png")
    plt.show()

    return dipole_correlation, leading_term_avg

#--------- Tested and confirmed 3/29/25
def numerical_derivative(signal, steps, step_size):
    #derivative = [(signal[i+1] - signal[i]) / (step_size) for i in range(steps - 1)] 
    #derivative = [(signal[i+2] - signal[i]) / (2 * step_size) for i in range(steps - 2)]

    derivative = savgol_filter(signal, window_length=5, polyorder=3, deriv=1)

    return derivative

#----------------FFT--- currently testing 3/29/25 somehow the real and imag are flipped

def perform_fft(signal, sampling_interval):

    num_samples = int(len(signal))

    # Perform the FFT
    fft_result = np.fft.fft(signal)
    fft_freqs = np.fft.fftfreq(num_samples, d=(sampling_interval))

    # Extract the real and imaginary parts 
    fft_real = np.real(fft_result)
    fft_imag = np.imag(fft_result)

    # Only keep the positive frequencies
    positive_freqs = fft_freqs[:num_samples // 2]
    real_fft = fft_real[:num_samples // 2]
    imag_fft = - fft_imag[:num_samples // 2] # Phase information 180*

    # Create a figure with two subplots (2 rows, 1 column)
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    # Plot the real part of the FFT on the first subplot
    axes[0].plot(positive_freqs, real_fft, label="Real Part")
    axes[0].set_title('Real Part of FFT')
    axes[0].set_xlabel('Frequency (Hz)')
    axes[0].set_ylabel('Amplitude')
    axes[0].set_xlim(0, max(positive_freqs))  # Limit x-axis to Nyquist frequency
    axes[0].grid(True)
    axes[0].legend()

    # Plot the imaginary part of the FFT on the second subplot
    axes[1].plot(positive_freqs, imag_fft, label="Imaginary Part", color='orange')
    axes[1].set_title('Imaginary Part of FFT')
    axes[1].set_xlabel('Frequency (Hz)')
    axes[1].set_ylabel('Amplitude')
    axes[1].set_xlim(0, max(positive_freqs))  # Limit x-axis to Nyquist frequency
    axes[1].grid(True)
    axes[1].legend()

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Display the plot
    plt.show()
    
    return real_fft, imag_fft, positive_freqs


##### We are here #####
### IDK we may need to rethink this equation

def apply_constant(real_fft, imag_fft, vol, leading_term_avg):
    # The dipole moment is calculated in olcao in atomic units.The dipole moment 
    #   in olcao is a_0e or Bohr radii * electroncharge. Boltzmann const in a.u.
    #   is 1 Hartree/K. Temperature is Kelvin. 
    
    # Temperature is Kelvin. 
    temp = 300

    # The volume that is pulled in from OLCAO is in Angstroms, it 
    #   needs to be converted into Bohr Radii for the calculation
    #   Convert from Bohr radii to angstrom. 1a_o = 0.529A.
    #   Since it is cubed, 1 a_0^3 = 0.1481847435 A^3
    vol = vol / 0.1481847435

    # Boltzmann in au Ha/K * eV/Ha
    k_B = 3.1668119646e-6 #* 27.2114

    # Boltzamm in Hz/K
    #k_B = 2.083661912e10

    #Boltzmann in eV/K
    #k_B = 8.61733262e-5

    # Leading Term in a.u. 
    term = 4 * np.pi / (vol * k_B * temp)

    print(f"Term: {term}")

    # Apply leading term
    dipole_angst_imag = (leading_term_avg * imag_fft)
    dipole_angst_real = (leading_term_avg * real_fft)

    # So this is without the leading term 
    #dipole_angst_imag = imag_fft
    #dipole_angst_real = real_fft

    # This is leaving me with weird units
    permittivity_imag = 1 + (term * dipole_angst_imag)
    permittivity_real = 1 + (term * dipole_angst_real)

    return permittivity_imag, permittivity_real


def plot_results(num_steps, dipole_correlation_result, PHI_derivative, negative_phi_derivative, freqs, real_part, imaginary_part, sampling_rate):
    num_steps = len(dipole_correlation_results)

    with open("permittivity.csv", "a", newline='') as csvfile:
        writer = csv.writer(csvfile)
        for i in range(len(freqs)):
            writer.writerow([freqs[i], real_part[i], imaginary_part[i]])

    # Plot 1: Dipole Correlation
    plt.subplot(3, 1, 1)
    plt.plot(dipole_correlation_result)
    plt.title("Dipole Correlation")
    plt.xlabel('Time Step')
    plt.ylabel('Correlation')

    # Plot 2: Phi Derivative
    plt.subplot(3, 1, 2)
    plt.plot(PHI_derivative)
    plt.title("Phi Derivative")
    plt.xlabel('Time Step')
    plt.ylabel('Derivative')

    # Plot 3: Negative Phi Derivative
    plt.subplot(3, 1, 3)
    plt.plot(negative_phi_derivative)
    plt.title("Neg. Phi Derivative")
    plt.xlabel('Time Step')
    plt.ylabel('Derivative')

    # Adjust layout to avoid overlap
    plt.tight_layout()
    
    # Convert frequencies from Hz to eV
    freqs_eV = freqs #* 4.1357e-15

    # Plot the real part of the Fourier transform
    plt.figure(figsize=(8, 6))

    # Real part plot
    plt.subplot(2, 1, 1)
    plt.plot(freqs_eV, real_part, label='Real Part', color='blue')
    plt.axhline(y=0, color='black', linestyle='-', linewidth=0.8)  # Solid line at y=0
    plt.xlim(left=0)  # Set x-axis limit to start at 0    
    plt.title('Real Part of Permittivity')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Amplitude')
    plt.legend()

    # Imaginary part plot
    plt.subplot(2, 1, 2)
    plt.plot(freqs_eV, imaginary_part, label='Imaginary Part', color='red')
    plt.axhline(y=0, color='black', linestyle='-', linewidth=0.8)  # Solid line at y=0
    plt.xlim(left=0)  # Set x-axis limit to start at 0
    plt.title('Imaginary Part of Permittivity')
    plt.xlabel('Energy (eV)')
    plt.ylabel('Amplitude')
    plt.legend()

    # Adjust layout and save the figure
    plt.xlim(left=0)  # Set x-axis limit to start at 0
    plt.tight_layout()
    plt.savefig('permitt.png', dpi=300, bbox_inches='tight')
    plt.show()


# Read in the information from the csv file 3/29/25
def load_vectors_from_csv(filename):
    try:
        data = pd.read_csv(filename, header=None)
    except filenotfounderror:
        print("error: file not found.")
        sys.exit(1)
        
    vectors = data.values.tolist()
    
    for vector in vectors:
        if len(vector) != 3:
            print("error: each row should have exactly three values representing a 3d vector")
            sys.exit(1)
   
    return vectors





if __name__ == "__main__":
    if "-f" not in sys.argv or len(sys.argv) < 5:
        print("python script.py -f <filename> -s <sampling_step")
        sys.exit(1)
        
    # Get the filename from the command line arguments
    filename_index = sys.argv.index("-f") + 1
    filename = sys.argv[filename_index]

    # Get the sampling frequency from the command line arguments
    sampling_step_index = sys.argv.index("-s") + 1
    sampling_step = int(sys.argv[sampling_step_index])
    print(f"Sampling Frequency: {sampling_step}")

    # Read in volume. Volume has units of Angstrom. (3/28/25)
    with open('vol', 'r') as file:
        volume_data = file.read().strip()
        vol = float(volume_data)
    print(f"Volume of Cell: {vol}")

    # Read in vector information (3/28/25)  
    vector_list = load_vectors_from_csv(filename)
    vector1 = vector_list[0]
    vectors = vector_list[0:]

    # Initialize the array
    dipole_correlation_results = []

    # "Normalize the vectors"
    norm_vectors = normalize_vectors(vector1,vectors)

    # Calculate the Auto Correlation Function
    dipole_correlation, leading_term_avg = dipole_correlation_function(norm_vectors)     
    dipole_correlation_result = dipole_correlation     

    # Information for the derivative
    aimd_time_step = 0.5e-15 #need to automate
    sampling_rate = sampling_step * aimd_time_step
    print(f"Sampling Rate: {sampling_rate}")

    time_steps = len(dipole_correlation_result)
    print(f"Num steps of derivative: {time_steps}")

    PHI_derivative = numerical_derivative(dipole_correlation_result, time_steps, sampling_rate)
    
    # Assuming PHI_derivative is a list
    PHI_derivative = np.array(PHI_derivative)

    # Now, you can safely negate it
    negative_PHI_derivative = - PHI_derivative

    # Information for the fft. The sampling rate is 1 / time between samples
    #   aimd time step is 2e-15 and currently we are taking a sample
    #   every 500 time steps.
    # Notes about the Nyquist frequency:
    #   The Nyquist frequency is half of the sampling rate and 
    #   represents the highest frequency that can be accurately 
    #   captured when sampling a continuous signal. Meaning that
    #   a sampling_step of 500 (1THz) will capture up to 500 GHz.
    #   If you want to sample up to 6 THz, you would nee a 
    #   sampling_step of 42. Nyquist frequency is 1/2 the 
    #   sampling_rate (sr). Frequency "f" you want in Hz convert to cycles/sec
    #   sr = 1/f then sr/aimd_time_step. Then /2.
    real_fft, imag_fft, freqs = perform_fft(negative_PHI_derivative, sampling_rate)

    # Apply the constants and calculate permittivity
    permittivity_imag, permittivity_real = apply_constant(real_fft, imag_fft, vol, leading_term_avg)

    plot_results(time_steps, dipole_correlation_result, PHI_derivative, negative_PHI_derivative, freqs, permittivity_real, permittivity_imag, sampling_rate)

