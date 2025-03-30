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
# LAST MODIFIED: March 29, 2025
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
from statsmodels.graphics.tsaplots import plot_acf
from statsmodels.tsa.stattools import acf
import scipy as sp
import csv

##################################################
# Define routines (?)
##################################################


#---------Tested and confirmed to work *Do NOT CHANGE 3/27/25
def normalize_vectors(vector1, vectors):

    # Normalization may not be the correct term for what this part
    #   of the script does. What this is doing is setting the first
    #   vector of the dipole moment as 0,0,0. Since we anticipate 
    #   using this on bulk materials that do not have an intrinsic
    #   dipole moment, this logically seems appropriate. If we were
    #   to apply this to materials that may naturally have a dipole 
    #   moment, this would need to be changed to accommodate that. 
    #   Another idea would be to get rid of this completely, but we 
    #   will wait and see if that is necessary. avw 3/27/35

    # Translation of the first vector to the origin (0,0,0)
    vector1_float = np.array(vector1, dtype=float)    
    vectors_float = np.array(vectors, dtype=float)

    # Set the first vector to the origin and adjust the others accordingly
    vector1_norm = vector1_float - vector1_float
    vectors_norm = vectors_float - vector1_float
    
    # Combine vector1_norm with vectors_norm to include it as the first line
    norm_vectors = [vector1_norm] + list(vectors_norm)

    # Write to CSV
    with open("dipole_norm.csv", "a", newline='') as csvfile:
        writer = csv.writer(csvfile)
        for vec in norm_vectors:
            writer.writerow([float(x) for x in vec])  # Ensure each element is written as float

    return norm_vectors



#--------- Tested ---- I think this works but I am not really
#   sure what exactly I am supposed to be looking for in the
#   data. 3/28/25
def dipole_correlation_function(vectors):
    
    # Initialize arrays
    acf_tot = []
    dipole_correlation = []

    # Ensure that the vectors are in array format
    vectors = np.array(vectors)
    
    # Extract data for x, y, and z dimensions
    x_values = vectors[:, 0]
    y_values = vectors[:, 1]
    z_values = vectors[:, 2]

    length = len(x_values) - 1 

    acf_x = acf(x_values, nlags=length, fft=True)
    acf_y = acf(y_values, nlags=length, fft=True)  
    acf_z = acf(z_values, nlags=length, fft=True) 
    # Combine ACF values into a single array
    
    for i in range(length + 1):
        acf_vector = [acf_x[i], acf_y[i], acf_z[i]]
        acf_tot.append(acf_vector)
        acf_avgs = np.mean(acf_vector) 
        dipole_correlation.append(acf_avgs) 

    # Save the average ACF values to a CSV file
    avg_acf_df = pd.DataFrame(dipole_correlation, columns=["ACF_tot"])
    avg_acf_df.to_csv("dipole_correlation_avg.csv", index=False)

    # Save xyz correlation data to CSV
    correlation_df = pd.DataFrame(acf_tot, columns=["ACF_x", "ACF_y", "ACF_z"])
    output_csv = "dipole_correlation_xyz.csv"
    correlation_df.to_csv(output_csv, index=False)

    # Plotting X component (input and correlation)
    plt.figure(figsize=(12, 8))

    plt.subplot(2, 2, 1)
    plt.plot(x_values, label="Input X", color="red", alpha=0.7)
    plt.plot(acf_x, label="Dipole Correlation X", color="darkred", linestyle="--")
    plt.xlabel("Index/Lag")
    plt.ylabel("Value/Correlation")
    plt.title("X Component - Input and Dipole Correlation")
    plt.legend()
    plt.grid(True)

    # Plotting Y component (input and correlation)
    plt.subplot(2, 2, 2)
    plt.plot(y_values, label="Input Y", color="green", alpha=0.7)
    plt.plot(acf_y, label="Dipole Correlation Y", color="darkgreen", linestyle="--")
    plt.xlabel("Index/Lag")
    plt.ylabel("Value/Correlation")
    plt.title("Y Component - Input and Dipole Correlation")
    plt.legend()
    plt.grid(True)

    # Plotting Z component (input and correlation)
    plt.subplot(2, 2, 3)
    plt.plot(z_values, label="Input Z", color="blue", alpha=0.7)
    plt.plot(acf_z, label="Dipole Correlation Z", color="darkblue", linestyle="--")
    plt.xlabel("Index/Lag")
    plt.ylabel("Value/Correlation")
    plt.title("Z Component - Input and Dipole Correlation")
    plt.legend()
    plt.grid(True)

    # Plotting Total Dipole Correlation
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

    return dipole_correlation

#--------- Tested and confirmed 3/29/25
def numerical_derivative(signal, steps):
    derivative = [(signal[i+1] - signal[i]) for i in range(steps - 1)]
    negative_derivative = [-d for d in derivative]
    return negative_derivative

#----------------FFT--- currently testing 3/29/25

def perform_fft(signal, sampling_rate):

    num_samples = int(len(signal))

    # Perform the FFT
    fft_result = np.fft.fft(signal)
    fft_freqs = np.fft.fftfreq(num_samples, d=(1/sampling_rate))

    # Extract the real and imaginary parts
    fft_real = np.real(fft_result)
    fft_imag = np.imag(fft_result)

    # Only keep the positive frequencies
    positive_freqs = fft_freqs[:num_samples // 2]
    real_fft = fft_real[:num_samples // 2]
    imag_fft = fft_imag[:num_samples // 2]

    # Create a figure with two subplots (2 rows, 1 column)
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    # Plot the real part of the FFT on the first subplot
    axes[0].plot(positive_freqs / 1e9, real_fft, 'o-', label="Real Part")
    axes[0].set_title('Real Part of FFT')
    axes[0].set_xlabel('Frequency (GHz)')
    axes[0].set_ylabel('Amplitude')
    axes[0].set_xlim(0, sampling_rate / 2 / 1e9)  # Limit x-axis to Nyquist frequency
    axes[0].grid(True)
    axes[0].legend()

    # Plot the imaginary part of the FFT on the second subplot
    axes[1].plot(positive_freqs / 1e9, imag_fft, 'o-', label="Imaginary Part", color='orange')
    axes[1].set_title('Imaginary Part of FFT')
    axes[1].set_xlabel('Frequency (GHz)')
    axes[1].set_ylabel('Amplitude')
    axes[1].set_xlim(0, sampling_rate / 2 / 1e9)  # Limit x-axis to Nyquist frequency
    axes[1].grid(True)
    axes[1].legend()

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Save the entire figure as a PNG file
    plt.savefig('fft_parts_combined.png', dpi=300, bbox_inches='tight')

    # Display the plot
    plt.show()
    
    return real_fft, imag_fft, positive_freqs








##### We are here #####
### IDK we may need to rethink this equation

def apply_constant(real_fft, imag_fft, vol):
    # c.g.s. units were originally used for this

    boltzmann_constant = const.k
    
    #Temperature is Kelvin. 
    temp = 300

    # Boltzmann in cgs erg/K
    boltzmann = 1.380649e-16

    # Change volume from Angstom to cm^3
    volume = vol * 1e-24

    # Define pi
    pi = 3.14159265358979323846

    # erg to Hz - Unsure if this is necessary.
    # erg_to_hz = 1.509190311676e+26

    permittivity_imag = ((4*pi*boltzmann)/volume)*(imag_fft)
    permittivity_real = ((4*pi*boltzmann)/volume)*(real_fft)

    return permittivity_imag, permittivity_real


def plot_results(num_steps, dipole_correlation_result, negative_phi_derivative, freqs, real_part, imaginary_part):
    num_steps = len(dipole_correlation_results)

    with open("permittivity.csv", "a", newline='') as csvfile:
        writer = csv.writer(csvfile)
        for i in range(len(freqs)):
            writer.writerow([freqs[i], real_part[i], imaginary_part[i]])

    plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(dipole_correlation_result)
    plt.title("Dipole Correlation")

    phi_deriv =  [-d for d in negative_phi_derivative]
    plt.subplot(2,2,4)
    plt.plot(phi_deriv)
    plt.title("Phi Derivative")

    plt.subplot(2, 2, 2)
    plt.plot(negative_phi_derivative)
    plt.title("Neg Phi Derivative")
    
    plt.savefig('derivs.png', dpi=300, bbox_inches='tight')



    # plot the real and imaginary parts of the fourier transform
    plt.figure()
    plt.plot(freqs, real_part, label='real part')
    plt.plot(freqs, imaginary_part, label='imaginary part')
    plt.legend()
    plt.title('real and imaginary parts of fourier transform')
    plt.xlabel('frequency')
    plt.ylabel('amplitude')
    
    # set x-axis limit to start at 0
    plt.xlim(left=0)
    plt.xlim(right=20e9)
    plt.show()

    plt.tight_layout()
    plt.savefig('permitt.png', dpi=300, bbox_inches='tight')







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

    print(f"Filename: {filename}")
    print(f"Sampling Frequency: {sampling_freq}")


    # Read in volume (3/28/25)
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
    dipole_correlation = dipole_correlation_function(norm_vectors)     
    dipole_correlation_result = dipole_correlation     

    # Information for the derivative
    time_steps = len(dipole_correlation_result)
    print(f"Num steps of derivative: {time_steps}")
    num_steps = np.arange(0, time_steps)
    PHI_derivative = numerical_derivative(dipole_correlation_result, len(num_steps))
    negative_PHI_derivative = -PHI_derivative

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
    #   sampling_rate.

    aimd_time_step = 2e-15
    sampling_rate = sampling_step * aimd_time_step
    print(f"Sampling Rate: {sampling_rate}")    
    real_fft, imag_fft, freqs = perform_fft(negative_PHI_derivative, sampling_rate)

    # Apply the constants and calculate permittivity
    permittivity_imag, permittivity_real = apply_constant(real_fft, imag_fft, boltzmann_constant, vol)

    plot_results(num_steps, dipole_correlation_result, negative_PHI_derivative, freqs, permittivity_real, permittivity_imag)

