import numpy as np
import matplotlib.pyplot as plt

# Function to calculate the volume of the unit cell
def calculate_unit_cell_volume(a, b, c):
    # Cross product of b and c
    cross_product = np.cross(b, c)
    # Dot product of a with the cross product of b and c
    volume = np.abs(np.dot(a, cross_product))
    return abs(volume)

# Function to read and process the XDATCAR file
def process_file(filename):
    volumes = []  # List to store volumes for each repetition
    with open(filename, 'r') as file:
        # Skip the first two lines
        
        next(file)
        next(file)
        
        while True:
            try:
                # Read the next three lines for lattice vectors
                a = np.array([float(x.strip()) for x in next(file).split()])
                b = np.array([float(x.strip()) for x in next(file).split()])
                c = np.array([float(x.strip()) for x in next(file).split()])
               
               # Skip line 6
                next(file)
                
                # Read line 7 and sum the values
                line_7_values = [float(x) for x in next(file).split()]
                sum_line_7 = sum(line_7_values)
                
                # Skip the number of lines based on the sum of line 7 + 1
                for _ in range(int(sum_line_7) + 1):
                    next(file)

                # Calculate the volume of the unit cell for the current repetition
                volume = calculate_unit_cell_volume(a, b, c)
                volumes.append(volume)

                next(file)
                next(file)

            except StopIteration:
                # End of file reached
                break
    
    return volumes

# Read and process the file
filename = 'XDATCAR'  # Change this to the path of your XDATCAR file
volumes = process_file(filename)

# Plot the volumes
plt.figure(figsize=(10, 6))
plt.plot(volumes, label='Unit Cell Volume', color='b', marker='o')
plt.xlabel('Repetition Number')
plt.ylabel('Volume (cubic units)')
plt.title('Unit Cell Volume for Each Repetition')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
