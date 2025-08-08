import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
HARMONIC_CONSTANT = 0.35  # Scaling factor for harmonics
MAX_ITERATIONS = 10

# Step 1: Initialize Lattice with Inner-Outward Growth
def initialize_lattice(binary_data, harmonic_constant=HARMONIC_CONSTANT):
    normalized_data = binary_data / 255.0  # Normalize binary data to [0, 1]
    
    # Calculate lattice size
    data_length = len(normalized_data)
    lattice_size = int(np.ceil(np.cbrt(data_length)))
    lattice = np.zeros((lattice_size, lattice_size, lattice_size), dtype=np.float64)

    # Prioritize inner regions for mapping
    center = lattice_size // 2
    offset = 0
    for idx, value in enumerate(normalized_data):
        x, y, z = (center + offset) % lattice_size, (center - offset) % lattice_size, (center + offset) % lattice_size
        lattice[x, y, z] += value * harmonic_constant
        offset += 1
    return lattice, data_length

# Step 2: Retrieve Data from Lattice
def retrieve_from_lattice(lattice, harmonic_constant=HARMONIC_CONSTANT, data_length=None):
    flattened_data = lattice.flatten() / harmonic_constant  # Scale back by harmonic constant
    binary_data = np.round(flattened_data[:data_length] * 255).astype(np.uint8)  # Crop to original size
    return binary_data

# Step 3: Apply Feedback Correction
def feedback_correction(lattice, binary_data, harmonic_constant=HARMONIC_CONSTANT):
    retrieved_data = retrieve_from_lattice(lattice, harmonic_constant, data_length=len(binary_data))
    error = (binary_data - retrieved_data) / 255.0  # Normalize the error
    
    for idx, value in enumerate(error):
        x, y, z = idx % lattice.shape[0], (idx // lattice.shape[0]) % lattice.shape[1], idx // (lattice.shape[0] ** 2)
        lattice[x, y, z] += value * harmonic_constant  # Correct the lattice harmonically
    return lattice

# Step 4: Apply Reflective Gain
def apply_reflective_gain(lattice, gain_factor=0.1):
    center = lattice.shape[0] // 2
    for x in range(lattice.shape[0]):
        for y in range(lattice.shape[1]):
            for z in range(lattice.shape[2]):
                distance = np.sqrt((x - center)**2 + (y - center)**2 + (z - center)**2)
                lattice[x, y, z] += gain_factor / (1 + distance)
    return lattice

# Visualization of the 3D Lattice
def visualize_lattice(lattice, iteration):
    x, y, z = np.nonzero(lattice)  # Get non-zero positions
    values = lattice[x, y, z]  # Corresponding harmonic values
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x, y, z, c=values, cmap='viridis', s=10)

    # Add labels and title
    ax.set_title(f"3D Lattice Visualization of Harmonics - Iteration {iteration}", fontsize=16)
    ax.set_xlabel("X-axis", fontsize=12)
    ax.set_ylabel("Y-axis", fontsize=12)
    ax.set_zlabel("Z-axis", fontsize=12)
    plt.colorbar(scatter, ax=ax, label="Harmonic Values")
    plt.show()

# Main Execution
if __name__ == "__main__":
    # Load binary bits from your file
    with open(r'd:\\colecovision.rom', 'rb') as file:
        binary_data = np.frombuffer(file.read(), dtype=np.uint8)  # Read binary as bytes

    # Initialize the lattice
    lattice, data_length = initialize_lattice(binary_data)

    for iteration in range(1, MAX_ITERATIONS + 1):
        print(f"Iteration {iteration}")
        
        # Apply feedback correction
        lattice = feedback_correction(lattice, binary_data)
        
        # Apply reflective gain
        lattice = apply_reflective_gain(lattice, gain_factor=0.05)
        
        # Retrieve data
        retrieved_data = retrieve_from_lattice(lattice, data_length=data_length)

        # Visualize the lattice
        visualize_lattice(lattice, iteration)

        # Outputs: Compare Original and Retrieved Data
        print("Lattice Shape:", lattice.shape)
        print("Original Data (First 10 Bytes):", binary_data[:10])
        print("Retrieved Data (First 10 Bytes):", retrieved_data[:10])

        # Validate if original and retrieved data match
        data_matches = np.array_equal(binary_data, retrieved_data)
        print("Data matches:", data_matches)
        if data_matches:
            print("Data successfully recovered!")
            break
        else:
            print("Differences detected in data.")
