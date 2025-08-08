import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
HARMONIC_CONSTANT = 0.35
TOLERANCE = 0.01
MAX_ITERATIONS = 10
REFLECTIVE_GAIN = 0.1

# Step 1: Encode Data into a 3D Lattice
def store_in_lattice(binary_data, harmonic_constant=HARMONIC_CONSTANT):
    normalized_data = binary_data / 255.0
    lattice_size = int(np.cbrt(len(normalized_data))) + 1
    lattice = np.zeros((lattice_size, lattice_size, lattice_size), dtype=np.float64)
    for idx, value in enumerate(normalized_data):
        x, y, z = idx % lattice_size, (idx // lattice_size) % lattice_size, idx // (lattice_size ** 2)
        lattice[x, y, z] += value * harmonic_constant
    return lattice

# Step 2: Retrieve Data from the 3D Lattice
def retrieve_from_lattice(lattice, harmonic_constant=HARMONIC_CONSTANT):
    flattened_data = lattice.flatten() / harmonic_constant
    return np.round(flattened_data * 255).astype(np.uint8)

# Step 3: Feedback Correction
def feedback_correction(lattice, binary_data, harmonic_constant=HARMONIC_CONSTANT):
    retrieved_data = retrieve_from_lattice(lattice, harmonic_constant)
    error = (binary_data - retrieved_data[:len(binary_data)]) / 255.0
    for idx, value in enumerate(error):
        x, y, z = idx % lattice.shape[0], (idx // lattice.shape[0]) % lattice.shape[1], idx // (lattice.shape[0] ** 2)
        lattice[x, y, z] += value * harmonic_constant
    return lattice

# Step 4: Apply Reflective Gain to Center Values
def apply_reflective_gain(lattice, gain_factor=REFLECTIVE_GAIN):
    center = np.array(lattice.shape) // 2
    for x in range(lattice.shape[0]):
        for y in range(lattice.shape[1]):
            for z in range(lattice.shape[2]):
                distance = np.linalg.norm(np.array([x, y, z]) - center)
                lattice[x, y, z] += gain_factor / (1 + distance)
    return lattice

# Visualization of the 3D Lattice
def visualize_lattice(lattice, iteration):
    x, y, z = np.nonzero(lattice)
    values = lattice[x, y, z]
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x, y, z, c=values, cmap='viridis', s=10)
    ax.set_title(f"3D Lattice Visualization of Harmonics - Iteration {iteration}", fontsize=16)
    ax.set_xlabel("X-axis", fontsize=12)
    ax.set_ylabel("Y-axis", fontsize=12)
    ax.set_zlabel("Z-axis", fontsize=12)
    plt.colorbar(scatter, ax=ax, label="Harmonic Values")
    plt.show()

# Main Execution
if __name__ == "__main__":
    # Load binary data from the file
    with open(r'd:\\colecovision.rom', 'rb') as file:
        binary_data = np.frombuffer(file.read(), dtype=np.uint8)
    
    # Store binary data in harmonic lattice
    lattice = store_in_lattice(binary_data)
    
    # Iterative feedback and visualization
    for iteration in range(MAX_ITERATIONS):
        print(f"Iteration {iteration + 1}")
        lattice = feedback_correction(lattice, binary_data)
        lattice = apply_reflective_gain(lattice)
        visualize_lattice(lattice, iteration + 1)
    
    # Retrieve data from the lattice
    retrieved_data = retrieve_from_lattice(lattice)
    
    # Outputs: Compare original and retrieved data
    print("Lattice Shape:", lattice.shape)
    print("Original Data (First 10 Bytes):", binary_data[:10])
    print("Retrieved Data (First 10 Bytes):", retrieved_data[:10])
    data_matches = np.array_equal(binary_data, retrieved_data[:len(binary_data)])
    print("Data matches:", data_matches)
    if not data_matches:
        print("Differences detected in data.")
