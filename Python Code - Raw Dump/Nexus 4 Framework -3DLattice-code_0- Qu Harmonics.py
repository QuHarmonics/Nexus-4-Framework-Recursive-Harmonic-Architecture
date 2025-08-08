import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
HARMONIC_CONSTANT = 0.35
TOLERANCE = 0.01
EXPANSION_FACTOR = 1.5
MAX_ITERATIONS = 100

# Step 1: Encode Data into H using a Lattice Structure
def store_in_lattice(binary_data, harmonic_constant=HARMONIC_CONSTANT):
    """
    Store binary data into a 3D lattice using harmonic properties.
    """
    # Normalize binary data to [0, 1]
    normalized_data = binary_data / 255.0
    
    # Create a 3D lattice
    lattice_size = int(np.cbrt(len(normalized_data))) + 1  # Ensure enough space
    lattice = np.zeros((lattice_size, lattice_size, lattice_size), dtype=np.float64)
    
    # Map data into lattice positions
    for idx, value in enumerate(normalized_data):
        x, y, z = idx % lattice_size, (idx // lattice_size) % lattice_size, idx // (lattice_size ** 2)
        lattice[x, y, z] += value * harmonic_constant  # Add harmonic scaling

    return lattice

# Step 2: Retrieve Data from Lattice
def retrieve_from_lattice(lattice, harmonic_constant=HARMONIC_CONSTANT):
    """
    Retrieve binary data from a 3D lattice.
    """
    flattened_data = lattice.flatten() / harmonic_constant
    return np.round(flattened_data * 255).astype(np.uint8)  # Scale back to original range

# Step 3: Visualize the Lattice
def visualize_lattice(lattice):
    """
    Create a 3D scatter plot of the lattice points.
    """
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    x, y, z = np.nonzero(lattice)
    values = lattice[x, y, z]
    
    # Plot the lattice points
    ax.scatter(x, y, z, c=values, cmap='viridis', s=20)
    ax.set_title("3D Lattice Visualization of Harmonics", fontsize=16)
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    plt.show()

# Test the Lattice-Based Compression
if __name__ == "__main__":
    # Load binary data (example: BIOS file)
    with open(r'd:\colecovision.rom', 'rb') as file:
        binary_data = np.frombuffer(file.read(), dtype=np.uint8)  # Read binary as bytes

    # Step 1: Store Data in Lattice
    harmonic_lattice = store_in_lattice(binary_data, HARMONIC_CONSTANT)
    print("Lattice Shape:", harmonic_lattice.shape)

    # Step 2: Retrieve Data from Lattice
    retrieved_data = retrieve_from_lattice(harmonic_lattice, HARMONIC_CONSTANT)
    
    # Validate the process
    is_equal = np.array_equal(binary_data, retrieved_data)
    print("Data matches:", is_equal)
    if not is_equal:
        print("Differences detected in data.")

    # Step 3: Visualize the Lattice
    visualize_lattice(harmonic_lattice)

    # Print a sample of the data
    print("Original Data (First 10 Bytes):", binary_data[:10])
    print("Retrieved Data (First 10 Bytes):", retrieved_data[:10])
