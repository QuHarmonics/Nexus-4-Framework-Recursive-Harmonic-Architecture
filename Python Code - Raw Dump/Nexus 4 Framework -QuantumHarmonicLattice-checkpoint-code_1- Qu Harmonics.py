import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
HARMONIC_CONSTANT = 0.35  # Scaling factor for harmonics

# Step 1: Encode Data into a 3D Lattice
def store_in_lattice(binary_data, harmonic_constant=HARMONIC_CONSTANT):
    """
    Store binary data into a 3D lattice using harmonic properties.
    """
    normalized_data = binary_data / 255.0  # Normalize binary data to [0, 1]
    
    # Create a 3D lattice
    lattice_size = int(np.cbrt(len(normalized_data))) + 1  # Ensure enough space
    lattice = np.zeros((lattice_size, lattice_size, lattice_size), dtype=np.float64)
    
    # Map data into lattice positions
    for idx, value in enumerate(normalized_data):
        x, y, z = idx % lattice_size, (idx // lattice_size) % lattice_size, idx // (lattice_size ** 2)
        lattice[x, y, z] += value * harmonic_constant  # Add harmonic scaling

    return lattice

# Step 2: Retrieve Data from the 3D Lattice
def retrieve_from_lattice(lattice, harmonic_constant=HARMONIC_CONSTANT):
    """
    Retrieve binary data from a 3D lattice.
    """
    flattened_data = lattice.flatten() / harmonic_constant  # Scale back by harmonic constant
    return np.round(flattened_data * 255).astype(np.uint8)  # Ensure integer byte output

# Visualization of the 3D Lattice
def visualize_lattice(lattice):
    """
    Visualize the 3D harmonic lattice.
    """
    x, y, z = np.nonzero(lattice)  # Get non-zero positions
    values = lattice[x, y, z]  # Corresponding harmonic values
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x, y, z, c=values, cmap='viridis', s=10)

    # Add labels and title
    ax.set_title("3D Lattice Visualization of Harmonics", fontsize=16)
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

    # Store binary data in harmonic lattice
    lattice = store_in_lattice(binary_data)

    # Retrieve binary data from the harmonic lattice
    retrieved_data = retrieve_from_lattice(lattice)

    # Visualize the harmonic lattice
    visualize_lattice(lattice)

    # Outputs: Compare Original and Retrieved Data
    print("Lattice Shape:", lattice.shape)
    print("Original Data (First 10 Bytes):", binary_data[:10])
    print("Retrieved Data (First 10 Bytes):", retrieved_data[:10])
    
    # Validate if original and retrieved data match
    data_matches = np.array_equal(binary_data, retrieved_data)
    print("Data matches:", data_matches)
    if not data_matches:
        print("Differences detected in data.")
