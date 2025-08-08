import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
HARMONIC_CONSTANT = 0.35  # Scaling factor for harmonics
MAX_ITERATIONS = 5

# Step 1: Encode Data into a 3D Lattice
def store_in_lattice(binary_data, harmonic_constant=HARMONIC_CONSTANT):
    normalized_data = binary_data / 255.0  # Normalize binary data to [0, 1]
    
    # Adjust lattice size to fit all data
    data_length = len(normalized_data)
    lattice_size = int(np.ceil(np.cbrt(data_length)))  # Exact size needed
    padded_length = lattice_size ** 3
    
    # Pad the normalized data with zeros if needed
    padded_data = np.zeros(padded_length, dtype=np.float64)
    padded_data[:data_length] = normalized_data

    # Create the 3D lattice
    lattice = padded_data.reshape((lattice_size, lattice_size, lattice_size)) * harmonic_constant
    return lattice, data_length

# Step 2: Retrieve Data from the 3D Lattice
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

# Visualization of the 3D Lattice
def visualize_lattice(lattice):
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

    # Step 1: Store binary data in harmonic lattice
    lattice, data_length = store_in_lattice(binary_data)

    for iteration in range(MAX_ITERATIONS):
        print(f"Iteration {iteration + 1}")
        lattice = feedback_correction(lattice, binary_data)
        retrieved_data = retrieve_from_lattice(lattice, data_length=data_length)

        # Step 2: Visualize the harmonic lattice
        visualize_lattice(lattice)

        # Step 3: Outputs: Compare Original and Retrieved Data
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
