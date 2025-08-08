import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants
HARMONIC_CONSTANT = 0.35  # Scaling factor for harmonics
MAX_ITERATIONS = 10

# Step 1: Encode Data into a 3D Lattice
def store_in_lattice(binary_data, harmonic_constant=HARMONIC_CONSTANT):
    normalized_data = binary_data / 255.0  # Normalize binary data to [0, 1]
    data_length = len(normalized_data)
    lattice_size = int(np.ceil(np.cbrt(data_length)))  # Calculate lattice size
    padded_length = lattice_size ** 3

    # Pad data to fit the lattice
    padded_data = np.zeros(padded_length, dtype=np.float64)
    padded_data[:data_length] = normalized_data

    # Create 3D lattice
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
    center = np.array(lattice.shape) // 2  # Find the center of the lattice

    for idx, value in enumerate(error):
        x, y, z = idx % lattice.shape[0], (idx // lattice.shape[0]) % lattice.shape[1], idx // (lattice.shape[0] ** 2)
        distance_from_center = np.linalg.norm(np.array([x, y, z]) - center)
        scaling_factor = np.exp(-distance_from_center / lattice.shape[0])  # Decay with distance
        lattice[x, y, z] += value * harmonic_constant * scaling_factor  # Adjust with scaling factor

    return lattice

# Step 4: Normalize the Lattice
def normalize_lattice(lattice):
    lattice -= np.mean(lattice)
    lattice /= np.max(np.abs(lattice))
    lattice *= HARMONIC_CONSTANT
    return lattice

# Visualization of the 3D Lattice
def visualize_lattice(lattice, iteration):
    x, y, z = np.nonzero(lattice)  # Get non-zero positions
    values = lattice[x, y, z]  # Corresponding harmonic values

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    scatter = ax.scatter(x, y, z, c=values, cmap='viridis', s=10)

    ax.set_title(f"3D Lattice Visualization of Harmonics - Iteration {iteration}", fontsize=16)
    ax.set_xlabel("X-axis", fontsize=12)
    ax.set_ylabel("Y-axis", fontsize=12)
    ax.set_zlabel("Z-axis", fontsize=12)
    plt.colorbar(scatter, ax=ax, label="Harmonic Values")
    plt.show()

# Visualize Byte-wise Differences
def visualize_byte_differences(original_data, retrieved_data):
    differences = np.abs(original_data - retrieved_data)
    plt.figure(figsize=(10, 6))
    plt.bar(range(len(differences)), differences, color='orange')
    plt.title("Byte-wise Differences Between Original and Retrieved Data")
    plt.xlabel("Byte Index")
    plt.ylabel("Difference")
    plt.show()

# Main Execution
if __name__ == "__main__":
    # Load binary bits from your file
    with open(r'd:\\colecovision.rom', 'rb') as file:
        binary_data = np.frombuffer(file.read(), dtype=np.uint8)  # Read binary as bytes

    # Step 1: Store binary data in harmonic lattice
    lattice, data_length = store_in_lattice(binary_data)

    for iteration in range(1, MAX_ITERATIONS + 1):
        print(f"Iteration {iteration}")
        
        # Step 3: Apply feedback correction
        lattice = feedback_correction(lattice, binary_data)
        
        # Step 4: Normalize lattice values
        lattice = normalize_lattice(lattice)
        
        # Step 5: Retrieve data from lattice
        retrieved_data = retrieve_from_lattice(lattice, data_length=data_length)

        # Step 6: Visualize lattice and differences
        visualize_lattice(lattice, iteration)
        visualize_byte_differences(binary_data[:100], retrieved_data[:100])  # Compare first 100 bytes

        # Validate if original and retrieved data match
        data_matches = np.array_equal(binary_data, retrieved_data)
        print("Lattice Shape:", lattice.shape)
        print("Original Data (First 10 Bytes):", binary_data[:10])
        print("Retrieved Data (First 10 Bytes):", retrieved_data[:10])
        print("Data matches:", data_matches)
        
        if data_matches:
            print("Data successfully recovered!")
            break
        else:
            print("Differences detected in data.")
