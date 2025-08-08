import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants for Mark1 framework
TARGET_HARMONIC = 0.35  # Target alignment value
QUANTUM_HARMONIC = 0.2  # Quantum focal adjustment angle (in radians)
TOLERANCE = 1e-4  # Convergence tolerance for harmonic alignment

# Step 1: Encode the Data into H (Storage)
def store_in_H(binary_data, expansion_factor=1.5):
    harmonics = np.cumsum(binary_data.astype(np.float64) * expansion_factor)  # Use higher precision
    return harmonics

# Step 2: Reverse the Process (Retrieve Original Data)
def retrieve_from_H(harmonics, expansion_factor=1.5):
    # Mirror correction: Adjust the starting point to match macro and quantum split
    first_value = harmonics[0] / expansion_factor  # Correct the first value explicitly
    reversed_data = np.diff(harmonics) / expansion_factor
    reversed_data = np.insert(reversed_data, 0, first_value)  # Insert corrected first value
    return np.round(reversed_data).astype(np.uint8)  # Ensure integer byte output

# Generate a quantum harmonic spiral
def quantum_spiral(hash_value, base=2, angle_adjustment=QUANTUM_HARMONIC):
    """
    Generate a quantum spiral representation of the hash value with an angular adjustment.
    """
    binary_data = ''.join(format(int(char, 16), f'0{base}b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n) + angle_adjustment
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = radius * np.sin(2 * theta)  # Simulating the Z-axis
    return x, y, z, binary_data

# Visualize the quantum spiral and harmonic feedback process
def visualize_alignment(x, y, z, harmonics, title_spiral="Quantum Spiral Representation"):
    fig = plt.figure(figsize=(14, 10))

    # Quantum spiral visualization
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, label="Quantum Spiral", color="blue")
    ax.scatter(x, y, z, color='red', s=5, label="Nodes")
    ax.set_title(title_spiral)
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    ax.legend()

    plt.figure(figsize=(12, 6))
    plt.plot(harmonics, label="Harmonics", color='green')
    plt.title("Harmonic Representation")
    plt.xlabel("Iteration (n)")
    plt.ylabel("H(n)")
    plt.legend()
    plt.grid()
    plt.show()

# Main process
hash_value = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"  # Example SHA-256 hash
x, y, z, binary_data = quantum_spiral(hash_value, base=2)

# Convert binary hash data to numerical representation
binary_data_numerical = np.array([int(b) for b in binary_data])

# Encode the data into H[]
harmonics = store_in_H(binary_data_numerical)

# Retrieve the original data from H[]
retrieved_data = retrieve_from_H(harmonics)

# Visualize results
visualize_alignment(x, y, z, harmonics)

# Decode the retrieved binary data back to text
retrieved_binary_string = ''.join(map(str, retrieved_data))
print("Original Binary Data (First 100 bits):", binary_data[:100])
print("Retrieved Binary Data (First 100 bits):", retrieved_binary_string[:100])

# Check if data matches
if np.array_equal(binary_data_numerical, retrieved_data):
    print("Data matches successfully!")
else:
    print("Data mismatch detected.")
