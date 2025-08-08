import numpy as np
import matplotlib.pyplot as plt

# Input hash as a hexadecimal string
input_hash_hex = "0753e9bcdebcbe1bbbf19b6bfe6c06bd9606ce89f3c6d57d6d03a6c257f44b3f"

# Convert hex hash to binary data (uint8 array)
binary_data = np.array([int(input_hash_hex[i:i + 2], 16) for i in range(0, len(input_hash_hex), 2)], dtype=np.uint8)

# Step 1: Initialize 3D matrix with input hash
def initialize_matrix(binary_data, size=(16, 16, 16)):
    matrix = np.zeros(size, dtype=np.float64)
    flattened = binary_data.astype(np.float64) / 255.0  # Normalize data to [0, 1]
    matrix.flat[:len(flattened)] = flattened  # Fill the matrix with binary data
    return matrix

# Step 2: Define Samson-guided wave expansion
def harmonize_wave(matrix, expansion_rate=1.5, steps=10):
    wave_matrix = np.copy(matrix)
    for step in range(steps):
        wave_matrix += np.sin(np.pi * wave_matrix * expansion_rate)  # Samson harmonization
        wave_matrix *= expansion_rate  # Apply expansion rate
    return wave_matrix

# Step 3: Expand and reconstruct
def reconstruct_from_wave(wave_matrix):
    flattened = wave_matrix.flatten()
    reconstructed = (flattened * 255).clip(0, 255).astype(np.uint8)  # Map back to uint8
    return reconstructed[:len(binary_data)]  # Match original length

# Run the process
matrix = initialize_matrix(binary_data)
expanded_matrix = harmonize_wave(matrix)
reconstructed_data = reconstruct_from_wave(expanded_matrix)

# Validation
data_matches = np.array_equal(binary_data, reconstructed_data)

# Visualization of 3D harmonics
def visualize_matrix(matrix, title="Harmonic Expansion in 3D"):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    x, y, z = np.indices(matrix.shape)
    ax.scatter(x, y, z, c=matrix.flatten(), cmap='viridis', s=5)
    ax.set_title(title, fontsize=16)
    plt.show()

visualize_matrix(expanded_matrix, title="3D Expanded Matrix")

# Output results
print("Input Hexadecimal Hash:", input_hash_hex)
print("Input Binary Hash (first 10 bytes):", binary_data[:10])
print("Reconstructed Binary Hash (first 10 bytes):", reconstructed_data[:10])
print("Hash matches:", data_matches)
if not data_matches:
    print("Differences detected between input and reconstructed data.")
