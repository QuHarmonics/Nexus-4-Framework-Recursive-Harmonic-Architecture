import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load binary data (example)
binary_data = np.random.randint(0, 256, size=1024, dtype=np.uint8)  # Example random binary data

# Step 1: Compress Data into H
def compress_to_H(binary_data, compression_factor=2):
    """
    Compress binary data into harmonics using a specified compression factor.

    Args:
        binary_data (np.array): Input binary data as bytes.
        compression_factor (int): Factor to reduce the data size.

    Returns:
        np.array: Compressed harmonics.
    """
    harmonics = np.cumsum(binary_data.astype(np.float64))  # Cumulative sum for harmonic encoding
    return harmonics[::compression_factor]  # Downsample based on the compression factor

# Step 2: Decompress Data from H
def decompress_from_H(compressed_harmonics, original_size, compression_factor=2):
    """
    Retrieve the original binary data from compressed harmonics.

    Args:
        compressed_harmonics (np.array): Compressed harmonic data.
        original_size (int): Original size of the binary data.
        compression_factor (int): Factor used during compression.

    Returns:
        np.array: Reconstructed binary data as bytes.
    """
    # Reconstruct cumulative harmonics
    reconstructed_harmonics = np.interp(
        np.linspace(0, len(compressed_harmonics) - 1, original_size),
        np.arange(len(compressed_harmonics)),
        compressed_harmonics
    )
    # Reverse cumulative sum to retrieve original data
    reversed_data = np.diff(np.insert(reconstructed_harmonics, 0, 0))  # Reverse cumulative sum
    return np.round(reversed_data).astype(np.uint8)  # Ensure integer byte output

# Compress the binary data
compression_factor = 2
compressed_harmonics = compress_to_H(binary_data, compression_factor=compression_factor)

# Decompress to retrieve the original binary data
retrieved_data = decompress_from_H(compressed_harmonics, len(binary_data), compression_factor=compression_factor)

# Visualization of the Compressed Harmonics
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Generate 3D coordinates for visualization
x = np.arange(len(compressed_harmonics))
y = compressed_harmonics
z = np.sin(x / 10.0)  # Arbitrary sine wave for depth

ax.plot(x, y, z, label="Compressed H(n)", color='g', lw=2)
ax.scatter(x, y, z, color='r', s=5, label="Compressed Nodes")

# Add labels and legend
ax.set_title("3D Visualization of Compressed H(n)", fontsize=16)
ax.set_xlabel("Iteration (n)", fontsize=12)
ax.set_ylabel("H(n)", fontsize=12)
ax.set_zlabel("Z-axis Wave", fontsize=12)
ax.legend()
plt.show()

# Outputs: Full Original and Decoded Data Comparison
print("=== Full Original Binary Data ===")
print(binary_data.tolist())  # Output full data as a list for readability

print("\n=== Full Compressed Harmonics ===")
print(compressed_harmonics.tolist())

print("\n=== Full Decoded Binary Data ===")
print(retrieved_data.tolist())

# Validate the process
is_equal = np.array_equal(binary_data, retrieved_data)
print("\nData matches:", is_equal)

if not is_equal:
    print("\nDifferences Detected:")
    mismatched_indices = np.where(binary_data != retrieved_data)[0]
    print("Indices where data mismatch occurs:", mismatched_indices.tolist())
    print("Original values at mismatched indices:", binary_data[mismatched_indices].tolist())
    print("Decoded values at mismatched indices:", retrieved_data[mismatched_indices].tolist())
