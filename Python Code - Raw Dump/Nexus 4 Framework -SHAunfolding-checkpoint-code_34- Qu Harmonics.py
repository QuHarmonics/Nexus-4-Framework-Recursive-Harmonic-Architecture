import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load binary data (example placeholder for file loading)
# Replace 'binary_data' with actual binary data for your use case
binary_data = np.random.randint(0, 256, size=1024, dtype=np.uint8)  # Example random binary data

# Step 1: Compress Data into H
def compress_to_H(binary_data, compression_factor=0.5):
    """
    Compress binary data into harmonics using a specified compression factor.

    Args:
        binary_data (np.array): Input binary data as bytes.
        compression_factor (float): Factor to reduce the data size.

    Returns:
        np.array: Compressed harmonics.
    """
    harmonics = np.cumsum(binary_data.astype(np.float64) * compression_factor)
    return harmonics[::int(1 / compression_factor)]  # Sample based on compression factor

# Step 2: Decompress Data from H
def decompress_from_H(compressed_harmonics, compression_factor=0.5):
    """
    Retrieve the original binary data from compressed harmonics.

    Args:
        compressed_harmonics (np.array): Compressed harmonic data.
        compression_factor (float): Factor used during compression.

    Returns:
        np.array: Reconstructed binary data as bytes.
    """
    reversed_data = np.diff(compressed_harmonics) / compression_factor
    first_value = compressed_harmonics[0] / compression_factor
    reversed_data = np.insert(reversed_data, 0, first_value)
    return np.round(reversed_data).astype(np.uint8)  # Ensure integer byte output

# Compress the binary data
compression_factor = 0.5
compressed_harmonics = compress_to_H(binary_data, compression_factor=compression_factor)

# Decompress the binary data
retrieved_data = decompress_from_H(compressed_harmonics, compression_factor=compression_factor)

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

# Outputs: Original and Decoded Data Comparison
print("Original Binary Data (First 10 bytes):", binary_data)
print("Compressed Harmonics (First 10):", compressed_harmonics)
print("Decoded Binary Data (First 10 bytes):", retrieved_data)

# Validate the process
is_equal = np.array_equal(binary_data, retrieved_data)
print("Data matches:", is_equal)
if not is_equal:
    print("Differences detected in data.")
