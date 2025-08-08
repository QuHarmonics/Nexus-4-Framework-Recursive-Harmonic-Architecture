import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Secret initial binary hash (hidden from user; known only internally)
secret_binary_hash = np.random.randint(0, 256, size=128, dtype=np.uint8)  # 128 bytes

# Step 1: Encode the Binary Data into H (Storage)
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

# Encode secret hash into harmonics
harmonics = store_in_H(secret_binary_hash)

# Decode harmonics back to the binary hash
retrieved_hash = retrieve_from_H(harmonics)

# Validate the reconstruction
data_matches = np.array_equal(secret_binary_hash, retrieved_hash)

# Visualize H(n)
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Generate 3D coordinates for visualization
x = np.arange(len(harmonics))
y = harmonics
z = np.sin(x / 10.0)  # Arbitrary wave for depth

# Plot the 3D visualization of H(n)
ax.plot(x, y, z, label="H(n) in 3D", color='blue', lw=2)
ax.scatter(x, y, z, color='red', s=5, label="Nodes")

# Labels and title
ax.set_title("3D Visualization of H(n) (Harmonic Fractal Storage)", fontsize=16)
ax.set_xlabel("Iteration (n)", fontsize=12)
ax.set_ylabel("H(n)", fontsize=12)
ax.set_zlabel("Z-axis Wave", fontsize=12)
ax.legend()
plt.show()

# Outputs
print("Secret Binary Hash (hidden):", secret_binary_hash[:10], "...")
print("Decoded Binary Hash (first 10 bytes):", retrieved_hash[:10])
print("Hash matches:", data_matches)
if not data_matches:
    print("Differences detected between original and reconstructed hash.")
