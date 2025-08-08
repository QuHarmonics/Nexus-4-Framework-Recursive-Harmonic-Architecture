import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load binary bits from your BIOS file
with open(r'd:\test.wav', 'rb') as file:
    binary_data = np.frombuffer(file.read(), dtype=np.uint8)  # Read binary as bytes

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

# Run the process
harmonics = store_in_H(binary_data)
retrieved_data = retrieve_from_H(harmonics)

# Visualize H(n)
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Generate 3D coordinates
x = np.arange(len(harmonics))
y = harmonics
z = np.sin(x / 10.0)  # Arbitrary wave to add depth

# Plot the 3D visualization
ax.plot(x, y, z, label="H(n) in 3D", color='b', lw=2)
ax.scatter(x, y, z, color='r', s=5, label="Nodes")

# Labels and title
ax.set_title("3D Visualization of H(n)", fontsize=16)
ax.set_xlabel("Iteration (n)", fontsize=12)
ax.set_ylabel("H(n)", fontsize=12)
ax.set_zlabel("Z-axis Wave", fontsize=12)
ax.legend()
plt.show()

# Outputs: Original and Decoded Data Comparison
print("Original Binary Data (First 10 bytes):", binary_data[:1000])
print("Decoded Binary Data (First 10 bytes):", retrieved_data[:1000])

# Validate the process
is_equal = np.array_equal(binary_data, retrieved_data)
print("Data matches:", is_equal)
if not is_equal:
    print("Differences detected in data.")