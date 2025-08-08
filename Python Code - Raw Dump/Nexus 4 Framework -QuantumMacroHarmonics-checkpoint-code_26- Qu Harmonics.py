import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Input text data
text_data = """
01101000011001010110110001101100011011111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000101000
"""  # Replace this with your text input

# Convert text to numerical representation
binary_data = np.frombuffer(text_data.encode('utf-8'), dtype=np.uint8)  # Encode text as bytes

# Step 1: Encode the Text Data into H (Storage)
def store_in_H(binary_data):
    # cumsum with int32
    cumsum_arr = np.cumsum(binary_data.astype(np.int32))
    return cumsum_arr

def retrieve_from_H(cumsum_arr):
    # difference in int32
    first_val = cumsum_arr[0]
    diffs = np.diff(cumsum_arr)
    diffs = np.insert(diffs, 0, first_val)
    return diffs.astype(np.uint8)

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
print("Original Text Data (First 100 characters):", text_data[:100])
print("Decoded Text Data (First 100 characters):", retrieved_data.tobytes().decode('utf-8', errors='replace')[:100])

# Validate the process
is_equal = text_data.encode('utf-8') == retrieved_data.tobytes()
print("Data matches:", is_equal)
if not is_equal:
    print("Differences detected in data.")
