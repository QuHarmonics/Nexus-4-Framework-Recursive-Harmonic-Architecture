import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Input text data
text_data = """
01010100011010000110010100100000011100010111010101101001011000110110101100100000011000100111001001101111011101110110111000100000011001100110111101111000001000000110101001110101011011010111000001110011001000000110111101110110011001010111001000100000011101000110100001100101001000000110110001100001011110100111100100100000011001000110111101100111100000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000101011000
"""  # Replace this with your text input

# Convert text to numerical representation
binary_data = np.frombuffer(text_data.encode('utf-8'), dtype=np.uint8)  # Encode text as bytes

# Step 1: Encode the Text Data into H (Storage)
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
#z = np.sin(x / 10.0)  # Arbitrary wave to add depth

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
print("Original Text Data:", text_data[:10000])
print("Decoded Text Data:", retrieved_data.tobytes().decode('utf-8', errors='replace')[:10000])
print(harmonics)
original_binary = ''.join(text_data.split())  # Remove whitespace
decoded_binary = ''.join(map(str, retrieved_data))
is_equal = original_binary == decoded_binary
# Validate the process
is_equal = text_data.encode('utf-8') == retrieved_data.tobytes()
print("Data matches:", is_equal)
if not is_equal:
    print("Differences detected in data.")
