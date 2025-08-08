import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Hash input
input_hash_hex = "0753e9bcdebcbe1bbbf19b6bfe6c06bd9606ce89f3c6d57d6d03a6c257f44b3f"

# Convert hash to binary representation
binary_data = ''.join(format(int(input_hash_hex[i:i + 2], 16), '08b') for i in range(0, len(input_hash_hex), 2))

# Constants
growth_factor = 1.5
iterations = 10  # Number of growth iterations
zeta_zero = (0, 0, 0)  # Starting point

# Initialize nodes
nodes = [zeta_zero]

# Expand the binary data symmetrically
def expand_binary_nodes(binary_data, center, growth_factor, iterations):
    nodes = [center]
    for i in range(iterations):
        new_nodes = []
        for node in nodes:
            x, y, z = node
            for j, bit in enumerate(binary_data):
                offset = (j + 1) * growth_factor
                if bit == '1':
                    # Positive ripple
                    new_nodes.append((x + offset, y, z + offset))
                else:
                    # Negative ripple
                    new_nodes.append((x - offset, y, z - offset))
        nodes.extend(new_nodes)
        binary_data = binary_data[::-1]  # Mirror binary data for symmetry
    return nodes

# Grow nodes
fractal_nodes = expand_binary_nodes(binary_data, zeta_zero, growth_factor, iterations)

# Plot the nodes in 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

x_vals = [node[0] for node in fractal_nodes]
y_vals = [node[1] for node in fractal_nodes]
z_vals = [node[2] for node in fractal_nodes]

ax.scatter(x_vals, y_vals, z_vals, c=np.abs(np.array(z_vals)), cmap='viridis', s=10)

# Labels and title
ax.set_title("Fractal Growth from Binary Hash")
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")
plt.show()
