import numpy as np
import matplotlib.pyplot as plt

# Convert hash to binary and reshape into a grid
def hash_to_matrix(hash_str):
    binary_representation = bin(int(hash_str, 16))[2:].zfill(256)  # Convert to binary (256 bits)
    matrix = np.array([int(b) for b in binary_representation]).reshape(16, 16)  # Reshape to 16x16 grid
    return matrix

# Create a lattice for 3D visualization
def create_lattice(matrix):
    x, y = np.meshgrid(range(matrix.shape[0]), range(matrix.shape[1]))
    z = matrix
    return x, y, z

# Plot the 3D lattice
def plot_3d_lattice(x, y, z, title):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, z, cmap='viridis', edgecolor='none')
    ax.set_title(title)
    plt.show()

# Overlay two 3D hashes
def overlay_hashes(x1, y1, z1, x2, y2, z2):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x1, y1, z1, cmap='viridis', alpha=0.7, edgecolor='none')
    ax.plot_surface(x2, y2, z2, cmap='plasma', alpha=0.5, edgecolor='none')
    ax.set_title("Overlay of Hash 1 and Hash 2")
    plt.show()

# Hashes
hash1 = "dffd6021bb2bd5b0af676290809ec3a53191dd81c7f70a4b28688a362182986f"  # Hello, World!
hash2 = "fab895797a2e77d11bd7fc6c54b73ed9d8b0bd59540f9bb1c72ca3eed0ff95bb"  # Hello, Quantum!

# Generate matrices and lattices
matrix1 = hash_to_matrix(hash1)
matrix2 = hash_to_matrix(hash2)

x1, y1, z1 = create_lattice(matrix1)
x2, y2, z2 = create_lattice(matrix2)

# Plot individual 3D lattices
plot_3d_lattice(x1, y1, z1, "3D View of Hash 1 (Hello, World!)")
plot_3d_lattice(x2, y2, z2, "3D View of Hash 2 (Hello, Quantum!)")

# Overlay and compare
overlay_hashes(x1, y1, z1, x2, y2, z2)
