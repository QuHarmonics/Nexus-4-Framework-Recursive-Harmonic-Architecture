import numpy as np
import matplotlib.pyplot as plt

# Parameters
spiral_size = 1024  # Initial size of the spiral
hash_size = 256  # Diameter of the hash (side view)
threshold = spiral_size  # Threshold for exponential compression
rings = []  # Storage for the rings of the spiral
compressed_hash = np.zeros(hash_size)  # Initial hash representation

# Generate the Spiral (Pre-Threshold)
def generate_spiral(data_length):
    rings = []
    for i in range(data_length):
        angle = 2 * np.pi * i / data_length
        x = np.cos(angle) * i / data_length
        y = np.sin(angle) * i / data_length
        rings.append({"x": x, "y": y, "data": np.random.randint(0, 2)})  # Random binary data
    return rings

# Lens Projection: Map the Spiral to the Hash
def project_to_hash(rings, hash_size):
    hash_projection = np.zeros(hash_size)
    for i, ring in enumerate(rings):
        position = int((i / len(rings)) * hash_size)
        hash_projection[position] += ring["data"]  # Summing perpendicular contributions
    return hash_projection

# Compression Phase: Exponential Growth Beyond Threshold
def compress_spiral(rings, hash_size, threshold):
    compressed_hash = np.zeros(hash_size)
    for i, ring in enumerate(rings):
        if i < threshold:
            position = int((i / threshold) * hash_size)
            compressed_hash[position] += ring["data"]  # Pre-threshold: Linear mapping
        else:
            # Post-threshold: Exponential compression
            compression_factor = (i - threshold + 1)
            position = int(hash_size / compression_factor)  # Exponential stacking
            compressed_hash[position] += ring["data"]
    return compressed_hash

# Visualization Function
def visualize_simulation(rings, hash_projection, compressed_hash, title):
    plt.figure(figsize=(12, 12))

    # Top View: Full Spiral
    plt.subplot(3, 1, 1)
    x_coords = [ring["x"] for ring in rings]
    y_coords = [ring["y"] for ring in rings]
    data_vals = [ring["data"] for ring in rings]
    plt.scatter(x_coords, y_coords, c=data_vals, cmap='coolwarm', s=5, alpha=0.7)
    plt.title(f"{title}: Top View (Spiral)")

    # Side View: Hash Projection
    plt.subplot(3, 1, 2)
    plt.bar(range(len(hash_projection)), hash_projection, color='orange')
    plt.title(f"{title}: Side View (Hash Projection)")

    # Side View: Compressed Hash
    plt.subplot(3, 1, 3)
    plt.bar(range(len(compressed_hash)), compressed_hash, color='green')
    plt.title(f"{title}: Side View (Compressed Hash)")

    plt.tight_layout()
    plt.show()

# Run the Simulation
rings = generate_spiral(spiral_size)
hash_projection = project_to_hash(rings, hash_size)
compressed_hash = compress_spiral(rings, hash_size, threshold)

# Visualize Results
visualize_simulation(rings, hash_projection, compressed_hash, "SHA-256 Spiral Compression")
