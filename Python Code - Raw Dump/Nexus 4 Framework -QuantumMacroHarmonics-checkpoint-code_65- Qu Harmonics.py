import numpy as np
import matplotlib.pyplot as plt

# Generate nested rings based on hash compression
def generate_rings(hash_hex, max_layers=10, resolution=1000):
    binary_data = np.array([int(bit) for bit in ''.join(format(int(c, 16), '04b') for c in hash_hex)])
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, resolution)

    rings = []
    for layer in range(max_layers):
        radius = 1 - (layer / max_layers)  # Compress inward
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        # Match data length to resolution by repeating or truncating binary_data
        data = np.tile(binary_data, resolution // n + 1)[:resolution]
        rings.append({"x": x, "y": y, "data": data})
    return rings

# Visualize nested rings
def plot_rings(rings, title):
    plt.figure(figsize=(8, 8))
    for ring in rings:
        plt.scatter(ring["x"], ring["y"], c=ring["data"], cmap='coolwarm', s=10, alpha=0.7)
    plt.title(title)
    plt.colorbar(label="Binary Value")
    plt.show()

# Input hash (SHA-256 for 'abc')
hash_hex = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Generate nested rings
rings = generate_rings(hash_hex)

# Visualize nested rings
plot_rings(rings, title="Nested Rings Representation of SHA-256 ('abc')")
