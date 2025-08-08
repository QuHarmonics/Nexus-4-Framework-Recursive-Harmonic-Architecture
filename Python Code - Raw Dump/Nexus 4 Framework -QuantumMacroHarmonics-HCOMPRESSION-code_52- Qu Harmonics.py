import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create a nested ring structure from the binary representation of the hash
def generate_rings_tilted(hash_hex, num_rings=32, tilt_angle=0):
    binary_data = np.array([int(bit) for bit in ''.join(format(int(c, 16), '04b') for c in hash_hex)])
    num_points_per_ring = len(binary_data) // num_rings

    rings = []
    for i in range(num_rings):
        ring_data = binary_data[i * num_points_per_ring: (i + 1) * num_points_per_ring]
        angle = np.linspace(0, 2 * np.pi, num_points_per_ring)
        radius = i / num_rings  # Spread rings evenly from center to edge
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        z = np.tan(np.radians(tilt_angle)) * radius  # Tilt based on angle
        rings.append({"x": x, "y": y, "z": z, "data": ring_data})
    return rings

# Plot the tilted rings
def plot_tilted_rings(rings, title="Tilted View of Nested Rings"):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    for ring in rings:
        ax.scatter(ring["x"], ring["y"], ring["z"], c=ring["data"], cmap="coolwarm", s=10, alpha=0.8)
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis (Tilted View)")
    plt.show()

# SHA-256 hash for 'abc'
hash_hex = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Generate tilted rings
tilt_angle = 90  # Adjust the angle to "tilt" the view into the plane of the hash
tilted_rings = generate_rings_tilted(hash_hex, tilt_angle=tilt_angle)

# Visualize the tilted rings
plot_tilted_rings(tilted_rings, title=f"Tilted View of Nested Rings at {tilt_angle}°")
