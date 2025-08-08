import numpy as np
import matplotlib.pyplot as plt

# Generate concentric rings with density contributions
def generate_density_rings(hash_binary, num_rings=16):
    angles = np.linspace(0, 2 * np.pi, len(hash_binary), endpoint=False)
    rings = []
    for r in range(1, num_rings + 1):
        # Each ring is proportional to its distance from the center
        contributions = np.array([int(bit) * (1 / r) for bit in hash_binary])
        x = r * np.cos(angles)
        y = r * np.sin(angles)
        rings.append({"x": x, "y": y, "contributions": contributions})
    return rings

# Calculate perpendicular summation at each node
def calculate_perpendicular_sums(rings, num_nodes):
    node_sums = np.zeros(num_nodes)
    for ring in rings:
        node_sums += ring["contributions"]
    return node_sums

# Visualize the density map
def plot_density_map(rings, node_sums, title):
    plt.figure(figsize=(8, 8))
    
    # Plot rings
    for ring in rings:
        plt.scatter(ring["x"], ring["y"], c=ring["contributions"], cmap='coolwarm', s=10, alpha=0.7)
    
    # Overlay perpendicular summations on circumference
    angles = np.linspace(0, 2 * np.pi, len(node_sums), endpoint=False)
    x = np.cos(angles)
    y = np.sin(angles)
    plt.scatter(x, y, c=node_sums, cmap='viridis', s=50, edgecolor='black', label="Perpendicular Sums")
    
    plt.title(title)
    plt.colorbar(label="Density Contribution")
    plt.legend()
    plt.show()

# Example SHA-256 binary hash
hash_hex = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
hash_binary = [int(bit) for bit in ''.join(format(int(c, 16), '04b') for c in hash_hex)]

# Generate density rings
rings = generate_density_rings(hash_binary, num_rings=32)

# Calculate perpendicular summations
node_sums = calculate_perpendicular_sums(rings, len(hash_binary))

# Visualize the result
plot_density_map(rings, node_sums, title="Perpendicular Summation Density Map of SHA-256")
