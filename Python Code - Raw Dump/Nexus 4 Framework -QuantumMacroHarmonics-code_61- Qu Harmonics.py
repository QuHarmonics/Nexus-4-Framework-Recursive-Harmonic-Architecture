import numpy as np
import matplotlib.pyplot as plt

# Generate spiral for a hash
def hash_spiral(hash_hex, layers=3, resolution=1000):
    binary_data = np.array([int(bit) for bit in ''.join(format(int(c, 16), '04b') for c in hash_hex)])
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi * layers, n)  # Spiral through layers
    r = np.linspace(1, 0, n)  # Radius decreases toward center
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, binary_data

# Visualize spiral
def plot_spiral(x, y, binary_data, title):
    plt.figure(figsize=(10, 8))
    plt.scatter(x, y, c=binary_data, cmap='coolwarm', s=10)
    plt.title(title)
    plt.colorbar(label="Binary Value")
    plt.show()

# Hash for 'abc'
hash_hex = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Generate spiral
x, y, binary_data = hash_spiral(hash_hex)

# Visualize the spiral
plot_spiral(x, y, binary_data, title="Hash Spiral for SHA-256 ('abc')")
