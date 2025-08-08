import numpy as np
import matplotlib.pyplot as plt

# Generate a spiral based on hash diameter
def generate_spiral(hash_hex, diameter=256, resolution=1000):
    binary_data = np.array([int(bit) for bit in ''.join(format(int(c, 16), '04b') for c in hash_hex)])
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, resolution)
    r = np.linspace(diameter / 2, 0, resolution)  # Spiral inward
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y, binary_data

# Visualize the spiral
def plot_spiral(x, y, binary_data, title):
    plt.figure(figsize=(10, 8))
    plt.scatter(x, y, c=binary_data, cmap='coolwarm', s=10)
    plt.title(title)
    plt.colorbar(label="Binary Value")
    plt.show()

# Unfold the spiral
def unfold_spiral(binary_data, resolution=1000):
    unfolded = np.zeros_like(binary_data)
    for i in range(len(binary_data)):
        # Simulate unfolding by reversing the spiral's inward flow
        unfolded[i] = binary_data[-(i + 1)]
    return unfolded

# Hash for 'abc'
hash_hex = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Generate the spiral
x, y, binary_data = generate_spiral(hash_hex)

# Visualize the spiral
plot_spiral(x, y, binary_data, title="Spiral Representation of SHA-256 ('abc')")

# Unfold the spiral
unfolded_data = unfold_spiral(binary_data)

# Visualize the unfolded data
plt.figure(figsize=(12, 6))
plt.plot(unfolded_data, label="Unfolded Spiral Data", color='green')
plt.title("Unfolded Spiral Data")
plt.legend()
plt.grid()
plt.show()
