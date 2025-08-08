import numpy as np
import matplotlib.pyplot as plt

# Simulate circular weaving of hash input
def circular_weaving(hash_hex, iterations=1000):
    binary_data = np.array([int(bit) for bit in ''.join(format(int(c, 16), '04b') for c in hash_hex)])
    n = len(binary_data)
    circle = np.zeros((n, n))  # Create a circular "fabric"

    # Weave data inward
    for i in range(iterations):
        idx = i % n  # Circular index
        circle[idx, :] += binary_data
        binary_data = np.roll(binary_data, -1)  # Twist the string inward

    return circle

# Visualize the circular weaving
def visualize_weaving(circle, title):
    plt.figure(figsize=(8, 8))
    plt.imshow(circle, cmap="viridis", aspect="auto")
    plt.title(title)
    plt.colorbar(label="Weaving Intensity")
    plt.show()

# Example: SHA-256 for 'abc'
hash_hex = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Simulate circular weaving
circle = circular_weaving(hash_hex)

# Visualize the weaving pattern
visualize_weaving(circle, title="Circular Weaving of SHA-256 ('abc')")
