import numpy as np
import matplotlib.pyplot as plt

# Generate a spiral representation of a hash
def hash_to_spiral(hash_hex, resolution=1000):
    binary_data = np.array([int(bit) for bit in ''.join(format(int(c, 16), '04b') for c in hash_hex)])
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi * (n / resolution), n)
    r = np.linspace(0, 1, n)  # Radius shrinking toward the center
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = binary_data - 0.5  # Map binary to waveform
    return x, y, z

# Recursive unwinding of the spiral
def unwind_spiral(x, y, z, iterations=1000, alpha=0.1):
    unwound_z = z.copy()
    for i in range(iterations):
        # Align outward from the center
        corrections = alpha * np.sin(np.linspace(0, 2 * np.pi, len(z)))
        unwound_z += corrections
        unwound_z = np.clip(unwound_z, -0.5, 0.5)  # Keep within bounds
    return x, y, unwound_z

# Visualize the spiral
def visualize_spiral(x, y, z, title):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, z, label=title)
    ax.set_title(title)
    ax.legend()
    plt.show()

# Example: SHA-256 for 'abc'
hash_hex = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Generate spiral
x, y, z = hash_to_spiral(hash_hex)

# Visualize original spiral
visualize_spiral(x, y, z, title="Original Hash Spiral")

# Unwind spiral
unwound_x, unwound_y, unwound_z = unwind_spiral(x, y, z)

# Visualize unwound spiral
visualize_spiral(unwound_x, unwound_y, unwound_z, title="Unwound Spiral")
