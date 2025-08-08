import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Original binary data (replace this with your hash data)
binary_data = np.array([7, 83, 233, 188, 222, 188, 190, 27, 187, 241], dtype=np.uint8)

# Step 1: Harmonic Expansion Function
def harmonic_expand(data, scale_factor=1.5):
    harmonics = np.cumsum(data.astype(np.float64) * scale_factor)
    return harmonics

# Harmonic expansion
harmonics = harmonic_expand(binary_data)

# Step 2: Reconstruct Wave for 3D Visualization
def reconstruct_wave(harmonics, depth=10):
    x = np.arange(len(harmonics))
    y = harmonics
    z = np.sin(x / depth) * np.cos(y / (2 * depth))
    return x, y, z

# Reconstruct the wave
x, y, z = reconstruct_wave(harmonics)

# 2D Plot
plt.figure(figsize=(10, 6))
plt.bar(np.arange(len(harmonics)), harmonics, label="Reconstructed Wave")
plt.title("Harmonic Expansion of Hash")
plt.xlabel("Iteration Step")
plt.ylabel("Amplitude")
plt.legend()
plt.show()

# 3D Plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x, y, z, label="Harmonic Lattice in 3D", color='blue')
ax.scatter(x, y, z, color='red', s=5, label="Nodes")
ax.set_title("3D Visualization of Harmonic Expansion", fontsize=16)
ax.set_xlabel("Iteration Step", fontsize=12)
ax.set_ylabel("Amplitude (H)", fontsize=12)
ax.set_zlabel("Depth Waveform", fontsize=12)
ax.legend()
plt.show()

# Output results
print("Original Binary Data (First 10 bytes):", binary_data[:10])
print("Reconstructed Wave Data (First 10 values):", harmonics[:10])
