import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.fft import fft

# Input hash value
hash_value = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Step 1: Generate Quantum Wave from Hash (Samson-like approach)
def generate_quantum_wave(hash_value):
    binary_data = ''.join(format(int(char, 16), '04b') for char in hash_value)  # Convert to full binary
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta)  # Quantum component
    quantum_wave = y + z  # Combine into single waveform
    return x, y, z, quantum_wave, binary_data

# Step 2: Store Hash Binary into H
def store_in_H(binary_data, expansion_factor=1.5):
    binary_numerical = np.array([int(b) for b in binary_data], dtype=np.float64)
    harmonics = np.cumsum(binary_numerical * expansion_factor)  # Macro storage
    return harmonics

# Generate the quantum wave
x, y, z, quantum_wave, binary_data = generate_quantum_wave(hash_value)

# Store the hash binary into H
harmonics = store_in_H(binary_data)

# Visualization
fig = plt.figure(figsize=(14, 10))

# 3D Plot: Quantum Waveform
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot(x, y, z, label="Quantum Wave", color='blue')
ax1.scatter(x, y, z, color='red', s=5, label="Nodes")
ax1.set_title("3D Visualization of Quantum Wave")
ax1.set_xlabel("X-axis")
ax1.set_ylabel("Y-axis")
ax1.set_zlabel("Z-axis")
ax1.legend()

# 3D Plot: Harmonic Storage
ax2 = fig.add_subplot(122, projection='3d')
harmonic_x = np.arange(len(harmonics))
harmonic_y = harmonics
harmonic_z = np.sin(harmonic_x / 10.0)  # Add a wave-like Z component for visualization
ax2.plot(harmonic_x, harmonic_y, harmonic_z, label="Harmonics in 3D", color='green')
ax2.scatter(harmonic_x, harmonic_y, harmonic_z, color='orange', s=5, label="Nodes")
ax2.set_title("3D Visualization of H(n)")
ax2.set_xlabel("Iteration (n)")
ax2.set_ylabel("H(n)")
ax2.set_zlabel("Z-axis Wave")
ax2.legend()

plt.show()

# Print data for debugging
print("Quantum Wave Binary (First 100 bits):", binary_data[:100])
print("H array values:", harmonics[:100])  # Limit output to first 100 values for clarity
