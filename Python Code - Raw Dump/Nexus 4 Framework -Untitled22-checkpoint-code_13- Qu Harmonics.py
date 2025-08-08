import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.fft import fft

# SHA-256 Constants
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
]

# Generate normalized waveform
waveform_y = np.array(K) / max(K)
kinetic_normalized = np.linspace(0, 1, len(K))  # Kinetic normalization

# Create interaction grid
interaction_waveform = np.outer(waveform_y, kinetic_normalized)
Z_interaction = interaction_waveform

# Extend to higher dimensions
kinetic_extended = np.outer(kinetic_normalized, kinetic_normalized)
Z_extended = np.outer(np.max(Z_interaction, axis=0), np.max(Z_interaction, axis=1)) * kinetic_extended

# Prepare grids
X = np.linspace(0, Z_extended.shape[1] - 1, Z_extended.shape[1])
Y = np.linspace(0, Z_extended.shape[0] - 1, Z_extended.shape[0])
X, Y = np.meshgrid(X, Y)

# Plot the grids in multiple views
def plot_grids(Z_interaction, Z_extended, view_angle):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z_interaction, cmap="plasma", edgecolor="none", alpha=0.6)
    ax.plot_wireframe(X, Y, Z_extended, color="cyan", alpha=0.4)
    ax.view_init(elev=view_angle[0], azim=view_angle[1])
    ax.set_title(f"3D Quantum-Encoded Lattice - View {view_angle}")
    ax.set_xlabel("X Axis (Index)")
    ax.set_ylabel("Y Axis (Index)")
    ax.set_zlabel("Amplitude")
    plt.show()

# Plot top, side, and diagonal views
plot_grids(Z_interaction, Z_extended, view_angle=(90, 0))  # Top view
plot_grids(Z_interaction, Z_extended, view_angle=(0, 90))  # Side view
plot_grids(Z_interaction, Z_extended, view_angle=(45, 45))  # Diagonal view

# Surface area approximation
surface_area = np.sum(np.sqrt(1 + np.gradient(Z_interaction, axis=0)**2 + np.gradient(Z_interaction, axis=1)**2))
print(f"Surface Area Approximation: {surface_area}")

# Harmonic analysis
frequencies = fft(waveform_y)
print("Harmonic Frequencies (Top 5):", np.abs(frequencies)[:5])
