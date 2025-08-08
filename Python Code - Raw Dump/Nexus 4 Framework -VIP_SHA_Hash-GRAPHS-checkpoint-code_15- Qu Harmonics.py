# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants from SHA-256
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

# Normalize constants for visualization
waveform_y = np.array(K) / max(K)
kinetic_normalized = np.linspace(0, 1, len(K))

# Generate the lattice structures
interaction_waveform = np.outer(waveform_y, kinetic_normalized)
Z_interaction = interaction_waveform

# Mirrored lattice
Z_mirrored = -Z_interaction

# Refined lattice interaction
Z_combined = Z_interaction + Z_mirrored

# Create meshgrid for 3D plots
X = np.linspace(0, len(K) - 1, len(K))
Y = np.linspace(0, len(K) - 1, len(K))
X, Y = np.meshgrid(X, Y)

# Plot Forward Kinetics Lattice
fig = plt.figure(figsize=(16, 10))

# Forward Kinetics
ax1 = fig.add_subplot(131, projection='3d')
ax1.plot_surface(X, Y, Z_interaction, cmap='viridis', edgecolor='none')
ax1.set_title("Forward Kinetics Lattice")
ax1.set_xlabel("X Axis")
ax1.set_ylabel("Y Axis")
ax1.set_zlabel("Amplitude")

# Mirrored Kinetics
ax2 = fig.add_subplot(132, projection='3d')
ax2.plot_surface(X, Y, Z_mirrored, cmap='plasma', edgecolor='none')
ax2.set_title("Mirrored Kinetics Lattice")
ax2.set_xlabel("X Axis")
ax2.set_ylabel("Y Axis")
ax2.set_zlabel("Amplitude")

# Combined Interaction Lattice
ax3 = fig.add_subplot(133, projection='3d')
ax3.plot_surface(X, Y, Z_combined, cmap='coolwarm', edgecolor='none')
ax3.set_title("Refined Interaction Lattice")
ax3.set_xlabel("X Axis")
ax3.set_ylabel("Y Axis")
ax3.set_zlabel("Amplitude")

plt.tight_layout()
plt.show()

# Surface Area Calculations
surface_area_forward = np.sum(np.sqrt(1 + np.gradient(Z_interaction, axis=0)**2 + np.gradient(Z_interaction, axis=1)**2))
surface_area_mirrored = np.sum(np.sqrt(1 + np.gradient(Z_mirrored, axis=0)**2 + np.gradient(Z_mirrored, axis=1)**2))
surface_area_combined = np.sum(np.sqrt(1 + np.gradient(Z_combined, axis=0)**2 + np.gradient(Z_combined, axis=1)**2))

results = {
    "Surface Area - Forward Kinetics": surface_area_forward,
    "Surface Area - Mirrored Kinetics": surface_area_mirrored,
    "Surface Area - Combined Interaction": surface_area_combined,
}

results
