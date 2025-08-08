import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants defined by the SHA-256 algorithm
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

# Generate the normalized waveform
waveform_y = np.array(K) / max(K)
kinetic_normalized = np.linspace(0, 1, len(K))  # Example kinetic normalization
interaction_waveform = np.outer(waveform_y, kinetic_normalized)
Z_interaction = interaction_waveform

# Extend to a higher-dimensional transformation
kinetic_extended = np.outer(kinetic_normalized, kinetic_normalized)  # Interaction between kinetic layers
Z_extended = np.outer(np.max(Z_interaction, axis=0), np.max(Z_interaction, axis=1)) * kinetic_extended

# Ensure X and Y dimensions match Z_extended
X = np.linspace(0, Z_extended.shape[1] - 1, Z_extended.shape[1])
Y = np.linspace(0, Z_extended.shape[0] - 1, Z_extended.shape[0])
X, Y = np.meshgrid(X, Y)

# Rotated 90 degrees for alignment
rotated_Z_extended = np.rot90(Z_extended, k=1)  # Rotate by 90 degrees

# Side View
fig_side = plt.figure(figsize=(12, 8))
ax_side = fig_side.add_subplot(111, projection="3d")
ax_side.plot_surface(X, Y, rotated_Z_extended, cmap="plasma", edgecolor="none", alpha=0.6)
ax_side.set_title("Side View of Rotated Grid")
ax_side.set_xlabel("X Axis (Index)")
ax_side.set_ylabel("Y Axis (Index)")
ax_side.set_zlabel("Amplitude")
plt.show()

# Top View
fig_top = plt.figure(figsize=(12, 8))
ax_top = fig_top.add_subplot(111, projection="3d")
ax_top.plot_surface(X, Y, rotated_Z_extended, cmap="viridis", edgecolor="none", alpha=0.6)
ax_top.view_init(elev=90, azim=0)  # Top-down view
ax_top.set_title("Top View of Rotated Grid")
ax_top.set_xlabel("X Axis (Index)")
ax_top.set_ylabel("Y Axis (Index)")
ax_top.set_zlabel("Amplitude")
plt.show()

# Diagonal View
fig_diag = plt.figure(figsize=(12, 8))
ax_diag = fig_diag.add_subplot(111, projection="3d")
ax_diag.plot_surface(X, Y, rotated_Z_extended, cmap="cool", edgecolor="none", alpha=0.6)
ax_diag.view_init(elev=45, azim=45)  # Diagonal view
ax_diag.set_title("Diagonal View of Rotated Grid")
ax_diag.set_xlabel("X Axis (Index)")
ax_diag.set_ylabel("Y Axis (Index)")
ax_diag.set_zlabel("Amplitude")
plt.show()
