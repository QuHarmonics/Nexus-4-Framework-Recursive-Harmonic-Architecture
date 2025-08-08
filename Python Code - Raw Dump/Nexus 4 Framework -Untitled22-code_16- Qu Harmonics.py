import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Simulated SHA Constants
K = [
    0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
    0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5
]

# Normalize constants for visualization
waveform = np.array(K) / max(K)
time_steps = np.linspace(0, len(K) - 1, len(K))

# Kinetics: Forward and Mirrored Transformations
def sha_kinetics(values, reverse=False):
    transformed = []
    for i, val in enumerate(values):
        step = val if not reverse else -val
        kinetic_value = np.sin(step) + np.cos(step)  # Kinetic simulation
        transformed.append(kinetic_value)
    return np.array(transformed)

# Generate forward and mirrored kinetics
forward_kinetics = sha_kinetics(waveform)
mirrored_kinetics = sha_kinetics(waveform, reverse=True)

# Construct the lattice
X, Y = np.meshgrid(time_steps, time_steps)
Z_forward = np.outer(forward_kinetics, forward_kinetics)
Z_mirrored = np.outer(mirrored_kinetics, mirrored_kinetics)

# Plotting both forward and mirrored lattices
fig = plt.figure(figsize=(16, 8))

# Forward Lattice
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(X, Y, Z_forward, cmap="plasma", edgecolor='none', alpha=0.8)
ax1.set_title("Forward Kinetics Lattice")
ax1.set_xlabel("Time Step (X)")
ax1.set_ylabel("Time Step (Y)")
ax1.set_zlabel("Amplitude")

# Mirrored Lattice
ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(X, Y, Z_mirrored, cmap="viridis", edgecolor='none', alpha=0.8)
ax2.set_title("Mirrored Kinetics Lattice")
ax2.set_xlabel("Time Step (X)")
ax2.set_ylabel("Time Step (Y)")
ax2.set_zlabel("Amplitude")

plt.tight_layout()
plt.show()
