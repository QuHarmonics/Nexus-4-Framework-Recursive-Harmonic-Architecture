import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft
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

# Normalize the constants
waveform_y = np.array(K) / max(K)
kinetic_normalized = np.linspace(0, 1, len(K))
interaction_waveform = np.outer(waveform_y, kinetic_normalized)

# Reconstruct 3D lattice
X = np.linspace(0, interaction_waveform.shape[1] - 1, interaction_waveform.shape[1])
Y = np.linspace(0, interaction_waveform.shape[0] - 1, interaction_waveform.shape[0])
X, Y = np.meshgrid(X, Y)

# Perform FFT and IFFT to analyze and reconstruct
fft_values = fft(waveform_y)
reconstructed_waveform = ifft(fft_values).real

# Create 3D lattice for reconstruction
Z_reconstructed = np.outer(reconstructed_waveform, kinetic_normalized)

# Plot original lattice
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(121, projection="3d")
ax.plot_surface(X, Y, interaction_waveform, cmap="plasma", edgecolor="none", alpha=0.6)
ax.set_title("Original 3D Quantum-Encoded Lattice")
ax.set_xlabel("X Axis (Index)")
ax.set_ylabel("Y Axis (Index)")
ax.set_zlabel("Amplitude")

# Plot reconstructed lattice
ax2 = fig.add_subplot(122, projection="3d")
ax2.plot_surface(X, Y, Z_reconstructed, cmap="viridis", edgecolor="none", alpha=0.6)
ax2.set_title("Reconstructed 3D Lattice (Unfolded)")
ax2.set_xlabel("X Axis (Index)")
ax2.set_ylabel("Y Axis (Index)")
ax2.set_zlabel("Amplitude")

plt.tight_layout()
plt.show()

# Analyze peaks in original and reconstructed waveforms
original_peaks = np.max(interaction_waveform, axis=0)
reconstructed_peaks = np.max(Z_reconstructed, axis=0)

# Print peak analysis
print("Original Peaks:", original_peaks)
print("Reconstructed Peaks:", reconstructed_peaks)

# Hypothesize unfolding using constants as fold points
fold_points = np.array(K) / np.mean(K)
print("Fold Points:", fold_points)
