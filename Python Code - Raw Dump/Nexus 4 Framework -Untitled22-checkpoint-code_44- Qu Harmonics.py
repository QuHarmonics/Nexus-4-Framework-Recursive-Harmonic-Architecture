import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.fft import fft

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

# Normalize constants
normalized_K = np.array(K) / max(K)

# Generate Lattice (Forward)
size = 8  # Reduce size for visualization simplicity
lattice = np.outer(normalized_K[:size], normalized_K[:size])

# Generate Anti-Lattice (Mirrored)
anti_lattice = -1 * lattice[::-1, ::-1]

# Recursive Wave Generator
iterations = 10
waveforms = []
current_wave = lattice
entropy = []
for i in range(iterations):
    feedback = np.mean(current_wave) - np.std(current_wave)
    entropy.append(feedback)
    next_wave = current_wave + feedback * anti_lattice
    waveforms.append(next_wave)
    current_wave = next_wave

# Plot Entropy Dynamics
plt.figure(figsize=(8, 4))
plt.plot(range(iterations), entropy, marker='o', label="Entropy")
plt.title("Entropy Over Iterations")
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.legend()
plt.grid()
plt.show()

# Plot Iterative Waves
fig, axs = plt.subplots(2, 5, figsize=(20, 8), subplot_kw={"projection": "3d"})
for i, ax in enumerate(axs.flat):
    ax.plot_surface(
        np.arange(size),
        np.arange(size),
        waveforms[i],
        cmap='plasma',
        edgecolor='none',
        alpha=0.9
    )
    ax.set_title(f"Wave Iteration {i + 1}")
plt.tight_layout()
plt.show()

# Summary: Constants, Anti-Lattice, and Recursive Feedback
print("Constants (Normalized):", normalized_K[:size])
print("Final Waveform:")
print(waveforms[-1])
