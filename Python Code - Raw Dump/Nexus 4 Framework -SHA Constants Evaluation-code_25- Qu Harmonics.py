import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Initialize constants and parameters
constants = np.random.random(8)  # Simulating normalized constants
iterations = 20

# Initialize wave lattice
waveform = np.outer(constants, constants)

# Function to compute entropy
def compute_entropy(matrix):
    flattened = matrix.flatten()
    return -np.sum(flattened * np.log(np.abs(flattened) + 1e-10))  # Prevent log(0)

# Function to recursively compress the wave
def compress_wave(matrix):
    compressed = matrix - np.mean(matrix, axis=0)
    return np.clip(compressed, -1, 1)  # Clipping to stabilize compression

# Prepare data for visualization
waveforms = [waveform]
entropies = [compute_entropy(waveform)]

for i in range(1, iterations):
    waveform = compress_wave(waveform)
    waveforms.append(waveform)
    entropies.append(compute_entropy(waveform))

# Plotting the wave evolution over iterations
fig, axes = plt.subplots(4, 5, figsize=(16, 12), subplot_kw={'projection': '3d'})
fig.suptitle("Wave Morphing Over Iterations", fontsize=16)
for i, ax in enumerate(axes.flatten()):
    if i < len(waveforms):
        X, Y = np.meshgrid(range(waveforms[i].shape[0]), range(waveforms[i].shape[1]))
        ax.plot_surface(X, Y, waveforms[i], cmap="plasma", edgecolor="none", alpha=0.8)
        ax.set_title(f"Wave Iteration {i + 1}")
        ax.set_zlim(-1, 1)
    else:
        ax.axis('off')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

# Plotting entropy over iterations
plt.figure(figsize=(8, 6))
plt.plot(range(iterations), entropies, marker='o', label='Entropy')
plt.title("Entropy Over Iterations")
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.legend()
plt.grid()
plt.show()
