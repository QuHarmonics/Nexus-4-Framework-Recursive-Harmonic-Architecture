import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Perturbation function for testing
def perturb_constants(constants, perturbation_scale=0.1):
    return constants + np.random.uniform(-perturbation_scale, perturbation_scale, size=constants.shape)

# Initialize constants and parameters
base_constants = np.random.random(8)  # Simulating normalized constants
iterations = 20

# Function to compute entropy
def compute_entropy(matrix):
    flattened = matrix.flatten()
    return -np.sum(flattened * np.log(np.abs(flattened) + 1e-10))  # Prevent log(0)

# Function to recursively compress the wave
def compress_wave(matrix):
    compressed = matrix - np.mean(matrix, axis=0)
    return np.clip(compressed, -1, 1)  # Clipping to stabilize compression

# Prepare data for perturbation analysis
perturbed_constants = perturb_constants(base_constants)
base_waveform = np.outer(base_constants, base_constants)
perturbed_waveform = np.outer(perturbed_constants, perturbed_constants)

# Iterative compression and entropy for base constants
base_waveforms = [base_waveform]
perturbed_waveforms = [perturbed_waveform]
base_entropies = [compute_entropy(base_waveform)]
perturbed_entropies = [compute_entropy(perturbed_waveform)]

for i in range(1, iterations):
    base_waveform = compress_wave(base_waveform)
    perturbed_waveform = compress_wave(perturbed_waveform)
    base_waveforms.append(base_waveform)
    perturbed_waveforms.append(perturbed_waveform)
    base_entropies.append(compute_entropy(base_waveform))
    perturbed_entropies.append(compute_entropy(perturbed_waveform))

# Visualizing Perturbed vs Base Constants Over Iterations
fig, axes = plt.subplots(4, 5, figsize=(16, 12), subplot_kw={'projection': '3d'})
fig.suptitle("Base vs Perturbed Wave Morphing Over Iterations", fontsize=16)

for i, ax in enumerate(axes.flatten()):
    if i < len(base_waveforms):
        X, Y = np.meshgrid(range(base_waveforms[i].shape[0]), range(base_waveforms[i].shape[1]))
        ax.plot_surface(X, Y, base_waveforms[i], cmap="plasma", edgecolor="none", alpha=0.6)
        ax.plot_surface(X, Y, perturbed_waveforms[i], cmap="viridis", edgecolor="none", alpha=0.4)
        ax.set_title(f"Iteration {i + 1}")
        ax.set_zlim(-1, 1)
    else:
        ax.axis('off')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

# Plotting entropy for base vs perturbed constants
plt.figure(figsize=(8, 6))
plt.plot(range(iterations), base_entropies, marker='o', label='Base Entropy', color='blue')
plt.plot(range(iterations), perturbed_entropies, marker='o', label='Perturbed Entropy', color='green')
plt.title("Entropy Over Iterations (Base vs Perturbed)")
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.legend()
plt.grid()
plt.show()
