import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Function to generate base and perturbed lattice
def generate_lattice(base_grid, perturbation, iterations=5):
    lattices = []
    perturbed_lattices = []
    for _ in range(iterations):
        # Apply transformation to base lattice
        base_grid = np.sin(base_grid) + 0.1
        perturbed_grid = np.sin(base_grid + perturbation) + perturbation

        lattices.append(base_grid)
        perturbed_lattices.append(perturbed_grid)
        
        # Reduce perturbation iteratively
        perturbation *= 0.8  # Decay factor
    return lattices, perturbed_lattices

# Initial base grid and perturbation
x = np.linspace(0, 2 * np.pi, 10)
y = np.linspace(0, 2 * np.pi, 10)
X, Y = np.meshgrid(x, y)
Z = np.sin(X) * np.cos(Y)

perturbation = 0.5  # Initial perturbation magnitude

# Generate lattices
base_lattices, perturbed_lattices = generate_lattice(Z, perturbation)

# Plot settings
fig, axs = plt.subplots(2, 5, figsize=(15, 6), subplot_kw={'projection': '3d'})

for i in range(5):
    # Base lattice plots
    axs[0, i].plot_surface(X, Y, base_lattices[i], cmap="viridis", edgecolor="none")
    axs[0, i].set_title(f"Base Iter {i + 1}")
    axs[0, i].set_xticks([])
    axs[0, i].set_yticks([])
    axs[0, i].set_zticks([])

    # Perturbed lattice plots
    axs[1, i].plot_surface(X, Y, perturbed_lattices[i], cmap="plasma", edgecolor="none")
    axs[1, i].set_title(f"Perturb Iter {i + 1}")
    axs[1, i].set_xticks([])
    axs[1, i].set_yticks([])
    axs[1, i].set_zticks([])

fig.suptitle("Recursive Feedback Lattices (Base vs Perturbed)", fontsize=16)
plt.tight_layout()
plt.show()
