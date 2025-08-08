# --- Full setup ---

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Parameters
lattice_size = 50
balanced_iterations = 30
entropy_threshold = 0.1  # Define your threshold

# Create synthetic final lattices (simulate results of previous stage)
final_base_lattice = np.sin(np.linspace(0, 2 * np.pi, lattice_size))
final_perturbed_lattice = np.sin(np.linspace(0, 2 * np.pi, lattice_size)) + 0.05 * np.random.randn(lattice_size)

# Storage for results
base_lattices_balanced = []
perturbed_lattices_balanced = []
base_entropy_balanced = []
perturbed_entropy_balanced = []

# Start with the final known lattices as initial conditions
current_base_lattice = final_base_lattice.copy()
current_perturbed_lattice = final_perturbed_lattice.copy()

# --- Balancing Loop ---

for _ in range(balanced_iterations):
    # Refine base lattice
    current_base_lattice = np.abs(np.sin(current_base_lattice + 0.05))
    base_lattices_balanced.append(current_base_lattice)
    base_entropy_balanced.append(np.std(current_base_lattice))
    
    # Refine perturbed lattice
    current_perturbed_lattice = np.abs(np.sin(current_perturbed_lattice + 0.07))
    perturbed_lattices_balanced.append(current_perturbed_lattice)
    perturbed_entropy_balanced.append(np.std(current_perturbed_lattice))

# --- Plot Entropy Trends ---

plt.figure(figsize=(10, 6))
plt.plot(range(1, balanced_iterations + 1), base_entropy_balanced, marker='o', label='Base Entropy (Balanced)')
plt.plot(range(1, balanced_iterations + 1), perturbed_entropy_balanced, marker='o', label='Perturbed Entropy (Balanced)')
plt.axhline(y=entropy_threshold, color='red', linestyle='--', label='Entropy Threshold')
plt.title("Entropy Trends During Balancing Phase")
plt.xlabel("Iteration")
plt.ylabel("Entropy (Standard Deviation)")
plt.legend()
plt.grid(True)
plt.show()

# --- Final Balanced Lattices ---

final_base_balanced = base_lattices_balanced[-1]
final_perturbed_balanced = perturbed_lattices_balanced[-1]

# --- 3D Visualization ---

fig, ax = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={"projection": "3d"})

# Create meshgrids
X, Y = np.meshgrid(range(lattice_size), range(lattice_size))

# For visualization, outer products create 2D surfaces
Z_base_balanced = np.outer(final_base_balanced, final_base_balanced)
Z_perturbed_balanced = np.outer(final_perturbed_balanced, final_perturbed_balanced)

# Plot base lattice
ax[0].plot_surface(X, Y, Z_base_balanced, cmap="viridis", edgecolor="none")
ax[0].set_title("Final Balanced Base Lattice")
ax[0].set_xticks([])
ax[0].set_yticks([])
ax[0].set_zticks([])

# Plot perturbed lattice
ax[1].plot_surface(X, Y, Z_perturbed_balanced, cmap="plasma", edgecolor="none")
ax[1].set_title("Final Balanced Perturbed Lattice")
ax[1].set_xticks([])
ax[1].set_yticks([])
ax[1].set_zticks([])

plt.tight_layout()
plt.show()

# --- Summary DataFrame ---

summary_df_balanced = pd.DataFrame({
    "Iteration": range(1, balanced_iterations + 1),
    "Base Entropy (Balanced)": base_entropy_balanced,
    "Perturbed Entropy (Balanced)": perturbed_entropy_balanced
})

# Display summary table
print("\nSummary of Entropy Over Balancing Iterations:")
print(summary_df_balanced)

# Optional: show last entropy values
print("\nFinal Base Entropy:", base_entropy_balanced[-1])
print("Final Perturbed Entropy:", perturbed_entropy_balanced[-1])
