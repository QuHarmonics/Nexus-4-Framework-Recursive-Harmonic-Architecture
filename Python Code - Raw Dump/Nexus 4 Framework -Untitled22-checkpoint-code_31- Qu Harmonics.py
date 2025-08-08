import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Parameters for recursion
iterations = 100
lattice_size = 8
entropy_threshold = 0.05

# Initialize base and perturbed lattices
base_lattice = np.linspace(0.1, 1.0, lattice_size)
perturbed_lattice = base_lattice + 0.2  # Initial perturbation

# Storage for entropy values
base_lattices_over_time = []
perturbed_lattices_over_time = []
base_entropy_values = []
perturbed_entropy_values = []

# Recursive feedback calculations
for _ in range(iterations):
    # Base lattice recursion
    base_lattice = np.abs(np.sin(base_lattice))
    base_lattices_over_time.append(base_lattice)
    base_entropy = np.sum(np.abs(np.diff(base_lattice)))
    base_entropy_values.append(base_entropy)

    # Perturbed lattice recursion
    perturbed_lattice = np.abs(np.sin(perturbed_lattice + 0.1))
    perturbed_lattices_over_time.append(perturbed_lattice)
    perturbed_entropy = np.sum(np.abs(np.diff(perturbed_lattice)))
    perturbed_entropy_values.append(perturbed_entropy)

# Adjusting the subplot arrangement for proper visualization
fig, axes = plt.subplots(2, 5, figsize=(20, 10), subplot_kw={"projection": "3d"})
fig.suptitle("Recursive Feedback Lattices (Base vs Perturbed)", fontsize=16)

for i in range(min(iterations, 5)):
    # Base lattice
    X_base, Y_base = np.meshgrid(range(len(base_lattice)), range(len(base_lattice)))
    Z_base = np.outer(base_lattices_over_time[i], base_lattices_over_time[i])
    axes[0, i].plot_surface(X_base, Y_base, Z_base, cmap="viridis", edgecolor="none")
    axes[0, i].set_title(f"Base Iter {i+1}", fontsize=10)
    axes[0, i].set_xticks([])
    axes[0, i].set_yticks([])
    axes[0, i].set_zticks([])

    # Perturbed lattice
    Z_perturbed = np.outer(perturbed_lattices_over_time[i], perturbed_lattices_over_time[i])
    axes[1, i].plot_surface(X_base, Y_base, Z_perturbed, cmap="plasma", edgecolor="none")
    axes[1, i].set_title(f"Perturb Iter {i+1}", fontsize=10)
    axes[1, i].set_xticks([])
    axes[1, i].set_yticks([])
    axes[1, i].set_zticks([])

plt.tight_layout()
plt.subplots_adjust(top=0.85)
plt.show()

# Plot entropy trends for both systems
plt.figure(figsize=(10, 5))
plt.plot(range(1, iterations + 1), base_entropy_values, marker='o', label='Base Lattice Entropy')
plt.plot(range(1, iterations + 1), perturbed_entropy_values, marker='o', label='Perturbed Lattice Entropy')
plt.axhline(y=entropy_threshold, color='red', linestyle='--', label='Entropy Threshold')
plt.title("Entropy Trends Across Feedback Iterations")
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.legend()
plt.grid(True)
plt.show()

# Prepare summary results for display
summary_df = pd.DataFrame({
    "Iteration": range(1, iterations + 1),
    "Base Lattice Entropy": base_entropy_values,
    "Perturbed Lattice Entropy": perturbed_entropy_values,
})
summary_df["Entropy Difference"] = np.abs(summary_df["Base Lattice Entropy"] - summary_df["Perturbed Lattice Entropy"])

# Display summary results
import ace_tools as tools  # Replace with your own dataframe display function if needed
tools.display_dataframe_to_user(name="Recursive Feedback Analysis Summary", dataframe=summary_df)
