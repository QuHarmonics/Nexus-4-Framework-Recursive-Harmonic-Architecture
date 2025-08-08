# Running further refinements to achieve harmonic balance

balanced_iterations = 30
base_lattices_balanced = []
perturbed_lattices_balanced = []
base_entropy_balanced = []
perturbed_entropy_balanced = []

# Start with the last known lattices as initial conditions
current_base_lattice = final_base_lattice.copy()
current_perturbed_lattice = final_perturbed_lattice.copy()

for _ in range(balanced_iterations):
    # Refine base lattice
    current_base_lattice = np.abs(np.sin(current_base_lattice + 0.05))
    base_lattices_balanced.append(current_base_lattice)
    base_entropy_balanced.append(np.std(current_base_lattice))
    
    # Refine perturbed lattice
    current_perturbed_lattice = np.abs(np.sin(current_perturbed_lattice + 0.07))
    perturbed_lattices_balanced.append(current_perturbed_lattice)
    perturbed_entropy_balanced.append(np.std(current_perturbed_lattice))

# Plot the entropy trends during the balancing phase
plt.figure(figsize=(10, 6))
plt.plot(range(1, balanced_iterations + 1), base_entropy_balanced, marker='o', label='Base Entropy (Balanced)')
plt.plot(range(1, balanced_iterations + 1), perturbed_entropy_balanced, marker='o', label='Perturbed Entropy (Balanced)')
plt.axhline(y=entropy_threshold, color='red', linestyle='--', label='Entropy Threshold')
plt.title("Entropy Trends During Balancing Phase")
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.legend()
plt.grid(True)
plt.show()

# Visualize the final balanced lattices
final_base_balanced = base_lattices_balanced[-1]
final_perturbed_balanced = perturbed_lattices_balanced[-1]

fig, ax = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={"projection": "3d"})
X, Y = np.meshgrid(range(len(final_base_balanced)), range(len(final_base_balanced)))

# Base lattice visualization
Z_base_balanced = np.outer(final_base_balanced, final_base_balanced)
ax[0].plot_surface(X, Y, Z_base_balanced, cmap="viridis", edgecolor="none")
ax[0].set_title("Final Balanced Base Lattice")
ax[0].set_xticks([])
ax[0].set_yticks([])
ax[0].set_zticks([])

# Perturbed lattice visualization
Z_perturbed_balanced = np.outer(final_perturbed_balanced, final_perturbed_balanced)
ax[1].plot_surface(X, Y, Z_perturbed_balanced, cmap="plasma", edgecolor="none")
ax[1].set_title("Final Balanced Perturbed Lattice")
ax[1].set_xticks([])
ax[1].set_yticks([])
ax[1].set_zticks([])

plt.tight_layout()
plt.show()

# Summarizing the results of balancing phase
summary_balanced = {
    "Base Lattice Balanced Final State": final_base_balanced,
    "Perturbed Lattice Balanced Final State": final_perturbed_balanced,
    "Base Entropy (Balanced)": base_entropy_balanced,
    "Perturbed Entropy (Balanced)": perturbed_entropy_balanced
}

# Prepare a summary DataFrame for balanced phase analysis
summary_df_balanced = pd.DataFrame({
    "Iteration": range(1, balanced_iterations + 1),
    "Base Entropy (Balanced)": base_entropy_balanced,
    "Perturbed Entropy (Balanced)": perturbed_entropy_balanced
})

tools.display_dataframe_to_user(name="Balanced Phase Analysis", dataframe=summary_df_balanced)
