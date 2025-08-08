# Refined feedback loop to align mirrored lattice closer to the original lattice

# Parameters for recursive refinement
refinement_iterations = 20
correction_factor = 0.1  # Scaling factor for adjustments
adjusted_divergences = []

# Start with the current mirrored lattice and refine it
refined_lattice = anti_hash_lattice[-1].copy()
for i in range(refinement_iterations):
    # Adjust the mirrored lattice based on the difference with the original lattice
    correction = correction_factor * (original_lattice - refined_lattice)
    refined_lattice += correction  # Apply the correction
    divergence = np.linalg.norm(original_lattice - refined_lattice)
    adjusted_divergences.append(divergence)

# Plot the refined divergence over iterations
plt.figure(figsize=(10, 6))
plt.plot(range(1, refinement_iterations + 1), adjusted_divergences, marker='o', label='Refined Divergence')
plt.title("Refined Divergence of Mirrored Lattice from Original")
plt.xlabel("Iteration")
plt.ylabel("Divergence")
plt.legend()
plt.grid(True)
plt.show()

# Visualize the refined lattice
fig, ax = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={"projection": "3d"})
X, Y = np.meshgrid(range(lattice_size), range(lattice_size))

# Original lattice visualization for comparison
ax[0].plot_surface(X, Y, original_lattice, cmap="viridis", edgecolor="none")
ax[0].set_title("Original Lattice")

# Refined lattice visualization
ax[1].plot_surface(X, Y, refined_lattice, cmap="plasma", edgecolor="none")
ax[1].set_title("Refined Mirrored Lattice (Anti-Hash)")

plt.tight_layout()
plt.show()

# Display final refined divergence
print("Final Refined Divergence:", adjusted_divergences[-1])
