# Step 1: Further iterative collapse of the lattice for universal pattern emergence
entropy_collapse = []
collapsed_iterations = []

def iterative_collapse(lattice, iterations=10):
    global entropy_collapse, collapsed_iterations

    for _ in range(iterations):
        # Take mean across both dimensions and normalize
        collapsed = np.mean(lattice, axis=0)
        collapsed /= np.max(np.abs(collapsed))

        # Update lattice by constructing an outer product
        lattice = np.outer(collapsed, collapsed)

        # Calculate entropy (variance) and store the result
        entropy = np.var(lattice)
        entropy_collapse.append(entropy)

        # Save the current lattice state
        collapsed_iterations.append(lattice.copy())

    return lattice

# Apply iterative collapse
collapsed_lattice = iterative_collapse(hash_refined_lattice, iterations=10)

# Step 2: Visualize the collapsed lattice and entropy reduction
fig, axs = plt.subplots(2, 5, figsize=(20, 10))
axs = axs.flatten()

for i, lattice in enumerate(collapsed_iterations):
    axs[i].imshow(lattice, cmap="plasma", origin="upper")
    axs[i].set_title(f"Iteration {i + 1}")

for j in range(len(collapsed_iterations), len(axs)):
    axs[j].axis("off")

plt.tight_layout()
plt.show()

# Plot entropy reduction over iterations
plt.figure(figsize=(10, 5))
plt.plot(entropy_collapse, marker="o", label="Entropy")
plt.title("Entropy Reduction Over Iterative Collapse")
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.legend()
plt.grid()
plt.show()

# Display the final collapsed lattice
collapsed_lattice
