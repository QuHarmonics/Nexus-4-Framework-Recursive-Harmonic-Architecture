import numpy as np
import matplotlib.pyplot as plt

# Initialize constants and setup iterative refinement process
constants = np.array([0.27264283, 0.46389402, 0.74472339, 0.9576116, 0.23449206, 0.36852961, 0.59942198, 0.70111437])

# Define initial lattice (waveform)
lattice = np.outer(constants, constants)

# Track entropy for each iteration
entropy_values = []
wave_iterations = []

# Recursive refinement function
def refine_waveform(lattice, iterations=10):
    global entropy_values, wave_iterations

    for i in range(iterations):
        # Compress lattice by taking mean and normalizing
        compressed = np.mean(lattice, axis=0)
        compressed /= np.max(np.abs(compressed))

        # Update lattice by outer product with the compressed state
        lattice = np.outer(compressed, compressed)

        # Calculate entropy as variance of lattice values (proxy for disorder)
        entropy = np.var(lattice)
        entropy_values.append(entropy)

        # Store wave for visualization
        wave_iterations.append(lattice.copy())

    return lattice


# Perform recursive refinement
refined_lattice = refine_waveform(lattice, iterations=10)

# Plot the results
fig, axs = plt.subplots(3, 4, figsize=(18, 12))
axs = axs.flatten()

# Plot each wave iteration
for i, wave in enumerate(wave_iterations):
    ax = axs[i]
    ax.imshow(wave, cmap="plasma", origin="upper")
    ax.set_title(f"Wave Iteration {i+1}")

# Hide any unused plots in the grid
for j in range(len(wave_iterations), len(axs)):
    axs[j].axis("off")

plt.tight_layout()
plt.show()

# Plot entropy reduction over iterations
plt.figure(figsize=(10, 5))
plt.plot(entropy_values, marker="o", label="Entropy")
plt.title("Entropy Over Iterations")
plt.xlabel("Iteration")
plt.ylabel("Entropy")
plt.legend()
plt.grid()
plt.show()

# Display the refined lattice
refined_lattice
