import numpy as np
import matplotlib.pyplot as plt

# Simulate anti-hash by creating mirrored lattice feedback and tracking divergence

# Initial settings
iterations = 20
lattice_size = 8
original_lattice = np.random.rand(lattice_size, lattice_size)  # Simulated original lattice
hash_constants = np.random.rand(lattice_size)  # Simulated hash constants

# Storage for feedback process
hash_lattice = []
anti_hash_lattice = []
divergences = []

# Forward hash lattice generation
current_lattice = original_lattice.copy()
for _ in range(iterations):
    current_lattice = np.abs(np.sin(current_lattice + 0.1))
    hash_lattice.append(current_lattice)

# Generate mirrored lattice for anti-hash
for step in range(iterations):
    mirror_lattice = -1 * np.flip(hash_lattice[step], axis=0)  # Mirrored along one axis
    anti_hash_lattice.append(mirror_lattice)
    divergence = np.linalg.norm(original_lattice - mirror_lattice)  # Track divergence
    divergences.append(divergence)

# Plot divergence over iterations
plt.figure(figsize=(10, 6))
plt.plot(range(1, iterations + 1), divergences, marker='o', label='Divergence (Anti-Hash)')
plt.title("Divergence of Mirrored Lattice from Original")
plt.xlabel("Iteration")
plt.ylabel("Divergence")
plt.legend()
plt.grid(True)
plt.show()

# Visualizing the forward and mirrored lattices at final iteration
fig, ax = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={"projection": "3d"})
X, Y = np.meshgrid(range(lattice_size), range(lattice_size))

# Original lattice visualization
Z_original = hash_lattice[-1]
ax[0].plot_surface(X, Y, Z_original, cmap="viridis", edgecolor="none")
ax[0].set_title("Final Forward Hash Lattice")

# Mirrored lattice visualization
Z_mirrored = anti_hash_lattice[-1]
ax[1].plot_surface(X, Y, Z_mirrored, cmap="plasma", edgecolor="none")
ax[1].set_title("Final Mirrored (Anti-Hash) Lattice")

plt.tight_layout()
plt.show()

# Display constants and final divergences for numeric comparison
hash_constants_normalized = hash_constants / np.max(hash_constants)
print("Hash Constants (Normalized):", hash_constants_normalized)
print("Final Divergence:", divergences[-1])
