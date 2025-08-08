# Decompose the final lattice states for emergent patterns and symmetry analysis

# Analyze the final base and perturbed waveforms
final_base_waveform = base_waveforms[-1]
final_perturbed_waveform = perturbed_waveforms[-1]

# Perform symmetry analysis
def analyze_symmetry(matrix):
    horizontal_symmetry = np.allclose(matrix, matrix[:, ::-1])
    vertical_symmetry = np.allclose(matrix, matrix[::-1, :])
    return horizontal_symmetry, vertical_symmetry

base_symmetry = analyze_symmetry(final_base_waveform)
perturbed_symmetry = analyze_symmetry(final_perturbed_waveform)

# Plot the final lattice states for visual comparison
fig, ax = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': '3d'})
X, Y = np.meshgrid(range(final_base_waveform.shape[0]), range(final_base_waveform.shape[1]))

# Base waveform
ax[0].plot_surface(X, Y, final_base_waveform, cmap="plasma", edgecolor="none", alpha=0.8)
ax[0].set_title("Final Base Waveform")
ax[0].set_zlim(-1, 1)

# Perturbed waveform
ax[1].plot_surface(X, Y, final_perturbed_waveform, cmap="viridis", edgecolor="none", alpha=0.8)
ax[1].set_title("Final Perturbed Waveform")
ax[1].set_zlim(-1, 1)

plt.tight_layout()
plt.show()

# Output symmetry results for final waveforms
base_symmetry, perturbed_symmetry
