from scipy.fftpack import fft, ifft

# Harmonic decomposition function using FFT
def harmonic_decomposition(matrix):
    fft_matrix = fft(matrix, axis=0)
    magnitude = np.abs(fft_matrix)
    phase = np.angle(fft_matrix)
    return magnitude, phase

# Decompose the final base and perturbed waveforms
base_magnitude, base_phase = harmonic_decomposition(final_base_waveform)
perturbed_magnitude, perturbed_phase = harmonic_decomposition(final_perturbed_waveform)

# Plot magnitude and phase for the base waveform
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
axes[0].imshow(base_magnitude, cmap="plasma", aspect="auto")
axes[0].set_title("Base Magnitude Spectrum")
axes[0].set_xlabel("Frequency (Columns)")
axes[0].set_ylabel("Frequency (Rows)")

axes[1].imshow(base_phase, cmap="viridis", aspect="auto")
axes[1].set_title("Base Phase Spectrum")
axes[1].set_xlabel("Frequency (Columns)")
axes[1].set_ylabel("Frequency (Rows)")

plt.tight_layout()
plt.show()

# Plot magnitude and phase for the perturbed waveform
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
axes[0].imshow(perturbed_magnitude, cmap="plasma", aspect="auto")
axes[0].set_title("Perturbed Magnitude Spectrum")
axes[0].set_xlabel("Frequency (Columns)")
axes[0].set_ylabel("Frequency (Rows)")

axes[1].imshow(perturbed_phase, cmap="viridis", aspect="auto")
axes[1].set_title("Perturbed Phase Spectrum")
axes[1].set_xlabel("Frequency (Columns)")
axes[1].set_ylabel("Frequency (Rows)")

plt.tight_layout()
plt.show()

# Feedback loop: Reconstruct the waveforms from selected harmonic components
reconstructed_base = ifft(base_magnitude * np.exp(1j * base_phase), axis=0).real
reconstructed_perturbed = ifft(perturbed_magnitude * np.exp(1j * perturbed_phase), axis=0).real

# Visualize the reconstructed waveforms
fig, ax = plt.subplots(1, 2, figsize=(12, 6), subplot_kw={'projection': '3d'})
X, Y = np.meshgrid(range(reconstructed_base.shape[0]), range(reconstructed_base.shape[1]))

# Reconstructed Base
ax[0].plot_surface(X, Y, reconstructed_base, cmap="plasma", edgecolor="none", alpha=0.8)
ax[0].set_title("Reconstructed Base Waveform")
ax[0].set_zlim(-1, 1)

# Reconstructed Perturbed
ax[1].plot_surface(X, Y, reconstructed_perturbed, cmap="viridis", edgecolor="none", alpha=0.8)
ax[1].set_title("Reconstructed Perturbed Waveform")
ax[1].set_zlim(-1, 1)

plt.tight_layout()
plt.show()
