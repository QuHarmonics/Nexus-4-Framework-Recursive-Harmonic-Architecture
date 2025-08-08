import numpy as np
import matplotlib.pyplot as plt

# Reinitialize parameters for the grid and FFT processing
grid_size = 500
iterations = 20

# Define two random input lattices to simulate the hashes
np.random.seed(42)  # For reproducibility
lattice1 = np.random.rand(grid_size, grid_size)
lattice2 = np.random.rand(grid_size, grid_size)

# Compute FFTs for both lattices
fft_lattice1 = np.fft.fftshift(np.fft.fft2(lattice1))
fft_lattice2 = np.fft.fftshift(np.fft.fft2(lattice2))

# Amplitude and phase components
amplitude1 = np.abs(fft_lattice1)
amplitude2 = np.abs(fft_lattice2)
phase1 = np.angle(fft_lattice1)
phase2 = np.angle(fft_lattice2)

# Calculate differences
amplitude_diff = np.abs(amplitude1 - amplitude2)
phase_diff = np.abs(phase1 - phase2)

# Backpropagation placeholder (needs iterative refinement in real process)
backprop_amplitude_diff = np.zeros_like(amplitude_diff)
backprop_phase_diff = np.zeros_like(phase_diff)

# Visualize the results
fig, axs = plt.subplots(3, 3, figsize=(15, 15))
axs = axs.ravel()

# Plot amplitude and phase differences
axs[0].imshow(amplitude1, cmap='viridis')
axs[0].set_title('Lattice 1 Amplitude (Frequency Space)')
axs[1].imshow(amplitude2, cmap='viridis')
axs[1].set_title('Lattice 2 Amplitude (Frequency Space)')
axs[2].imshow(amplitude_diff, cmap='viridis')
axs[2].set_title('Amplitude Difference (Frequency Space)')

axs[3].imshow(phase1, cmap='twilight')
axs[3].set_title('Lattice 1 Phase (Frequency Space)')
axs[4].imshow(phase2, cmap='twilight')
axs[4].set_title('Lattice 2 Phase (Frequency Space)')
axs[5].imshow(phase_diff, cmap='twilight')
axs[5].set_title('Phase Difference (Frequency Space)')

# Backpropagation plots
axs[6].imshow(backprop_amplitude_diff, cmap='viridis')
axs[6].set_title('Backpropagated Amplitude Difference')
axs[7].imshow(backprop_phase_diff, cmap='twilight')
axs[7].set_title('Backpropagated Phase Difference')

# Empty placeholder for comparison
axs[8].text(0.5, 0.5, "Comparison Placeholder", ha='center', va='center', fontsize=12)
axs[8].set_title('Placeholder')
axs[8].axis('off')

plt.tight_layout()
plt.show()
