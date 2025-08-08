import numpy as np
import matplotlib.pyplot as plt

# Initialize larger lattice (512x512 grid)
grid_size = 512

# Create two random lattices to simulate hash differences
lattice1 = np.random.rand(grid_size, grid_size)
lattice2 = np.random.rand(grid_size, grid_size)

# Fourier transforms to get frequency space representations
fft1 = np.fft.fftshift(np.fft.fft2(lattice1))
fft2 = np.fft.fftshift(np.fft.fft2(lattice2))

# Amplitude and phase in frequency space
amplitude1 = np.abs(fft1)
phase1 = np.angle(fft1)

amplitude2 = np.abs(fft2)
phase2 = np.angle(fft2)

# Compute amplitude and phase differences
amplitude_diff = np.abs(amplitude1 - amplitude2)
phase_diff = np.abs(phase1 - phase2)

# Backpropagation (inverse FFT) to analyze differences in the spatial domain
backprop_amplitude_diff = np.fft.ifft2(np.fft.ifftshift(amplitude_diff))
backprop_phase_diff = np.fft.ifft2(np.fft.ifftshift(phase_diff))

# Plot results for analysis
fig, axes = plt.subplots(3, 3, figsize=(12, 12))

# Plot amplitude grids
axes[0, 0].imshow(amplitude1, cmap='viridis', extent=(-grid_size//2, grid_size//2, -grid_size//2, grid_size//2))
axes[0, 0].set_title('Lattice 1 Amplitude (Frequency Space)')
axes[0, 1].imshow(amplitude2, cmap='viridis', extent=(-grid_size//2, grid_size//2, -grid_size//2, grid_size//2))
axes[0, 1].set_title('Lattice 2 Amplitude (Frequency Space)')
axes[0, 2].imshow(amplitude_diff, cmap='viridis', extent=(-grid_size//2, grid_size//2, -grid_size//2, grid_size//2))
axes[0, 2].set_title('Amplitude Difference (Frequency Space)')

# Plot phase grids
axes[1, 0].imshow(phase1, cmap='twilight', extent=(-grid_size//2, grid_size//2, -grid_size//2, grid_size//2))
axes[1, 0].set_title('Lattice 1 Phase (Frequency Space)')
axes[1, 1].imshow(phase2, cmap='twilight', extent=(-grid_size//2, grid_size//2, -grid_size//2, grid_size//2))
axes[1, 1].set_title('Lattice 2 Phase (Frequency Space)')
axes[1, 2].imshow(phase_diff, cmap='twilight', extent=(-grid_size//2, grid_size//2, -grid_size//2, grid_size//2))
axes[1, 2].set_title('Phase Difference (Frequency Space)')

# Plot backpropagated differences
axes[2, 0].imshow(np.abs(backprop_amplitude_diff), cmap='viridis', extent=(0, grid_size, 0, grid_size))
axes[2, 0].set_title('Backpropagated Amplitude Difference')
axes[2, 1].imshow(np.abs(backprop_phase_diff), cmap='viridis', extent=(0, grid_size, 0, grid_size))
axes[2, 1].set_title('Backpropagated Phase Difference')

# Leave the last subplot as a placeholder
axes[2, 2].axis('off')
axes[2, 2].set_title('Placeholder for Future Analysis')

plt.tight_layout()
plt.show()
