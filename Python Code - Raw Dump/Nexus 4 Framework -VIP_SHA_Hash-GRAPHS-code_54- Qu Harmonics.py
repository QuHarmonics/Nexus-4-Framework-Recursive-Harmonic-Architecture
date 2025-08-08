import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft2, fftshift

# Define lattice size and grid
grid_size = 512
lattice1 = np.random.rand(grid_size, grid_size)  # Placeholder for real hash lattice 1
lattice2 = np.random.rand(grid_size, grid_size)  # Placeholder for real hash lattice 2

# Compute difference lattice
difference_lattice = lattice2 - lattice1

# Apply Fourier Transform
fft_lattice1 = fftshift(fft2(lattice1))
fft_lattice2 = fftshift(fft2(lattice2))
fft_difference = fftshift(fft2(difference_lattice))

# Compute amplitude and phase
amplitude1 = np.abs(fft_lattice1)
phase1 = np.angle(fft_lattice1)
amplitude2 = np.abs(fft_lattice2)
phase2 = np.angle(fft_lattice2)
amplitude_diff = np.abs(fft_difference)
phase_diff = np.angle(fft_difference)

# Visualization: Amplitude and Phase in Frequency Space
fig, axs = plt.subplots(2, 3, figsize=(18, 12))

# Lattice 1: Amplitude and Phase
axs[0, 0].imshow(amplitude1, cmap='viridis', extent=[-grid_size//2, grid_size//2, -grid_size//2, grid_size//2])
axs[0, 0].set_title("Lattice 1 Amplitude (Frequency Space)")
axs[0, 0].set_xlabel("Frequency X")
axs[0, 0].set_ylabel("Frequency Y")

axs[1, 0].imshow(phase1, cmap='twilight', extent=[-grid_size//2, grid_size//2, -grid_size//2, grid_size//2])
axs[1, 0].set_title("Lattice 1 Phase (Frequency Space)")
axs[1, 0].set_xlabel("Frequency X")
axs[1, 0].set_ylabel("Frequency Y")

# Lattice 2: Amplitude and Phase
axs[0, 1].imshow(amplitude2, cmap='viridis', extent=[-grid_size//2, grid_size//2, -grid_size//2, grid_size//2])
axs[0, 1].set_title("Lattice 2 Amplitude (Frequency Space)")
axs[0, 1].set_xlabel("Frequency X")
axs[0, 1].set_ylabel("Frequency Y")

axs[1, 1].imshow(phase2, cmap='twilight', extent=[-grid_size//2, grid_size//2, -grid_size//2, grid_size//2])
axs[1, 1].set_title("Lattice 2 Phase (Frequency Space)")
axs[1, 1].set_xlabel("Frequency X")
axs[1, 1].set_ylabel("Frequency Y")

# Difference Lattice: Amplitude and Phase
axs[0, 2].imshow(amplitude_diff, cmap='viridis', extent=[-grid_size//2, grid_size//2, -grid_size//2, grid_size//2])
axs[0, 2].set_title("Difference Lattice Amplitude (Frequency Space)")
axs[0, 2].set_xlabel("Frequency X")
axs[0, 2].set_ylabel("Frequency Y")

axs[1, 2].imshow(phase_diff, cmap='twilight', extent=[-grid_size//2, grid_size//2, -grid_size//2, grid_size//2])
axs[1, 2].set_title("Difference Lattice Phase (Frequency Space)")
axs[1, 2].set_xlabel("Frequency X")
axs[1, 2].set_ylabel("Frequency Y")

plt.tight_layout()
plt.show()

# Return results for further analysis
amplitude1.max(), amplitude2.max(), amplitude_diff.max()
