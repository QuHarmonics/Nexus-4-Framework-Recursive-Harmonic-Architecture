import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft2, ifft2, fftshift

# Generate larger grid for expanded lattice
grid_size = 512
x = np.linspace(0, 1, grid_size)
y = np.linspace(0, 1, grid_size)
X, Y = np.meshgrid(x, y)

# Define example data for two SHA256 hashes
# These are placeholders to represent differences in lattice due to hash inputs
hash_1_waveform = np.sin(4 * np.pi * X) * np.cos(4 * np.pi * Y)
hash_2_waveform = np.sin(6 * np.pi * X) * np.cos(6 * np.pi * Y)

# Compute the difference lattice
difference_waveform = hash_1_waveform - hash_2_waveform

# Apply 2D Fourier Transforms
fft_hash_1 = fft2(hash_1_waveform)
fft_hash_2 = fft2(hash_2_waveform)
fft_diff = fft2(difference_waveform)  # <<<< (You missed this line!)

# Shift zero frequency to the center
fft_hash_1_shifted = fftshift(fft_hash_1)
fft_hash_2_shifted = fftshift(fft_hash_2)
fft_diff_shifted = fftshift(fft_diff)

# --- Visualization of FFT Magnitudes ---

plt.figure(figsize=(18, 6))

# FFT of first hash
plt.subplot(1, 3, 1)
plt.imshow(np.log(np.abs(fft_hash_1_shifted) + 1), extent=[-0.5, 0.5, -0.5, 0.5], cmap='viridis')
plt.title('FFT Magnitude - Hash 1')
plt.colorbar()

# FFT of second hash
plt.subplot(1, 3, 2)
plt.imshow(np.log(np.abs(fft_hash_2_shifted) + 1), extent=[-0.5, 0.5, -0.5, 0.5], cmap='plasma')
plt.title('FFT Magnitude - Hash 2')
plt.colorbar()

# FFT of the difference
plt.subplot(1, 3, 3)
plt.imshow(np.log(np.abs(fft_diff_shifted) + 1), extent=[-0.5, 0.5, -0.5, 0.5], cmap='magma')
plt.title('FFT Magnitude - Difference')
plt.colorbar()

plt.tight_layout()
plt.show()
