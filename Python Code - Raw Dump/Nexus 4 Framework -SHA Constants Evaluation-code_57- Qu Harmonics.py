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

# Apply Fourier Transform
fft_hash_1 = fft2(hash_1_waveform)
fft_hash_2 = fft2(hash_2_waveform)
fft_diff
