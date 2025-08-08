import numpy as np
import matplotlib.pyplot as plt

# Step 1: Generate a larger lattice (512x512)
lattice1 = np.random.rand(512, 512)
lattice2 = np.random.rand(512, 512)

# Step 2: Compute Fourier Transforms
fft_lattice1 = np.fft.fftshift(np.fft.fft2(lattice1))
fft_lattice2 = np.fft.fftshift(np.fft.fft2(lattice2))

# Step 3: Compute Amplitude and Phase
amplitude1 = np.abs(fft_lattice1)
amplitude2 = np.abs(fft_lattice2)
phase1 = np.angle(fft_lattice1)
phase2 = np.angle(fft_lattice2)

# Compute Differences
amplitude_diff = np.abs(amplitude1 - amplitude2)
phase_diff = np.abs(phase1 - phase2)

# Step 4: Backpropagate Differences
backprop_amplitude_diff = np.fft.ifft2(np.fft.ifftshift(amplitude_diff))
backprop_phase_diff = np.fft.ifft2(np.fft.ifftshift(phase_diff))

# Plot results
fig, axs = plt.subplots(3, 3, figsize=(15, 15))

# Top Row: Amplitude
axs[0, 0].imshow(amplitude1, cmap='viridis')
axs[0, 0].set_title("Lattice 1 Amplitude (Frequency Space)")
axs[0, 1].imshow(amplitude2, cmap='viridis')
axs[0, 1].set_title("Lattice 2 Amplitude (Frequency Space)")
axs[0, 2].imshow(amplitude_diff, cmap='viridis')
axs[0, 2].set_title("Amplitude Difference (Frequency Space)")

# Middle Row: Phase
axs[1, 0].imshow(phase1, cmap='twilight_shifted')
axs[1, 0].set_title("Lattice 1 Phase (Frequency Space)")
axs[1, 1].imshow(phase2, cmap='twilight_shifted')
axs[1, 1].set_title("Lattice 2 Phase (Frequency Space)")
axs[1, 2].imshow(phase_diff, cmap='twilight_shifted')
axs[1, 2].set_title("Phase Difference (Frequency Space)")

# Bottom Row: Backpropagated Differences
axs[2, 0].imshow(np.abs(backprop_amplitude_diff), cmap='plasma')
axs[2, 0].set_title("Backpropagated Amplitude Difference")
axs[2, 1].imshow(np.abs(backprop_phase_diff), cmap='plasma')
axs[2, 1].set_title("Backpropagated Phase Difference")
axs[2, 2].axis('off')  # Leave one blank

plt.tight_layout()
plt.show()
