import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

# Original hash and expanded data from Byte1 (example values from the visualization)
original_hash = [3, 155, 243, 11, 227, 123, 211, 235, 195, 91, 179, 203, 163, 59, 147, 171, 131, 27, 115, 139]
expanded_data = [3, 155, 243, 11, 227, 123, 211, 235, 195, 91, 179, 203, 163, 59, 147, 171, 131, 27, 115, 139, 147, 131, 227, 59, 195, 91, 179, 235, 115, 147, 59, 27, 243, 211, 195, 123, 171, 91, 139, 59, 131, 195, 27, 203, 115, 11, 139, 179, 59, 195, 131, 171, 115, 203, 27, 91, 59, 11, 195, 179, 139, 59, 155, 243]

# Truncate expanded data to match original hash length for comparison
expanded_data_truncated = expanded_data[:len(original_hash)]

# Byte-by-byte differences
differences = np.array(original_hash) - np.array(expanded_data_truncated)

# Compute statistical measures
mean_abs_diff = np.mean(np.abs(differences))
correlation_coeff = np.corrcoef(original_hash, expanded_data_truncated)[0, 1]

# FFT comparison
fft_original = np.abs(fft(original_hash))
fft_expanded = np.abs(fft(expanded_data_truncated))

# Visualization
fig, axes = plt.subplots(3, 1, figsize=(12, 10))
fig.suptitle("Hash vs. Expanded Byte1 Analysis")

# Byte-by-byte differences
axes[0].bar(range(len(differences)), differences, color="blue", alpha=0.7)
axes[0].set_title("Byte-by-Byte Differences (Original - Expanded)")
axes[0].set_xlabel("Byte Index")
axes[0].set_ylabel("Difference")

# FFT spectra comparison
axes[1].plot(fft_original, label="FFT of Original Hash", color="green", alpha=0.7)
axes[1].plot(fft_expanded, label="FFT of Expanded Data", color="red", alpha=0.7, linestyle="--")
axes[1].set_title("FFT Spectra Comparison")
axes[1].set_xlabel("Frequency Index")
axes[1].set_ylabel("Magnitude")
axes[1].legend()

# Scatter plot of original vs expanded data
axes[2].scatter(original_hash, expanded_data_truncated, color="purple", alpha=0.7)
axes[2].set_title("Scatter Plot: Original vs. Expanded Data")
axes[2].set_xlabel("Original Hash Values")
axes[2].set_ylabel("Expanded Data Values")

plt.tight_layout()
plt.show()

# Output statistical measures
mean_abs_diff, correlation_coeff
