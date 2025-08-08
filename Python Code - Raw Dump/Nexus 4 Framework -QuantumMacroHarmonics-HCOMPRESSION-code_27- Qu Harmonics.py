import numpy as np
import matplotlib.pyplot as plt
from mpmath import zetazero

# Generate zeta zeros
def generate_zeta_zeros(num_zeros):
    return [zetazero(n).real for n in range(1, num_zeros + 1)]

# Transform binary data into harmonic space
def binary_to_waveform(binary_data):
    n = len(binary_data)
    x = np.linspace(0, 2 * np.pi, n)
    return np.sin(x) * (binary_data - 0.5)

# Overlay zeta zeros on harmonic space
def filter_with_zeta_zeros(waveform, zeta_zeros, harmonic_length):
    filtered_waveform = waveform.copy()
    for zero in zeta_zeros:
        zero_index = int((zero / max(zeta_zeros)) * harmonic_length)
        if 0 <= zero_index < len(filtered_waveform):
            filtered_waveform[zero_index] = 0  # Cancel predicted noise at zeta zero
    return filtered_waveform

# Convert harmonic space back to binary
def waveform_to_binary(waveform, threshold=0.0):
    return np.array([1 if h > threshold else 0 for h in waveform], dtype=int)

# Input string and hash
input_string = "abc"
hashed_binary = np.array([int(b) for b in hash_to_binary(input_string)])

# Transform the hash into harmonic space
hashed_waveform = binary_to_waveform(hashed_binary)

# Generate zeta zeros
num_zeros = 100  # Number of zeta zeros to use
zeta_zeros = generate_zeta_zeros(num_zeros)

# Filter the hashed waveform using zeta zeros
filtered_waveform = filter_with_zeta_zeros(hashed_waveform, zeta_zeros, len(hashed_waveform))

# Convert the filtered waveform back to binary
reconstructed_binary = waveform_to_binary(filtered_waveform)

# Display results
print("\nOriginal Hashed Binary (First 128 Bits):")
print(hashed_binary[:128])
print("\nReconstructed Binary After Zeta Filtering (First 128 Bits):")
print(reconstructed_binary[:128])

# Plot the comparison
plt.figure(figsize=(12, 6))
plt.subplot(2, 1, 1)
plt.plot(hashed_waveform[:128], label="Original Hashed Waveform", color='orange')
plt.title("Original Hashed Waveform (First 128 Points)")
plt.grid()
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(filtered_waveform[:128], label="Filtered Waveform (Zeta Zeros Applied)", color='green')
plt.title("Filtered Waveform After Zeta Zeros (First 128 Points)")
plt.grid()
plt.legend()

plt.tight_layout()
plt.show()
