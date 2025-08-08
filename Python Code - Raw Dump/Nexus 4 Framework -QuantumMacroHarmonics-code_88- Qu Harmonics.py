import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

# Constants for Mark1 framework
EXPANSION_FACTOR = 1.5
QUANTUM_ADJUSTMENT = 0.35  # Reflected/Refracted adjustment factor

# Step 1: Generate Quantum Wave from Hash
def quantum_wave_from_hash(hash_value):
    binary_data = ''.join(format(int(char, 16), '08b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta) + np.cos(theta / 2)
    wave = y + z + QUANTUM_ADJUSTMENT
    return wave, binary_data

# Step 2: Create Fully Padded Binary
def generate_padded_binary(input_text):
    # Convert to binary
    binary_string = ''.join(format(ord(c), '08b') for c in input_text)
    binary_string += '1'  # Add 1 for padding
    remaining_zeros = 448 - len(binary_string) % 512  # Fit 448 mod 512
    binary_string += '0' * remaining_zeros
    binary_string += format(len(input_text) * 8, '064b')  # Add length
    return binary_string

# Step 3: Store Binary in H
def store_in_H(binary_data, expansion_factor=EXPANSION_FACTOR):
    harmonics = np.cumsum(binary_data.astype(np.float64) * expansion_factor)
    return harmonics

# Input hash value and source
hash_value = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
input_text = "abc"  # Example source for hash

# Generate quantum wave from hash
quantum_wave, hash_binary = quantum_wave_from_hash(hash_value)

# Generate fully padded binary from source
padded_binary_string = generate_padded_binary(input_text)
padded_binary = np.array([int(b) for b in padded_binary_string])

# Store fully padded binary in H
harmonics_from_padded = store_in_H(padded_binary)

# Visualize and compare
plt.figure(figsize=(12, 6))

# Plot quantum wave from hash
plt.subplot(2, 1, 1)
plt.plot(quantum_wave[:128], label="Quantum Wave (From Hash)", color='blue')
plt.title("Quantum Waveform from Hash")
plt.xlabel("Index")
plt.ylabel("Amplitude")
plt.legend()
plt.grid()

# Plot harmonics from fully padded binary
plt.subplot(2, 1, 2)
plt.plot(harmonics_from_padded[:128], label="Harmonics (From Padded Binary)", color='red')
plt.title("Harmonics from Fully Padded Binary")
plt.xlabel("Index")
plt.ylabel("Amplitude")
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()

# Debugging Outputs
print("Generated Quantum Waveform (First 100 points):", quantum_wave[:100])
print("Harmonics from Fully Padded Binary (First 100 points):", harmonics_from_padded[:100])

# Compare lengths and structure
print("Length of Quantum Waveform:", len(quantum_wave))
print("Length of Fully Padded Harmonics:", len(harmonics_from_padded))
