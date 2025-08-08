import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft

# Constants
EXPANSION_FACTOR = 1.5
TOLERANCE = 1e-4

# Step 1: Generate Quantum Wave
def generate_quantum_wave(hash_value):
    """
    Generate a quantum wave from hash data.
    """
    binary_data = ''.join(format(int(char, 16), '04b') for char in hash_value)
    wave = np.array([int(b) for b in binary_data], dtype=np.float64)
    wave = np.sin(wave * np.pi)  # Create a quantum-like sine wave
    return wave, binary_data

# Step 2: Apply Fourier Transform to Extract Harmonics
def extract_harmonics_from_wave(quantum_wave):
    """
    Decompose quantum wave into harmonics using Fourier Transform.
    """
    harmonics = np.abs(fft(quantum_wave))  # Magnitudes of frequency components
    return harmonics

# Step 3: Store Harmonics in H-Array
def store_in_H(harmonics):
    """
    Map the frequency components into H array.
    """
    cumulative_harmonics = np.cumsum(harmonics) * EXPANSION_FACTOR
    return cumulative_harmonics

# Step 4: Retrieve Binary Data from H
def retrieve_from_H(harmonics):
    """
    Retrieve binary data from the H array.
    """
    first_value = harmonics[0] / EXPANSION_FACTOR
    reversed_data = np.diff(harmonics) / EXPANSION_FACTOR
    reversed_data = np.insert(reversed_data, 0, first_value)
    return np.round(reversed_data).astype(np.uint8)

# Step 5: Convert H Back to Quantum Wave
def reconstruct_wave_from_H(harmonics):
    """
    Use the inverse FFT to reconstruct the original quantum wave from H array.
    """
    reconstructed_wave = ifft(harmonics).real  # Real part of inverse FFT
    return reconstructed_wave

# Hash input
hash_value = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Generate the quantum wave
quantum_wave, binary_data = generate_quantum_wave(hash_value)

# Extract harmonics
harmonics = extract_harmonics_from_wave(quantum_wave)

# Store in H
H_array = store_in_H(harmonics)

# Retrieve binary from H
retrieved_binary = retrieve_from_H(H_array)

# Reconstruct the wave
reconstructed_wave = reconstruct_wave_from_H(harmonics)

# Outputs
print("Original Binary Data:", binary_data[:100])
print("Retrieved Binary Data:", ''.join(map(str, retrieved_binary[:100])))
print("H array values:", H_array[:10])

# Visualization
plt.figure(figsize=(12, 6))
plt.plot(H_array, label="H(n) Harmonics", color="blue")
plt.title("Harmonic Storage Representation")
plt.xlabel("Index")
plt.ylabel("H(n)")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(12, 6))
plt.plot(quantum_wave, label="Original Quantum Wave", color="green")
plt.plot(reconstructed_wave, label="Reconstructed Quantum Wave", linestyle='dashed', color="red")
plt.title("Quantum Wave Comparison")
plt.legend()
plt.grid()
plt.show()
