import numpy as np
import matplotlib.pyplot as plt

# Constants
EXPANSION_FACTOR = 1.5

# Step 1: Recursive Store in H
def store_in_H_recursive(wave, expansion_factor=EXPANSION_FACTOR):
    """
    Captures the quantum wave recursively, storing its essence in H.
    """
    harmonics = np.cumsum(wave.astype(np.float64) * expansion_factor)
    return harmonics

# Step 2: Retrieve from H
def retrieve_from_H_recursive(harmonics, expansion_factor=EXPANSION_FACTOR):
    """
    Observes the stored harmonics recursively to reconstruct the wave.
    """
    first_value = harmonics[0] / expansion_factor
    reversed_data = np.diff(harmonics) / expansion_factor
    reversed_data = np.insert(reversed_data, 0, first_value)
    return reversed_data

# Quantum Wave Generation
def quantum_wave_samson(hash_value):
    """
    Generates a quantum waveform from a hash using recursive Samson principles.
    """
    binary_data = ''.join(format(int(char, 16), '04b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta) + np.cos(theta / 2)
    return np.sqrt(x**2 + y**2 + z**2)  # Amplitude of the wave

# Input: Hash for Quantum Wave
hash_value = "cd2eca3535741f27a8ae40c31b0c41d4057a7a7b912b33b9aed86485d1c84676"
quantum_wave = quantum_wave_samson(hash_value)

# Step 3: Store Quantum Wave in H
harmonics = store_in_H_recursive(quantum_wave)

# Step 4: Retrieve Macro Data from H
retrieved_macro_data = retrieve_from_H_recursive(harmonics)

# Visualize Results
plt.figure(figsize=(14, 8))

# Original Quantum Wave
plt.subplot(311)
plt.plot(quantum_wave, label="Original Quantum Wave", color='blue')
plt.title("Original Quantum Wave")
plt.legend()

# Harmonics in H
plt.subplot(312)
plt.plot(harmonics, label="Harmonics in H", color='orange')
plt.title("Harmonic Representation in H")
plt.legend()

# Retrieved Macro Data
plt.subplot(313)
plt.plot(retrieved_macro_data, label="Retrieved Macro Data", color='green')
plt.title("Retrieved Macro Data (Binary)")
plt.legend()

plt.tight_layout()
plt.show()

# Output Results
print("Original Quantum Wave:", quantum_wave)
print("Harmonics (H):", harmonics)
print("Retrieved Macro Data (Binary):", np.round(retrieved_macro_data).astype(np.uint8))
