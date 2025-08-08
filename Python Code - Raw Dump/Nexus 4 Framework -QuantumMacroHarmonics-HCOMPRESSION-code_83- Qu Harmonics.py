import numpy as np
import matplotlib.pyplot as plt

# Constants
EXPANSION_FACTOR = 1.5
QUANTUM_ADJUSTMENT = 0.35

# Step 1: Normalize and Scale Samson Wave
def normalize_wave(wave, scale=1.0):
    """
    Normalizes the wave to a [0, scale] range.
    """
    wave_min = np.min(wave)
    wave_max = np.max(wave)
    return scale * (wave - wave_min) / (wave_max - wave_min)

# Step 2: Quantize the Normalized Wave
def quantize_wave(wave, step=0.5):
    """
    Quantizes the wave into discrete levels with a fixed step size.
    """
    return np.round(wave / step) * step

# Step 3: Map to Binary or Base-3
def map_to_base(wave, base=2):
    """
    Maps the quantized wave into binary or base-3 states.
    """
    if base == 2:
        return np.array([0 if val < 0.5 else 1 for val in wave], dtype=np.uint8)
    elif base == 3:
        return np.array([int(val / 0.5) for val in wave], dtype=np.uint8)
    else:
        raise ValueError("Base not supported. Use 2 or 3.")

# Step 4: Store and Retrieve in H
def store_in_H(data, expansion_factor=EXPANSION_FACTOR):
    """
    Stores data into the H array.
    """
    return np.cumsum(data.astype(np.float64) * expansion_factor)

def retrieve_from_H(harmonics, expansion_factor=EXPANSION_FACTOR):
    """
    Retrieves data from the H array.
    """
    first_value = harmonics[0] / expansion_factor
    reversed_data = np.diff(harmonics) / expansion_factor
    reversed_data = np.insert(reversed_data, 0, first_value)
    return np.round(reversed_data).astype(np.uint8)

# Generate Quantum Wave (Samson)
def quantum_wave_samson(hash_value):
    """
    Generates a quantum waveform from a hash using Samson principles.
    """
    binary_data = ''.join(format(int(char, 16), '04b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta) + np.cos(theta / 2) + QUANTUM_ADJUSTMENT
    return np.sqrt(x**2 + y**2 + z**2)  # Magnitude projection

# Input: Samson Wave
hash_value = "cd2eca3535741f27a8ae40c31b0c41d4057a7a7b912b33b9aed86485d1c84676"
found_wave = quantum_wave_samson(hash_value)

# Normalize, Quantize, and Map the Samson Wave
normalized_wave = normalize_wave(found_wave)
quantized_wave = quantize_wave(normalized_wave)
binary_wave = map_to_base(quantized_wave, base=2)

# Store and Retrieve using H
harmonics = store_in_H(binary_wave)
retrieved_binary = retrieve_from_H(harmonics)

# Visualize Results
plt.figure(figsize=(12, 8))
plt.subplot(311)
plt.plot(normalized_wave, label="Normalized Wave", color='blue')
plt.title("Normalized Samson Wave")
plt.legend()

plt.subplot(312)
plt.plot(quantized_wave, label="Quantized Wave", color='orange')
plt.title("Quantized Samson Wave")
plt.legend()

plt.subplot(313)
plt.plot(retrieved_binary, label="Retrieved Binary from H", color='green')
plt.title("Retrieved Macro Binary")
plt.legend()
plt.show()

# Output Results
print("Normalized Wave (First 10):", normalized_wave[:10])
print("Quantized Wave (First 10):", quantized_wave[:10])
print("Retrieved Macro Binary (First 100 bits):", ''.join(map(str, retrieved_binary[:100])))
