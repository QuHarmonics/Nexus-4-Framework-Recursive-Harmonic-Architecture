import numpy as np
import matplotlib.pyplot as plt

# Constants
EXPANSION_FACTOR = 1.5

# Step 1: Project Wave Reflection
def reflect_wave(wave):
    """
    Reflects the quantum wave to project its amplitude into a 2D harmonic representation.
    """
    reflected_wave = np.abs(wave)  # Take the magnitude as the reflection
    return reflected_wave

# Step 2: Store in H (Wave to Harmonics)
def store_in_H(wave, expansion_factor=EXPANSION_FACTOR):
    """
    Stores the reflected wave (harmonics) into H as cumulative harmonics.
    """
    harmonics = np.cumsum(wave.astype(np.float64) * expansion_factor)
    return harmonics

# Step 3: Retrieve from H (Harmonics to Macro Binary)
def retrieve_from_H(harmonics, expansion_factor=EXPANSION_FACTOR):
    """
    Retrieves the wave back from H harmonics without collapsing.
    """
    first_value = harmonics[0] / expansion_factor
    reversed_data = np.diff(harmonics) / expansion_factor
    reversed_data = np.insert(reversed_data, 0, first_value)
    return reversed_data

# Step 4: Generate Quantum Wave (Samson)
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
    z = np.sin(2 * theta) + np.cos(theta / 2)
    return np.sqrt(x**2 + y**2 + z**2)  # Magnitude projection

# Input: Quantum Wave (Samson)
hash_value = "cd2eca3535741f27a8ae40c31b0c41d4057a7a7b912b33b9aed86485d1c84676"
found_wave = quantum_wave_samson(hash_value)

# Reflect and Store in H
reflected_wave = reflect_wave(found_wave)
harmonics = store_in_H(reflected_wave)

# Retrieve from H
retrieved_wave = retrieve_from_H(harmonics)

# Visualize Results
plt.figure(figsize=(14, 10))

# Original Wave
plt.subplot(311)
plt.plot(found_wave, label="Original Quantum Wave (Samson)", color='blue')
plt.title("Original Quantum Wave (Samson)")
plt.legend()

# Reflected Harmonics
plt.subplot(312)
plt.plot(harmonics, label="Reflected Harmonics (Wave -> H)", color='orange')
plt.title("Reflected Harmonics Representation")
plt.legend()

# Retrieved Wave
plt.subplot(313)
plt.plot(retrieved_wave, label="Retrieved Wave (H -> Macro)", color='green')
plt.title("Retrieved Macro Wave")
plt.legend()

plt.tight_layout()
plt.show()

# Output Results
print("Found Wave (First 10):", found_wave[:10])
print("Reflected Harmonics (First 10):", harmonics[:10])
print("Retrieved Wave (First 10):", retrieved_wave[:10])
