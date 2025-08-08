import numpy as np
import matplotlib.pyplot as plt

# Constants
EXPANSION_FACTOR = 1.5
PHASE_SCALE = 0.1  # To ensure phase stays within manageable bounds

# Step 1: Generate Samson's Wave
def quantum_wave_samson(hash_value):
    """
    Generates a quantum waveform from a hash using amplitude and phase encoding.
    """
    binary_data = ''.join(format(int(char, 16), '04b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    amplitude = np.sin(2 * theta) + np.cos(theta / 2)
    phase = np.angle(np.exp(1j * theta))  # Phase of the wave
    return amplitude, phase

# Step 2: Store Waveform in H
def store_waveform_in_H(amplitude, phase, expansion_factor=EXPANSION_FACTOR):
    """
    Encodes the amplitude and phase directly into H as harmonic components.
    """
    H_amplitude = np.cumsum(amplitude * expansion_factor)
    H_phase = np.cumsum(phase * PHASE_SCALE)
    return H_amplitude, H_phase

# Step 3: Retrieve Waveform from H
def retrieve_waveform_from_H(H_amplitude, H_phase, expansion_factor=EXPANSION_FACTOR):
    """
    Decodes the harmonic components back into amplitude and phase.
    """
    amplitude = np.diff(H_amplitude, prepend=0) / expansion_factor
    phase = np.diff(H_phase, prepend=0) / PHASE_SCALE
    return amplitude, phase

# Input: Example Hash Value
hash_value = "cd2eca3535741f27a8ae40c31b0c41d4057a7a7b912b33b9aed86485d1c84676"

# Generate Waveform
amplitude, phase = quantum_wave_samson(hash_value)

# Encode into H
H_amplitude, H_phase = store_waveform_in_H(amplitude, phase)

# Retrieve Waveform
retrieved_amplitude, retrieved_phase = retrieve_waveform_from_H(H_amplitude, H_phase)

# Visualize Results
plt.figure(figsize=(12, 8))

# Original Amplitude
plt.subplot(311)
plt.plot(amplitude, label="Original Amplitude", color='blue')
plt.title("Original Quantum Amplitude")
plt.legend()

# Retrieved Amplitude
plt.subplot(312)
plt.plot(retrieved_amplitude, label="Retrieved Amplitude", color='green', linestyle='--')
plt.title("Retrieved Quantum Amplitude")
plt.legend()

# Phase Comparison
plt.subplot(313)
plt.plot(phase, label="Original Phase", color='blue')
plt.plot(retrieved_phase, label="Retrieved Phase", color='green', linestyle='--')
plt.title("Phase Comparison")
plt.legend()

plt.tight_layout()
plt.show()

# Output Results
print("Original Amplitude (First 10):", amplitude[:10])
print("Retrieved Amplitude (First 10):", retrieved_amplitude[:10])
print("Original Phase (First 10):", phase[:10])
print("Retrieved Phase (First 10):", retrieved_phase[:10])
