import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.fft import fft

# Constants for Testing and H Framework
EXPANSION_FACTOR = 1.5

# Step 1: Simulate a Quantum Waveform from the Hash (SHA-like mechanism)
def generate_quantum_wave(hash_value):
    """
    Generate a quantum-like waveform from the input hash value.
    """
    binary_data = ''.join(format(int(char, 16), '04b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta) + np.cos(theta / 2)
    wave = y + z
    return x, y, z, wave, binary_data

# Step 2: Store Quantum Data in H
def store_in_H(binary_data, expansion_factor=EXPANSION_FACTOR):
    harmonics = np.cumsum(binary_data.astype(np.float64) * expansion_factor)
    return harmonics

# Step 3: Convert Macro Hash into a Quantum Waveform
def macro_to_quantum_wave(macro_binary_data):
    """
    Transform a macro binary hash into its quantum representation.
    """
    n = len(macro_binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta) + np.cos(theta / 3)
    wave = y + z
    return x, y, z, wave

# Step 4: Visualize Combined Quantum and Macro Waves
def visualize_combined(x1, y1, z1, wave1, x2, y2, z2, wave2, title="Combined Visualization"):
    fig = plt.figure(figsize=(14, 8))
    
    # Quantum Wave
    ax1 = fig.add_subplot(121, projection='3d')
    ax1.plot(x1, y1, z1, label="Quantum Wave", color="blue")
    ax1.scatter(x1, y1, z1, color="red", s=5)
    ax1.set_title("Quantum Waveform")
    ax1.legend()
    
    # Macro Wave
    ax2 = fig.add_subplot(122, projection='3d')
    ax2.plot(x2, y2, z2, label="Macro Wave", color="green")
    ax2.scatter(x2, y2, z2, color="orange", s=5)
    ax2.set_title("Macro Waveform")
    ax2.legend()
    
    plt.suptitle(title, fontsize=16)
    plt.show()

# Hash Value and Binary Setup
hash_value = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Generate Quantum Waveform from Hash
x1, y1, z1, quantum_wave, binary_data = generate_quantum_wave(hash_value)

# Convert Macro Hash into Quantum Waveform
binary_data_array = np.array([int(b) for b in binary_data])
harmonics = store_in_H(binary_data_array)
x2, y2, z2, macro_wave = macro_to_quantum_wave(binary_data_array)

# Visualize the Combined Result
visualize_combined(x1, y1, z1, quantum_wave, x2, y2, z2, macro_wave, title="Quantum and Macro Wave Comparison")

# Print Harmonic Array and Binary Output for Validation
print("H array values:", harmonics[:100])  # Show only the first 100 for brevity
print("Quantum Wave Binary (First 100 bits):", binary_data[:100])
