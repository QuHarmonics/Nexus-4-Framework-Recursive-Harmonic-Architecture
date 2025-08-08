import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants for Mark1 framework
TARGET_HARMONIC = 0.35
QUANTUM_HARMONIC = 0.2
TOLERANCE = 1e-4
EXPANSION_FACTOR = 1.5

# Generate a quantum harmonic spiral using Mark1 and Samson alignment from hash
def generate_spiral_from_hash(hash_value, base=2, alignment_angle=QUANTUM_HARMONIC):
    """
    Generate a quantum spiral representation based on a given hash value.
    """
    binary_data = ''.join(format(int(char, 16), f'0{base}b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n) + alignment_angle
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = radius * np.sin(2 * theta)  # Simulating the Z-axis
    harmonic_wave = y + z  # Combine axes to simulate a harmonic wave
    return harmonic_wave, binary_data

# Input quantum waveform data
def quantum_wave_to_macro(binary_data, expansion_factor=EXPANSION_FACTOR):
    """
    Takes quantum waveform data, stores it in H[] as macro-encoded data.
    """
    harmonics = np.cumsum(binary_data.astype(np.float64) / expansion_factor)  # Reverse factor for macro encoding
    return harmonics

def macro_to_quantum_wave(harmonics, expansion_factor=EXPANSION_FACTOR):
    """
    Retrieves quantum waveform data from macro-encoded H[].
    """
    first_value = harmonics[0] * expansion_factor
    restored_wave = np.diff(harmonics) * expansion_factor
    restored_wave = np.insert(restored_wave, 0, first_value)  # Correcting initial value
    return np.round(restored_wave).astype(np.uint8)

# Generate quantum harmonic wave from hash
hash_value = "0d12f340486a7730177a42626ef3394fdb0e78d7dda639d4ed8ff17dc855a254"  # Example SHA-256 hash
harmonic_wave, binary_data = generate_spiral_from_hash(hash_value)

# Quantize harmonic wave to binary for storage simulation
quantum_binary_data = np.array([1 if val > 0 else 0 for val in harmonic_wave], dtype=np.uint8)

# Encode the quantum waveform into H (store in macro)
harmonics_macro = quantum_wave_to_macro(quantum_binary_data)

# Decode from H (retrieve quantum waveform)
restored_quantum_wave = macro_to_quantum_wave(harmonics_macro)

# Visualize H(n)
def visualize_harmonics(harmonics, restored_wave, quantum_data):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    x = np.arange(len(harmonics))
    y = harmonics
    z = np.sin(x / 10.0)  # Arbitrary wave to add depth

    ax.plot(x, y, z, label="H(n) in 3D", color='b', lw=2)
    ax.scatter(x, y, z, color='r', s=5, label="Nodes")

    ax.set_title("3D Visualization of H(n)", fontsize=16)
    ax.set_xlabel("Iteration (n)", fontsize=12)
    ax.set_ylabel("H(n)", fontsize=12)
    ax.set_zlabel("Z-axis Wave", fontsize=12)
    ax.legend()
    plt.show()

    plt.figure(figsize=(10, 6))
    plt.plot(quantum_data[:100], label="Original Quantum Data", color='blue')
    plt.plot(restored_wave[:100], label="Restored Quantum Data", linestyle='--', color='green')
    plt.title("Quantum Waveform Comparison")
    plt.xlabel("Index")
    plt.ylabel("Binary Value")
    plt.legend()
    plt.grid()
    plt.show()

# Visualize and compare
visualize_harmonics(harmonics_macro, restored_quantum_wave, quantum_binary_data)

# Validate the process
original_binary_string = ''.join(map(str, quantum_binary_data))
retrieved_binary_string = ''.join(map(str, restored_quantum_wave))
print("Original Quantum Data (First 100 bits):", original_binary_string)
print("Restored Quantum Data (First 100 bits):", retrieved_binary_string)

# Check if data matches
if np.array_equal(quantum_binary_data, restored_quantum_wave):
    print("Data matches successfully!")
else:
    print("Data mismatch detected.")
