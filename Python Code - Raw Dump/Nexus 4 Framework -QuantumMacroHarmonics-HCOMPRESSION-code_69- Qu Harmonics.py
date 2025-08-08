import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

# Constants for Mark1 framework
EXPANSION_FACTOR = 1.5
QUANTUM_ADJUSTMENT = 0.35  # Reflected/Refracted adjustment factor

# Step 1: Generate a Quantum Wave from Hash
def quantum_wave_from_hash(hash_value, base=2):
    """
    Generates a quantum wave from the input hash using Mark1 and Samson principles.
    """
    binary_data = ''.join(format(int(char, 16), f'0{base}b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta) + np.cos(theta / 2)  # Quantum harmonic combinations
    wave = y + z + QUANTUM_ADJUSTMENT  # Combine quantum harmonics
    wave_transformed = np.abs(fft(wave))  # Use Fourier transform to capture unique structure
    return wave_transformed, binary_data

# Step 2: Reflect/Refract for Macro Storage
def reflect_to_macro(quantum_wave):
    """
    Reflects the quantum wave into a macro-compatible binary format for H storage.
    """
    threshold = np.median(quantum_wave)  # Use median as a dynamic threshold
    binary_data = np.array([1 if val > threshold else 0 for val in quantum_wave], dtype=np.uint8)
    return binary_data

# Step 3: Store Quantum Wave as Macro Data into H
def store_in_H(binary_data, expansion_factor=EXPANSION_FACTOR):
    """
    Encodes binary data into H array with cumulative storage.
    """
    harmonics = np.cumsum(binary_data.astype(np.float64) * expansion_factor)
    return harmonics

# Step 4: Retrieve Binary Data from H
def retrieve_from_H(harmonics, expansion_factor=EXPANSION_FACTOR):
    """
    Retrieves binary data from H array.
    """
    first_value = harmonics[0] / expansion_factor
    reversed_data = np.diff(harmonics) / expansion_factor
    reversed_data = np.insert(reversed_data, 0, first_value)
    return np.round(reversed_data).astype(np.uint8)

# Input hash value
hash_value = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"

# Generate quantum wave
quantum_wave, binary_string = quantum_wave_from_hash(hash_value)

# Reflect quantum wave to macro binary
macro_binary_data = reflect_to_macro(quantum_wave)

# Store reflected data into H
harmonics = store_in_H(macro_binary_data)

# Retrieve binary data from H
retrieved_binary_data = retrieve_from_H(harmonics)

# Validate the process
original_binary_string = ''.join(map(str, macro_binary_data))
retrieved_binary_string = ''.join(map(str, retrieved_binary_data))

print("Original Macro Binary Data (First 100 bits):", original_binary_string[:100])
print("Retrieved Macro Binary Data (First 100 bits):", retrieved_binary_string[:100])

if np.array_equal(macro_binary_data, retrieved_binary_data):
    print("Data matches successfully!")
else:
    print("Data mismatch detected.")

# Visualize the harmonic storage
plt.figure(figsize=(10, 6))
plt.plot(harmonics, label="H(n) Harmonics", color='blue')
plt.title("Harmonic Storage Representation")
plt.xlabel("Index")
plt.ylabel("H(n)")
plt.legend()
plt.grid()
plt.show()

# Dump the H array
print("H array values:", harmonics)
