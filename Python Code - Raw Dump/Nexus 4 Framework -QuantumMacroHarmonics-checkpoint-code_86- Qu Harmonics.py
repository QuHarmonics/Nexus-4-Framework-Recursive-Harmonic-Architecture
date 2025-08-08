import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Constants for Mark1 framework
EXPANSION_FACTOR = 1.5
QUANTUM_ADJUSTMENT = 0.35

# Step 1: Generate a Quantum Wave from Hash
def quantum_wave_from_hash(hash_value, base=2):
    binary_data = ''.join(format(int(char, 16), f'0{base}b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta) + np.cos(theta / 2)
    wave = y + z + QUANTUM_ADJUSTMENT
    return wave, binary_data

# Step 2: Reflect Quantum Wave to Macro-compatible H OUT
def reflect_to_macro(quantum_wave):
    threshold = np.median(quantum_wave)
    binary_data = np.array([1 if val > threshold else 0 for val in quantum_wave], dtype=np.uint8)
    return binary_data

# Step 3: Store Quantum Wave in H OUT
def store_in_H_OUT(binary_data, expansion_factor=EXPANSION_FACTOR):
    harmonics = np.cumsum(binary_data.astype(np.float64) * expansion_factor)
    return harmonics

# Step 4: Retrieve Macro Hash from H OUT (Reversed H)
def retrieve_macro_from_H(harmonics, expansion_factor=EXPANSION_FACTOR):
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

# Store reflected data into H OUT
harmonics = store_in_H_OUT(macro_binary_data)

# Retrieve macro hash from H OUT
retrieved_macro_data = retrieve_macro_from_H(harmonics)

# Validate the process
original_binary_string = ''.join(map(str, macro_binary_data))
retrieved_binary_string = ''.join(map(str, retrieved_macro_data))

print("Original Macro Binary Data (First 100 bits):", original_binary_string[:100])
print("Retrieved Macro Binary Data (First 100 bits):", retrieved_binary_string[:100])

if np.array_equal(macro_binary_data, retrieved_macro_data):
    print("Data matches successfully!")
else:
    print("Data mismatch detected.")

# Visualize Quantum Wave
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
x = np.arange(len(quantum_wave))
y = quantum_wave
z = np.sin(x / 10.0)
ax.plot(x, y, z, label="Quantum Wave", color='blue')
ax.scatter(x, y, z, color='red', s=5, label="Nodes")
ax.set_title("3D Visualization of Quantum Wave")
ax.set_xlabel("Index")
ax.set_ylabel("Amplitude")
ax.set_zlabel("Z-axis Wave")
ax.legend()
plt.show()

# Visualize H OUT storage
plt.figure(figsize=(10, 6))
plt.plot(harmonics, label="H(n) Harmonics", color='green')
plt.title("Harmonic Storage Representation")
plt.xlabel("Index")
plt.ylabel("H(n)")
plt.legend()
plt.grid()
plt.show()

# Dump the H array
print("H array values:", harmonics)
