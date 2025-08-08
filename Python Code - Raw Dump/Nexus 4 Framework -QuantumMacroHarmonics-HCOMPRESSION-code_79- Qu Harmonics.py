import numpy as np
import matplotlib.pyplot as plt

# Constants
EXPANSION_FACTOR = 1.5
QUANTUM_ADJUSTMENT = 0.35

# Step 1: H Mechanics (Macro to Quantum and Back)
def store_in_H(binary_data, expansion_factor=EXPANSION_FACTOR):
    """
    Encodes macro binary data into the H array as a quantum-compatible representation.
    """
    harmonics = np.cumsum(binary_data.astype(np.float64) * expansion_factor)
    return harmonics

def retrieve_from_H(harmonics, expansion_factor=EXPANSION_FACTOR):
    """
    Decodes the H array back into macro binary data.
    """
    first_value = harmonics[0] / expansion_factor
    reversed_data = np.diff(harmonics) / expansion_factor
    reversed_data = np.insert(reversed_data, 0, first_value)
    return np.round(reversed_data).astype(np.uint8)

# Step 2: Samson Wave Generation (Quantum Creation)
def quantum_wave_samson(hash_value):
    """
    Generates a quantum waveform from a given hash using Samson principles.
    """
    binary_data = ''.join(format(int(char, 16), '04b') for char in hash_value)
    n = len(binary_data)
    theta = np.linspace(0, 2 * np.pi, n)
    radius = np.linspace(0, 1, n)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    z = np.sin(2 * theta) + np.cos(theta / 2) + QUANTUM_ADJUSTMENT
    return x, y, z

# Step 3: Collapse Found Wave (Quantum to Macro)
def found_wave_to_binary(found_wave):
    """
    Collapses a found quantum waveform into binary data for macro representation.
    """
    threshold = np.mean(found_wave)
    binary_data = np.array([1 if point > threshold else 0 for point in found_wave], dtype=np.uint8)
    return binary_data

# Step 4: Validation of H Mechanics
def validate_H_mechanics(seed_binary):
    """
    Validates the 1:1 storage and retrieval mechanics of H using known macro data.
    """
    # Store in H
    harmonics = store_in_H(seed_binary)
    # Retrieve from H
    retrieved_binary = retrieve_from_H(harmonics)
    # Check fidelity
    data_matches = np.array_equal(seed_binary, retrieved_binary)
    return harmonics, retrieved_binary, data_matches

# Testing with Known Data
# Input hash value for quantum wave generation
hash_value = "cd2eca3535741f27a8ae40c31b0c41d4057a7a7b912b33b9aed86485d1c84676"  # Corrected SHA-256 hash for "hello"

# Generate Quantum Wave using Samson
x_q, y_q, z_q = quantum_wave_samson(hash_value)

# Input binary data from seed + padding
def sha_padding(seed):
    """
    Generates the binary representation of the seed with SHA-like padding.
    """
    seed_binary = ''.join(format(ord(c), '08b') for c in seed)
    seed_binary += '1'
    while len(seed_binary) % 512 != 448:
        seed_binary += '0'
    seed_length = len(seed) * 8
    seed_binary += format(seed_length, '064b')
    return np.array([int(bit) for bit in seed_binary], dtype=np.uint8)

seed = "hellodave"
seed_binary = sha_padding(seed)

# Validate H Mechanics
harmonics, retrieved_binary, data_matches = validate_H_mechanics(seed_binary)

# Generate Macro Wave for Seed + Padding Binary
x_m, y_m, z_m = quantum_wave_samson(''.join(map(str, seed_binary)))

# Visualize the Waves and Compare
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot Quantum Wave (Samson)
ax.plot(x_q, y_q, z_q, label="Quantum Wave (Samson)", color='blue')
ax.scatter(x_q, y_q, z_q, color='red', s=5)

# Plot Macro Wave (Seed Binary)
ax.plot(x_m, y_m, z_m, label="Macro Wave (Binary)", color='green')
ax.scatter(x_m, y_m, z_m, color='yellow', s=5)

# Labels and Visualization
ax.set_title("Comparison of Quantum and Macro Waves")
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")
ax.legend()
plt.show()

# Output Results
print("Seed Binary (First 100 bits):", ''.join(map(str, seed_binary[:100])))
print("Retrieved Binary (First 100 bits):", ''.join(map(str, retrieved_binary[:100])))
print("Data Matches:", data_matches)
