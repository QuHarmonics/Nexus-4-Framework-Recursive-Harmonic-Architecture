import numpy as np
import matplotlib.pyplot as plt

# Constants
EXPANSION_FACTOR = 1.5
QUANTUM_ADJUSTMENT = 0.35

# Step 1: H Mechanics (Macro to Quantum and Back)
def store_in_H(binary_data, expansion_factor=EXPANSION_FACTOR):
    harmonics = np.cumsum(binary_data.astype(np.float64) * expansion_factor)
    return harmonics

def retrieve_from_H(harmonics, expansion_factor=EXPANSION_FACTOR):
    first_value = harmonics[0] / expansion_factor
    reversed_data = np.diff(harmonics) / expansion_factor
    reversed_data = np.insert(reversed_data, 0, first_value)
    return np.round(reversed_data).astype(np.uint8)

# Step 2: Samson Wave Generation (Quantum Creation)
def quantum_wave_samson(hash_value):
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
    threshold = np.mean(found_wave)
    binary_data = np.array([1 if point > threshold else 0 for point in found_wave], dtype=np.uint8)
    return binary_data

# Step 4: Validation of H Mechanics
def validate_H_mechanics(seed_binary):
    harmonics = store_in_H(seed_binary)
    retrieved_binary = retrieve_from_H(harmonics)
    data_matches = np.array_equal(seed_binary, retrieved_binary)

    # Debugging: Print actual data
    print("Seed Binary (First 100 bits):", ''.join(map(str, seed_binary[:100])))
    print("Retrieved Binary (First 100 bits):", ''.join(map(str, retrieved_binary[:100])))
    return harmonics, retrieved_binary, data_matches

# Input hash value for quantum wave generation
hash_value = "cd2eca3535741f27a8ae40c31b0c41d4057a7a7b912b33b9aed86485d1c84676"  # SHA-256 for "hello"

# Generate Quantum Wave using Samson
x_q, y_q, z_q = quantum_wave_samson(hash_value)

# Input binary data from seed + padding
def sha_padding(seed):
    seed_binary = ''.join(format(ord(c), '08b') for c in seed)
    seed_binary += '1'
    while len(seed_binary) % 512 != 448:
        seed_binary += '0'
    seed_length = len(seed) * 8
    seed_binary += format(seed_length, '064b')
    return np.array([int(bit) for bit in seed_binary], dtype=np.uint8)

seed = "hello"
seed_binary = sha_padding(seed)

# Change the seed to test independence
seed_alt = "world"  # New test seed
seed_binary_alt = sha_padding(seed_alt)

# Validate H Mechanics for both seeds
print("Testing Original Seed...")
harmonics, retrieved_binary, data_matches = validate_H_mechanics(seed_binary)

print("\nTesting Alternative Seed...")
harmonics_alt, retrieved_binary_alt, data_matches_alt = validate_H_mechanics(seed_binary_alt)

# Visualize the Waves and Compare
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Plot Quantum Wave (Samson)
ax.plot(x_q, y_q, z_q, label="Quantum Wave (Samson)", color='blue')
ax.scatter(x_q, y_q, z_q, color='red', s=5)

# Labels and Visualization
ax.set_title("Comparison of Quantum and Macro Waves")
ax.set_xlabel("X-axis")
ax.set_ylabel("Y-axis")
ax.set_zlabel("Z-axis")
ax.legend()
plt.show()

# Results
print("Original Seed Binary Matches:", data_matches)
print("Alternative Seed Binary Matches:", data_matches_alt)
